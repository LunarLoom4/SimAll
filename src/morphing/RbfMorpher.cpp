// =============================================================================
// SimAll Beta - Morphing Subsystem
// File   : src/morphing/RbfMorpher.cpp
// =============================================================================
#include "morphing/RbfMorpher.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace simall::morphing
{

namespace
{

inline double wendland_c2(double r, double R)
{
    if (r >= R)
        return 0.0;
    const double t = 1.0 - r / R;
    const double t2 = t * t;
    return t2 * t2 * (4.0 * r / R + 1.0);
}

inline double dist(const util::Vec3d& a, const util::Vec3d& b)
{
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/// Dense LDL^T factorisation in place on A (n×n symmetric).
/// Stores D on the diagonal, L (unit-diagonal) strictly below.
bool ldlt(std::vector<double>& A, int n)
{
    for (int j = 0; j < n; ++j) {
        double D = A[j * n + j];
        for (int k = 0; k < j; ++k)
            D -= A[j * n + k] * A[j * n + k] * A[k * n + k];
        if (std::abs(D) < 1e-30)
            return false;
        A[j * n + j] = D;
        for (int i = j + 1; i < n; ++i) {
            double s = A[i * n + j];
            for (int k = 0; k < j; ++k)
                s -= A[i * n + k] * A[j * n + k] * A[k * n + k];
            A[i * n + j] = s / D;
        }
    }
    return true;
}

void ldlt_solve(const std::vector<double>& A, int n, std::vector<double>& b)
{
    // L y = b
    for (int i = 0; i < n; ++i) {
        double s = b[i];
        for (int k = 0; k < i; ++k)
            s -= A[i * n + k] * b[k];
        b[i] = s;
    }
    // D z = y
    for (int i = 0; i < n; ++i)
        b[i] /= A[i * n + i];
    // L^T x = z
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int k = i + 1; k < n; ++k)
            s -= A[k * n + i] * b[k];
        b[i] = s;
    }
}

} // namespace

void RbfMorpher::set_controls(const std::vector<util::Vec3d>& p, const std::vector<util::Vec3d>& d)
{
    if (p.size() != d.size())
        throw std::invalid_argument("RbfMorpher: controls/displacements size mismatch");
    ctrl_ = p;
    disp_ = d;
    trained_ = false;
}

void RbfMorpher::train(RbfOptions opt)
{
    const int n = static_cast<int>(ctrl_.size());
    havePoly_ = opt.addPolynomial;
    const int m = havePoly_ ? 4 : 0;
    const int N = n + m;
    if (n == 0) {
        trained_ = true;
        return;
    }

    // Auto-set support radius from bbox if user gave 0.
    if (opt.supportRadius <= 0.0) {
        util::Vec3d mn = ctrl_[0], mx = ctrl_[0];
        for (const auto& c : ctrl_) {
            mn.x = std::min(mn.x, c.x);
            mx.x = std::max(mx.x, c.x);
            mn.y = std::min(mn.y, c.y);
            mx.y = std::max(mx.y, c.y);
            mn.z = std::min(mn.z, c.z);
            mx.z = std::max(mx.z, c.z);
        }
        const double dx = mx.x - mn.x, dy = mx.y - mn.y, dz = mx.z - mn.z;
        R_ = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (R_ <= 0.0)
            R_ = 1.0;
    } else
        R_ = opt.supportRadius;

    // Assemble dense saddle-point matrix:
    //  [ Φ + λI   P ]
    //  [   P^T    0 ]
    std::vector<double> A(N * N, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i * N + j] = wendland_c2(dist(ctrl_[i], ctrl_[j]), R_);
        }
        A[i * N + i] += opt.regularization;
        if (havePoly_) {
            A[i * N + n + 0] = 1.0;
            A[i * N + n + 1] = ctrl_[i].x;
            A[i * N + n + 2] = ctrl_[i].y;
            A[i * N + n + 3] = ctrl_[i].z;
            A[(n + 0) * N + i] = 1.0;
            A[(n + 1) * N + i] = ctrl_[i].x;
            A[(n + 2) * N + i] = ctrl_[i].y;
            A[(n + 3) * N + i] = ctrl_[i].z;
        }
    }
    // Tiny ridge on the polynomial block to keep LDL stable
    if (havePoly_)
        for (int k = 0; k < 4; ++k)
            A[(n + k) * N + n + k] += opt.regularization;

    if (!ldlt(A, N))
        throw std::runtime_error("RbfMorpher: LDL^T factorisation failed");

    auto solveComp = [&](int comp, std::vector<double>& alpha, double poly[4]) {
        std::vector<double> rhs(N, 0.0);
        for (int i = 0; i < n; ++i)
            rhs[i] = (comp == 0 ? disp_[i].x : comp == 1 ? disp_[i].y : disp_[i].z);
        ldlt_solve(A, N, rhs);
        alpha.assign(rhs.begin(), rhs.begin() + n);
        if (havePoly_)
            for (int k = 0; k < 4; ++k)
                poly[k] = rhs[n + k];
    };
    solveComp(0, ax_, bx_);
    solveComp(1, ay_, by_);
    solveComp(2, az_, bz_);
    trained_ = true;
    SIMALL_LOG_INFO("Morph", "RBF trained: ", n, " ctrl pts, R=", R_);
}

util::Vec3d RbfMorpher::evaluate(const util::Vec3d& x) const
{
    util::Vec3d s{0, 0, 0};
    if (!trained_)
        return s;
    for (std::size_t i = 0; i < ctrl_.size(); ++i) {
        const double w = wendland_c2(dist(x, ctrl_[i]), R_);
        s.x += w * ax_[i];
        s.y += w * ay_[i];
        s.z += w * az_[i];
    }
    if (havePoly_) {
        s.x += bx_[0] + bx_[1] * x.x + bx_[2] * x.y + bx_[3] * x.z;
        s.y += by_[0] + by_[1] * x.x + by_[2] * x.y + by_[3] * x.z;
        s.z += bz_[0] + bz_[1] * x.x + bz_[2] * x.y + bz_[3] * x.z;
    }
    return s;
}

void RbfMorpher::apply(meshing::Mesh& mesh) const
{
    if (!trained_)
        return;
    auto& N = mesh.nodes();
    const std::size_t nN = N.size();
    for (std::size_t i = 0; i < nN; ++i) {
        const util::Vec3d s = evaluate({N.x[i], N.y[i], N.z[i]});
        N.x[i] += s.x;
        N.y[i] += s.y;
        N.z[i] += s.z;
    }
    mesh.compute_geometry();
}

} // namespace simall::morphing
