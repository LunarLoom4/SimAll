// =============================================================================
// SimAll Beta - Morphing Subsystem
// File   : src/morphing/FfdBox.cpp
// =============================================================================
#include "morphing/FfdBox.hpp"

#include <algorithm>
#include <cmath>

namespace simall::morphing
{

namespace
{

// Bernstein polynomial B_i^n(t) = C(n,i) · t^i · (1-t)^(n-i).
double bernstein(std::uint32_t i, std::uint32_t n, double t)
{
    if (n == 0)
        return 1.0;
    // Compute binomial coefficient iteratively (n ≤ ~8 in typical FFD).
    double C = 1.0;
    for (std::uint32_t k = 0; k < i; ++k)
        C = C * double(n - k) / double(k + 1);
    return C * std::pow(t, double(i)) * std::pow(1.0 - t, double(n - i));
}

} // namespace

FfdBox::FfdBox(
    const Vec3& origin, const Vec3& size, std::uint32_t nu, std::uint32_t nv, std::uint32_t nw)
    : origin_(origin)
    , size_(size)
    , nu_(std::max(nu, 2u))
    , nv_(std::max(nv, 2u))
    , nw_(std::max(nw, 2u))
    , P_(std::size_t(nu_) * nv_ * nw_)
{
    for (std::uint32_t k = 0; k < nw_; ++k) {
        for (std::uint32_t j = 0; j < nv_; ++j) {
            for (std::uint32_t i = 0; i < nu_; ++i) {
                Vec3 p;
                p.x = origin_.x + (double(i) / double(nu_ - 1)) * size_.x;
                p.y = origin_.y + (double(j) / double(nv_ - 1)) * size_.y;
                p.z = origin_.z + (double(k) / double(nw_ - 1)) * size_.z;
                P_[idx(i, j, k)] = p;
            }
        }
    }
    P0_ = P_;
}

const Vec3& FfdBox::control_point(std::uint32_t i, std::uint32_t j, std::uint32_t k) const
{
    return P_[idx(i, j, k)];
}
void FfdBox::set_control_point(std::uint32_t i, std::uint32_t j, std::uint32_t k, const Vec3& p)
{
    P_[idx(i, j, k)] = p;
}
void FfdBox::displace_control_point(std::uint32_t i,
                                    std::uint32_t j,
                                    std::uint32_t k,
                                    const Vec3& d)
{
    auto& p = P_[idx(i, j, k)];
    p.x += d.x;
    p.y += d.y;
    p.z += d.z;
}
void FfdBox::reset()
{
    P_ = P0_;
}

bool FfdBox::world_to_lattice(const Vec3& w, Vec3& stu) const
{
    const double s = (w.x - origin_.x) / std::max(size_.x, 1e-30);
    const double t = (w.y - origin_.y) / std::max(size_.y, 1e-30);
    const double u = (w.z - origin_.z) / std::max(size_.z, 1e-30);
    stu = {s, t, u};
    return s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0;
}

Vec3 FfdBox::evaluate(double s, double t, double u) const
{
    // Tensor product Bernstein basis.  Precompute 1-D bases.
    std::vector<double> Bs(nu_), Bt(nv_), Bu(nw_);
    for (std::uint32_t i = 0; i < nu_; ++i)
        Bs[i] = bernstein(i, nu_ - 1, s);
    for (std::uint32_t j = 0; j < nv_; ++j)
        Bt[j] = bernstein(j, nv_ - 1, t);
    for (std::uint32_t k = 0; k < nw_; ++k)
        Bu[k] = bernstein(k, nw_ - 1, u);
    Vec3 out{};
    for (std::uint32_t k = 0; k < nw_; ++k) {
        for (std::uint32_t j = 0; j < nv_; ++j) {
            for (std::uint32_t i = 0; i < nu_; ++i) {
                const double w = Bs[i] * Bt[j] * Bu[k];
                const auto& P = P_[idx(i, j, k)];
                out.x += w * P.x;
                out.y += w * P.y;
                out.z += w * P.z;
            }
        }
    }
    return out;
}

void FfdBox::morph(const std::vector<Vec3>& in, std::vector<Vec3>& out) const
{
    if (&in != &out)
        out.assign(in.size(), Vec3{});
    for (std::size_t p = 0; p < in.size(); ++p) {
        Vec3 stu;
        if (world_to_lattice(in[p], stu)) {
            out[p] = evaluate(stu.x, stu.y, stu.z);
        } else if (&in != &out) {
            out[p] = in[p];
        }
    }
}

} // namespace simall::morphing
