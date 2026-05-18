// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/LinearSolvers.cpp
//
// Krylov methods + ILU(0) / Jacobi preconditioners. Implementations follow
// Saad ("Iterative Methods for Sparse Linear Systems", 2nd ed).
// =============================================================================
#include "solver/LinearSolvers.hpp"
#include "solver/AMGPreconditioner.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace simall::solver {

namespace {

using Vec = util::aligned_vector<double>;

double dot(const Vec& a, const Vec& b) {
    double s = 0; const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}
double norm2(const Vec& a) { return std::sqrt(dot(a, a)); }

void axpy(double alpha, const Vec& x, Vec& y) {
    const std::size_t n = x.size();
    for (std::size_t i = 0; i < n; ++i) y[i] += alpha * x[i];
}
void scal(double a, Vec& x) {
    for (double& v : x) v *= a;
}

// ----------------------------------------------------------------- Jacobi
class JacobiPreconditioner final : public IPreconditioner {
public:
    void setup(const CSRMatrix& A) override {
        const std::size_t n = A.rows();
        invDiag_.assign(n, 1.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (int k = A.rowPtr[i]; k < A.rowPtr[i+1]; ++k)
                if (A.colIdx[k] == static_cast<int>(i)) {
                    invDiag_[i] = (A.values[k] != 0.0) ? 1.0 / A.values[k] : 1.0;
                    break;
                }
        }
    }
    void apply(const Vec& r, Vec& z) const override {
        z.resize(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) z[i] = invDiag_[i] * r[i];
    }
private:
    Vec invDiag_;
};

// ----------------------------------------------------------------- ILU(0)
class ILU0Preconditioner final : public IPreconditioner {
public:
    void setup(const CSRMatrix& A) override {
        // Zero-fill incomplete LU: factorize A in place on its own sparsity
        // pattern. L is unit-lower, U is upper; both stored in `lu_`.
        rowPtr_ = A.rowPtr;
        colIdx_ = A.colIdx;
        lu_     = A.values;
        const std::size_t n = A.rows();

        std::vector<int> diagPtr(n, -1);
        for (std::size_t i = 0; i < n; ++i)
            for (int k = rowPtr_[i]; k < rowPtr_[i+1]; ++k)
                if (colIdx_[k] == static_cast<int>(i)) { diagPtr[i] = k; break; }

        for (std::size_t i = 1; i < n; ++i) {
            for (int kp = rowPtr_[i]; kp < rowPtr_[i+1]; ++kp) {
                const int k = colIdx_[kp];
                if (k >= static_cast<int>(i)) break;
                if (diagPtr[k] < 0) continue;
                const double piv = lu_[diagPtr[k]];
                if (piv == 0.0) continue;
                const double m = lu_[kp] / piv;
                lu_[kp] = m;
                // Update row i ← row i - m * row k (only on existing pattern).
                for (int jp = diagPtr[k] + 1; jp < rowPtr_[k+1]; ++jp) {
                    const int j = colIdx_[jp];
                    // find (i, j) in row i
                    for (int qp = kp + 1; qp < rowPtr_[i+1]; ++qp) {
                        if (colIdx_[qp] == j) { lu_[qp] -= m * lu_[jp]; break; }
                        if (colIdx_[qp] > j) break;
                    }
                }
            }
        }
        diagPtr_ = std::move(diagPtr);
    }
    void apply(const Vec& r, Vec& z) const override {
        const std::size_t n = r.size();
        Vec y(n, 0.0);
        // Forward: L y = r  (L unit-lower)
        for (std::size_t i = 0; i < n; ++i) {
            double s = r[i];
            for (int kp = rowPtr_[i]; kp < rowPtr_[i+1]; ++kp) {
                const int j = colIdx_[kp];
                if (j >= static_cast<int>(i)) break;
                s -= lu_[kp] * y[j];
            }
            y[i] = s;
        }
        // Backward: U z = y
        z.assign(n, 0.0);
        for (std::ptrdiff_t i = n - 1; i >= 0; --i) {
            double s = y[i];
            const int kd = diagPtr_[i];
            for (int kp = kd + 1; kp < rowPtr_[i+1]; ++kp)
                s -= lu_[kp] * z[colIdx_[kp]];
            z[i] = (kd >= 0 && lu_[kd] != 0.0) ? s / lu_[kd] : s;
        }
    }
private:
    std::vector<int>    rowPtr_, colIdx_, diagPtr_;
    Vec                 lu_;
};

// ----------------------------------------------------------------- SSOR
// Symmetric Successive Over-Relaxation:
//      M = (1/(ω(2-ω))) (D + ωL) D^-1 (D + ωU)
// Applying M^-1 r consists of a forward sweep + scale + backward sweep.
// ω = 1 reduces to Symmetric Gauss-Seidel.
class SsorPreconditioner final : public IPreconditioner {
public:
    explicit SsorPreconditioner(double omega = 1.0) : omega_(omega) {}

    void setup(const CSRMatrix& A) override {
        rowPtr_ = A.rowPtr;
        colIdx_ = A.colIdx;
        values_ = A.values;
        const std::size_t n = A.rows();
        diag_   .assign(n, 1.0);
        invDiag_.assign(n, 1.0);
        diagPtr_.assign(n, -1);
        for (std::size_t i = 0; i < n; ++i) {
            for (int k = rowPtr_[i]; k < rowPtr_[i+1]; ++k) {
                if (colIdx_[k] == static_cast<int>(i)) {
                    diagPtr_[i] = k;
                    diag_[i]    = values_[k];
                    invDiag_[i] = (values_[k] != 0.0) ? 1.0 / values_[k] : 1.0;
                    break;
                }
            }
        }
    }
    void apply(const Vec& r, Vec& z) const override {
        const std::size_t n = r.size();
        Vec y(n, 0.0);
        // Forward sweep: solve (D + ωL) y = ω r  → y_i = ω * (r_i - Σ_{j<i} ω L_ij y_j) / D_i
        // L is the strict lower triangle.
        // Here we use (D + ωL) y = ω r equivalently y_i = (ω r_i - Σ_{j<i} ω a_ij y_j) / a_ii
        for (std::size_t i = 0; i < n; ++i) {
            double s = omega_ * r[i];
            for (int kp = rowPtr_[i]; kp < rowPtr_[i+1]; ++kp) {
                const int j = colIdx_[kp];
                if (j >= static_cast<int>(i)) break;
                s -= omega_ * values_[kp] * y[j];
            }
            y[i] = s * invDiag_[i];
        }
        // Scale by D / ω(2-ω)   →   ŷ_i = (2-ω) * D_i * y_i / ω · (something)
        // Standard formulation: SSOR(ω) preconditioner M = (1/(2-ω)) (D/ω + L) (D/ω)^-1 (D/ω + U).
        // → applying M^-1 r:  intermediate y solves (D/ω + L) y = r, then v = D y,
        //                     then solve (D/ω + U) z = v, then scale by (2-ω).
        // Equivalent formula: y_forward → multiply by D / ω → backward sweep.
        // Re-scale:
        Vec v(n, 0.0);
        const double scale = (2.0 - omega_);
        for (std::size_t i = 0; i < n; ++i) v[i] = diag_[i] * y[i] / omega_;
        // Backward sweep: (D/ω + U) z = v  →  z_i = ω * (v_i - Σ_{j>i} a_ij z_j) / D_i
        z.assign(n, 0.0);
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 1; i >= 0; --i) {
            double s = v[i];
            for (int kp = rowPtr_[i]; kp < rowPtr_[i+1]; ++kp) {
                const int j = colIdx_[kp];
                if (j <= static_cast<int>(i)) continue;
                s -= values_[kp] * z[j];
            }
            z[i] = s * omega_ * invDiag_[i];
        }
        // Final scale: M^-1 ↔ multiply by (2-ω).
        for (double& zv : z) zv *= scale;
    }

private:
    double               omega_ = 1.0;
    std::vector<int>     rowPtr_, colIdx_, diagPtr_;
    Vec                  values_, diag_, invDiag_;
};

}  // namespace

std::unique_ptr<IPreconditioner> make_preconditioner(PreconditionerKind k) {
    switch (k) {
        case PreconditionerKind::Jacobi: return std::make_unique<JacobiPreconditioner>();
        case PreconditionerKind::ILU:    return std::make_unique<ILU0Preconditioner>();
        case PreconditionerKind::SSOR:   return std::make_unique<SsorPreconditioner>(1.0);
        case PreconditionerKind::AMG:    return std::make_unique<AMGPreconditioner>();
        default:                          return std::make_unique<JacobiPreconditioner>();
    }
}

// ============================================================ GMRES(m)
namespace {

class GMRESmSolver final : public ILinearSolver {
public:
    explicit GMRESmSolver(LinearSolverConfig c)
        : cfg_(c), M_(make_preconditioner(c.preconditioner)) {}

    int solve(const CSRMatrix& A, const Vec& b, Vec& x) override {
        M_->setup(A);
        const std::size_t n = A.rows();
        x.assign(n, 0.0);
        const int m = std::max(1, cfg_.restart);
        Vec r(n), w(n), z(n);
        std::vector<Vec> V(m + 1, Vec(n, 0.0));
        std::vector<std::vector<double>> H(m + 1, std::vector<double>(m, 0.0));
        std::vector<double> cs(m, 0.0), sn(m, 0.0), s(m + 1, 0.0);

        const double bnorm = std::max(norm2(b), 1e-30);
        int iter = 0;
        for (int outer = 0; outer < cfg_.maxIterations; ++outer) {
            A.spmv(x, w);
            for (std::size_t i = 0; i < n; ++i) r[i] = b[i] - w[i];
            M_->apply(r, z);
            double beta = norm2(z);
            if (beta < cfg_.tolerance * bnorm) { last_ = beta / bnorm; return iter; }
            V[0] = z; scal(1.0 / beta, V[0]);
            std::fill(s.begin(), s.end(), 0.0);
            s[0] = beta;

            int j = 0;
            for (j = 0; j < m && iter < cfg_.maxIterations; ++j, ++iter) {
                A.spmv(V[j], w);
                M_->apply(w, z);                       // z = M^-1 A v_j
                for (int i = 0; i <= j; ++i) {
                    H[i][j] = dot(V[i], z);
                    axpy(-H[i][j], V[i], z);
                }
                H[j+1][j] = norm2(z);
                if (H[j+1][j] != 0.0) {
                    V[j+1] = z; scal(1.0 / H[j+1][j], V[j+1]);
                }
                // Givens rotations
                for (int i = 0; i < j; ++i) {
                    const double t = cs[i]*H[i][j] + sn[i]*H[i+1][j];
                    H[i+1][j] = -sn[i]*H[i][j] + cs[i]*H[i+1][j];
                    H[i][j]   = t;
                }
                const double dd = std::hypot(H[j][j], H[j+1][j]);
                cs[j] = H[j][j] / dd; sn[j] = H[j+1][j] / dd;
                H[j][j]   = dd; H[j+1][j] = 0.0;
                const double t = cs[j]*s[j];
                s[j+1] = -sn[j]*s[j]; s[j] = t;
                last_  = std::abs(s[j+1]) / bnorm;
                if (last_ < cfg_.tolerance) { ++j; break; }
            }
            // Solve upper-triangular H y = s, then x += V y
            std::vector<double> y(j, 0.0);
            for (int i = j - 1; i >= 0; --i) {
                double sum = s[i];
                for (int k = i + 1; k < j; ++k) sum -= H[i][k] * y[k];
                y[i] = sum / H[i][i];
            }
            for (int i = 0; i < j; ++i) axpy(y[i], V[i], x);
            if (last_ < cfg_.tolerance) return iter;
        }
        return iter;
    }
    double last_residual() const override { return last_; }
private:
    LinearSolverConfig                cfg_;
    std::unique_ptr<IPreconditioner>  M_;
    double                            last_ = 0.0;
};

// ============================================================ BiCGSTAB
class BiCGSTABSolver final : public ILinearSolver {
public:
    explicit BiCGSTABSolver(LinearSolverConfig c)
        : cfg_(c), M_(make_preconditioner(c.preconditioner)) {}

    int solve(const CSRMatrix& A, const Vec& b, Vec& x) override {
        M_->setup(A);
        const std::size_t n = A.rows();
        x.assign(n, 0.0);
        Vec r(n), rh(n), p(n, 0.0), v(n, 0.0), s(n), t(n), y(n), z(n);
        A.spmv(x, v);
        for (std::size_t i = 0; i < n; ++i) r[i] = b[i] - v[i];
        rh = r;
        double rho = 1.0, alpha = 1.0, omega = 1.0;
        const double bnorm = std::max(norm2(b), 1e-30);
        std::fill(v.begin(), v.end(), 0.0);
        for (int it = 0; it < cfg_.maxIterations; ++it) {
            double rhoNew = dot(rh, r);
            if (rhoNew == 0.0) { last_ = norm2(r)/bnorm; return it; }
            const double beta = (rhoNew/rho) * (alpha/omega);
            for (std::size_t i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega*v[i]);
            M_->apply(p, y); A.spmv(y, v);
            alpha = rhoNew / dot(rh, v);
            for (std::size_t i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
            if (norm2(s) < cfg_.tolerance * bnorm) {
                axpy(alpha, y, x); last_ = norm2(s)/bnorm; return it+1;
            }
            M_->apply(s, z); A.spmv(z, t);
            omega = dot(t, s) / std::max(dot(t, t), 1e-30);
            for (std::size_t i = 0; i < n; ++i) {
                x[i] += alpha*y[i] + omega*z[i];
                r[i] = s[i] - omega*t[i];
            }
            rho = rhoNew;
            last_ = norm2(r)/bnorm;
            if (last_ < cfg_.tolerance) return it+1;
        }
        return cfg_.maxIterations;
    }
    double last_residual() const override { return last_; }
private:
    LinearSolverConfig                cfg_;
    std::unique_ptr<IPreconditioner>  M_;
    double                            last_ = 0.0;
};

// ============================================================ CG (SPD)
class CGSolver final : public ILinearSolver {
public:
    explicit CGSolver(LinearSolverConfig c)
        : cfg_(c), M_(make_preconditioner(c.preconditioner)) {}

    int solve(const CSRMatrix& A, const Vec& b, Vec& x) override {
        M_->setup(A);
        const std::size_t n = A.rows();
        x.assign(n, 0.0);
        Vec r = b, z(n), p(n), Ap(n);
        M_->apply(r, z); p = z;
        double rz = dot(r, z);
        const double bnorm = std::max(norm2(b), 1e-30);
        for (int it = 0; it < cfg_.maxIterations; ++it) {
            A.spmv(p, Ap);
            const double pAp = dot(p, Ap);
            if (pAp == 0) { last_ = std::sqrt(rz) / bnorm; return it; }
            const double alpha = rz / pAp;
            axpy(alpha, p, x);
            axpy(-alpha, Ap, r);
            const double rn = norm2(r);
            last_ = rn / bnorm;
            if (last_ < cfg_.tolerance) return it + 1;
            M_->apply(r, z);
            const double rzNew = dot(r, z);
            const double beta  = rzNew / rz;
            for (std::size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
            rz = rzNew;
        }
        return cfg_.maxIterations;
    }
    double last_residual() const override { return last_; }
private:
    LinearSolverConfig                cfg_;
    std::unique_ptr<IPreconditioner>  M_;
    double                            last_ = 0.0;
};

}  // namespace

std::unique_ptr<ILinearSolver> make_gmres   (LinearSolverConfig cfg) { return std::make_unique<GMRESmSolver>(cfg); }
std::unique_ptr<ILinearSolver> make_bicgstab(LinearSolverConfig cfg) { return std::make_unique<BiCGSTABSolver>(cfg); }
std::unique_ptr<ILinearSolver> make_cg      (LinearSolverConfig cfg) { return std::make_unique<CGSolver>(cfg); }

// ============================================================ TFQMR
// Freund (1993) "A Transpose-Free Quasi-Minimal Residual Algorithm for
// Non-Hermitian Linear Systems", SIAM J. Sci. Comput. 14, 470-482.
//
// TFQMR is a CGS-derivative that avoids the irregular convergence of
// BiCGSTAB by applying a quasi-minimal residual smoothing.  It does NOT
// need A^T (unlike QMR / BiCG).  Per iteration cost: 2 SpMVs, 2 M^-1
// applications, comparable to BiCGSTAB but with smoother residual history
// on indefinite / strongly non-symmetric systems.
//
// Implementation follows Freund's Algorithm 5.1.  All names match the
// reference for traceability.
namespace {
class TfqmrSolver final : public ILinearSolver {
public:
    explicit TfqmrSolver(LinearSolverConfig c)
        : cfg_(c), M_(make_preconditioner(c.preconditioner)) {}

    int solve(const CSRMatrix& A, const Vec& b, Vec& x) override {
        M_->setup(A);
        const std::size_t n = A.rows();
        x.assign(n, 0.0);

        Vec r(n), rh(n), w(n), d(n, 0.0);
        Vec u1(n), u2(n), v(n), Au1(n), Au2(n);
        Vec z(n);                  // preconditioner output scratch

        // r0 = b - A x0   (x0 = 0)
        r = b;
        rh = r;                     // shadow residual
        w  = r;
        M_->apply(r, u1);           // u1 = M^-1 r
        A.spmv(u1, v);              // v  = A u1
        Au1 = v;                    // store A * u1
        double tau   = norm2(r);
        const double bnorm = std::max(norm2(b), 1.0e-30);
        if (tau / bnorm < cfg_.tolerance) { last_ = tau / bnorm; return 0; }
        double rho   = dot(rh, r);
        double theta = 0.0;
        double eta   = 0.0;
        double alpha = 0.0;

        for (int it = 0; it < cfg_.maxIterations; ++it) {
            const double sigma = dot(rh, v);
            if (std::abs(sigma) < 1.0e-30) { last_ = tau / bnorm; return it; }
            alpha = rho / sigma;
            // u2 = u1 - α M^-1 v = (M^-1 r)_{k+1/2}
            // Need: M^-1 (r - α v) but Freund's algorithm operates on
            // preconditioned vectors directly. We approximate:
            //   u2 = u1 - α * (M^-1 v)
            M_->apply(v, z);                  // z = M^-1 v
            for (std::size_t i = 0; i < n; ++i) u2[i] = u1[i] - alpha * z[i];
            A.spmv(u2, Au2);                  // Au2 = A u2

            // Two-step inner loop: m = 2k, 2k+1  (Freund's "double-step")
            // Step 1: w ← w - α A u1
            for (std::size_t i = 0; i < n; ++i) w[i] -= alpha * Au1[i];
            double theta_new = norm2(w) / std::max(tau, 1.0e-30);
            double c         = 1.0 / std::sqrt(1.0 + theta_new * theta_new);
            double tau_new   = tau * theta_new * c;
            double eta_new   = c * c * alpha;
            for (std::size_t i = 0; i < n; ++i) {
                d[i] = u1[i] + (theta * theta * eta / std::max(alpha, 1.0e-30)) * d[i];
                x[i] += eta_new * d[i];
            }
            tau   = tau_new;
            theta = theta_new;
            eta   = eta_new;
            // Step 2: w ← w - α A u2
            for (std::size_t i = 0; i < n; ++i) w[i] -= alpha * Au2[i];
            theta_new = norm2(w) / std::max(tau, 1.0e-30);
            c         = 1.0 / std::sqrt(1.0 + theta_new * theta_new);
            tau_new   = tau * theta_new * c;
            eta_new   = c * c * alpha;
            for (std::size_t i = 0; i < n; ++i) {
                d[i] = u2[i] + (theta * theta * eta / std::max(alpha, 1.0e-30)) * d[i];
                x[i] += eta_new * d[i];
            }
            tau   = tau_new;
            theta = theta_new;
            eta   = eta_new;

            // Approximate residual bound: ||r_k|| ≤ √(2k+1) τ_k
            const double resBound = tau * std::sqrt(2.0 * (it + 1) + 1.0);
            last_ = resBound / bnorm;
            if (last_ < cfg_.tolerance) {
                // Final true residual check
                A.spmv(x, z);
                for (std::size_t i = 0; i < n; ++i) z[i] = b[i] - z[i];
                last_ = norm2(z) / bnorm;
                if (last_ < cfg_.tolerance) return it + 1;
            }

            // Update rho, β, and u1, v for next outer iteration.
            // r_{k+1} = r_k - α A (u1 + u2)   (CGS identity)
            for (std::size_t i = 0; i < n; ++i)
                r[i] -= alpha * (Au1[i] + Au2[i]);
            const double rhoNew = dot(rh, r);
            if (std::abs(rho) < 1.0e-30) { last_ = tau / bnorm; return it; }
            const double beta = rhoNew / rho;
            // u1 = M^-1 r + β (u2 + β u1)    (CGS direction update,
            // preconditioned variant)
            M_->apply(r, z);
            for (std::size_t i = 0; i < n; ++i)
                u1[i] = z[i] + beta * (u2[i] + beta * u1[i]);
            A.spmv(u1, Au1);
            // v = A u1 + β (A u2 + β v)
            for (std::size_t i = 0; i < n; ++i)
                v[i] = Au1[i] + beta * (Au2[i] + beta * v[i]);
            rho = rhoNew;
        }
        return cfg_.maxIterations;
    }
    double last_residual() const override { return last_; }
private:
    LinearSolverConfig                cfg_;
    std::unique_ptr<IPreconditioner>  M_;
    double                            last_ = 0.0;
};
}  // namespace

std::unique_ptr<ILinearSolver> make_tfqmr(LinearSolverConfig cfg) {
    return std::make_unique<TfqmrSolver>(cfg);
}

}  // namespace simall::solver
