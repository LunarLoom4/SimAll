// =============================================================================
// SimAll Beta - Adjoint Subsystem
// File   : src/adjoint/DiscreteAdjoint.cpp
// =============================================================================
#include "adjoint/DiscreteAdjoint.hpp"

#include <algorithm>
#include <cmath>

namespace simall::adjoint
{

void spmv_transpose(const CsrMatrix& A, const double* x, double* y) noexcept
{
    std::fill(y, y + A.n, 0.0);
    for (std::size_t i = 0; i < A.n; ++i) {
        const double xi = x[i];
        for (std::size_t k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k) {
            y[A.colIdx[k]] += A.values[k] * xi;
        }
    }
}

namespace
{

double dot(const double* a, const double* b, std::size_t n) noexcept
{
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        s += a[i] * b[i];
    return s;
}
double norm2(const double* a, std::size_t n) noexcept
{
    return std::sqrt(dot(a, a, n));
}

void jacobi_inv_diag(const CsrMatrix& A, std::vector<double>& invD)
{
    invD.assign(A.n, 1.0);
    for (std::size_t i = 0; i < A.n; ++i) {
        for (std::size_t k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k) {
            if (A.colIdx[k] == i) {
                const double d = A.values[k];
                invD[i] = (std::abs(d) > 1e-300) ? 1.0 / d : 1.0;
                break;
            }
        }
    }
}

// Restarted GMRES on Aᵀ with diagonal-Jacobi preconditioning.
bool gmres_transpose(const CsrMatrix& A,
                     const std::vector<double>& b,
                     std::vector<double>& x,
                     std::size_t maxIter,
                     std::size_t mRestart,
                     double tol,
                     std::size_t& itersOut,
                     double& resOut)
{
    const std::size_t n = A.n;
    x.assign(n, 0.0);
    std::vector<double> invD;
    jacobi_inv_diag(A, invD);

    std::vector<double> r(n), w(n);
    const double bnorm = std::max(norm2(b.data(), n), 1e-300);
    itersOut = 0;
    while (itersOut < maxIter) {
        // r = b - Aᵀ x   (preconditioned by left-Jacobi)
        spmv_transpose(A, x.data(), r.data());
        for (std::size_t i = 0; i < n; ++i)
            r[i] = invD[i] * (b[i] - r[i]);
        const double beta = norm2(r.data(), n);
        resOut = beta;
        if (beta / bnorm < tol)
            return true;

        const std::size_t m = mRestart;
        std::vector<std::vector<double>> V(m + 1, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> H(m + 1, std::vector<double>(m, 0.0));
        std::vector<double> g(m + 1, 0.0);
        std::vector<double> cs(m, 0.0), sn(m, 0.0);
        g[0] = beta;
        for (std::size_t i = 0; i < n; ++i)
            V[0][i] = r[i] / beta;

        std::size_t inner = 0;
        for (std::size_t j = 0; j < m; ++j) {
            // w = M⁻¹ Aᵀ V[j]
            spmv_transpose(A, V[j].data(), w.data());
            for (std::size_t i = 0; i < n; ++i)
                w[i] *= invD[i];
            // Modified Gram-Schmidt
            for (std::size_t i = 0; i <= j; ++i) {
                H[i][j] = dot(w.data(), V[i].data(), n);
                for (std::size_t k = 0; k < n; ++k)
                    w[k] -= H[i][j] * V[i][k];
            }
            H[j + 1][j] = norm2(w.data(), n);
            if (H[j + 1][j] > 1e-300) {
                for (std::size_t k = 0; k < n; ++k)
                    V[j + 1][k] = w[k] / H[j + 1][j];
            }
            // Apply previous Givens
            for (std::size_t i = 0; i < j; ++i) {
                const double t = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                H[i][j] = t;
            }
            // New Givens
            const double denom = std::hypot(H[j][j], H[j + 1][j]);
            cs[j] = H[j][j] / denom;
            sn[j] = H[j + 1][j] / denom;
            H[j][j] = denom;
            H[j + 1][j] = 0.0;
            const double t = cs[j] * g[j];
            g[j + 1] = -sn[j] * g[j];
            g[j] = t;
            inner = j + 1;
            ++itersOut;
            if (std::abs(g[j + 1]) / bnorm < tol)
                break;
            if (itersOut >= maxIter)
                break;
        }
        // Back-substitute y = H⁻¹ g
        std::vector<double> y(inner, 0.0);
        for (std::ptrdiff_t i = std::ptrdiff_t(inner) - 1; i >= 0; --i) {
            double s = g[i];
            for (std::size_t k = i + 1; k < inner; ++k)
                s -= H[i][k] * y[k];
            y[i] = (std::abs(H[i][i]) > 1e-300) ? s / H[i][i] : 0.0;
        }
        for (std::size_t k = 0; k < inner; ++k)
            for (std::size_t i = 0; i < n; ++i)
                x[i] += y[k] * V[k][i];

        resOut = std::abs(g[inner]);
        if (resOut / bnorm < tol)
            return true;
    }
    return resOut / std::max(bnorm, 1e-300) < tol;
}

} // namespace

DiscreteAdjointResult solve_discrete_adjoint(const DiscreteAdjointInputs& in,
                                             std::size_t nDV,
                                             DiscreteAdjointOptions opt)
{
    DiscreteAdjointResult r;
    if (!in.A) {
        r.error = "Null Jacobian";
        return r;
    }
    if (in.dJdQ.size() != in.A->n) {
        r.error = "dJdQ size mismatch";
        return r;
    }
    r.psi.assign(in.A->n, 0.0);
    if (!gmres_transpose(
            *in.A, in.dJdQ, r.psi, opt.maxIter, opt.restart, opt.tol, r.iterations, r.residual)) {
        r.error = "Adjoint GMRES did not converge to tol";
        // Continue: caller often wants the best-effort gradient anyway.
    }

    // gradient = dJdAlpha - (dR/dα)ᵀ ψ
    r.gradient.assign(nDV, 0.0);
    if (in.dRdAlphaTransposeTimes) {
        std::vector<double> dRtPsi(nDV, 0.0);
        in.dRdAlphaTransposeTimes(dRtPsi);
        for (std::size_t i = 0; i < nDV; ++i)
            r.gradient[i] = -dRtPsi[i];
    }
    if (!in.dJdAlpha.empty()) {
        const std::size_t k = std::min(nDV, in.dJdAlpha.size());
        for (std::size_t i = 0; i < k; ++i)
            r.gradient[i] += in.dJdAlpha[i];
    }
    r.ok = true;
    return r;
}

} // namespace simall::adjoint
