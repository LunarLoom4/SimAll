// =============================================================================
// SimAll Beta - ROM Subsystem
// File   : src/rom/OperatorInference.cpp
// =============================================================================
#include "rom/OperatorInference.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace simall::rom {

namespace {

// Tiny symmetric eigendecomposition (Jacobi).  Returns eigenpairs sorted
// by eigenvalue descending.
void sym_eig(std::vector<std::vector<double>>& A,
              std::vector<std::vector<double>>& V,
              std::vector<double>& w) {
    const int n = int(A.size());
    V.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) V[i][i] = 1.0;
    for (int sweep = 0; sweep < 120; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < n; ++p) for (int q = p + 1; q < n; ++q) off += std::abs(A[p][q]);
        if (off < 1e-14) break;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                if (std::abs(A[p][q]) < 1e-20) continue;
                const double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
                const double t = (theta >= 0)
                    ?  1.0 / (theta + std::sqrt(1.0 + theta*theta))
                    : -1.0 / (-theta + std::sqrt(1.0 + theta*theta));
                const double c = 1.0 / std::sqrt(1.0 + t*t);
                const double s = t * c;
                const double app = A[p][p], aqq = A[q][q], apq = A[p][q];
                A[p][p] = app - t * apq;
                A[q][q] = aqq + t * apq;
                A[p][q] = 0.0; A[q][p] = 0.0;
                for (int i = 0; i < n; ++i) {
                    if (i != p && i != q) {
                        const double aip = A[i][p], aiq = A[i][q];
                        A[i][p] = c * aip - s * aiq; A[p][i] = A[i][p];
                        A[i][q] = s * aip + c * aiq; A[q][i] = A[i][q];
                    }
                    const double vip = V[i][p], viq = V[i][q];
                    V[i][p] = c * vip - s * viq;
                    V[i][q] = s * vip + c * viq;
                }
            }
        }
    }
    w.assign(n, 0.0);
    for (int i = 0; i < n; ++i) w[i] = A[i][i];
    std::vector<int> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return w[a] > w[b]; });
    auto Vc = V; auto wc = w;
    for (int i = 0; i < n; ++i) {
        w[i] = wc[ord[i]];
        for (int r = 0; r < n; ++r) V[r][i] = Vc[r][ord[i]];
    }
}

// Solve  Aᵀ A β = Aᵀ b  with Tikhonov reg  λI.  A is (m × p), b is (m).
std::vector<double> ridge_solve(const std::vector<std::vector<double>>& A,
                                 const std::vector<double>& b,
                                 double lambda) {
    const std::size_t m = A.size();
    const std::size_t p = m ? A[0].size() : 0;
    if (p == 0) return {};
    std::vector<std::vector<double>> AtA(p, std::vector<double>(p, 0.0));
    std::vector<double>              Atb(p, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < p; ++j) {
            Atb[j] += A[i][j] * b[i];
            for (std::size_t k = 0; k < p; ++k) AtA[j][k] += A[i][j] * A[i][k];
        }
    }
    for (std::size_t j = 0; j < p; ++j) AtA[j][j] += lambda;
    // Gaussian elimination.
    std::vector<std::vector<double>> M = AtA;
    std::vector<double> rhs = Atb;
    for (std::size_t i = 0; i < p; ++i) {
        std::size_t piv = i;
        for (std::size_t r = i + 1; r < p; ++r) if (std::abs(M[r][i]) > std::abs(M[piv][i])) piv = r;
        std::swap(M[i], M[piv]); std::swap(rhs[i], rhs[piv]);
        const double d = M[i][i];
        if (std::abs(d) < 1e-30) continue;
        for (std::size_t r = i + 1; r < p; ++r) {
            const double f = M[r][i] / d;
            for (std::size_t c = i; c < p; ++c) M[r][c] -= f * M[i][c];
            rhs[r] -= f * rhs[i];
        }
    }
    std::vector<double> x(p, 0.0);
    for (std::ptrdiff_t i = std::ptrdiff_t(p) - 1; i >= 0; --i) {
        double s = rhs[i];
        for (std::size_t c = i + 1; c < p; ++c) s -= M[i][c] * x[c];
        x[i] = (std::abs(M[i][i]) > 1e-30) ? s / M[i][i] : 0.0;
    }
    return x;
}

}  // namespace

InferredOperators fit_operator_inference(
        const std::vector<std::vector<double>>& X,
        const std::vector<std::vector<double>>& dXdt,
        const std::vector<std::vector<double>>& U,
        OperatorInferenceOptions opt) {
    InferredOperators op;
    if (X.empty() || X[0].empty()) return op;
    const std::size_t nx = X.size();
    const std::size_t nt = X[0].size();
    const std::size_t nu = U.empty() ? 0 : U.size();

    // POD via XᵀX eigendecomposition (small n_t typical for snapshots).
    std::vector<std::vector<double>> XtX(nt, std::vector<double>(nt, 0.0));
    for (std::size_t i = 0; i < nt; ++i) {
        for (std::size_t j = i; j < nt; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < nx; ++k) s += X[k][i] * X[k][j];
            XtX[i][j] = s; XtX[j][i] = s;
        }
    }
    std::vector<std::vector<double>> evec;
    std::vector<double>              eval;
    sym_eig(XtX, evec, eval);
    double totEnergy = 0.0;
    for (auto v : eval) totEnergy += std::max(v, 0.0);
    std::size_t r = opt.reducedDim;
    if (r == 0) {
        double cum = 0.0;
        for (std::size_t i = 0; i < eval.size(); ++i) {
            cum += std::max(eval[i], 0.0);
            if (cum / std::max(totEnergy, 1e-30) >= opt.energyThreshold) { r = i + 1; break; }
        }
        if (r == 0) r = std::min<std::size_t>(eval.size(), 4);
    }
    r = std::min(r, nt);
    op.r = r;
    // V_{ki} = (X · evec_i) / sqrt(eval_i)
    op.V.assign(nx, std::vector<double>(r, 0.0));
    for (std::size_t i = 0; i < r; ++i) {
        const double sv = std::sqrt(std::max(eval[i], 1e-30));
        for (std::size_t k = 0; k < nx; ++k) {
            double s = 0.0;
            for (std::size_t j = 0; j < nt; ++j) s += X[k][j] * evec[j][i];
            op.V[k][i] = s / sv;
        }
    }
    // Project: xHat_{i,t} = Σ_k V_{k,i} X_{k,t};  dxHat_{i,t} likewise.
    std::vector<std::vector<double>> xHat (r,  std::vector<double>(nt, 0.0));
    std::vector<std::vector<double>> dxHat(r,  std::vector<double>(nt, 0.0));
    for (std::size_t i = 0; i < r; ++i) {
        for (std::size_t t = 0; t < nt; ++t) {
            double s1 = 0.0, s2 = 0.0;
            for (std::size_t k = 0; k < nx; ++k) {
                s1 += op.V[k][i] * X[k][t];
                s2 += op.V[k][i] * dXdt[k][t];
            }
            xHat[i][t]  = s1;
            dxHat[i][t] = s2;
        }
    }
    // Build regression matrix D (nt × p) with columns [x̂, u, (x̂⊗x̂)?].
    const std::size_t qsize = opt.includeQuadratic ? r * r : 0;
    const std::size_t p = r + nu + qsize;
    std::vector<std::vector<double>> D(nt, std::vector<double>(p, 0.0));
    for (std::size_t t = 0; t < nt; ++t) {
        for (std::size_t i = 0; i < r; ++i) D[t][i] = xHat[i][t];
        for (std::size_t j = 0; j < nu; ++j) D[t][r + j] = U[j][t];
        if (opt.includeQuadratic) {
            for (std::size_t i = 0; i < r; ++i)
                for (std::size_t j = 0; j < r; ++j)
                    D[t][r + nu + i * r + j] = xHat[i][t] * xHat[j][t];
        }
    }
    // Fit each reduced row independently.
    op.A_hat.assign(r, std::vector<double>(r, 0.0));
    op.B_hat.assign(r, std::vector<double>(nu, 0.0));
    op.H_hat.assign(r, std::vector<double>(qsize, 0.0));
    for (std::size_t i = 0; i < r; ++i) {
        std::vector<double> b(nt);
        for (std::size_t t = 0; t < nt; ++t) b[t] = dxHat[i][t];
        auto beta = ridge_solve(D, b, opt.regularization);
        for (std::size_t j = 0; j < r;  ++j) op.A_hat[i][j] = beta[j];
        for (std::size_t j = 0; j < nu; ++j) op.B_hat[i][j] = beta[r + j];
        for (std::size_t j = 0; j < qsize; ++j) op.H_hat[i][j] = beta[r + nu + j];
    }
    return op;
}

std::vector<double> rom_rhs(const InferredOperators& op,
                              const std::vector<double>& xHat,
                              const std::vector<double>& u) {
    const std::size_t r = op.r;
    std::vector<double> out(r, 0.0);
    if (xHat.size() != r) return out;
    for (std::size_t i = 0; i < r; ++i) {
        for (std::size_t j = 0; j < r; ++j) out[i] += op.A_hat[i][j] * xHat[j];
        for (std::size_t j = 0; j < op.B_hat[i].size() && j < u.size(); ++j)
            out[i] += op.B_hat[i][j] * u[j];
        if (!op.H_hat[i].empty()) {
            for (std::size_t a = 0; a < r; ++a)
                for (std::size_t b = 0; b < r; ++b)
                    out[i] += op.H_hat[i][a * r + b] * xHat[a] * xHat[b];
        }
    }
    return out;
}

}  // namespace simall::rom
