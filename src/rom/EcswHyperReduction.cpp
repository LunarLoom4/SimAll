// =============================================================================
// SimAll Beta - ROM Subsystem
// File   : src/rom/EcswHyperReduction.cpp
// =============================================================================
#include "rom/EcswHyperReduction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace simall::rom {

namespace {

// Lawson-Hanson NNLS: min ‖G w − b‖ s.t. w ≥ 0.
//   G is (m × n), b is (m), result w is (n).
std::vector<double> nnls(const std::vector<std::vector<double>>& G,
                          const std::vector<double>& b,
                          std::size_t maxIter) {
    const std::size_t m = G.size();
    const std::size_t n = m ? G[0].size() : 0;
    std::vector<double> w(n, 0.0);                  // gradient
    std::vector<double> x(n, 0.0);                  // solution
    std::vector<bool>   P(n, false);                // passive set membership
    if (n == 0) return x;

    auto dot_col_resid = [&](std::size_t j, const std::vector<double>& r) {
        double s = 0.0;
        for (std::size_t i = 0; i < m; ++i) s += G[i][j] * r[i];
        return s;
    };

    std::vector<double> r = b;                      // residual
    std::size_t iter = 0;
    const std::size_t hardCap = maxIter == 0 ? 3 * n : maxIter;

    while (iter++ < hardCap) {
        // w = Gᵀ r
        for (std::size_t j = 0; j < n; ++j)
            w[j] = P[j] ? -std::numeric_limits<double>::infinity()
                         : dot_col_resid(j, r);
        // find argmax over inactive set
        std::size_t jStar = n;
        double wMax = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (!P[j] && w[j] > wMax) { wMax = w[j]; jStar = j; }
        }
        if (jStar == n || wMax < 1e-12) break;
        P[jStar] = true;

        // Solve unconstrained LS over passive set.
        while (true) {
            std::vector<std::size_t> pSet;
            for (std::size_t j = 0; j < n; ++j) if (P[j]) pSet.push_back(j);
            const std::size_t p = pSet.size();
            // Normal equations  (Gpᵀ Gp) s = Gpᵀ b
            std::vector<std::vector<double>> AtA(p, std::vector<double>(p, 0.0));
            std::vector<double>              Atb(p, 0.0);
            for (std::size_t i = 0; i < m; ++i) {
                for (std::size_t a = 0; a < p; ++a) {
                    const double gia = G[i][pSet[a]];
                    Atb[a] += gia * b[i];
                    for (std::size_t bIdx = 0; bIdx < p; ++bIdx)
                        AtA[a][bIdx] += gia * G[i][pSet[bIdx]];
                }
            }
            // Tiny regularisation.
            for (std::size_t a = 0; a < p; ++a) AtA[a][a] += 1e-12;
            // Gauss elimination
            std::vector<std::vector<double>> M = AtA;
            std::vector<double> rhs = Atb;
            for (std::size_t i = 0; i < p; ++i) {
                std::size_t piv = i;
                for (std::size_t rr = i + 1; rr < p; ++rr)
                    if (std::abs(M[rr][i]) > std::abs(M[piv][i])) piv = rr;
                std::swap(M[i], M[piv]); std::swap(rhs[i], rhs[piv]);
                const double d = M[i][i];
                if (std::abs(d) < 1e-30) continue;
                for (std::size_t rr = i + 1; rr < p; ++rr) {
                    const double f = M[rr][i] / d;
                    for (std::size_t c = i; c < p; ++c) M[rr][c] -= f * M[i][c];
                    rhs[rr] -= f * rhs[i];
                }
            }
            std::vector<double> s(p, 0.0);
            for (std::ptrdiff_t i = std::ptrdiff_t(p) - 1; i >= 0; --i) {
                double sv = rhs[i];
                for (std::size_t c = i + 1; c < p; ++c) sv -= M[i][c] * s[c];
                s[i] = (std::abs(M[i][i]) > 1e-30) ? sv / M[i][i] : 0.0;
            }

            // Check non-negativity.
            double alpha = 1.0;
            std::size_t jBlock = p;
            for (std::size_t a = 0; a < p; ++a) {
                if (s[a] <= 0.0) {
                    const double xj = x[pSet[a]];
                    const double t  = xj / (xj - s[a] + 1e-30);
                    if (t < alpha) { alpha = t; jBlock = a; }
                }
            }
            if (jBlock == p) {
                // Accept this passive-set solution.
                for (std::size_t a = 0; a < p; ++a) x[pSet[a]] = s[a];
                break;
            }
            // Move along the constrained direction by alpha, drop the worst.
            for (std::size_t a = 0; a < p; ++a)
                x[pSet[a]] = x[pSet[a]] + alpha * (s[a] - x[pSet[a]]);
            for (std::size_t a = 0; a < p; ++a)
                if (x[pSet[a]] <= 1e-14) {
                    x[pSet[a]] = 0.0;
                    P[pSet[a]] = false;
                }
        }
        // Recompute residual.
        for (std::size_t i = 0; i < m; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < n; ++j) s += G[i][j] * x[j];
            r[i] = b[i] - s;
        }
        double rn = 0.0;
        for (auto v : r) rn += v * v;
        if (rn < 1e-24) break;
    }
    return x;
}

}  // namespace

EcswResult ecsw_hyper_reduction(
        const std::vector<std::vector<std::vector<double>>>& er,
        double tau, std::size_t maxIter) {
    EcswResult out;
    if (er.empty() || er[0].empty()) return out;
    const std::size_t nSnap = er.size();
    const std::size_t nElem = er[0].size();
    const std::size_t r     = er[0][0].size();

    // Build G (m × n_e) and b (m) where m = r · n_s.
    const std::size_t m = r * nSnap;
    std::vector<std::vector<double>> G(m, std::vector<double>(nElem, 0.0));
    std::vector<double>              b(m, 0.0);
    for (std::size_t s = 0; s < nSnap; ++s) {
        for (std::size_t e = 0; e < nElem; ++e) {
            for (std::size_t i = 0; i < r; ++i) {
                G[s * r + i][e] = er[s][e][i];
                b[s * r + i]   += er[s][e][i];
            }
        }
    }

    auto w = nnls(G, b, maxIter);
    double rn = 0.0, bn = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < nElem; ++j) s += G[i][j] * w[j];
        rn += (b[i] - s) * (b[i] - s);
        bn += b[i] * b[i];
    }
    out.relativeError = std::sqrt(rn / std::max(bn, 1e-30));
    (void)tau;  // tau is informational here — NNLS minimises ‖Gw - b‖ outright

    for (std::size_t j = 0; j < nElem; ++j) {
        if (w[j] > 1e-10) {
            out.indices.push_back(j);
            out.weights.push_back(w[j]);
        }
    }
    return out;
}

}  // namespace simall::rom
