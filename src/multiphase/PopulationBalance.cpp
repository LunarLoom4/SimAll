// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/PopulationBalance.cpp
// =============================================================================
#include "multiphase/PopulationBalance.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace simall::multiphase {

void PopulationBalance::initialize(const meshing::Mesh& m, solver::FieldRegistry& f,
                                   PbmProps p) {
    mesh_ = &m; fields_ = &f; props_ = std::move(p);
    const std::size_t nC = m.cells().size();
    for (int k = 0; k < props_.numMoments; ++k)
        f.scalar(moment_name(k), nC);
}

std::string PopulationBalance::moment_name(int k) const {
    return "pbm_m" + std::to_string(k);
}

bool PopulationBalance::pd_algorithm(const std::vector<double>& m,
                                     std::vector<double>& w,
                                     std::vector<double>& L)
{
    const int twoN = static_cast<int>(m.size());
    if (twoN < 2 || (twoN & 1)) return false;
    const int N = twoN / 2;
    if (m[0] <= 0.0) return false;

    // Build product-difference table (Gordon 1968).
    std::vector<std::vector<double>> P(twoN + 1,
                                       std::vector<double>(twoN + 1, 0.0));
    P[0][0] = 1.0;
    for (int i = 0; i <= twoN - 1; ++i) P[i][1] = (i == 0 ? 1.0 : 0.0);
    for (int i = 0; i < twoN; ++i) P[i][2] = std::pow(-1.0, i) * m[i];
    for (int j = 3; j <= twoN; ++j) {
        for (int i = 0; i <= twoN - j + 1; ++i) {
            P[i][j] = P[0][j-1] * P[i+1][j-2] - P[0][j-2] * P[i+1][j-1];
        }
    }
    std::vector<double> alpha(twoN, 0.0);
    for (int j = 1; j < twoN; ++j) {
        const double denom = P[0][j] * P[0][j+1];
        alpha[j] = (std::abs(denom) > 1e-30) ? P[0][j+2] / denom : 0.0;
    }
    // Tridiagonal symmetric Jacobi matrix from continued-fraction coefs.
    std::vector<double> a(N, 0.0), b(N, 0.0);
    for (int i = 0; i < N; ++i) {
        a[i] = alpha[2*i+1] + (i > 0 ? alpha[2*i] : 0.0);
        if (i < N-1) {
            const double prod = alpha[2*i+1] * alpha[2*i+2];
            b[i] = (prod > 0) ? std::sqrt(prod) : 0.0;
        }
    }
    // QL algorithm on the symmetric tridiagonal Jacobi matrix.
    std::vector<double> d = a, e(N, 0.0);
    for (int i = 0; i < N-1; ++i) e[i] = b[i];
    e[N-1] = 0.0;
    std::vector<std::vector<double>> Z(N, std::vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) Z[i][i] = 1.0;

    for (int l = 0; l < N; ++l) {
        for (int iter = 0; iter < 60; ++iter) {
            int mIdx = l;
            for (; mIdx < N-1; ++mIdx) {
                const double dd = std::abs(d[mIdx]) + std::abs(d[mIdx+1]);
                if (std::abs(e[mIdx]) + dd == dd) break;
            }
            if (mIdx == l) break;
            double g = (d[l+1] - d[l]) / (2.0 * e[l]);
            double r = std::hypot(g, 1.0);
            g = d[mIdx] - d[l] + e[l] / (g + std::copysign(r, g));
            double s = 1.0, c = 1.0, p = 0.0;
            for (int i = mIdx - 1; i >= l; --i) {
                const double f = s * e[i];
                const double b1 = c * e[i];
                r = std::hypot(f, g);
                e[i+1] = r;
                if (r == 0.0) { d[i+1] -= p; e[mIdx] = 0.0; break; }
                s = f / r; c = g / r;
                g = d[i+1] - p;
                const double t = (d[i] - g) * s + 2.0 * c * b1;
                p = s * t; d[i+1] = g + p; g = c * t - b1;
                for (int k = 0; k < N; ++k) {
                    const double zz = Z[k][i+1];
                    Z[k][i+1] = s * Z[k][i] + c * zz;
                    Z[k][i]   = c * Z[k][i] - s * zz;
                }
            }
            if (r == 0.0 && mIdx >= l) continue;
            d[l] -= p; e[l] = g; e[mIdx] = 0.0;
        }
    }
    // Eigenvalues = abscissae; weights = m0 * (first eigenvector component)².
    L.assign(N, 0.0); w.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        L[i] = d[i];
        w[i] = m[0] * Z[0][i] * Z[0][i];
    }
    return true;
}

bool PopulationBalance::quadrature_nodes(std::size_t c,
                                         std::vector<double>& w,
                                         std::vector<double>& L) const
{
    if (!fields_) return false;
    std::vector<double> m(props_.numMoments);
    for (int k = 0; k < props_.numMoments; ++k) {
        const auto* f = fields_->find_scalar(moment_name(k));
        if (!f || c >= f->size()) return false;
        m[k] = (*f)[c];
    }
    return pd_algorithm(m, w, L);
}

void PopulationBalance::integrate_sources(double dt) {
    if (!fields_ || !mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    const int K = props_.numMoments;
    const int N = K / 2;
    std::vector<double> w, L;
    std::vector<solver::ScalarField*> mom(K);
    for (int k = 0; k < K; ++k) mom[k] = fields_->find_scalar(moment_name(k));

    for (std::size_t c = 0; c < nC; ++c) {
        if (!quadrature_nodes(c, w, L)) continue;
        for (int k = 0; k < K; ++k) {
            double S = 0.0;
            // Growth: G_k = k * G * Σ w_i L_i^{k-1}
            if (k > 0 && props_.growthRateG != 0.0) {
                for (int i = 0; i < N; ++i)
                    S += w[i] * std::pow(L[i], k-1);
                S *= k * props_.growthRateG;
            }
            // Aggregation: B_k - D_k = 0.5 ΣΣ w_i w_j β(L_i,L_j)((L_i^3+L_j^3)^{k/3} - L_i^k - L_j^k)
            if (props_.aggregationKernel) {
                double agg = 0.0;
                for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    const double b = props_.aggregationKernel(L[i], L[j]);
                    const double Lnew = std::cbrt(std::pow(L[i],3.0)+std::pow(L[j],3.0));
                    agg += w[i]*w[j]*b * (0.5*std::pow(Lnew,k) - 0.5*std::pow(L[i],k));
                }
                S += agg;
            }
            // Breakage: B_k - D_k = Σ w_i g(L_i) [∫ L^k β(L,L_i) dL - L_i^k]
            // Approximate the fragment integral by symmetric binary split.
            if (props_.breakageFrequency) {
                double brk = 0.0;
                for (int i = 0; i < N; ++i) {
                    const double g = props_.breakageFrequency(L[i]);
                    const double Ld = std::cbrt(0.5) * L[i];   // equal-volume split
                    brk += w[i] * g * (2.0*std::pow(Ld,k) - std::pow(L[i],k));
                }
                S += brk;
            }
            if (mom[k]) (*mom[k])[c] += dt * S;
        }
    }
}

}  // namespace simall::multiphase
