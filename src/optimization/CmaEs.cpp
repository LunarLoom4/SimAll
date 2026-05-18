// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/CmaEs.cpp
// =============================================================================
#include "optimization/CmaEs.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>

namespace simall::optimization {

namespace {

// Symmetric eigendecomposition (Jacobi rotations).  Small n only.
void jacobi_eig(std::vector<std::vector<double>>& A,
                 std::vector<std::vector<double>>& V,
                 std::vector<double>& w) {
    const int n = int(A.size());
    V.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) V[i][i] = 1.0;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q)
                off += std::abs(A[p][q]);
        if (off < 1e-14) break;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const double apq = A[p][q];
                if (std::abs(apq) < 1e-20) continue;
                const double theta = (A[q][q] - A[p][p]) / (2.0 * apq);
                const double t = (theta >= 0)
                    ?  1.0 / (theta + std::sqrt(1.0 + theta*theta))
                    : -1.0 / (-theta + std::sqrt(1.0 + theta*theta));
                const double c = 1.0 / std::sqrt(1.0 + t*t);
                const double s = t * c;
                const double app = A[p][p];
                const double aqq = A[q][q];
                A[p][p] = app - t * apq;
                A[q][q] = aqq + t * apq;
                A[p][q] = 0.0; A[q][p] = 0.0;
                for (int i = 0; i < n; ++i) {
                    if (i != p && i != q) {
                        const double aip = A[i][p];
                        const double aiq = A[i][q];
                        A[i][p] = c * aip - s * aiq; A[p][i] = A[i][p];
                        A[i][q] = s * aip + c * aiq; A[q][i] = A[i][q];
                    }
                    const double vip = V[i][p];
                    const double viq = V[i][q];
                    V[i][p] = c * vip - s * viq;
                    V[i][q] = s * vip + c * viq;
                }
            }
        }
    }
    w.assign(n, 0.0);
    for (int i = 0; i < n; ++i) w[i] = std::max(A[i][i], 0.0);
}

}  // namespace

CmaEsResult run_cmaes(CmaEsOptions opt) {
    CmaEsResult r;
    const std::size_t n = opt.dim;
    if (n == 0 || !opt.evaluate) { r.error = "Bad CMA-ES options"; return r; }
    if (opt.lambda == 0) opt.lambda = std::size_t(4 + std::floor(3.0 * std::log(double(n))));
    const std::size_t lambda = opt.lambda;
    const std::size_t mu     = lambda / 2;

    std::vector<double> weights(mu, 0.0);
    double sumW = 0.0;
    for (std::size_t i = 0; i < mu; ++i) {
        weights[i] = std::log(double(mu) + 0.5) - std::log(double(i + 1));
        sumW += weights[i];
    }
    for (auto& w : weights) w /= sumW;
    double muEff = 0.0;
    for (auto w : weights) muEff += w * w;
    muEff = 1.0 / muEff;

    const double cs    = (muEff + 2.0) / (double(n) + muEff + 5.0);
    const double ds    = 1.0 + 2.0 * std::max(0.0, std::sqrt((muEff - 1.0) / (double(n) + 1.0)) - 1.0) + cs;
    const double cc    = (4.0 + muEff / double(n)) / (double(n) + 4.0 + 2.0 * muEff / double(n));
    const double c1    = 2.0 / (std::pow(double(n) + 1.3, 2.0) + muEff);
    const double cmu   = std::min(1.0 - c1,
                                   2.0 * (muEff - 2.0 + 1.0 / muEff) /
                                   (std::pow(double(n) + 2.0, 2.0) + muEff));
    const double chiN  = std::sqrt(double(n)) * (1.0 - 1.0 / (4.0 * double(n))
                                  + 1.0 / (21.0 * double(n) * double(n)));

    std::vector<double> mean = opt.mean0.empty() ? std::vector<double>(n, 0.0) : opt.mean0;
    if (mean.size() != n) mean.assign(n, 0.0);
    const std::vector<double> xLow = opt.xLow.empty()
        ? std::vector<double>(n, -std::numeric_limits<double>::infinity())
        : opt.xLow;
    const std::vector<double> xUp  = opt.xUp.empty()
        ? std::vector<double>(n,  std::numeric_limits<double>::infinity())
        : opt.xUp;

    double sigma = opt.sigma0;
    std::vector<std::vector<double>> C(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) C[i][i] = 1.0;
    std::vector<double> pSigma(n, 0.0), pC(n, 0.0);

    std::mt19937_64 rng(opt.seed);
    std::normal_distribution<double> N01(0.0, 1.0);

    r.bestF = std::numeric_limits<double>::infinity();
    r.bestX = mean;

    for (std::size_t gen = 0; gen < opt.maxGen; ++gen) {
        // Eigendecompose C = B D² Bᵀ.
        std::vector<std::vector<double>> Ccopy = C;
        std::vector<std::vector<double>> B;
        std::vector<double>              d2;
        jacobi_eig(Ccopy, B, d2);
        std::vector<double> D(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) D[i] = std::sqrt(std::max(d2[i], 1e-30));

        // Sample λ offspring.
        std::vector<std::vector<double>> X(lambda, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> Z(lambda, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> Y(lambda, std::vector<double>(n, 0.0));
        std::vector<double>              fvals(lambda, 0.0);
        for (std::size_t k = 0; k < lambda; ++k) {
            for (std::size_t i = 0; i < n; ++i) Z[k][i] = N01(rng);
            for (std::size_t i = 0; i < n; ++i) {
                double s = 0.0;
                for (std::size_t j = 0; j < n; ++j) s += B[i][j] * D[j] * Z[k][j];
                Y[k][i] = s;
                X[k][i] = std::min(std::max(mean[i] + sigma * s, xLow[i]), xUp[i]);
            }
            fvals[k] = opt.evaluate(X[k]);
        }

        // Sort by fitness ascending.
        std::vector<std::size_t> idx(lambda);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b){ return fvals[a] < fvals[b]; });

        if (fvals[idx[0]] < r.bestF) {
            r.bestF = fvals[idx[0]];
            r.bestX = X[idx[0]];
        }
        r.history.push_back(r.bestF);

        // Recombination of best μ.
        std::vector<double> newMean(n, 0.0);
        std::vector<double> yMean(n, 0.0);
        for (std::size_t k = 0; k < mu; ++k) {
            for (std::size_t i = 0; i < n; ++i) {
                newMean[i] += weights[k] * X[idx[k]][i];
                yMean[i]   += weights[k] * Y[idx[k]][i];
            }
        }
        // pSigma  = (1-cs) pSigma + sqrt(cs(2-cs)muEff) * C^{-1/2} yMean
        // C^{-1/2} y = B D^{-1} Bᵀ y
        std::vector<double> BtY(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) BtY[i] += B[j][i] * yMean[j];
        std::vector<double> invCSqrtY(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                invCSqrtY[i] += B[i][j] * (BtY[j] / D[j]);

        for (std::size_t i = 0; i < n; ++i)
            pSigma[i] = (1.0 - cs) * pSigma[i]
                       + std::sqrt(cs * (2.0 - cs) * muEff) * invCSqrtY[i];

        const double pSigmaNorm = std::sqrt(std::inner_product(pSigma.begin(), pSigma.end(), pSigma.begin(), 0.0));
        const bool   hSig = pSigmaNorm / std::sqrt(1.0 - std::pow(1.0 - cs, 2.0 * (gen + 1)))
                            < (1.4 + 2.0 / (double(n) + 1.0)) * chiN;

        for (std::size_t i = 0; i < n; ++i)
            pC[i] = (1.0 - cc) * pC[i]
                   + (hSig ? std::sqrt(cc * (2.0 - cc) * muEff) * yMean[i] : 0.0);

        // Rank-1 + rank-mu update of C.
        const double delta = (1.0 - (hSig ? 1.0 : 0.0)) * cc * (2.0 - cc);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                double rankMu = 0.0;
                for (std::size_t k = 0; k < mu; ++k)
                    rankMu += weights[k] * Y[idx[k]][i] * Y[idx[k]][j];
                C[i][j] = (1.0 - c1 - cmu) * C[i][j]
                         + c1 * (pC[i] * pC[j] + delta * C[i][j])
                         + cmu * rankMu;
            }
        }
        // Sigma update.
        sigma *= std::exp((cs / ds) * (pSigmaNorm / chiN - 1.0));
        sigma = std::min(sigma, 1e10);

        mean = newMean;
        r.generations = gen + 1;
        if (sigma < opt.tolSigma) break;
        if (r.history.size() > 5
            && std::abs(r.history.back() - r.history[r.history.size() - 5]) < opt.tolFitness)
            break;
    }
    r.ok = true;
    return r;
}

}  // namespace simall::optimization
