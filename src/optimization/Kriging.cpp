// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/Kriging.cpp
// =============================================================================
#include "optimization/Kriging.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace simall::optimization {

namespace {

double kernel(const std::vector<double>& a, const std::vector<double>& b,
              const std::vector<double>& theta) {
    double s = 0.0;
    const std::size_t d = a.size();
    for (std::size_t i = 0; i < d; ++i) {
        const double diff = a[i] - b[i];
        s += theta[i] * diff * diff;
    }
    return std::exp(-s);
}

bool cholesky(std::vector<std::vector<double>>& A) {
    const std::size_t n = A.size();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double s = A[i][j];
            for (std::size_t k = 0; k < j; ++k) s -= A[i][k] * A[j][k];
            if (i == j) {
                if (s <= 0.0) return false;
                A[i][i] = std::sqrt(s);
            } else {
                A[i][j] = s / A[j][j];
            }
        }
        for (std::size_t j = i + 1; j < n; ++j) A[i][j] = 0.0;
    }
    return true;
}

void chol_solve(const std::vector<std::vector<double>>& L,
                 const std::vector<double>& b,
                 std::vector<double>& x) {
    const std::size_t n = L.size();
    std::vector<double> z(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double s = b[i];
        for (std::size_t k = 0; k < i; ++k) s -= L[i][k] * z[k];
        z[i] = s / L[i][i];
    }
    x.assign(n, 0.0);
    for (std::ptrdiff_t i = std::ptrdiff_t(n) - 1; i >= 0; --i) {
        double s = z[i];
        for (std::size_t k = i + 1; k < n; ++k) s -= L[k][i] * x[k];
        x[i] = s / L[i][i];
    }
}

double log_likelihood(const std::vector<std::vector<double>>& X,
                       const std::vector<double>& y,
                       const std::vector<double>& theta,
                       double& sigma2Out,
                       std::vector<std::vector<double>>& Lout,
                       std::vector<double>& alphaOut,
                       double& muOut) {
    const std::size_t n = X.size();
    const double nug = 1e-8;
    std::vector<std::vector<double>> K(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const double k = kernel(X[i], X[j], theta);
            K[i][j] = (i == j ? k + nug : k);
            K[j][i] = K[i][j];
        }
    }
    Lout = K;
    if (!cholesky(Lout)) return -1e18;

    // Universal Kriging: μ = (1ᵀ K⁻¹ 1)⁻¹ 1ᵀ K⁻¹ y
    std::vector<double> ones(n, 1.0);
    std::vector<double> Kinv1, KinvY;
    chol_solve(Lout, ones, Kinv1);
    chol_solve(Lout, y,    KinvY);
    const double oneKinv1 = std::inner_product(ones.begin(), ones.end(), Kinv1.begin(), 0.0);
    const double oneKinvY = std::inner_product(ones.begin(), ones.end(), KinvY.begin(), 0.0);
    muOut = oneKinvY / std::max(oneKinv1, 1e-300);

    std::vector<double> yc(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) yc[i] = y[i] - muOut;
    chol_solve(Lout, yc, alphaOut);
    sigma2Out = std::inner_product(yc.begin(), yc.end(), alphaOut.begin(), 0.0) / double(n);

    // log L = -0.5 * (n log σ² + 2 Σ log L_ii + n)
    double logDet = 0.0;
    for (std::size_t i = 0; i < n; ++i) logDet += std::log(Lout[i][i]);
    return -0.5 * (double(n) * std::log(std::max(sigma2Out, 1e-300)) + 2.0 * logDet);
}

}  // namespace

bool KrigingModel::fit(const std::vector<std::vector<double>>& X,
                        const std::vector<double>& y) {
    if (X.empty() || X.size() != y.size()) return false;
    const std::size_t d = X[0].size();
    X_     = X;
    y_     = y;
    theta_.assign(d, 1.0);

    // Coordinate-descent search on log θ.
    std::vector<std::vector<double>> bestL;
    std::vector<double>              bestAlpha;
    double bestSigma2 = 1.0, bestMu = 0.0;
    double bestLL = log_likelihood(X_, y_, theta_, bestSigma2, bestL, bestAlpha, bestMu);

    for (int sweep = 0; sweep < 4; ++sweep) {
        for (std::size_t i = 0; i < d; ++i) {
            double bestT = theta_[i];
            for (double lg : {-3.0, -1.5, -0.5, 0.0, 0.5, 1.5, 3.0}) {
                std::vector<double> t = theta_;
                t[i] = std::pow(10.0, lg);
                double s2 = 1.0, mu = 0.0;
                std::vector<std::vector<double>> L;
                std::vector<double> a;
                const double ll = log_likelihood(X_, y_, t, s2, L, a, mu);
                if (ll > bestLL) {
                    bestLL = ll; bestT = t[i];
                    bestL = L; bestAlpha = a; bestSigma2 = s2; bestMu = mu;
                }
            }
            theta_[i] = bestT;
        }
    }
    L_      = bestL;
    alpha_  = bestAlpha;
    sigma2_ = bestSigma2;
    mu_     = bestMu;
    return !L_.empty();
}

double KrigingModel::predict(const std::vector<double>& x) const {
    if (X_.empty()) return mu_;
    double s = 0.0;
    for (std::size_t i = 0; i < X_.size(); ++i)
        s += kernel(X_[i], x, theta_) * alpha_[i];
    return mu_ + s;
}

void KrigingModel::predict(const std::vector<double>& x, double& m, double& v) const {
    m = predict(x);
    if (X_.empty()) { v = sigma2_; return; }
    const std::size_t n = X_.size();
    std::vector<double> kxX(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) kxX[i] = kernel(X_[i], x, theta_);
    std::vector<double> KinvKxX;
    chol_solve(L_, kxX, KinvKxX);
    double dot = 0.0;
    for (std::size_t i = 0; i < n; ++i) dot += kxX[i] * KinvKxX[i];
    v = std::max(sigma2_ * (1.0 - dot), 0.0);
}

}  // namespace simall::optimization
