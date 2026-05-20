// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/Kriging.hpp
// Week   : 18
//
// Universal Kriging surrogate (Gaussian process regression) with an
// anisotropic squared-exponential kernel
//
//     k(x, x') = σ² · exp( -Σ_i θ_i · (x_i - x'_i)² )
//
// Trained on a DoE-style dataset {(x_i, y_i)}.  Hyperparameters (θ, σ²)
// are estimated by maximum likelihood via a coordinate-descent search on
// log θ.  Predictions return both mean and Bayesian variance.
//
// Usage pattern (typical CFD design loop):
//   1.  DOE/LHS → high-fidelity CFD evaluations
//   2.  KrigingModel.fit(X, y)
//   3.  CmaEs minimises (or MMA samples) the Kriging surrogate's
//       lower-confidence-bound  μ(x) - κ σ(x).
//   4.  Add the best surrogate optimum to the training set, re-fit.
// =============================================================================
#pragma once

#include <cstddef>
#include <vector>

namespace simall::optimization
{

class KrigingModel
{
public:
    KrigingModel() = default;

    /// Train on `X` (n rows, d cols) and outputs `y` (size n).
    /// Returns true if Cholesky succeeded.
    bool fit(const std::vector<std::vector<double>>& X, const std::vector<double>& y);

    /// Predict mean at `x`.
    [[nodiscard]] double predict(const std::vector<double>& x) const;

    /// Predict mean and Bayesian variance at `x`.
    void predict(const std::vector<double>& x, double& mean, double& var) const;

    [[nodiscard]] const std::vector<double>& theta() const noexcept { return theta_; }
    [[nodiscard]] double sigma2() const noexcept { return sigma2_; }
    [[nodiscard]] std::size_t n_train() const noexcept { return X_.size(); }

private:
    std::vector<std::vector<double>> X_;
    std::vector<double> y_;
    std::vector<double> theta_;
    std::vector<std::vector<double>> L_; // lower Cholesky factor
    std::vector<double> alpha_;          // K⁻¹ (y - μ)
    double mu_ = 0.0;
    double sigma2_ = 1.0;
};

} // namespace simall::optimization
