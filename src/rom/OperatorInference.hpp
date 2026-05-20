// =============================================================================
// SimAll Beta - ROM Subsystem
// File   : src/rom/OperatorInference.hpp
// Week   : 18
//
// Operator Inference (Peherstorfer & Willcox, 2016).  Builds a reduced-
// order model
//
//     dx̂/dt = Â x̂ + B̂ u(t) + Ĥ (x̂ ⊗ x̂)
//
// from *snapshot* trajectories of the full-order state, without ever
// touching the full-order operator code.  Given:
//
//     X    n_x  × n_t   state snapshots
//     dXdt n_x  × n_t   time-derivative snapshots (or finite-difference)
//     U    n_u  × n_t   input snapshots          (optional)
//
//   * Compute a POD basis V (n_x × r)  by truncated SVD of X.
//   * Project:  x̂ = Vᵀ x.
//   * Fit  Â, B̂, Ĥ by solving a least-squares problem
//
//        min  Σ_t ‖ Â x̂_t + B̂ u_t + Ĥ (x̂_t ⊗ x̂_t)  −  Vᵀ dXdt_t ‖²
//
//     via normal equations (Tikhonov-regularised SPD system).
//
// The quadratic Hessian Ĥ is stored as a flattened r × r² matrix.
// =============================================================================
#pragma once

#include <cstddef>
#include <vector>

namespace simall::rom
{

struct OperatorInferenceOptions
{
    std::size_t reducedDim = 0; // 0 → choose by 99% energy
    double energyThreshold = 0.99;
    double regularization = 1e-8;
    bool includeQuadratic = false;
};

struct InferredOperators
{
    std::size_t r = 0;
    std::vector<std::vector<double>> V;     // basis (n_x × r)
    std::vector<std::vector<double>> A_hat; // r × r
    std::vector<std::vector<double>> B_hat; // r × n_u
    std::vector<std::vector<double>> H_hat; // r × r² (row-major Kron)
};

[[nodiscard]] InferredOperators fit_operator_inference(
    const std::vector<std::vector<double>>& X,    // n_x × n_t
    const std::vector<std::vector<double>>& dXdt, // n_x × n_t
    const std::vector<std::vector<double>>& U,    // n_u × n_t (may be empty)
    OperatorInferenceOptions opt);

/// Evaluate dx̂/dt = Â x̂ + B̂ u + Ĥ (x̂ ⊗ x̂)
[[nodiscard]] std::vector<double> rom_rhs(const InferredOperators& op,
                                          const std::vector<double>& xHat,
                                          const std::vector<double>& u);

} // namespace simall::rom
