// =============================================================================
// SimAll Beta - ROM Subsystem
// File   : src/rom/EcswHyperReduction.hpp
// Week   : 18
//
// Energy-Conserving Sampling and Weighting (ECSW), Farhat et al. 2014.
// Builds a sparse cubature  {e_k, w_k}  such that for every snapshot
// the projected residual  Σ_k w_k Vᵀ r_{e_k}  approximates the full
// Galerkin projection  Σ_e Vᵀ r_e  to within an L²-tolerance τ.
//
// Solves the non-negative least-squares (NNLS) problem
//
//        min ‖ G w − b ‖₂   s.t.  w ≥ 0,
//
// with  G ∈ R^{(r·n_s) × n_e}   and  b = column sums of G.  We use the
// Lawson-Hanson active-set NNLS, which is the textbook choice for ECSW
// and converges in a small number of iterations even for thousands of
// element columns.
//
// Inputs (caller-supplied):
//   * `elementResiduals[s][e][i]` — reduced element residual r̂_e^s  ∈ R^r
//   * `tau`                        — relative L² tolerance, e.g. 1e-3.
// Outputs:
//   * `indices`  — chosen element ids  (size  N_select)
//   * `weights`  — positive cubature weights  (size N_select)
// =============================================================================
#pragma once

#include <cstddef>
#include <vector>

namespace simall::rom
{

struct EcswResult
{
    std::vector<std::size_t> indices;
    std::vector<double> weights;
    double relativeError = 0.0;
};

[[nodiscard]] EcswResult ecsw_hyper_reduction(
    const std::vector<std::vector<std::vector<double>>>& elementResiduals,
    double tau = 1e-3,
    std::size_t maxIter = 0);

} // namespace simall::rom
