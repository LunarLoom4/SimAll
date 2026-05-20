// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/FluxCorrectedTransport.hpp
// Phase  : 6.5 — Flux-Corrected Transport (Boris-Book / Zalesak)
//
// FCT blends a low-order monotone flux F_L (typically upwind) and a high-
// order non-monotone flux F_H (typically central-difference or QUICK) so
// that the corrected flux:
//
//      F_f = F_L_f + α_f · (F_H_f - F_L_f),    α_f ∈ [0,1]
//
// preserves the local monotonicity (no new extrema, no overshoots /
// undershoots) of the low-order solution.  α_f is chosen via Zalesak's
// limiter using cell-wise sum of incoming/outgoing antidiffusive fluxes.
//
// Algorithm (Zalesak 1979):
//   1)  φ^low_i = φ^n_i - Δt/V_i · Σ_f F_L_f                        (advance)
//   2)  A_f      = F_H_f - F_L_f                                     (antidiffusion)
//   3)  P^+_i    = Σ_{f: A_f flows in}  |A_f|,   P^-_i = Σ outgoing  |A_f|
//   4)  Q^+_i    = (φ^max_i - φ^low_i) · V_i / Δt
//       Q^-_i    = (φ^low_i - φ^min_i) · V_i / Δt
//   5)  R^+_i    = min(1, Q^+_i / P^+_i),   R^-_i  = min(1, Q^-_i / P^-_i)
//   6)  α_f      = min(R^+_recv, R^-_send)
//
// References:
//   Boris & Book, J. Comp. Phys. 11, 38-69 (1973).
//   Zalesak, J. Comp. Phys. 31, 335-362 (1979).
//   Kuzmin, "A guide to numerical methods for transport equations" (2010).
// =============================================================================
#pragma once

#include "FieldRegistry.hpp"

#include "meshing/MeshStorage.hpp"
#include "solver/Solver.hpp"

#include <vector>

namespace simall::solver
{

class FluxCorrectedTransport
{
public:
    /// Construct against a mesh + BC table. Field & mass-flux passed per call.
    FluxCorrectedTransport(const meshing::Mesh& mesh, const std::vector<BoundarySpec>& bcs)
        : mesh_(mesh), bcs_(bcs)
    {
    }

    /// Single-step explicit FCT advance of a passive scalar:
    ///     ∂φ/∂t + ∇·(U φ) = 0
    /// Inputs:
    ///   phi        cell-centred field at time n;  modified in-place to n+1
    ///   massFlux   ρ U · A at each face (positive = owner→neighbour)
    ///   dt         time step
    /// The low-order flux is first-order upwind (TVD by construction).
    /// The high-order flux is centred (second-order, non-monotone).
    void advance(util::aligned_vector<double>& phi,
                 const util::aligned_vector<double>& massFlux,
                 double dt);

    /// Compute Zalesak limiter coefficients α_f without advancing the field
    /// (useful for hybrid implicit schemes that need only the limiter values).
    /// alpha is resized to faces().size().
    void compute_limiter(const util::aligned_vector<double>& phi,
                         const util::aligned_vector<double>& massFlux,
                         double dt,
                         util::aligned_vector<double>& alpha);

private:
    const meshing::Mesh& mesh_;
    const std::vector<BoundarySpec>& bcs_;
};

} // namespace simall::solver
