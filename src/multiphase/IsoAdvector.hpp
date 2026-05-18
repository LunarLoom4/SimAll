// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/IsoAdvector.hpp
// Phase  : 12.10 — Geometric VOF advection (Roenby, Bredmose, Jasak 2016).
//
// "Iso" stands for iso-surface — the algorithm reconstructs the cell
// interface as a plane through the α-iso-surface and time-integrates the
// flux of liquid across each face by a sub-time-step swept-volume
// procedure.  Compared with MULES, IsoAdvector preserves sharp interfaces
// without algebraic compression and is bounded by construction.
//
//   1. Reconstruct interface plane per cell:
//        n_c    = -∇α / |∇α|     (cell normal)
//        x*     : isovalue α_iso ≈ 0.5 along the interpolated face values
//        cut-volume α_c maps uniquely to a plane offset δ_c.
//   2. For each face: estimate Δα_f over [t, t+Δt] by:
//        - Computing the swept region from the face's owner cell.
//        - Integrating its overlap with the interface half-space.
//        - Δα_f = ∫_{t}^{t+Δt} (φ_f · α_f^*) dt
//   3. Update α_c by mass-conservative summation of face contributions:
//        α_c^{n+1} = α_c^n - Δt/V_c · Σ_f Δα_f · sign(owner,c)
//
// This module wires the reconstruction (delegated to PlicReconstruction)
// with the geometric-flux integration.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "multiphase/PlicReconstruction.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <vector>

namespace simall::multiphase {

struct IsoAdvectorParams {
    double alpha_iso     = 0.5;       // iso-value for interface position
    int    subSteps      = 1;         // sub-time-steps within each Δt
    double tolBisection  = 1.0e-6;    // tolerance on plane-offset bisection
    int    maxBisection  = 32;        // maximum iterations
    bool   boundedness   = true;      // clip α ∈ [0,1] after update
};

class IsoAdvector {
public:
    void initialize(const meshing::Mesh& mesh, IsoAdvectorParams params = {});

    /// Advance α by Δt on the mesh using the face mass fluxes (kg/s).
    /// Reads "alpha" and "massFlux" (face-centred); writes back "alpha".
    /// Returns L1 bound violation observed before clipping (diagnostic).
    double step(double dt, solver::FieldRegistry& fields);

private:
    /// Compute swept-region "α-flux" through face f over the sub-step.
    double face_flux_alpha(meshing::FaceId f,
                           const solver::ScalarField& alpha,
                           const PlicPatch& patch,
                           double dt,
                           double massFlux) const;

    const meshing::Mesh* mesh_ = nullptr;
    IsoAdvectorParams    p_{};
    PlicReconstruction   plic_;
};

}  // namespace simall::multiphase
