// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/Rosseland.hpp
// Phase  : 11.4 — Rosseland (optically-thick) diffusion approximation.
//
// Valid for optical thickness τ = (κ + σ_s) L_ref ≫ 1.  Radiative transport
// degenerates to a non-linear diffusion of the gas temperature:
//
//   q_rad = -(16 σ n² T³) / (3 (κ + σ_s)) · ∇T = -k_R ∇T
//   ∇·q_rad → added to "S_rad" energy-equation source (W/m³)
//
// References:
//   Modest, "Radiative Heat Transfer", 3 ed., §16.3 (Rosseland limit).
//   Siegel & Howell, "Thermal Radiation Heat Transfer", 6 ed., §13.4.
//
// The implementation reads cell-centred "T" from FieldRegistry, computes
// the temperature-dependent radiative conductivity k_R per cell, then
// assembles ∇·(k_R ∇T) by face-flux divergence using the centroid-to-
// centroid orthogonal distance.  Result is **added** to "S_rad".
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

namespace simall::radiation {

struct RosselandProps {
    double absorption   = 0.5;    // κ   [1/m]
    double scattering   = 0.0;    // σ_s [1/m]
    double refractiveN  = 1.0;    // n   refractive index
    double T_min        = 50.0;   // safety clamp
    double T_max        = 5000.0;
    bool   accumulate   = true;   // add to S_rad; false = overwrite
};

class Rosseland {
public:
    bool initialize(const meshing::Mesh& mesh,
                    solver::FieldRegistry& fields,
                    const RosselandProps& props);

    /// Evaluates ∇·(k_R ∇T) and writes/adds to "S_rad".  Returns the
    /// L∞ norm of the source added this call.
    double apply();

    const RosselandProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh*   mesh_ = nullptr;
    solver::FieldRegistry* F_    = nullptr;
    RosselandProps         p_{};
};

}  // namespace simall::radiation
