// =============================================================================
// SimAll Beta - Heat-Transfer Subsystem
// File   : src/heat_transfer/Solidification.hpp
// Phase  : 17 — Enthalpy-porosity solidification / melting model
// (Voller & Prakash 1987, Brent et al. 1988).
//
// A single mushy-zone formulation:
//
//   ρ Du/Dt = -∇p + ∇·τ + ρ g - A_mush u (1 - f_l)² / (f_l³ + ε)
//   ρ DH/Dt = ∇·(k ∇T) + S_h        H = h + f_l L
//
// where f_l ∈ [0,1] is liquid fraction (linear between T_solidus and
// T_liquidus), A_mush ≈ 10^5–10^7 the Darcy permeability coefficient,
// L the latent heat of fusion, and ε ≈ 1e-3 a small number.
//
// Outputs written into FieldRegistry:
//   - "f_l"         (scalar) liquid fraction
//   - "S_DarcyMom"  (vector) Darcy momentum sink
//   - "S_LatentEn"  (scalar) latent heat source for energy equation
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

namespace simall::heat_transfer {

struct SolidificationProps {
    double T_solidus    = 1700.0;     // K
    double T_liquidus   = 1730.0;
    double latentHeat   = 2.7e5;      // J/kg
    double rho          = 7800.0;     // kg/m³
    double cp           = 750.0;      // J/(kg·K)
    double A_mush       = 1.0e6;      // Darcy coefficient (Pa·s/m²)
    double eps          = 1.0e-3;
};

class Solidification {
public:
    void initialize(const meshing::Mesh& mesh, solver::FieldRegistry& fields,
                    SolidificationProps props);

    /// Update liquid fraction, Darcy momentum sink, and latent-heat source
    /// from current T and U fields. Should be called every outer iteration
    /// before the SIMPLE assembly.
    void update(double dt, solver::FieldRegistry& fields);

private:
    const meshing::Mesh*  mesh_   = nullptr;
    SolidificationProps   p_{};
    // Previous-step liquid fraction for ∂f_l/∂t latent source.
    util::aligned_vector<double> fl_prev_;
};

}  // namespace simall::heat_transfer
