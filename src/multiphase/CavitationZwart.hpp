// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/CavitationZwart.hpp
// Phase  : 12.7 — Zwart-Gerber-Belamri (2004) cavitation model.
//
// Inception controlled by the nuclei volume fraction r_nuc rather than by
// the Schnerr-Sauer per-unit-volume nuclei count.
//
//   ṁ⁻ (vaporisation, p < p_v):
//       ṁ⁻ = F_vap · 3 ρ_v r_nuc (1 − α_v) / R_B · √( (2/3)(p_v − p)/ρ_l )
//
//   ṁ⁺ (condensation, p > p_v):
//       ṁ⁺ = F_cond · 3 ρ_v α_v / R_B · √( (2/3)(p − p_v)/ρ_l )
//
// Default constants per Zwart et al. (2004):
//   R_B    = 10⁻⁶ m   (representative bubble radius)
//   r_nuc  = 5×10⁻⁴   (volume fraction of nucleation sites)
//   F_vap  = 50,   F_cond = 0.01.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

namespace simall::multiphase
{

struct ZwartProps
{
    double rhoLiquid = 998.2;
    double rhoVapor = 0.02308;
    double pSaturation = 2339.0;
    double R_bubble = 1.0e-6;     // m
    double r_nucleation = 5.0e-4; // [-]
    double F_vap = 50.0;
    double F_cond = 0.01;
};

class CavitationZwart
{
public:
    void initialize(const meshing::Mesh& mesh, ZwartProps props);
    double apply(solver::FieldRegistry& fields);
    const ZwartProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    ZwartProps p_{};
};

} // namespace simall::multiphase
