// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/CavitationSchnerrSauer.hpp
// Phase  : 12.6 — Schnerr & Sauer (2001) homogeneous-mixture cavitation
// model for VOF / mixture flows.
//
// Mass-transfer between liquid (α_l = 1 - α_v) and vapor (α_v) is driven
// by the local pressure relative to saturation pressure p_v:
//
//   ṁ⁺ (condensation,  p > p_v):
//       ṁ⁺ = (3 ρ_v ρ_l / ρ_m) (α_v (1-α_v) / R_B) √(  (2/3)(p-p_v)/ρ_l )
//
//   ṁ⁻ (vaporization,   p < p_v):
//       ṁ⁻ = (3 ρ_v ρ_l / ρ_m) (α_v (1-α_v) / R_B) √( -(2/3)(p-p_v)/ρ_l )
//
//   R_B = ( α_v / ((1-α_v) (4/3) π n₀) )^(1/3)         bubble radius
//   ρ_m = α_v ρ_v + (1 - α_v) ρ_l                       mixture density
//
// The source ṁ (with sign) is written as
//   - "S_alpha_v" : volumetric source for the vapor volume-fraction equation
//   - "S_mass"    : continuity mass source coupling for the pressure equation
//
// n₀ is the per-unit-volume bubble nucleation density (default 1e13 1/m³).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"

namespace simall::multiphase
{

struct CavitationProps
{
    double rhoLiquid = 998.2;      // kg/m³
    double rhoVapor = 0.02308;     // kg/m³
    double pSaturation = 2339.0;   // Pa (water @ 20°C)
    double nucleiDensity = 1.0e13; // 1/m³
    double evapCoeff = 1.0;        // empirical multiplier on ṁ⁻
    double condCoeff = 1.0;        // empirical multiplier on ṁ⁺
};

class CavitationSchnerrSauer
{
public:
    void initialize(const meshing::Mesh& mesh, CavitationProps props);

    /// Reads fields["p"] and fields["alpha_v"], writes
    ///   fields["S_alpha_v"]   [1/s]   ← ṁ / ρ_v
    ///   fields["S_mass"]      [kg/m³/s] ← ṁ
    /// Returns the total mass-transfer integral over the mesh (diagnostic).
    double apply(solver::FieldRegistry& fields);

    const CavitationProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    CavitationProps p_{};
};

} // namespace simall::multiphase
