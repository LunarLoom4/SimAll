// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/LeeMassTransfer.hpp
// Phase  : 12.7 — Lee (1979) volumetric phase-change source for VOF /
// mixture flows with a single saturation temperature.
//
// Volumetric mass-transfer rate (kg / m³ / s) per cell:
//
//        evaporation  (T > T_sat):
//           ṁ = + c_evap · α_l ρ_l (T − T_sat) / T_sat
//
//        condensation (T < T_sat):
//           ṁ = − c_cond · α_v ρ_v (T_sat − T) / T_sat
//
// The latent-heat sink/source added to the energy equation:
//
//        S_h = − ṁ · L_v       [W/m³]
//
// Writes into FieldRegistry:
//
//        "S_alpha_v"  =  ṁ / ρ_v      [1/s]   (vapor-fraction equation)
//        "S_mass"     =  ṁ            [kg/m³/s] (continuity coupling)
//        "S_energy"   =  S_h          [W/m³]  (energy equation sink/source)
//
// Reads "T", "alpha_v" (or "alpha_l" = 1 - alpha_v).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

namespace simall::multiphase
{

struct LeeProps
{
    double T_sat = 373.15;     // K (water @ 1 atm)
    double L_vap = 2.257e6;    // J/kg latent heat
    double rho_liquid = 958.0; // kg/m³
    double rho_vapor = 0.598;  // kg/m³
    double c_evap = 0.1;       // Lee evaporation constant [1/s]
    double c_cond = 0.1;       // Lee condensation constant [1/s]
};

class LeeMassTransfer
{
public:
    void initialize(const meshing::Mesh& mesh, const LeeProps& props);
    /// Returns total integrated phase change (∫ṁ dV) [kg/s] (diagnostic).
    double apply(solver::FieldRegistry& fields);

    const LeeProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    LeeProps p_{};
};

} // namespace simall::multiphase
