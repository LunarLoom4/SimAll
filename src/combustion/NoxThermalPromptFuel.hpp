// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/NoxThermalPromptFuel.hpp
// Phase  : 11.17 — NO_x post-processor (thermal + prompt + fuel-N pathways).
//
// Three formation routes feed a single transported scalar NO:
//
//   1. THERMAL (extended Zeldovich, 1946; Lavoie-Heywood-Keck 1970):
//        N2  + O  ⇌ NO + N      k1f = 1.8e8 exp(-318/RT) m³/(mol·s)
//        N   + O2 ⇌ NO + O      k2f = 1.8e4  T  exp(-4680/T)
//        N   + OH ⇌ NO + H      k3f = 7.1e7 exp(-450/T)
//        (Quasi-steady N atom assumption → algebraic d[NO]/dt.)
//
//   2. PROMPT (Fenimore 1971, De Soete fit 1975):
//        d[NO]/dt = k_pr [N2] [O2]^a [Fuel]^b exp(-E_pr/RT)
//        with a=1, b=0, k_pr=6.4e6 m³/(mol·s), E_pr/R=36510 K.
//
//   3. FUEL (Bose-Wendt 1988, De Soete 1975 simplified):
//        d[NO]/dt = k_f [Fuel-N] [O2] exp(-E_f/RT)
//        with k_f=1.0e10 1/s, E_f/R=33700 K.
//
// Implementation: per-cell source ω̇_NO summed over the three contributions,
// written into "S_NO" for an external ScalarTransport on Y_NO.
//
// Reference:
//   Bowman, "Kinetics of Pollutant Formation and Destruction in Combustion",
//   Prog. Energy Combust. Sci. 1, 33–45 (1975).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

namespace simall::combustion
{

struct NoxProps
{
    bool includeThermal = true;
    bool includePrompt = true;
    bool includeFuel = false;
    double rho_default = 1.0; // [kg/m³]  (when no rho field present)
    double Mw_NO = 30.006e-3;
    double Mw_N2 = 28.014e-3;
    double Mw_O2 = 31.998e-3;
    double Mw_fuelN = 17.031e-3; // NH3 default
};

class NoxThermalPromptFuel
{
public:
    void initialize(const meshing::Mesh& mesh, NoxProps props = {});

    /// Reads "T", and any of Y_N2 / Y_O2 / Y_OH / Y_fuel / Y_fuelN / Y_NO
    /// (zero-default if missing); writes additive contribution to "S_NO"
    /// in [kg/(m³·s)].  Returns total integrated source.
    double apply(solver::FieldRegistry& fields);

private:
    const meshing::Mesh* mesh_ = nullptr;
    NoxProps p_{};
};

} // namespace simall::combustion
