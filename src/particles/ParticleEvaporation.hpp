// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ParticleEvaporation.hpp
// Phase  : 13.13 — Droplet evaporation (Spalding 1953 + Abramzon-Sirignano 1989).
//
// For each Lagrangian parcel:
//
//   dm_p/dt = -π d_p ρ_g D_v · Sh* · ln(1 + B_M)
//
//   B_M = (Y_s - Y_∞) / (1 - Y_s)             — Spalding mass transfer #
//   Y_s = (Mw_v · p_s) / (Mw_v · p_s + Mw_g · (p - p_s))   (Raoult/Clausius-
//         Clapeyron at droplet surface)
//   Sh* = 2 + (Sh_0 - 2) / F_M(B_M),  Sh_0 = 2 + 0.6 Re^{1/2} Sc^{1/3}
//   F_M = (1 + B_M)^{0.7} · ln(1 + B_M) / B_M     (Abramzon-Sirignano film)
//
// Vapor mass added to the gas via the per-cell source field "S_vapor"
// (kg/(m³·s)).  Parcel mass decreases until d_p hits the user-specified
// d_min (parcel deactivated at that point).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "particles/LagrangianTracker.hpp"
#include "solver/FieldRegistry.hpp"

namespace simall::particles {

struct EvaporationProps {
    double Mw_vapor   = 18.015e-3;   // kg/mol water default
    double Mw_gas     = 28.97e-3;    // air default
    double L_vap      = 2.45e6;      // latent heat of vaporisation [J/kg]
    double Tboil      = 373.15;
    double T_ref      = 298.15;
    double Dv_ref     = 2.6e-5;      // [m²/s] water-in-air
    double rho_l      = 998.2;
    double cp_l       = 4186.0;
    double mu_g       = 1.81e-5;
    double rho_g_ref  = 1.225;
    double Sc_t       = 0.7;
    double d_min      = 1.0e-7;      // parcel disappears below
    bool   abramzon   = true;        // use film-corrected Sherwood
};

class ParticleEvaporation {
public:
    void initialize(const meshing::Mesh& mesh, EvaporationProps props);

    /// Reads "T" (gas), "p" (gas), "Y_vapor" (gas mass fraction).
    /// Writes additive "S_vapor", "S_energy" (J/(m³·s)).
    /// Updates parcel mass, sets active=false when d_p < d_min.
    double apply(double dt, LagrangianTracker& tracker,
                 solver::FieldRegistry& fields,
                 double parcelTemperature = 320.0);

    const EvaporationProps& props() const noexcept { return p_; }

private:
    /// Clausius-Clapeyron saturation pressure (Pa) at temperature T (K).
    double p_sat(double T) const;

    const meshing::Mesh* mesh_ = nullptr;
    EvaporationProps     p_{};
};

}  // namespace simall::particles
