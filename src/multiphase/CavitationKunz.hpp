// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/CavitationKunz.hpp
// Phase  : 12.7 — Kunz et al. (2000) cavitation model.
//
// Two distinct empirical rate coefficients drive vaporisation (m⁻) and
// condensation (m⁺) processes:
//
//   ṁ⁻ = -C_dest · ρ_v · α_l · min(0, p − p_v) / (½ ρ_l U_∞² t_∞)
//   ṁ⁺ = +C_prod · ρ_v · α_v² (1 − α_v) / t_∞
//
// where t_∞ = L_ref / U_∞ is a reference time-scale.  The model is widely
// used for sheet/cloud cavitation around hydrofoils and pumps.  Sign
// convention: ṁ > 0 → liquid → vapor.
//
// FieldRegistry output:
//   "S_alpha_v"  =  ṁ / ρ_v       [1/s]
//   "S_mass"     =  ṁ              [kg/(m³·s)]
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

namespace simall::multiphase {

struct KunzProps {
    double rhoLiquid    = 998.2;       // kg/m³
    double rhoVapor     = 0.02308;     // kg/m³
    double pSaturation  = 2339.0;      // Pa
    double Uref         = 1.0;         // m/s reference velocity
    double Lref         = 1.0;         // m   reference length
    double C_dest       = 100.0;       // empirical (Kunz 2000)
    double C_prod       = 100.0;       // empirical (Kunz 2000)
};

class CavitationKunz {
public:
    void initialize(const meshing::Mesh& mesh, KunzProps props);
    double apply(solver::FieldRegistry& fields);
    const KunzProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    KunzProps            p_{};
};

}  // namespace simall::multiphase
