// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/RpiWallBoiling.hpp
// Phase  : 12.11 — Rensselaer Polytechnic Institute (RPI) wall-boiling
// model (Kurul & Podowski 1990; Tolubinsky & Kostanchuk 1970).
//
// Heat-flux partitioning at heated walls in subcooled / saturated boiling:
//
//     q_wall  =  q_conv  +  q_quench  +  q_evap
//
// where:
//     A_b      = nucleate-boiling area fraction  (Del Valle-Kenning 1985)
//     A_b      = min( π d_w² N_w / 4, 1 )
//     d_w      = bubble departure diameter (Tolubinsky-Kostanchuk):
//                d_w = d_ref · exp( -ΔT_sub/ΔT_0 )
//     N_w      = nucleation site density (Lemmert-Chawla):
//                N_w = ( 185 ΔT_sup )^1.805
//     f_w      = bubble departure frequency (Cole 1960):
//                f_w = √( 4 g (ρ_l - ρ_v) / (3 d_w ρ_l) )
//
//     q_conv   = (1 - A_b) h_c (T_w - T_l)
//     q_quench = 2 A_b √(λ_l ρ_l c_p,l f_w / π) (T_w - T_l)
//     q_evap   = π/6 · d_w³ · N_w · f_w · ρ_v · h_fg
//
// Bubble departure feeds the per-cell mass-transfer source ṁ_evap used by
// the vapour α-equation and momentum equation.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "solver/Solver.hpp"

#include <vector>

namespace simall::multiphase {

struct RpiProps {
    double T_sat       = 373.15;       // K (water @ 1 atm)
    double rho_liquid  = 958.0;
    double rho_vapor   = 0.598;
    double cp_liquid   = 4216.0;       // J/(kg K)
    double lambda_l    = 0.679;        // W/(m K)
    double h_fg        = 2.257e6;      // J/kg latent heat
    double d_ref       = 6.0e-4;       // m   Tolubinsky reference dia
    double dT_ref      = 45.0;         // K
    double g           = 9.81;
    double h_conv      = 1000.0;       // W/(m²K) single-phase HTC (fallback)
};

class RpiWallBoiling {
public:
    void initialize(const meshing::Mesh& mesh,
                    const std::vector<solver::BoundarySpec>& bcs,
                    RpiProps props);

    /// Per-cell evaluation:  reads "T" (cell T), "T_wall_zone" lookup,
    /// "alpha_v"; writes "S_alpha_v" + "S_mass" + "S_energy" augmenting
    /// any existing fields (additive).  Returns total integrated q_wall.
    double apply(solver::FieldRegistry& fields);

    const RpiProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh*              mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    RpiProps                          p_{};
};

}  // namespace simall::multiphase
