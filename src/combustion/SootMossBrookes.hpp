// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/SootMossBrookes.hpp
// Phase  : 11.18 — Moss-Brookes soot model (Brookes & Moss 1999).
//
// Two transport equations for soot particle number density N and soot
// mass fraction Y_s:
//
//   ∂(ρ Ñ)/∂t + ∇·(ρ Ũ Ñ) = ∇·(D_N ∇Ñ) + α_N - γ_N
//   ∂(ρ Ỹ_s)/∂t + ∇·(ρ Ũ Ỹ_s) = ∇·(D_s ∇Ỹ_s) + α_M + γ_M - β_M
//
// Sources (Brookes-Moss 1999 for methane/kerosene):
//   α_N   = C_α · N_A · [Fuel] · exp(-T_α/T)
//             - C_β · ((24 R T) / (ρ_s · N_A))^{1/2} · d_p^{1/2} · N²
//   γ_N   = oxidation contribution to number  (small; included via β_M)
//   α_M   = M_p C_α  [Fuel] exp(-T_α/T)              [g/(cm³·s)]
//   γ_M   = C_γ A_s (P_fuel) exp(-T_γ/T)             [surface growth]
//   β_M   = C_oxid η_coll A_s P_OH (T)^{1/2}         [OH-oxidation]
//
// d_p = (6 ρ Y_s / (π ρ_s N))^{1/3}, A_s = π d_p² N.
//
// Reference:
//   Brookes & Moss, "Predictions of soot and thermal radiation properties in
//   confined turbulent jet diffusion flames", Combust. Flame 116, 486-503
//   (1999); Lindstedt 1994 for nucleation kinetics.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

namespace simall::combustion
{

struct MossBrookesProps
{
    double C_alpha = 54.0;              // nucleation pre-exp     [1/s]
    double T_alpha = 21000.0;           // nucleation activation T [K]
    double C_beta = 1.0e-15;            // coagulation
    double C_gamma = 11700.0;           // surface growth pre-exp
    double T_gamma = 12100.0;           // surface growth Ea/R
    double C_oxid = 0.015;              // OH-oxidation collision
    double rho_soot = 1800.0;           // soot density [kg/m³]
    double Mw_fuel = 16.04e-3;          // CH4 default
    double C_atoms_per_p = 12.0 * 60.0; // PAH-style nucleation cluster
};

class SootMossBrookes
{
public:
    void initialize(const meshing::Mesh& mesh, MossBrookesProps props = {});

    /// Reads "T", "rho_mix", "Y_fuel", "Y_OH", "N_soot", "Y_soot";
    /// writes "S_N_soot" and "S_Y_soot" sources (kg or 1/m³ per second
    /// per unit volume).  Returns volume-integrated mass production rate
    /// (positive = formation dominates, negative = oxidation dominates).
    double apply(solver::FieldRegistry& fields);

    const MossBrookesProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    MossBrookesProps p_{};
};

} // namespace simall::combustion
