// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/EddyDissipationConcept.hpp
// Phase  : 10.4 — Turbulence-chemistry interaction via two related models:
//
//   1. EDM (Eddy-Dissipation Model, Magnussen & Hjertager 1976)
//      ω̇_k = ν_k W_k · A · (ρ ε / k) · min_R (Y_R / (ν_R W_R))
//
//   2. EDC (Eddy Dissipation Concept, Magnussen 1981)
//      γ_λ = C_γ · (ν ε / k²)^{1/4}     (fine-structure length fraction)
//      τ*  = C_τ · sqrt(ν / ε)          (fine-structure time scale)
//      ω̇_k = ρ · γ_λ² · (Y*_k - Y_k) / (τ* · (1 - γ_λ³))
//
//      Y*_k is obtained by a constant-pressure plug-flow reactor (PFR)
//      integration of the Arrhenius laminar source from Y_k over τ*,
//      using a numerically robust BDF1 (implicit Euler) step on the
//      species mass-fraction vector.
//
// Writes per-cell, per-species reaction-rate field "wdot_<species>" and
// the combined "S_combustion" energy source. Reads the upstream
// laminar reaction set provided via `set_reactions()`.
// =============================================================================
#pragma once

#include "combustion/LaminarFiniteRate.hpp"
#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <string>
#include <vector>

namespace simall::combustion
{

enum class TciMode
{
    EDM,
    EDC
};

struct EdcProps
{
    TciMode mode = TciMode::EDC;
    double A_EDM = 4.0;      // EDM A constant (mixing rate scaling)
    double B_EDM = 0.5;      // EDM B constant (oxidiser side)
    double C_gamma = 2.1377; // EDC γ_λ coefficient
    double C_tau = 0.4082;   // EDC τ*    coefficient
    double nu = 1.5e-5;      // kinematic viscosity (m²/s)
};

class EddyDissipationConcept
{
public:
    void initialize(const meshing::Mesh& mesh,
                    solver::FieldRegistry& fields,
                    std::vector<Species> species,
                    std::vector<Reaction> reactions,
                    EdcProps props);

    /// Apply TCI: compute ω̇_k for every cell, write into "wdot_<name>"
    /// fields, write combined heat-release into "S_combustion".
    void apply(solver::FieldRegistry& fields, double rho_default = 1.2);

private:
    /// Single-cell PFR step for EDC fine structure:
    ///   dY/dt = ω̇_lam(Y) / ρ
    /// Integrated by BDF1 from Y_in over time τ.
    void psr_step(const std::vector<double>& Yin,
                  double rho,
                  double T,
                  double tau,
                  std::vector<double>& Yout) const;

    /// Evaluate laminar Arrhenius rates ω̇_k for a single cell state.
    void laminar_rates(const std::vector<double>& Y,
                       double rho,
                       double T,
                       std::vector<double>& wdot) const;

    const meshing::Mesh* mesh_ = nullptr;
    EdcProps p_{};
    std::vector<Species> species_;
    std::vector<Reaction> reactions_;
};

} // namespace simall::combustion
