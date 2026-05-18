// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/LaminarFiniteRate.hpp
// Phase  : 10 — Laminar finite-rate Arrhenius combustion (multi-step).
//
// Solves species transport:
//
//   ∂(ρ Y_k)/∂t + ∇·(ρ U Y_k) = ∇·(ρ D_k ∇Y_k) + ω̇_k
//
// where ω̇_k = ν'_k · q for each elementary reaction with rate
//
//   q = A T^β exp(-E_a / RT) · Π_j [X_j]^{ν'_j}
//
// The energy equation receives an extra heat-release source
//
//   S_T = -Σ_k h_f0,k · ω̇_k          (LHV-based, neglects sensible enthalpy
//                                      variation across species).
//
// Coupling: this class owns one ScalarTransport per species and writes a
// scalar "S_combustion" field consumed by the energy equation.
//
// Fully 3-D polyhedral; thermodynamics use ideal-gas / constant cp per
// species (full NASA-7 polynomial fits arrive in Phase 10.3 hand-off).
// =============================================================================
#pragma once

#include "solver/Solver.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/FieldRegistry.hpp"
#include "meshing/MeshStorage.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <string>
#include <vector>

namespace simall::combustion {

struct Species {
    std::string name;
    double      molarMass;            // [kg/mol]
    double      formationEnthalpy;    // [J/kg]   h_f0
    double      diffusivity;          // [m²/s]   D_k
};

struct Reaction {
    std::vector<double> nuP;          // products  ν'   (size = nSpecies)
    std::vector<double> nuR;          // reactants ν''  (size = nSpecies)
    std::vector<double> order;        // forward-rate order per species
    double A   = 1.0e10;              // [SI cgs-equivalent units, user-supplied]
    double beta = 0.0;
    double Ea  = 0.0;                 // [J/mol]
};

class LaminarFiniteRate {
public:
    LaminarFiniteRate(meshing::Mesh& mesh,
                      solver::FieldRegistry& fields,
                      solver::ILinearSolver& linear);

    void add_species(Species s)                 { species_.push_back(std::move(s)); }
    void add_reaction(Reaction r)               { reactions_.push_back(std::move(r)); }
    void set_density(double rho)                { rho_ = rho; }
    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }

    /// Build per-species transport equations (call after add_species).
    void initialize();

    /// One outer iteration: compute reaction rates, update sources, transport
    /// each Y_k, store combined energy source in "S_combustion".
    void step();

private:
    meshing::Mesh&                            mesh_;
    solver::FieldRegistry&                    F_;
    solver::ILinearSolver&                    lin_;
    double                                    rho_ = 1.0;
    std::vector<Species>                      species_;
    std::vector<Reaction>                     reactions_;
    std::vector<solver::BoundarySpec>         bcs_;
    std::vector<std::unique_ptr<solver::ScalarTransport>> Yeqs_;
    std::vector<util::aligned_vector<double>> srcY_;     // ω̇_k storage
};

}  // namespace simall::combustion
