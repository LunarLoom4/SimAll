// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/SpalartAllmaras.hpp
// Phase  : 8.4 — Spalart-Allmaras one-equation RANS model.
//
// Reference: Spalart & Allmaras, "A one-equation turbulence model for
// aerodynamic flows", AIAA-92-0439.
//
//   ∂(ρν̃)/∂t + ∇·(ρUν̃) = (1/σ)∇·((ν+ν̃)∇ν̃) + (cb2/σ) ρ |∇ν̃|²
//                       + cb1 ρ Ŝ ν̃ - cw1 ρ fw (ν̃/d)²
//
//   μ_t = ρ ν̃ fv1,   fv1 = χ³/(χ³+cv1³),   χ = ν̃/ν
//
// Fully 3-D on arbitrary polyhedra. Owns its own LinearSolver and
// WallDistance computation (same architectural decision as KOmegaSST).
// =============================================================================
#pragma once

#include "turbulence/ITurbulenceModel.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/WallDistance.hpp"
#include "solver/Solver.hpp"
#include "solver/LinearSolvers.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <vector>

namespace simall::turbulence {

class SpalartAllmaras_Full final : public ITurbulenceModel {
public:
    std::string name() const override { return "SpalartAllmaras"; }
    void   initialize(meshing::Mesh& m, solver::FieldRegistry& f) override;
    void   solve(double dt, solver::FieldRegistry& f) override;
    double turbulent_viscosity(std::size_t c) const override
        { return c < mut_.size() ? mut_[c] : 0.0; }

    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }
    void set_density(double rho)   { rho_ = rho; }
    void set_viscosity(double mu)  { mu_  = mu;  }

private:
    meshing::Mesh*               mesh_   = nullptr;
    solver::FieldRegistry*       fields_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3;

    std::unique_ptr<solver::ILinearSolver>      lin_;
    std::unique_ptr<solver::IWallDistance>      wallDist_;
    std::unique_ptr<solver::ScalarTransport>    nuTildeEq_;

    util::aligned_vector<double> mut_;
    util::aligned_vector<double> Stilde_;
    util::aligned_vector<double> src_;     // explicit production - destruction
};

}  // namespace simall::turbulence
