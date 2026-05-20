// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/VolumeOfFluid.hpp
// Phase  : 12.1 — Volume of Fluid (VOF) for two immiscible fluids.
//
// Solves an advection equation for the volume fraction α ∈ [0,1]:
//   ∂α/∂t + ∇·(U α) + ∇·(U_r α (1-α)) = 0
//
// The artificial-compression term (∇·(U_r α(1-α))) is the MULES-style sharp-
// interface treatment used by OpenFOAM. Mixture density and viscosity used by
// the momentum equation:
//   ρ = α ρ_1 + (1-α) ρ_2
//   μ = α μ_1 + (1-α) μ_2
//
// Fully 3-D on arbitrary polyhedra. Honours the same BC table as the rest
// of the solver subsystem.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/Solver.hpp"

#include <memory>
#include <vector>

namespace simall::multiphase
{

struct VOFFluid
{
    double density;
    double viscosity;
};

class VolumeOfFluid
{
public:
    VolumeOfFluid(meshing::Mesh& mesh,
                  solver::FieldRegistry& fields,
                  solver::ILinearSolver& linear);

    void configure(VOFFluid primary,
                   VOFFluid secondary,
                   const std::vector<solver::BoundarySpec>& bcs);

    /// Enable Brackbill 1992 continuum surface force model.
    void set_surface_tension(double sigma) { sigma_ = sigma; }

    /// Advance α by one outer iteration. Updates ρ, μ in the registry,
    /// recomputes the body-force field "f_surfTension" (per unit volume).
    void step();

private:
    meshing::Mesh& mesh_;
    solver::FieldRegistry& F_;
    solver::ILinearSolver& lin_;
    VOFFluid f1_{1000.0, 1.0e-3};
    VOFFluid f2_{1.225, 1.8e-5};
    double sigma_ = 0.0; // [N/m] surface tension
    std::vector<solver::BoundarySpec> bcs_;
    std::unique_ptr<solver::ScalarTransport> alphaEq_;
};

} // namespace simall::multiphase
