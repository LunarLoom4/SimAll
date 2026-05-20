// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/EnergyEquation.hpp
// Phase  : 9 — Temperature (incompressible thermal) transport.
//
//   ∂(ρ c_p T)/∂t + ∇·(ρ c_p U T) = ∇·(k_eff ∇T) + S_T
//   k_eff = k + c_p μ_t / Pr_t      (turbulent Reynolds analogy)
//
// Implemented as a thin orchestrator on top of ScalarTransport. Reads the
// turbulent viscosity (mut) field if present and the molecular conductivity
// from the material database; supports Dirichlet, Neumann (heat flux), and
// Robin (convective film) wall BCs.
//
// Fully 3-D polyhedral; no symmetry assumptions.
// =============================================================================
#pragma once

#include "solver/ScalarTransport.hpp"
#include "solver/Solver.hpp"

#include <memory>
#include <vector>

namespace simall::solver
{

struct EnergyOptions
{
    double rho = 1.0;
    double cp = 1006.0; // [J/(kg·K)] air default
    double k = 0.0262;  // [W/(m·K)]
    double PrT = 0.85;  // turbulent Prandtl number
    double urf = 0.9;
};

class EnergyEquation
{
public:
    EnergyEquation(meshing::Mesh& mesh,
                   FieldRegistry& fields,
                   ILinearSolver& linear,
                   const std::vector<BoundarySpec>& bcs,
                   EnergyOptions opt);

    /// One outer iteration. Returns L2 residual of the discrete system.
    double iterate();

    /// Override an interface zone's BC at runtime (used by CHT coupler).
    void set_zone_bc(ScalarBC bc) { Teq_->set_zone_bc(bc); }
    /// Mean temperature on a given boundary zone, computed from face owner
    /// cells of that zone (Dirichlet handoff for partitioned coupling).
    double mean_zone_temperature(meshing::ZoneId z) const;
    /// Mean heat flux (W/m²) flowing OUT of the fluid through the given
    /// zone (sign convention: positive = fluid loses heat). Uses k_eff.
    double mean_zone_heat_flux(meshing::ZoneId z) const;

private:
    meshing::Mesh& mesh_;
    FieldRegistry& F_;
    EnergyOptions opt_;
    std::unique_ptr<ScalarTransport> Teq_;
};

} // namespace simall::solver
