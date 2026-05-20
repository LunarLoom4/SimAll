// =============================================================================
// SimAll Beta - Heat Transfer Subsystem
// File   : src/heat_transfer/ConjugateHeatTransfer.hpp
// Phase  : 9.4 — Conjugate Heat Transfer (fluid ↔ solid coupling).
//
// Each "solid region" carries its own independent mesh and its own
// transient conduction equation:
//
//     ρ_s c_s ∂T_s/∂t = ∇·(k_s ∇T_s) + Q_v
//
// At each fluid-solid interface zone-pair, a partitioned Dirichlet-Neumann
// coupling is applied (Verstraete & Vilmin, 2014):
//   Fluid side : Dirichlet  T_f|Γ = T_s|Γ           (most recent solid T)
//   Solid side : Neumann   q_s|Γ = -k_f ∂T_f/∂n     (heat flux from fluid)
//
// The coupler iterates fluid + solids until interface temperature & flux
// converge to a user-specified tolerance, providing fully implicit CHT for
// industrial use (turbine cooling, electronics, manifolds, motors…).
//
// Pure 3-D polyhedral. No symmetry / 1-D assumption.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/EnergyEquation.hpp"
#include "solver/Solver.hpp"

#include <memory>
#include <string>
#include <vector>

namespace simall::heat
{

struct SolidMaterial
{
    double rho = 7850.0;
    double cp = 460.0;
    double k = 50.0; // [W/(m·K)]   (e.g. steel)
};

struct InterfacePair
{
    meshing::ZoneId fluidZone;
    meshing::ZoneId solidZone;
};

class SolidRegion
{
public:
    SolidRegion(meshing::Mesh mesh, SolidMaterial mat, std::string name);
    void initialize(double T0);
    /// One implicit step (size dt). Wall BCs default Neumann (adiabatic);
    /// the CHT coupler overrides interface zone fluxes before calling.
    double step(double dt);

    meshing::Mesh& mesh() { return mesh_; }
    solver::FieldRegistry& fields() { return F_; }
    const std::string& name() const { return name_; }
    double conductivity() const { return mat_.k; }
    /// Set Neumann flux (W/m²) on every face of a given boundary zone.
    void set_interface_flux(meshing::ZoneId z, double q);
    /// Mean interface temperature (used by fluid-side Dirichlet handoff).
    double interface_temperature(meshing::ZoneId z) const;

private:
    meshing::Mesh mesh_;
    solver::FieldRegistry F_;
    SolidMaterial mat_;
    std::string name_;
    std::vector<solver::BoundarySpec> bcs_;
    std::unique_ptr<solver::ILinearSolver> lin_;
    std::unique_ptr<solver::EnergyEquation> Teq_;
};

class ConjugateHeatTransfer
{
public:
    void add_solid(std::unique_ptr<SolidRegion> r) { solids_.push_back(std::move(r)); }
    void add_interface(InterfacePair p) { ifaces_.push_back(p); }

    /// Drive a partitioned Dirichlet-Neumann iteration to convergence.
    /// `fluidEnergy` is the already-configured fluid-side energy solver.
    void step(double dt,
              solver::EnergyEquation& fluidEnergy,
              meshing::Mesh& fluidMesh,
              solver::FieldRegistry& fluidFields,
              double tol = 1e-4,
              int maxOuter = 20);

private:
    std::vector<std::unique_ptr<SolidRegion>> solids_;
    std::vector<InterfacePair> ifaces_;
};

} // namespace simall::heat
