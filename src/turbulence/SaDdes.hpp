// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/SaDdes.hpp
// Phase  : 8.6 — Delayed Detached-Eddy Simulation based on Spalart-Allmaras
// (SA-DDES, Spalart et al. 2006).
//
// Modifies the destruction-term length scale of SA from the wall distance
// d to a hybrid length scale d̃ blending RANS near walls with LES away:
//
//   d̃ = d - f_d · max(0, d - C_DES Δ)
//   f_d = 1 - tanh((8 r_d)³),
//   r_d = (ν_t + ν) / (sqrt(U_{i,j} U_{i,j}) κ² d²)
//
// where Δ is the local cubic-root cell-volume scale and C_DES ≈ 0.65.
//
// This module replaces the SA ν̃ destruction term while reusing the
// existing SpalartAllmaras_Full transport equation and source machinery.
// =============================================================================
#pragma once

#include "turbulence/ITurbulenceModel.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/WallDistance.hpp"
#include "solver/LinearSolvers.hpp"
#include "solver/Solver.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <vector>

namespace simall::turbulence {

class SaDdes_Full final : public ITurbulenceModel {
public:
    std::string name() const override { return "SA-DDES"; }
    void   initialize(meshing::Mesh& m, solver::FieldRegistry& f) override;
    void   solve(double dt, solver::FieldRegistry& f) override;
    double turbulent_viscosity(std::size_t c) const override
        { return c < mut_.size() ? mut_[c] : 0.0; }

    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }
    void set_density(double rho)   { rho_ = rho; }
    void set_viscosity(double mu)  { mu_  = mu;  }
    void set_C_DES(double c)       { C_DES_ = c; }

private:
    void compute_length_scales(const solver::FieldRegistry& f);
    void compute_sources(const solver::FieldRegistry& f);

    meshing::Mesh*               mesh_   = nullptr;
    solver::FieldRegistry*       fields_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3, C_DES_ = 0.65;

    std::unique_ptr<solver::ILinearSolver>   lin_;
    std::unique_ptr<solver::IWallDistance>   wallDist_;
    std::unique_ptr<solver::ScalarTransport> nuTildeEq_;

    util::aligned_vector<double> dHybrid_;
    util::aligned_vector<double> delta_;
    util::aligned_vector<double> mut_;
    util::aligned_vector<double> src_;
    util::aligned_vector<double> Smag_;     // strain magnitude
};

}  // namespace simall::turbulence
