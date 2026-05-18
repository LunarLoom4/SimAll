// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KOmegaSST.hpp
// Phase  : 8.2 — Menter (1994) k-ω SST turbulence model.
//
// Two transport equations:
//   ∂(ρk)/∂t + ∇·(ρUk)  = ∇·[(μ + σ_k μ_t) ∇k]  + P_k  - β* ρ k ω
//   ∂(ρω)/∂t + ∇·(ρUω) = ∇·[(μ + σ_ω μ_t) ∇ω]
//                       + α ρ S²  - β ρ ω²  +  2(1-F1) ρ σ_ω2 / ω ∇k·∇ω
//
// μ_t = ρ a1 k / max(a1 ω, S F2)
//
// All constants per Menter, Kuntz & Langtry (2003): SST 2003 revision.
// Fully 3-D (uses 3-D velocity gradient, 3-D wall distance, 3-D mesh).
// =============================================================================
#pragma once

#include "turbulence/ITurbulenceModel.hpp"
#include "solver/Solver.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/WallDistance.hpp"
#include "solver/LinearSolvers.hpp"

#include <memory>
#include <vector>

namespace simall::turbulence {

class KOmegaSST_Full final : public ITurbulenceModel {
public:
    KOmegaSST_Full();
    std::string name() const override { return "kOmegaSST"; }

    void  initialize(meshing::Mesh& mesh, solver::FieldRegistry& f) override;
    void  solve(double dt, solver::FieldRegistry& f) override;
    double turbulent_viscosity(std::size_t cellId) const override;

    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }
    void set_density(double rho)   { rho_ = rho; }
    void set_viscosity(double mu)  { mu_  = mu; }

private:
    void compute_strain_rate(const solver::FieldRegistry& f);
    void compute_blending();
    void update_mut(const solver::FieldRegistry& f);

    meshing::Mesh*                  mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    std::unique_ptr<solver::ILinearSolver> lin_;
    std::unique_ptr<solver::IWallDistance> wallDist_;
    std::unique_ptr<solver::ScalarTransport> kEq_;
    std::unique_ptr<solver::ScalarTransport> wEq_;

    util::aligned_vector<double> S_;       // strain-rate magnitude
    util::aligned_vector<double> F1_;      // SST blending function
    util::aligned_vector<double> F2_;
    util::aligned_vector<double> mut_;     // turbulent viscosity (μ_t / ρ → ν_t? we store μ_t)
    util::aligned_vector<double> Pk_;      // production source
    util::aligned_vector<double> CDk_;     // cross-diffusion term for ω eq

    double rho_ = 1.0;
    double mu_  = 1.8e-5;
};

}  // namespace simall::turbulence
