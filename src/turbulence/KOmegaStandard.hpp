// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KOmegaStandard.hpp
// Phase  : 8.7 — Wilcox 1988/2006 standard k-ω.
//
// Two transport equations, no F1/F2 blending, no cross-diffusion:
//   ∂(ρk)/∂t + ∇·(ρUk) = ∇·((μ+σ_k μ_t)∇k) + P_k − β* ρ k ω
//   ∂(ρω)/∂t + ∇·(ρUω) = ∇·((μ+σ_ω μ_t)∇ω) + α(ω/k) P_k − β ρ ω²
//
//   μ_t = ρ k / ω
//   α = 5/9,   β = 3/40,   β* = 9/100,   σ_k = σ_ω = 0.5
//
// 2006 stress-limited variant: ω̄ = max(ω, C_lim √(2 S_ij S_ij / β*));
//                              μ_t = ρ k / ω̄.
// =============================================================================
#pragma once

#include "turbulence/ITurbulenceModel.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/Solver.hpp"
#include "solver/LinearSolvers.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <vector>

namespace simall::turbulence {

class KOmegaStandard_Full final : public ITurbulenceModel {
public:
    std::string name() const override { return "kOmegaStandard"; }
    void   initialize(meshing::Mesh& m, solver::FieldRegistry& f) override;
    void   solve(double dt, solver::FieldRegistry& f) override;
    double turbulent_viscosity(std::size_t c) const override
        { return c < mut_.size() ? mut_[c] : 0.0; }

    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }
    void set_density(double rho)   { rho_ = rho; }
    void set_viscosity(double mu)  { mu_  = mu;  }

private:
    meshing::Mesh*                    mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3;

    std::unique_ptr<solver::ILinearSolver>   lin_;
    std::unique_ptr<solver::ScalarTransport> kEq_;
    std::unique_ptr<solver::ScalarTransport> wEq_;

    util::aligned_vector<double> mut_;
};

}  // namespace simall::turbulence
