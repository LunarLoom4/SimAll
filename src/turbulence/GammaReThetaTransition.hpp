// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/GammaReThetaTransition.hpp
// Phase  : 8.5 — Langtry & Menter (2009) γ-Reθ_t correlation-based
// transition model. Two additional transport equations:
//
//   ∂(ρ γ)/∂t      + ∇·(ρ U γ)     = ∇·[(μ + μ_t/σ_f) ∇γ] + P_γ - E_γ
//   ∂(ρ Re̅θ_t)/∂t + ∇·(ρ U Re̅θ_t) = ∇·[σ_θt (μ+μ_t)  ∇Re̅θ_t] + P_θt
//
// Couples to k-ω SST by producing an effective intermittency γ_eff written
// as scalar field "gamma_eff", which the SST production term must multiply.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/WallDistance.hpp"
#include "solver/LinearSolvers.hpp"

#include <memory>
#include <vector>

namespace simall::turbulence {

class GammaReThetaTransition {
public:
    GammaReThetaTransition();

    void initialize(meshing::Mesh& mesh, solver::FieldRegistry& fields,
                    double density, double viscosity);
    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }

    /// One outer iteration: assemble & solve γ and Reθ_t equations, update γ_eff.
    /// Returns max(γ_residual, Reθ_residual).
    double solve_iteration(double dt, solver::FieldRegistry& fields);

private:
    void compute_strain_and_vorticity(const solver::FieldRegistry& f);
    void compute_sources(const solver::FieldRegistry& f);
    void update_gamma_eff(const solver::FieldRegistry& f);
    /// Empirical correlation Reθ_t = f(Tu, λ_θ) for free-stream condition.
    static double correlation_ReThetaT(double Tu, double lambda);
    /// Onset function used in P_γ.
    static double F_onset(double rho, double mu, double mut, double S,
                          double dwall, double ReThetaT);

    meshing::Mesh* mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    std::unique_ptr<solver::ILinearSolver>   lin_;
    std::unique_ptr<solver::ScalarTransport> gammaEq_;
    std::unique_ptr<solver::ScalarTransport> ReThetaEq_;
    std::unique_ptr<solver::IWallDistance>   wd_;

    util::aligned_vector<double> S_, W_;      // strain and vorticity magnitudes
    util::aligned_vector<double> Pgamma_, ReThetaT_;
    double rho_ = 1.0, mu_ = 1.8e-5;
};

}  // namespace simall::turbulence
