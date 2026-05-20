// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KEpsilonStandard.hpp
// Phase  : 8.2 — Standard k-ε model (Launder & Spalding 1974).
//
//   ∂(ρk)/∂t + ∇·(ρUk) = ∇·((μ + μ_t/σ_k)∇k) + P_k - ρε
//   ∂(ρε)/∂t + ∇·(ρUε) = ∇·((μ + μ_t/σ_ε)∇ε) + ε/k (C1 P_k - C2 ρε)
//
//   μ_t = ρ C_μ k²/ε
//   P_k = μ_t |S|²,  |S|² = 2 S_ij S_ij
//
// Constants: Cμ=0.09, σ_k=1.0, σ_ε=1.3, C1=1.44, C2=1.92.
//
// Wall treatment: standard wall functions (when wallDistance is small the
// near-wall ε is computed from k^{3/2} / (κ y_p) — handled by Dirichlet BC
// on ε on near-wall faces). Fully 3-D polyhedral.
// =============================================================================
#pragma once

#include "solver/LinearSolvers.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/Solver.hpp"
#include "turbulence/ITurbulenceModel.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <vector>

namespace simall::turbulence
{

class KEpsilonStandard_Full final : public ITurbulenceModel
{
public:
    std::string name() const override { return "kEpsilonStandard"; }
    void initialize(meshing::Mesh& m, solver::FieldRegistry& f) override;
    void solve(double dt, solver::FieldRegistry& f) override;
    double turbulent_viscosity(std::size_t c) const override
    {
        return c < mut_.size() ? mut_[c] : 0.0;
    }

    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }
    void set_density(double rho) { rho_ = rho; }
    void set_viscosity(double mu) { mu_ = mu; }

private:
    meshing::Mesh* mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3;

    std::unique_ptr<solver::ILinearSolver> lin_;
    std::unique_ptr<solver::ScalarTransport> kEq_;
    std::unique_ptr<solver::ScalarTransport> eEq_;

    util::aligned_vector<double> mut_;
    util::aligned_vector<double> Pk_;
    util::aligned_vector<double> srcK_;
    util::aligned_vector<double> srcE_;
};

} // namespace simall::turbulence
