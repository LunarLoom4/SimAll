// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/LRR_Reynolds_Stress.hpp
// Phase  : 8.7 — Launder-Reece-Rodi (LRR-IP) Reynolds Stress Model (1975).
//
// Solves seven transport equations: six independent components of the
// symmetric Reynolds-stress tensor R_ij = ⟨u'_i u'_j⟩ (R11, R22, R33, R12,
// R13, R23) plus the turbulent dissipation rate ε.
//
//   ∂(ρR_ij)/∂t + ∇·(ρU R_ij) = ∇·(D_ij,ν+D_ij,t) + P_ij + Φ_ij − ε_ij
//
//   Production P_ij = − R_ik ∂U_j/∂x_k − R_jk ∂U_i/∂x_k
//   Pressure-strain  Φ_ij = Φ_ij,1 + Φ_ij,2
//     Φ_ij,1 = −C1 ρ ε/k (R_ij − ⅔ δ_ij k)         (Rotta return-to-isotropy)
//     Φ_ij,2 = −C2 (P_ij − ⅔ δ_ij P)               (IP — Isotropisation of Production)
//   Dissipation ε_ij = ⅔ δ_ij ε                     (isotropic-ε assumption)
//
//   Constants:  C1 = 1.8,  C2 = 0.6,  Cμ = 0.09,   σ_k = 1.0,  σ_ε = 1.3,
//               C1ε = 1.44,  C2ε = 1.92.
//
//   Effective viscosity for the momentum equation can use the trace
//   μ_t = ρ Cμ k²/ε with k = ½ trace(R_ij).
// =============================================================================
#pragma once

#include "solver/LinearSolvers.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/Solver.hpp"
#include "turbulence/ITurbulenceModel.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <array>
#include <memory>
#include <vector>

namespace simall::turbulence
{

class LRR_Reynolds_Stress_Full final : public ITurbulenceModel
{
public:
    std::string name() const override { return "LRR-RSM"; }
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
    static constexpr std::size_t NR = 6; // R11,R22,R33,R12,R13,R23

    meshing::Mesh* mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3;

    std::unique_ptr<solver::ILinearSolver> lin_;
    std::array<std::unique_ptr<solver::ScalarTransport>, NR> Req_;
    std::unique_ptr<solver::ScalarTransport> eEq_;

    util::aligned_vector<double> mut_;
};

} // namespace simall::turbulence
