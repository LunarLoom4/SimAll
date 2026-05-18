// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KKLOmegaTransition.hpp
// Phase  : 8.7 — Walters-Cokljat (2008) k_T - k_L - ω three-equation
// transition-sensitive turbulence model (a.k.a. k-kL-ω).
//
//   k_T  : turbulent kinetic energy (resolved fluctuations)
//   k_L  : laminar kinetic energy   (pre-transitional fluctuations)
//   ω    : specific dissipation rate
//
// Energy is transferred from k_L to k_T through bypass (R_BP) and natural
// (R_NAT) transition production terms.  The model uses a small-/large-
// scale energy split, a fully algebraic eddy viscosity μ_t = μ_t,s + μ_t,l,
// and engages an additional dissipation length scale λ_eff.
//
//   ∂(ρk_T)/∂t + ∇·(ρU k_T) = ∇·((ν + α_T/σ_k) ∇k_T) + P_kT + R_BP + R_NAT − ω k_T − D_T
//   ∂(ρk_L)/∂t + ∇·(ρU k_L) = ∇·(ν ∇k_L)             + P_kL − R_BP − R_NAT − D_L
//   ∂(ρω)/∂t  + ∇·(ρU ω)    = ∇·((ν + α_T/σ_ω) ∇ω)   + C_ω1 ω/k_T P_kT
//                            + (C_ωR/f_W − 1) ω (R_BP+R_NAT)/k_T − C_ω2 ω²
//                            + C_ω3 f_ω α_T f_W² √k_T / d³
//
// Default constants are from Walters & Cokljat 2008 (J. Fluids Eng. 130 121401).
// =============================================================================
#pragma once

#include "turbulence/ITurbulenceModel.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/Solver.hpp"
#include "solver/WallDistance.hpp"
#include "solver/LinearSolvers.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <vector>

namespace simall::turbulence {

class KKLOmegaTransition_Full final : public ITurbulenceModel {
public:
    std::string name() const override { return "kKLOmegaTransition"; }
    void   initialize(meshing::Mesh& m, solver::FieldRegistry& f) override;
    void   solve(double dt, solver::FieldRegistry& f) override;
    double turbulent_viscosity(std::size_t c) const override
        { return c < mut_.size() ? mut_[c] : 0.0; }

    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }
    void set_density(double rho)   { rho_ = rho; }
    void set_viscosity(double mu)  { mu_  = mu;  }

private:
    meshing::Mesh*                           mesh_ = nullptr;
    std::vector<solver::BoundarySpec>        bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3;

    std::unique_ptr<solver::ILinearSolver>   lin_;
    std::unique_ptr<solver::IWallDistance>   wallDist_;
    std::unique_ptr<solver::ScalarTransport> kTEq_;
    std::unique_ptr<solver::ScalarTransport> kLEq_;
    std::unique_ptr<solver::ScalarTransport> wEq_;

    util::aligned_vector<double> mut_;
};

}  // namespace simall::turbulence
