// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KEpsilonRng.hpp
// Phase  : 8.7 — RNG k-ε (Yakhot, Orszag, Thangam, Gatski, Speziale 1992).
//
// Derived from Renormalisation Group theory. Differs from standard k-ε in
// the ε source-term coefficient C2*  =  C2 + Cμ η³ (1 - η/η₀) / (1 + β η³)
// where η = S k/ε is a dimensionless strain rate.  RNG handles rapidly
// strained flows, swirling flows, and stagnation regions better than
// standard k-ε.  Also uses lower σ_k, σ_ε and slightly different Cμ.
//
// Transport equations (incompressible, constant density):
//   ∂(ρk)/∂t + ∇·(ρUk) = ∇·((μ+μ_t/σ_k)∇k) + P_k − ρε
//   ∂(ρε)/∂t + ∇·(ρUε) = ∇·((μ+μ_t/σ_ε)∇ε) + (ε/k)(C1 P_k − C2* ρε)
//
// Constants (RNG):
//   Cμ=0.0845, σ_k=σ_ε=0.7194, C1=1.42, C2=1.68, η₀=4.38, β=0.012.
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

class KEpsilonRng_Full final : public ITurbulenceModel {
public:
    std::string name() const override { return "kEpsilonRNG"; }
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
    std::unique_ptr<solver::ScalarTransport> eEq_;

    util::aligned_vector<double> mut_;
    util::aligned_vector<double> Pk_;
};

}  // namespace simall::turbulence
