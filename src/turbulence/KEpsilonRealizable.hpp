// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KEpsilonRealizable.hpp
// Phase  : 8.7 — Realizable k-ε (Shih, Liou, Shabbir, Yang, Zhu 1995).
//
// "Realisable" because it satisfies certain mathematical constraints on the
// Reynolds stresses consistent with the physics of turbulent flows (a
// positive ⟨u'u'⟩ is guaranteed when Cμ is computed from local strain & rotation).
//
//   Cμ = 1 / (A_0 + A_s · k U* / ε),       U* = √(S_ij S_ij + Ω̃_ij Ω̃_ij)
//   A_0 = 4.04,    A_s = √6 cos(φ),
//   φ   = (1/3) arccos( √6 · S_ij S_jk S_ki / (S_ij S_ij)^{3/2} )
//
// ε equation:
//   ∂(ρε)/∂t + ∇·(ρUε) = ∇·((μ+μ_t/σ_ε)∇ε)
//                       + ρ C1 S ε  − ρ C2 ε² / (k + √(νε))
//   C1 = max(0.43, η/(η+5)),   η = S k / ε.
//
// Constants: σ_k = 1.0, σ_ε = 1.2, C2 = 1.9.
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

class KEpsilonRealizable_Full final : public ITurbulenceModel
{
public:
    std::string name() const override { return "kEpsilonRealizable"; }
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
};

} // namespace simall::turbulence
