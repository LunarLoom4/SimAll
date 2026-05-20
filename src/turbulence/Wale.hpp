// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/Wale.hpp
// Phase  : 8.6 — WALE (Wall-Adapting Local Eddy-viscosity) LES model
//          (Nicoud & Ducros, 1999, Flow Turb. Combust., 62, 183-200).
//
//   μ_sgs = ρ (C_w Δ)² · (S^d:S^d)^{3/2}
//                       --------------------------
//                       (S:S)^{5/2} + (S^d:S^d)^{5/4}
//
//   g_ij  = ∂u_i/∂x_j
//   S_ij  = ½(g_ij + g_ji)                 strain rate
//   S^d_ij= ½(g_ij² + g_ji²) - ⅓ δ_ij g_kk²    traceless symmetric part of g²
//
// WALE recovers the correct cubic near-wall scaling y³ for μ_sgs without
// needing van-Driest damping or wall distance.
// =============================================================================
#pragma once

#include "solver/Solver.hpp"
#include "turbulence/ITurbulenceModel.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <vector>

namespace simall::turbulence
{

class WALE_LES final : public ITurbulenceModel
{
public:
    std::string name() const override { return "WALE"; }
    void initialize(meshing::Mesh&, solver::FieldRegistry&) override;
    void solve(double dt, solver::FieldRegistry&) override;
    double turbulent_viscosity(std::size_t c) const override
    {
        return c < mut_.size() ? mut_[c] : 0.0;
    }

    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }
    void set_density(double rho) { rho_ = rho; }
    void set_viscosity(double mu) { mu_ = mu; }
    void set_wale_constant(double cw) { Cw_ = cw; }

private:
    meshing::Mesh* mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3, Cw_ = 0.325;
    util::aligned_vector<double> mut_;
    util::aligned_vector<double> delta_;
};

} // namespace simall::turbulence
