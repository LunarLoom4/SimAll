// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/SasSst.hpp
// Phase  : 8.7 — Scale-Adaptive Simulation built on k-ω SST (Menter-Egorov 2010).
//
// Adds a Q_SAS source term to the ω-equation derived from the von-Kármán
// length scale L_vk = κ |S| / |∇²U|.  In unsteady regions where the
// resolved-scale length scale L = √k / (C_μ^{1/4} ω) becomes large relative
// to L_vk, Q_SAS activates and reduces μ_t, allowing LES-like resolution
// without an explicit grid-based switch.
//
//   Q_SAS  =  ρ max( ζ₂ κ S² (L/L_vk)²
//                    - C ρ k / σ_φ · max(|∇k|²/k², |∇ω|²/ω²) ,  0 )
//
//   ζ₂ = 3.51,  σ_φ = 2/3,  C = 2.
// =============================================================================
#pragma once

#include "solver/LinearSolvers.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/Solver.hpp"
#include "solver/WallDistance.hpp"
#include "turbulence/ITurbulenceModel.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <vector>

namespace simall::turbulence
{

class SasSst_Full final : public ITurbulenceModel
{
public:
    std::string name() const override { return "SAS-SST"; }
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
    double rho_ = 1.0, mu_ = 1.8e-5;

    std::unique_ptr<solver::ILinearSolver> lin_;
    std::unique_ptr<solver::IWallDistance> wallDist_;
    std::unique_ptr<solver::ScalarTransport> kEq_;
    std::unique_ptr<solver::ScalarTransport> wEq_;

    util::aligned_vector<double> mut_;
    util::aligned_vector<double> F1_, F2_;
};

} // namespace simall::turbulence
