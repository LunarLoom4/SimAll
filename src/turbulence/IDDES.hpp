// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/IDDES.hpp
// Phase  : 8.7 — Improved Delayed Detached-Eddy Simulation (Shur, Spalart,
// Strelets, Travin 2008).  SA-IDDES with shielding + WMLES branch.
//
//   l_IDDES = f̃_d (1 + f_e) l_RANS + (1 − f̃_d) l_LES
//   f̃_d    = max( (1 − f_dt) , f_d )
//   f_d    = 1 − tanh( (8 r_d)³ )                         (DDES shielding)
//   f_dt   = 1 − tanh( (8 r_dt)³ )                        (WMLES branch)
//   f_e    = max( (f_e1 − 1), 0 ) · ψ · f_e2              (WMLES elevation)
//
//   Δ_IDDES = min( max( C_w d_w, C_w h_max, h_wn ), h_max )
//   l_LES   = C_DES Δ_IDDES,    l_RANS = d_w
//
// Built on Spalart-Allmaras transport equation infrastructure (same as SaDdes).
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

class IDDES_Full final : public ITurbulenceModel
{
public:
    std::string name() const override { return "SA-IDDES"; }
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
    void compute_hmax();

    meshing::Mesh* mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3;

    std::unique_ptr<solver::ILinearSolver> lin_;
    std::unique_ptr<solver::IWallDistance> wallDist_;
    std::unique_ptr<solver::ScalarTransport> nuTildeEq_;

    util::aligned_vector<double> hmax_;  // max edge length per cell
    util::aligned_vector<double> delta_; // cubic-root volume per cell
    util::aligned_vector<double> mut_;
};

} // namespace simall::turbulence
