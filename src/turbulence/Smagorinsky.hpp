// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/Smagorinsky.hpp
// Phase  : 8.5 — Classical Smagorinsky LES sub-grid model.
//
//   μ_sgs = ρ (C_s Δ)² |S|,    |S| = sqrt(2 S_ij S_ij)
//   Δ     = V_cell^{1/3}                          (filter width)
//
// Optional van-Driest near-wall damping:
//   D = 1 - exp(-y+/A+),  A+ = 26
// Algebraic — no transport equations.
// =============================================================================
#pragma once

#include "solver/Solver.hpp"
#include "solver/WallDistance.hpp"
#include "turbulence/ITurbulenceModel.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <vector>

namespace simall::turbulence
{

class Smagorinsky_LES final : public ITurbulenceModel
{
public:
    std::string name() const override { return "smagorinsky"; }
    void initialize(meshing::Mesh&, solver::FieldRegistry&) override;
    void solve(double dt, solver::FieldRegistry&) override;
    double turbulent_viscosity(std::size_t c) const override
    {
        return c < mut_.size() ? mut_[c] : 0.0;
    }

    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }
    void set_density(double rho) { rho_ = rho; }
    void set_viscosity(double mu) { mu_ = mu; }
    void set_smagorinsky_constant(double cs) { Cs_ = cs; }
    void set_van_driest(bool enable) { vanDriest_ = enable; }

private:
    meshing::Mesh* mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3;
    double Cs_ = 0.17;
    bool vanDriest_ = true;

    std::unique_ptr<solver::IWallDistance> wallDist_;
    util::aligned_vector<double> mut_;
    util::aligned_vector<double> delta_;
};

} // namespace simall::turbulence
