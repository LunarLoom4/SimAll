// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/DynamicSmagorinsky.hpp
// Phase  : 8.7 — Germano-Lilly Dynamic Smagorinsky LES.
//
// Computes the Smagorinsky coefficient Cs²(x,t) locally from the resolved
// scales via a test-filter at 2Δ.  Using the Germano identity:
//
//   L_ij  =  ũ_i u_j  −  ũ_i ũ_j               (resolved-stress)
//   M_ij  =  2 Δ²  |S| S_ij  −  2 Δ̂²  |S̃| S̃_ij  (with Δ̂ = 2Δ)
//   C_s²  =  ⟨L_ij M_ij⟩ / ⟨M_ij M_ij⟩         (Lilly LS minimisation)
//
// The test filter is a face-neighbour box average.  C_s² is clipped to
// [0, Cs_max²].  An optional spatial average over neighbours is applied
// to stabilise the procedure.
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

class DynamicSmagorinsky_LES final : public ITurbulenceModel
{
public:
    std::string name() const override { return "dynamicSmagorinsky"; }
    void initialize(meshing::Mesh& m, solver::FieldRegistry& f) override;
    void solve(double dt, solver::FieldRegistry& f) override;
    double turbulent_viscosity(std::size_t c) const override
    {
        return c < mut_.size() ? mut_[c] : 0.0;
    }

    void set_boundaries(const std::vector<solver::BoundarySpec>& bcs) { bcs_ = bcs; }
    void set_density(double rho) { rho_ = rho; }
    void set_viscosity(double mu) { mu_ = mu; }
    void set_cs_max(double v) { csMax_ = v; }

private:
    meshing::Mesh* mesh_ = nullptr;
    std::vector<solver::BoundarySpec> bcs_;
    double rho_ = 1.0, mu_ = 1.0e-3;
    double csMax_ = 0.23; // physical upper bound

    util::aligned_vector<double> mut_;
    util::aligned_vector<double> Cs2_;
    util::aligned_vector<double> delta_; // grid filter width
};

} // namespace simall::turbulence
