// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/OversetBc.hpp
// Phase  : 16.9 — Overset boundary (Chimera interpolation receptor).
//
// On overset receptor faces the field value is interpolated from a donor
// mesh via meshing::OversetInterpolation.  At assembly time the BC simply
// imposes a Dirichlet condition with the donor-interpolated value.
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"

#include <vector>

namespace simall::solver::bc {

class OversetBc : public IBoundaryCondition {
public:
    OversetBc() = default;

    BcKind kind() const noexcept override { return BcKind::Overset; }
    const char* name() const noexcept override { return "Overset"; }

    /// Setter: provide the per-receptor-cell donor value for the named
    /// variable. The solver / OversetCouplingDriver populates this map
    /// each outer iteration.
    void setDonorValues(std::vector<std::pair<meshing::CellId,double>> v) {
        donor_ = std::move(v);
    }

    std::size_t apply(BcContext& ctx) override;

    std::unique_ptr<IBoundaryCondition> clone() const override {
        auto c = std::make_unique<OversetBc>();
        c->setZone(zone()); c->donor_ = donor_; return c;
    }

private:
    std::vector<std::pair<meshing::CellId, double>> donor_;
};

}  // namespace simall::solver::bc
