// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/SymmetryBc.hpp
// Phase  : 16.4 — Symmetry boundary (mirror plane).
//   Normal velocity   : zero (Dirichlet)
//   Tangential vel    : zero gradient
//   All scalars       : zero gradient
// Implemented as a pure no-op for diagonal-dominant scalars and as a
// projected-velocity Dirichlet for the normal component.
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"

namespace simall::solver::bc
{

class SymmetryBc : public IBoundaryCondition
{
public:
    BcKind kind() const noexcept override { return BcKind::Symmetry; }
    const char* name() const noexcept override { return "Symmetry"; }

    std::size_t apply(BcContext& ctx) override;
    std::unique_ptr<IBoundaryCondition> clone() const override
    {
        auto c = std::make_unique<SymmetryBc>();
        c->setZone(zone());
        return c;
    }
};

} // namespace simall::solver::bc
