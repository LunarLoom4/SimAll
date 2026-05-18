// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/FanBc.hpp
// Phase  : 16.8 — Fan / lumped pressure-rise BC.
//
//   Δp = f(Q)         polynomial pressure-rise vs. volumetric flow rate
//
// Implemented as a porous-jump variant where the jump is a fitted polynomial
// rather than the physical Darcy + Forchheimer formula.
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"
#include "materials/PolynomialFit.hpp"

namespace simall::solver::bc {

class FanBc : public IBoundaryCondition {
public:
    explicit FanBc(materials::PolynomialFit curve)
        : curve_(std::move(curve)) {}

    BcKind kind() const noexcept override { return BcKind::Fan; }
    const char* name() const noexcept override { return "Fan"; }

    std::size_t apply(BcContext& ctx) override;
    std::unique_ptr<IBoundaryCondition> clone() const override {
        auto c = std::make_unique<FanBc>(curve_); c->setZone(zone()); return c;
    }

    const materials::PolynomialFit& curve() const noexcept { return curve_; }

private:
    materials::PolynomialFit curve_;
};

}  // namespace simall::solver::bc
