// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/AxisymmetricBc.hpp
// Phase  : 16.5 — Axisymmetric axis condition.
//
// 2-D axisymmetric problems: the axis (typically y = 0) requires:
//   * Radial velocity   : v_r = 0     (Dirichlet)
//   * Axial   velocity  : ∂v_z/∂n = 0 (Neumann homogeneous)
//   * All scalars       : ∂φ/∂n = 0
//
// The axis direction is configurable.
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"

namespace simall::solver::bc
{

enum class AxisAlignment
{
    X,
    Y,
    Z
};

struct AxisymmetricParams
{
    AxisAlignment axis = AxisAlignment::X;
};

class AxisymmetricBc : public IBoundaryCondition
{
public:
    explicit AxisymmetricBc(AxisymmetricParams p = {}) : p_(p) {}

    BcKind kind() const noexcept override { return BcKind::Axisymmetric; }
    const char* name() const noexcept override { return "Axisymmetric"; }

    std::size_t apply(BcContext& ctx) override;
    std::unique_ptr<IBoundaryCondition> clone() const override
    {
        auto c = std::make_unique<AxisymmetricBc>(p_);
        c->setZone(zone());
        return c;
    }

    AxisymmetricParams& params() noexcept { return p_; }
    const AxisymmetricParams& params() const noexcept { return p_; }

private:
    AxisymmetricParams p_;
};

} // namespace simall::solver::bc
