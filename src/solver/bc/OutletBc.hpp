// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/OutletBc.hpp
// Phase  : 16.3 — Pressure outlet / Outflow (zero-gradient).
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"

namespace simall::solver::bc
{

enum class OutletKind
{
    PressureStatic,
    Outflow
};

struct OutletParams
{
    OutletKind kind = OutletKind::PressureStatic;
    double staticPressure = 0.0; // Pa
    double backflowT = 300.0;    // K (only used if reverse flow)
    double backflowK = 1.0e-3;
    double backflowOmega = 1.0;
};

class OutletBc : public IBoundaryCondition
{
public:
    explicit OutletBc(OutletParams p = {}) : p_(p) {}

    BcKind kind() const noexcept override { return BcKind::Outlet; }
    const char* name() const noexcept override { return "Outlet"; }

    std::size_t apply(BcContext& ctx) override;
    std::unique_ptr<IBoundaryCondition> clone() const override
    {
        auto c = std::make_unique<OutletBc>(p_);
        c->setZone(zone());
        return c;
    }

    OutletParams& params() noexcept { return p_; }
    const OutletParams& params() const noexcept { return p_; }

private:
    OutletParams p_;
};

} // namespace simall::solver::bc
