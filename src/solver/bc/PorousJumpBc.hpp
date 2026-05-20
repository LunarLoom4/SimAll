// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/PorousJumpBc.hpp
// Phase  : 16.7 — Porous-jump (thin perforated plate) Δp = (μ/α + ½ C2 ρ|v|) v Δm.
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"

namespace simall::solver::bc
{

struct PorousJumpParams
{
    double permeability = 1.0e-7;     // alpha [m^2]    (Darcy term)
    double inertialCoefficient = 0.0; // C2     [1/m]   (Forchheimer)
    double thickness = 1.0e-3;        // Δm     [m]
    double porosity = 1.0;            // 0 < ε ≤ 1
};

class PorousJumpBc : public IBoundaryCondition
{
public:
    explicit PorousJumpBc(PorousJumpParams p = {}) : p_(p) {}

    BcKind kind() const noexcept override { return BcKind::PorousJump; }
    const char* name() const noexcept override { return "PorousJump"; }

    std::size_t apply(BcContext& ctx) override;
    std::unique_ptr<IBoundaryCondition> clone() const override
    {
        auto c = std::make_unique<PorousJumpBc>(p_);
        c->setZone(zone());
        return c;
    }

    PorousJumpParams& params() noexcept { return p_; }
    const PorousJumpParams& params() const noexcept { return p_; }

private:
    PorousJumpParams p_;
};

} // namespace simall::solver::bc
