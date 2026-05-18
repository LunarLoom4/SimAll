// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/WallBc.hpp
// Phase  : 16.1 — Wall boundary condition (no-slip / slip / moving / heat).
//
//   Velocity         : no-slip (default), free-slip, or specified U_wall
//   Temperature      : isothermal, adiabatic, heat flux, or h(T - T_inf)
//   Turbulence       : standard wall functions (k, omega, epsilon)
//
// The BC is multiplexed: it inspects `ctx.variable` and dispatches to the
// appropriate sub-routine.  This avoids carrying ten separate Wall classes.
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"

namespace simall::solver::bc {

enum class WallVelocityMode  { NoSlip, FreeSlip, MovingWall };
enum class WallThermalMode   { Adiabatic, Isothermal, HeatFlux, Convection };

struct WallParams {
    WallVelocityMode  velocityMode = WallVelocityMode::NoSlip;
    double            U_wall[3]    {0,0,0};       // moving wall velocity
    WallThermalMode   thermalMode  = WallThermalMode::Adiabatic;
    double            T_wall       = 300.0;       // K (isothermal)
    double            heatFlux     = 0.0;         // W/m^2 (positive = into fluid)
    double            h_conv       = 0.0;         // W/(m^2 K)
    double            T_infinity   = 300.0;       // K
    double            roughness    = 0.0;         // m (Cebeci-Bradshaw)
    bool              useWallFunctions = true;
};

class WallBc : public IBoundaryCondition {
public:
    explicit WallBc(WallParams p = {}) : p_(p) {}

    BcKind kind() const noexcept override { return BcKind::Wall; }
    const char* name() const noexcept override { return "Wall"; }

    std::size_t apply(BcContext& ctx) override;
    std::unique_ptr<IBoundaryCondition> clone() const override {
        auto c = std::make_unique<WallBc>(p_); c->setZone(zone()); return c;
    }

    WallParams&       params()       noexcept { return p_; }
    const WallParams& params() const noexcept { return p_; }

private:
    WallParams p_;
};

}  // namespace simall::solver::bc
