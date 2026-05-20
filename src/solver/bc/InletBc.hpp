// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/InletBc.hpp
// Phase  : 16.2 — Inlet boundary condition (velocity / pressure / mass-flow).
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"

#include <string>

namespace simall::solver::bc
{

enum class InletKind
{
    Velocity,
    PressureTotal,
    MassFlow
};

struct InletParams
{
    InletKind kind = InletKind::Velocity;
    double velocity[3]{1, 0, 0};      // m/s for VelocityInlet
    double totalPressure = 101325.0;  // Pa  for PressureInlet
    double staticTemperature = 300.0; // K
    double massFlowRate = 1.0;        // kg/s   for MassFlowInlet
    double direction[3]{1, 0, 0};     // unit vector (mass-flow)
    double turbulence_k = 1.0e-3;     // m^2/s^2
    double turbulence_omega = 1.0;    // 1/s
    double turbulence_intensity = 0.05;
    double hydraulicDiameter = 0.0; // m  (turbulence specification)
    std::string speciesProfile;     // optional name
};

class InletBc : public IBoundaryCondition
{
public:
    explicit InletBc(InletParams p = {}) : p_(p) {}

    BcKind kind() const noexcept override { return BcKind::Inlet; }
    const char* name() const noexcept override { return "Inlet"; }

    std::size_t apply(BcContext& ctx) override;
    std::unique_ptr<IBoundaryCondition> clone() const override
    {
        auto c = std::make_unique<InletBc>(p_);
        c->setZone(zone());
        return c;
    }

    InletParams& params() noexcept { return p_; }
    const InletParams& params() const noexcept { return p_; }

private:
    InletParams p_;
};

} // namespace simall::solver::bc
