// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/RosinRammlerInjector.hpp
// Phase  : 13.10 — Rosin-Rammler injection (Rosin-Rammler-Sperling-Bennett).
//
// Samples N parcels with diameters drawn from
//
//   F(d) = 1 - exp[-(d / X)^n]
//
// i.e. PDF f(d) = (n/X) (d/X)^{n-1} exp[-(d/X)^n] for d ∈ [d_min, d_max].
// Inverse-CDF sampling using a uniform U:
//
//   d = X · (-ln(1 - U))^{1/n}
//
// followed by rejection against [d_min, d_max].  Position is sampled
// uniformly over a circular orifice (radius = nozzleRadius) with parcel
// velocity = magnitudeVelocity · axis_unit.
//
// Mass conservation: total injected mass per second = ṁ; the parcel mass
// = (4/3)π (d/2)³ ρ_l, and the number of physical droplets represented by
// the parcel is set so that Σ N_phys · m_p = ṁ · dt.
// =============================================================================
#pragma once

#include "particles/LagrangianTracker.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <random>

namespace simall::particles
{

struct RosinRammlerProps
{
    double X = 5.0e-5; // characteristic diameter [m]
    double n = 3.5;    // spread (shape) parameter
    double d_min = 1.0e-6;
    double d_max = 5.0e-4;
    double rho_l = 998.2;          // liquid density [kg/m³]
    double mass_flowrate = 1.0e-3; // ṁ [kg/s]
    util::Vec3d origin{0, 0, 0};
    util::Vec3d axis{1, 0, 0}; // injection direction (unit-normalised)
    double nozzleRadius = 1.0e-4;
    double speed = 50.0;
    std::uint64_t rngSeed = 0xDEADBEEF;
};

class RosinRammlerInjector
{
public:
    void initialize(RosinRammlerProps props);

    /// Inject `parcelsThisStep` parcels into the tracker over an integration
    /// interval dt; sets diameter + mass + position + velocity per parcel.
    /// Returns the actual number injected.
    std::size_t inject(double dt, std::size_t parcelsThisStep, LagrangianTracker& tracker);

    /// Sample a single diameter via inverse-CDF Rosin-Rammler.
    double sample_diameter();

    const RosinRammlerProps& props() const noexcept { return p_; }

private:
    RosinRammlerProps p_{};
    std::mt19937_64 rng_;
};

} // namespace simall::particles
