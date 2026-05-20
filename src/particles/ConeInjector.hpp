// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ConeInjector.hpp
// Phase  : 13.12 — Solid / hollow cone injector with swirl.
//
// Each parcel is launched at an angle drawn from a half-cone with
//
//   solid:   θ ∈ [0, θ_full/2]        uniform in cos θ
//   hollow:  θ ∈ [θ_inner, θ_outer]   uniform in cos θ
//
// azimuth uniform in [0, 2π], optionally biased by a swirl angle φ_sw that
// rotates the velocity tangentially around the cone axis:
//
//   v = U · [cos θ · ê_axis + sin θ · (cos(φ + φ_sw) ê_a + sin(φ + φ_sw) ê_b)]
//
// Diameter sourced from a constant or pluggable distribution (e.g. via a
// caller-supplied lambda from RosinRammlerInjector::sample_diameter).
// =============================================================================
#pragma once

#include "particles/LagrangianTracker.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <functional>
#include <random>

namespace simall::particles
{

enum class ConeKind
{
    Solid,
    Hollow
};

struct ConeInjectorProps
{
    ConeKind kind = ConeKind::Solid;
    double theta_outer = 0.35; // [rad] half-angle of outer cone
    double theta_inner = 0.0;  // hollow-cone inner half-angle
    double swirl_angle = 0.0;  // [rad] tangential bias
    double speed = 50.0;       // [m/s]
    util::Vec3d origin{0, 0, 0};
    util::Vec3d axis{1, 0, 0};
    double diameter_const = 5.0e-5; // fallback when no sampler supplied
    double rho_l = 998.2;
    std::uint64_t rngSeed = 0xABBA'FACE;
};

class ConeInjector
{
public:
    /// Diameter sampler (optional). If null, p_.diameter_const is used.
    using DiameterSampler = std::function<double()>;

    void initialize(ConeInjectorProps props, DiameterSampler sampler = nullptr);

    std::size_t inject(double dt, std::size_t Nparcels, LagrangianTracker& tracker);

    const ConeInjectorProps& props() const noexcept { return p_; }

private:
    ConeInjectorProps p_{};
    DiameterSampler diaSampler_;
    std::mt19937_64 rng_;
};

} // namespace simall::particles
