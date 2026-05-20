// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/SprayBreakup.hpp
// Phase  : 13.4 — Kelvin-Helmholtz / Rayleigh-Taylor (KH-RT) hybrid
// secondary atomisation model (Reitz 1987; Patterson & Reitz 1998).
//
// Operates on an existing LagrangianTracker droplet population. For each
// droplet the model computes:
//
//   KH wave (aerodynamic stripping):
//     Λ_KH = 9.02 r (1 + 0.45 Z^0.5)(1 + 0.4 T^0.7) / (1 + 0.87 We_g^1.67)^0.6
//     Ω_KH = (0.34 + 0.385 We_g^1.5) / ((1 + Z)(1 + 1.4 T^0.6))
//              · √(σ / (ρ_l r³))
//     τ_KH = 3.726 B1 r / (Λ_KH Ω_KH)        (B1 ≈ 1.73 default)
//
//   RT wave (deceleration-driven, only when a_d > 0):
//     Ω_RT² = (2/(3√3)) (|a_d|(ρ_l-ρ_g) / σ)^(3/2) · √(σ/(ρ_l+ρ_g))
//     Λ_RT  = C_RT · √(σ / (a_d(ρ_l - ρ_g)))   (C_RT ≈ 2π/√3 default)
//
//   Whichever timer reaches the cell residence time first triggers a
//   diameter reduction r ← r_child. RT trumps KH when both are mature.
//
// Dimensionless groups (with relative velocity u_r = ‖u - v_p‖):
//   We_g = ρ_g u_r² r / σ            We_l = ρ_l u_r² r / σ
//   Re_l = ρ_l u_r r / μ_l           Z    = √(We_l)/Re_l   (Ohnesorge)
//   T    = Z √(We_g)
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "particles/LagrangianTracker.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

#include <vector>

namespace simall::particles
{

struct SprayBreakupProps
{
    double sigma = 0.072; // surface tension liquid-gas [N/m] (water/air)
    double mu_l = 1.0e-3; // dynamic viscosity of liquid [Pa·s]
    double rho_l = 998.2; // density of liquid [kg/m³]
    double rho_g = 1.225; // density of gas    [kg/m³]
    double B0 = 0.61;     // KH child-radius ratio: r_KH = B0 · Λ_KH
    double B1 = 1.73;     // KH timescale coefficient
    double C_RT = 0.1;    // RT-wave size limiter (× critical Λ_RT)
    double Ctau_RT = 1.0; // RT timescale scaling
    double WeCrit = 6.0;  // KH only active if We_g > WeCrit
};

class SprayBreakup
{
public:
    void initialize(const SprayBreakupProps& props)
    {
        p_ = props;
        rt_age_.clear();
    }

    /// Advance the breakup state by dt for the given tracker / fluid.
    /// `accel` is the per-particle deceleration magnitude (zero ⇒ no RT).
    void apply(double dt, LagrangianTracker& tracker, const solver::FieldRegistry& fields);

    const SprayBreakupProps& props() const noexcept { return p_; }

private:
    SprayBreakupProps p_{};
    std::vector<double> rt_age_; ///< accumulated RT clock per particle
};

} // namespace simall::particles
