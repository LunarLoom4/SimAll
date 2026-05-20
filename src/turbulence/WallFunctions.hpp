// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/WallFunctions.hpp
// Phase  : 9.7 — Wall-function library for RANS / wall-modelled LES.
//
// Implements the universal velocity-temperature profile across the inner
// boundary layer using:
//
//   - Spalding (1961) blended law-of-the-wall (single expression valid
//     0 ≤ y⁺ < ~300, no kink at y⁺=11.6):
//
//        y⁺ = u⁺ + e^{-κB} ( e^{κu⁺} − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6 )
//
//   - Reichardt (1951) alternative for energy/Prandtl scaling.
//   - Werner-Wengle (1991) two-layer power law for LES wall stress.
//   - Kader (1981) thermal wall function:
//
//        T⁺ = Pr · y⁺ · e^{-Γ} + [2.12 ln(y⁺) + β(Pr)] · e^{-1/Γ}
//        β(Pr) = (3.85 Pr^{1/3} − 1.3)² + 2.12 ln(Pr)
//        Γ     = 0.01 (Pr y⁺)⁴ / (1 + 5 Pr³ y⁺)
//
// All functions are scalar in/out and stateless; safe for OpenMP/SIMD use
// inside boundary-condition loops.
// =============================================================================
#pragma once

namespace simall::turbulence::wallfn
{

inline constexpr double kKappa = 0.41;
inline constexpr double kB = 5.5;

/// Solve Spalding's law for y⁺ given u⁺ (closed form).
double spalding_yplus(double uplus);

/// Solve Spalding's law for u⁺ given y⁺ via Newton iteration. Robust on
/// 0 ≤ y⁺ ≤ 1e6; converges in ≤ 8 iterations for default tolerance.
double spalding_uplus(double yplus, double tol = 1e-9, int maxIter = 32);

/// Werner-Wengle (LES): returns wall shear stress τ_w given near-wall
/// tangential velocity u_p at distance y_p and fluid (μ, ρ).
double werner_wengle_tau(double u_p, double y_p, double rho, double mu);

/// Reichardt's law u⁺(y⁺).
double reichardt_uplus(double yplus);

/// Kader's blended thermal wall function: T⁺ as a function of y⁺ and Pr.
double kader_tplus(double yplus, double Pr);

/// Friction velocity u_τ from the parallel velocity u_p at wall distance y_p
/// (Newton-Raphson on Spalding law). Returns 0 when u_p == 0.
double friction_velocity(
    double u_p, double y_p, double rho, double mu, double tol = 1e-9, int maxIter = 32);

} // namespace simall::turbulence::wallfn
