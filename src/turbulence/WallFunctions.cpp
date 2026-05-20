// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/WallFunctions.cpp
// =============================================================================
#include "turbulence/WallFunctions.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence::wallfn
{

double spalding_yplus(double uplus)
{
    const double ku = kKappa * uplus;
    const double t = std::exp(ku) - 1.0 - ku - 0.5 * ku * ku - (ku * ku * ku) / 6.0;
    return uplus + std::exp(-kKappa * kB) * t;
}

double spalding_uplus(double yplus, double tol, int maxIter)
{
    if (yplus <= 0.0)
        return 0.0;
    // Initial guess: log law for large y⁺, linear for small.
    double u = (yplus < 11.6) ? yplus : (std::log(std::max(1e-12, yplus)) / kKappa + kB);
    for (int it = 0; it < maxIter; ++it) {
        const double ku = kKappa * u;
        const double eku = std::exp(ku);
        const double f =
            u + std::exp(-kKappa * kB) * (eku - 1.0 - ku - 0.5 * ku * ku - (ku * ku * ku) / 6.0)
            - yplus;
        const double fp = 1.0 + std::exp(-kKappa * kB) * kKappa * (eku - 1.0 - ku - 0.5 * ku * ku);
        const double du = f / std::max(1e-30, fp);
        u -= du;
        if (std::abs(du) < tol)
            break;
    }
    return std::max(0.0, u);
}

double reichardt_uplus(double yplus)
{
    if (yplus <= 0.0)
        return 0.0;
    const double a = std::log(1.0 + 0.4 * yplus) / kKappa;
    const double b =
        7.8 * (1.0 - std::exp(-yplus / 11.0) - (yplus / 11.0) * std::exp(-yplus / 3.0));
    return a + b;
}

double werner_wengle_tau(double u_p, double y_p, double rho, double mu)
{
    // Two-layer power law:
    //   |u_p| ≤ μ/(2ρy_p) (A^(1/(1-B)))² : viscous
    //   else                              : power
    constexpr double A = 8.3, B = 1.0 / 7.0;
    const double nu = mu / std::max(1e-30, rho);
    const double up = std::abs(u_p);
    const double thresh = 0.5 * nu / std::max(1e-30, y_p) * std::pow(A, 2.0 / (1.0 - B));
    double tau;
    if (up <= thresh) {
        tau = 2.0 * mu * up / std::max(1e-30, y_p);
    } else {
        const double term1 = (1.0 - B) / 2.0 * std::pow(A, (1.0 + B) / (1.0 - B))
                             * std::pow(nu / std::max(1e-30, y_p), 1.0 + B);
        const double term2 = (1.0 + B) / A * std::pow(nu / std::max(1e-30, y_p), B) * up;
        tau = rho * std::pow(term1 + term2, 2.0 / (1.0 + B));
    }
    return (u_p >= 0.0) ? tau : -tau;
}

double kader_tplus(double yplus, double Pr)
{
    if (yplus <= 0.0)
        return 0.0;
    const double beta = std::pow(3.85 * std::cbrt(Pr) - 1.3, 2) + 2.12 * std::log(Pr);
    const double Gamma = (0.01 * std::pow(Pr * yplus, 4)) / (1.0 + 5.0 * std::pow(Pr, 3) * yplus);
    return Pr * yplus * std::exp(-Gamma)
           + (2.12 * std::log(yplus) + beta) * std::exp(-1.0 / std::max(1e-30, Gamma));
}

double friction_velocity(double u_p, double y_p, double rho, double mu, double tol, int maxIter)
{
    if (u_p == 0.0 || y_p <= 0.0)
        return 0.0;
    const double nu = mu / std::max(1e-30, rho);
    const double sign = (u_p < 0) ? -1.0 : 1.0;
    const double up_abs = std::abs(u_p);
    // Initial guess from log law assuming y⁺ > 30.
    double utau = std::max(
        1e-12, up_abs * kKappa / std::log(std::max(1.1, up_abs * y_p * kKappa / nu) + 1e-9));
    for (int it = 0; it < maxIter; ++it) {
        const double yp = utau * y_p / nu;
        const double upGuess = spalding_uplus(yp, 1e-12, 64);
        const double f = utau * upGuess - up_abs;
        // df/du_τ ≈ u⁺ + u_τ · ∂u⁺/∂u_τ ; use finite difference.
        const double du = 1e-6 * utau;
        const double yp2 = (utau + du) * y_p / nu;
        const double up2 = spalding_uplus(yp2, 1e-12, 64);
        const double fp = ((utau + du) * up2 - up_abs - f) / du;
        const double step = f / std::max(1e-30, std::abs(fp)) * (fp < 0 ? -1.0 : 1.0);
        utau -= step;
        if (utau < 1e-12)
            utau = 1e-12;
        if (std::abs(step) < tol * utau)
            break;
    }
    return sign * utau;
}

} // namespace simall::turbulence::wallfn
