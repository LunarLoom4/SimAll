// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/RealGasEos.hpp
// Phase  : 15.2 — Real-gas equation of state.
//
//   * IdealGas     : p = ρ R_s T          (R_s = R / M_w)
//   * PengRobinson : cubic EOS with Soave α(T, ω):
//                    p = RT/(v - b) - aα(T)/(v² + 2 b v - b²)
//                    a  = 0.45724 R² T_c² / p_c
//                    b  = 0.07780 R T_c / p_c
//                    α  = [1 + κ (1 - √(T/T_c))]²
//                    κ  = 0.37464 + 1.54226 ω - 0.26992 ω²
//
//   * NASA7Cp      : NASA-7 polynomial for c_p(T) (two temperature ranges):
//                    c_p / R = a1 + a2 T + a3 T² + a4 T³ + a5 T⁴
//                    h / R T = a1 + a2 T/2 + a3 T²/3 + a4 T³/4 + a5 T⁴/5 + a6/T
//
// All routines are header-only structs with constexpr-friendly methods; the
// MaterialDatabase can attach them as PropertyFunctions via lambdas.
// =============================================================================
#pragma once

#include <array>
#include <cmath>

namespace simall::materials {

inline constexpr double kRgasUniversal = 8.314462618;   // J/(mol·K)

struct IdealGasEos {
    double molecularWeight = 0.02897;     // kg/mol (default = air)
    double Rspecific() const noexcept     { return kRgasUniversal / molecularWeight; }
    double density(double T, double p) const noexcept {
        return p / (Rspecific() * T);
    }
    double pressure(double T, double rho) const noexcept {
        return rho * Rspecific() * T;
    }
    double soundSpeed(double T, double gamma) const noexcept {
        return std::sqrt(gamma * Rspecific() * T);
    }
};

struct PengRobinsonEos {
    double molecularWeight = 0.044;    // kg/mol (CO2 example)
    double Tcrit           = 304.13;   // K
    double pcrit           = 7.377e6;  // Pa
    double acentric        = 0.225;

    double Rspecific() const noexcept { return kRgasUniversal / molecularWeight; }

    /// Newton iterate on molar volume v from (T, p). Returns ρ in kg/m³.
    double density(double T, double p) const {
        const double R = Rspecific() * molecularWeight; // J/(mol·K)
        const double a = 0.45724 * R * R * Tcrit * Tcrit / pcrit;
        const double b = 0.07780 * R * Tcrit / pcrit;
        const double kappa = 0.37464 + 1.54226 * acentric
                           - 0.26992 * acentric * acentric;
        const double sqrtTr = std::sqrt(T / Tcrit);
        const double alpha = (1.0 + kappa * (1.0 - sqrtTr));
        const double aalpha = a * alpha * alpha;

        // Initial guess: ideal-gas v = R T / p (m³/mol)
        double v = R * T / std::max(p, 1.0);
        for (int it = 0; it < 50; ++it) {
            const double denom1 = v - b;
            const double denom2 = v * v + 2.0 * b * v - b * b;
            const double f  = R * T / denom1 - aalpha / denom2 - p;
            const double df = -R * T / (denom1 * denom1)
                            + 2.0 * aalpha * (v + b) / (denom2 * denom2);
            const double dv = -f / df;
            v += dv;
            if (std::abs(dv) < 1e-14 * std::abs(v)) break;
            if (v < b * 1.01) v = b * 1.01;
        }
        const double rhoMolar = 1.0 / v;             // mol/m³
        return rhoMolar * molecularWeight;           // kg/m³
    }
    double pressure(double T, double rho) const {
        const double R = Rspecific() * molecularWeight;
        const double a = 0.45724 * R * R * Tcrit * Tcrit / pcrit;
        const double b = 0.07780 * R * Tcrit / pcrit;
        const double kappa = 0.37464 + 1.54226 * acentric
                           - 0.26992 * acentric * acentric;
        const double sqrtTr = std::sqrt(T / Tcrit);
        const double alpha = (1.0 + kappa * (1.0 - sqrtTr));
        const double aalpha = a * alpha * alpha;
        const double v = molecularWeight / std::max(rho, 1e-30);
        return R * T / (v - b) - aalpha / (v * v + 2.0 * b * v - b * b);
    }
};

struct NASA7Cp {
    /// 7-coefficient set for low (≤ Tmid) and high temperature ranges.
    std::array<double, 7> coeffsLow {{0,0,0,0,0,0,0}};
    std::array<double, 7> coeffsHigh{{0,0,0,0,0,0,0}};
    double Tlow  = 200.0;
    double Tmid  = 1000.0;
    double Thigh = 5000.0;
    double molecularWeight = 0.02897;

    /// c_p [J/(kg·K)] from polynomial.
    double cp(double T) const noexcept {
        const auto& a = (T < Tmid) ? coeffsLow : coeffsHigh;
        const double t = std::clamp(T, Tlow, Thigh);
        const double R = kRgasUniversal / molecularWeight;
        return R * (a[0] + a[1]*t + a[2]*t*t + a[3]*t*t*t + a[4]*t*t*t*t);
    }
    /// Sensible enthalpy h [J/kg] referenced so that h(0 K) cancels the
    /// a[5] constant of integration.
    double enthalpy(double T) const noexcept {
        const auto& a = (T < Tmid) ? coeffsLow : coeffsHigh;
        const double t = std::clamp(T, Tlow, Thigh);
        const double R = kRgasUniversal / molecularWeight;
        return R * t * (a[0] + 0.5*a[1]*t + (1.0/3.0)*a[2]*t*t
                              + 0.25*a[3]*t*t*t + 0.2*a[4]*t*t*t*t)
             + R * a[5];
    }
    /// Standard-state entropy s°(T) [J/(kg·K)].
    double entropy(double T) const noexcept {
        const auto& a = (T < Tmid) ? coeffsLow : coeffsHigh;
        const double t = std::clamp(T, Tlow, Thigh);
        const double R = kRgasUniversal / molecularWeight;
        return R * (a[0]*std::log(t) + a[1]*t + 0.5*a[2]*t*t
                   + (1.0/3.0)*a[3]*t*t*t + 0.25*a[4]*t*t*t*t + a[6]);
    }
};

}  // namespace simall::materials
