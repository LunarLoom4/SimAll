// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/TransportProperties.hpp
// Phase  : 11 — Pure-species and multi-component transport properties.
//
//   Sutherland viscosity     : μ(T) = μ_ref (T/T_ref)^{3/2} (T_ref+S)/(T+S)
//   Sutherland conductivity  : k(T) = k_ref (T/T_ref)^{3/2} (T_ref+S_k)/(T+S_k)
//   Eucken mass diffusivity  : D_im = (k / (ρ cp)) / Le      (constant Lewis)
//
//   Wilke mixture rule (J. Chem. Phys. 1950, 18, 517):
//   μ_mix = Σ_i (x_i μ_i) / Σ_j (x_j φ_ij)
//   with   φ_ij = (1 + (μ_i/μ_j)^{1/2} (M_j/M_i)^{1/4})²
//                 / sqrt(8 (1 + M_i/M_j))
//
//   Mathur-Saxena conductivity mixture (analogous form).
//
// Header-only; no .cpp needed. All routines noexcept.
// =============================================================================
#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace simall::materials {

struct SutherlandViscosity {
    double mu_ref = 1.716e-5;   // Pa·s (air @ 273.15K)
    double T_ref  = 273.15;     // K
    double S      = 110.4;      // K

    double mu(double T) const noexcept {
        const double r = T / T_ref;
        return mu_ref * r * std::sqrt(r) * (T_ref + S) / (T + S);
    }
};

struct SutherlandConductivity {
    double k_ref = 0.0241;      // W/(m·K) air @ 273.15K
    double T_ref = 273.15;
    double S_k   = 194.0;

    double k(double T) const noexcept {
        const double r = T / T_ref;
        return k_ref * r * std::sqrt(r) * (T_ref + S_k) / (T + S_k);
    }
};

struct SpeciesTransport {
    double molecularWeight = 0.02897;   // kg/mol
    SutherlandViscosity    visc;
    SutherlandConductivity cond;
};

/// Wilke mixture viscosity. xMole are mole fractions (Σ = 1), one per species.
inline double wilke_mixture_viscosity(
        double T,
        const std::vector<SpeciesTransport>& sp,
        const std::vector<double>& xMole) noexcept
{
    const std::size_t n = sp.size();
    if (n == 0 || xMole.size() != n) return 0.0;
    std::vector<double> mu(n);
    for (std::size_t i = 0; i < n; ++i) mu[i] = sp[i].visc.mu(T);
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double denom = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            const double mr = mu[i] / mu[j];
            const double Mr = sp[j].molecularWeight / sp[i].molecularWeight;
            const double num = 1.0 + std::sqrt(mr) * std::pow(Mr, 0.25);
            const double phi = (num*num) /
                std::sqrt(8.0 * (1.0 + sp[i].molecularWeight/sp[j].molecularWeight));
            denom += xMole[j] * phi;
        }
        if (denom > 1e-30) total += xMole[i] * mu[i] / denom;
    }
    return total;
}

/// Mass-fraction to mole-fraction conversion.
inline std::vector<double> mass_to_mole(
        const std::vector<double>& Y,
        const std::vector<SpeciesTransport>& sp) noexcept
{
    const std::size_t n = Y.size();
    std::vector<double> x(n, 0.0);
    if (n == 0 || sp.size() != n) return x;
    double mWmean_inv = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        mWmean_inv += Y[i] / sp[i].molecularWeight;
    const double mWmean = (mWmean_inv > 0.0) ? 1.0 / mWmean_inv : 0.0;
    for (std::size_t i = 0; i < n; ++i)
        x[i] = (sp[i].molecularWeight > 0)
             ? Y[i] * mWmean / sp[i].molecularWeight : 0.0;
    return x;
}

/// Mathur-Saxena mixture conductivity (same combination rule structure as Wilke).
inline double mathur_saxena_conductivity(
        double T,
        const std::vector<SpeciesTransport>& sp,
        const std::vector<double>& xMole) noexcept
{
    const std::size_t n = sp.size();
    if (n == 0 || xMole.size() != n) return 0.0;
    double sum1 = 0.0, sum2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double ki = sp[i].cond.k(T);
        sum1 += xMole[i] * ki;
        sum2 += (ki > 0) ? xMole[i] / ki : 0.0;
    }
    return 0.5 * (sum1 + (sum2 > 1e-30 ? 1.0/sum2 : 0.0));
}

}  // namespace simall::materials
