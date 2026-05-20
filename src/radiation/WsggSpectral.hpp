// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/WsggSpectral.hpp
// Phase  : 11.8 — Weighted Sum of Gray Gases (WSGG) non-grey model.
//
// Approximates the non-grey emissivity of an H₂O–CO₂ gas mixture as the
// weighted sum of N grey gases plus one transparent window:
//
//   ε(T, p_w·L) = Σ_{k=0..N} a_k(T) · (1 - exp(-κ_k · p_w · L))
//
// with a_k(T) polynomial coefficients (Smith-Shen-Friedman 1982 or
// Bordbar-Wecel-Hyppänen 2014 for varied molar ratios).  Coefficients for
// the canonical Smith 4-grey-gas model (κ_0 = 0 = transparent window) are
// hardcoded below for p_w·L pressure-path-length range 0.001 – 10 atm·m.
//
// API:
//   compute_kappa(T, p_w, L, k) → κ_k for grey-band k (1/m)
//   compute_weight(T, k)        → a_k(T) dimensionless
//   compute_emissivity(T, pwL)  → ε for diagnostic use
//   apply(FieldRegistry)        → writes per-cell "kappa_eff" using
//                                 the cell-mean Planck-mean expression
//                                 κ_eff = -ln(1-ε)/L_path with L_path
//                                 supplied at construction (mean beam
//                                 length, e.g. 3.6 V_total/A_total).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <array>
#include <vector>

namespace simall::radiation
{

constexpr int kWsggBands = 4; // 1 transparent + 3 grey

struct WsggProps
{
    double p_H2O = 0.2; // partial pressure [atm]
    double p_CO2 = 0.1;
    double meanBeamLength = 0.5; // L_m [m]
    bool useBordbar = false;     // true ⇒ Bordbar 2014 (else Smith 1982)
};

class WsggSpectral
{
public:
    bool initialize(const meshing::Mesh& mesh,
                    solver::FieldRegistry& fields,
                    const WsggProps& props);

    /// Per-cell evaluation of κ_eff(T) written into "kappa_eff".
    double apply();

    double compute_emissivity(double T, double pwL) const;
    double compute_weight(double T, int k) const;
    double compute_kappa(int k) const; // [1/(atm·m)]; multiply by (p_w+p_c)

    const WsggProps& props() const noexcept { return p_; }

private:
    void load_smith_coefficients();
    void load_bordbar_coefficients();

    const meshing::Mesh* mesh_ = nullptr;
    solver::FieldRegistry* F_ = nullptr;
    WsggProps p_{};
    // a_k(T) = b_{k,0} + b_{k,1} T + b_{k,2} T² + b_{k,3} T³
    std::array<std::array<double, 4>, kWsggBands> b_{};
    std::array<double, kWsggBands> kappa_{};
};

} // namespace simall::radiation
