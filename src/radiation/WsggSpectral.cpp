// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/WsggSpectral.cpp
//
// Coefficient sources:
//   * Smith, T.F., Shen, Z.F., Friedman, J.N. "Evaluation of Coefficients
//     for the Weighted Sum of Gray Gases Model", J. Heat Transfer 104
//     (1982) 602-608.    p_w / p_c = 2 mixture, 3 grey + 1 window.
//   * Bordbar, M.H., Wecel, G., Hyppänen, T. "A line by line based weighted
//     sum of gray gases model for inhomogeneous CO2-H2O mixture in oxy-fired
//     combustion", Combust. Flame 161 (2014) 2435-2445.
// =============================================================================
#include "radiation/WsggSpectral.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::radiation {

bool WsggSpectral::initialize(const meshing::Mesh& mesh,
                              solver::FieldRegistry& F,
                              const WsggProps& props) {
    mesh_ = &mesh; F_ = &F; p_ = props;
    if (p_.useBordbar) load_bordbar_coefficients();
    else               load_smith_coefficients();
    F.scalar("T",         mesh.cells().size());
    F.scalar("kappa_eff", mesh.cells().size());
    SIMALL_LOG_INFO("Radiation",
        "WSGG init: p_w=", p_.p_H2O, " p_c=", p_.p_CO2,
        " L_m=", p_.meanBeamLength,
        " coeffs=", (p_.useBordbar ? "Bordbar2014" : "Smith1982"));
    return true;
}

void WsggSpectral::load_smith_coefficients() {
    // Smith 1982, Table 2  (p_w / p_c = 2, p_w + p_c = 1 atm).
    // k=0 transparent window; κ_0 = 0.
    kappa_ = { 0.0, 0.4303, 7.055, 178.1 };       // [1/(atm·m)]
    // b_{k, 0..3}; T scaled in K.  Note Smith reports in 1/1000 K^j powers;
    // the values below already include those scalings.
    b_ = {{
        // k=0 transparent (a_0 = 1 - Σ a_k>0)
        { 1.0, 0.0, 0.0, 0.0 },
        { 5.150e-1, -2.303e-4,  9.779e-8, -1.494e-11 },
        { 7.749e-2,  3.399e-4, -2.297e-7,  3.770e-11 },
        { 1.907e-1, -1.824e-4,  5.608e-8, -5.122e-12 }
    }};
}

void WsggSpectral::load_bordbar_coefficients() {
    // Bordbar 2014, Table 2 (p_w/p_c = 2, oxy-fired).
    kappa_ = { 0.0, 0.192, 1.719, 11.370 };
    b_ = {{
        { 1.0, 0.0, 0.0, 0.0 },
        { 5.617e-1, -7.840e-4,  6.939e-7, -1.892e-10 },
        { 1.390e-1,  1.095e-3, -8.305e-7,  1.907e-10 },
        { 4.018e-2,  4.260e-4, -3.289e-7,  7.685e-11 }
    }};
}

double WsggSpectral::compute_weight(double T, int k) const {
    if (k < 0 || k >= kWsggBands) return 0.0;
    if (k == 0) {
        double s = 0.0;
        for (int j = 1; j < kWsggBands; ++j) s += compute_weight(T, j);
        return std::clamp(1.0 - s, 0.0, 1.0);
    }
    const auto& c = b_[k];
    return std::clamp(c[0] + c[1]*T + c[2]*T*T + c[3]*T*T*T, 0.0, 1.0);
}

double WsggSpectral::compute_kappa(int k) const {
    return (k >= 0 && k < kWsggBands) ? kappa_[k] : 0.0;
}

double WsggSpectral::compute_emissivity(double T, double pwL) const {
    double eps = 0.0;
    for (int k = 1; k < kWsggBands; ++k) {
        const double ak = compute_weight(T, k);
        eps += ak * (1.0 - std::exp(-kappa_[k] * pwL));
    }
    return std::clamp(eps, 0.0, 1.0);
}

double WsggSpectral::apply() {
    if (!mesh_ || !F_) return 0.0;
    const auto* T = F_->find_scalar("T");
    auto*       K = F_->find_scalar("kappa_eff");
    if (!T || !K) return 0.0;
    const double pwL = (p_.p_H2O + p_.p_CO2) * p_.meanBeamLength;
    double kMax = 0.0;
    for (std::size_t c = 0; c < T->size(); ++c) {
        const double Tc  = std::max((*T)[c], 200.0);
        const double eps = compute_emissivity(Tc, pwL);
        // Equivalent grey absorption (Bouguer inversion):
        //   κ_eff = -ln(1 - ε) / L_m
        const double k_eff = -std::log(std::max(1.0 - eps, 1e-12))
                           / std::max(p_.meanBeamLength, 1e-9);
        (*K)[c] = k_eff;
        kMax = std::max(kMax, k_eff);
    }
    return kMax;
}

}  // namespace simall::radiation
