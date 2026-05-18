// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/SprayBreakup.cpp
// =============================================================================
#include "particles/SprayBreakup.hpp"

#include <algorithm>
#include <cmath>

namespace simall::particles {

namespace { constexpr double kPi = 3.14159265358979323846; }

void SprayBreakup::apply(double dt,
                         LagrangianTracker& tracker,
                         const solver::FieldRegistry& F)
{
    auto& parts = tracker.mutable_particles();
    if (parts.empty()) return;
    rt_age_.resize(parts.size(), 0.0);

    const auto* U = F.find_vector("U");
    if (!U) return;

    for (std::size_t i = 0; i < parts.size(); ++i) {
        auto& pt = parts[i];
        if (!pt.active || pt.cell == static_cast<meshing::CellId>(-1)) continue;

        const std::size_t c = static_cast<std::size_t>(pt.cell);
        if (c >= U->x.size()) continue;
        const util::Vec3d u{U->x[c], U->y[c], U->z[c]};
        const util::Vec3d rel = u - pt.v;
        const double ur = rel.norm();
        if (ur < 1e-9) { rt_age_[i] = 0.0; continue; }

        const double r = 0.5 * std::cbrt(6.0 * pt.mass / (kPi * p_.rho_l));
        if (r < 1e-9) continue;

        const double We_g = p_.rho_g * ur*ur * r / p_.sigma;
        const double We_l = p_.rho_l * ur*ur * r / p_.sigma;
        const double Re_l = p_.rho_l * ur    * r / std::max(1e-30, p_.mu_l);
        const double Z    = std::sqrt(We_l) / std::max(1e-30, Re_l);
        const double Tn   = Z * std::sqrt(We_g);

        // KH wave.
        const double Lkh_num = 9.02 * r * (1.0 + 0.45*std::sqrt(Z))
                                       * (1.0 + 0.4 *std::pow(Tn, 0.7));
        const double Lkh_den = std::pow(1.0 + 0.87*std::pow(We_g, 1.67), 0.6);
        const double Lkh = Lkh_num / Lkh_den;
        const double Okh_num = 0.34 + 0.385*std::pow(We_g, 1.5);
        const double Okh_den = (1.0 + Z)*(1.0 + 1.4*std::pow(Tn, 0.6));
        const double Okh = (Okh_num / Okh_den)
                         * std::sqrt(p_.sigma / std::max(1e-30, p_.rho_l*r*r*r));
        const double tau_kh = 3.726 * p_.B1 * r / std::max(1e-30, Lkh * Okh);

        // Deceleration magnitude (estimate from drag |a_d| ≈ |F_drag|/m_p).
        // We approximate a_d via instantaneous relative-velocity-squared scaling.
        // Sufficient to drive RT activation; tracker provides exact accel
        // through future hooks.
        const double Cd = (Re_l > 0 && Re_l < 1000.0)
                          ? (24.0/Re_l)*(1.0 + 0.15*std::pow(Re_l,0.687))
                          : 0.44;
        const double a_d = 0.75 * Cd * p_.rho_g * ur*ur / (p_.rho_l * (2.0*r));

        bool rt_triggered = false;
        if (a_d > 0.0) {
            const double dRho = p_.rho_l - p_.rho_g;
            const double O2   = (2.0/(3.0*std::sqrt(3.0)))
                              * std::pow(a_d*dRho/p_.sigma, 1.5)
                              * std::sqrt(p_.sigma/(p_.rho_l + p_.rho_g));
            const double O_rt = std::sqrt(std::max(0.0, O2));
            const double L_rt = (2.0*kPi/std::sqrt(3.0))
                              * std::sqrt(p_.sigma / std::max(1e-30, a_d*dRho));
            const double tau_rt = p_.Ctau_RT / std::max(1e-30, O_rt);

            rt_age_[i] += dt;
            if (2.0*r > p_.C_RT * L_rt && rt_age_[i] > tau_rt) {
                rt_triggered = true;
                // RT child radius = C_RT · L_rt / 2.
                const double r_new = 0.5 * p_.C_RT * L_rt;
                const double m_new = (4.0/3.0)*kPi*r_new*r_new*r_new * p_.rho_l;
                pt.mass = m_new;
                rt_age_[i] = 0.0;
            }
        }
        if (!rt_triggered && We_g > p_.WeCrit) {
            // KH stripping: r → B0·Lkh over τ_kh.
            const double r_eq = p_.B0 * Lkh;
            if (r_eq < r) {
                const double r_new = r + (r_eq - r) * (dt / std::max(dt, tau_kh));
                const double rc = std::max(r_new, 1e-9);
                pt.mass = (4.0/3.0)*kPi*rc*rc*rc * p_.rho_l;
            }
        }
    }
}

}  // namespace simall::particles
