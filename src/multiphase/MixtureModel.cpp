// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/MixtureModel.cpp
// =============================================================================
#include "multiphase/MixtureModel.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase {

void MixtureModel::initialize(const meshing::Mesh& m,
                              std::vector<MixturePhase> phases,
                              double gx, double gy, double gz) {
    mesh_   = &m;
    phases_ = std::move(phases);
    g_[0]   = gx; g_[1] = gy; g_[2] = gz;
    SIMALL_LOG_INFO("Multiphase",
        "MixtureModel initialised with ", phases_.size(),
        " phases (", m.cells().size(), " cells)");
}

double MixtureModel::apply(solver::FieldRegistry& F) {
    if (!mesh_ || phases_.empty()) return 0.0;
    const std::size_t nC = mesh_->cells().size();
    const auto* Umix = F.find_vector("U_mix");
    if (!Umix) return 0.0;

    auto& rhoM = F.scalar("rho_mix", nC);
    auto& muM  = F.scalar("mu_mix",  nC);
    std::fill(rhoM.begin(), rhoM.end(), 0.0);
    std::fill(muM .begin(), muM .end(), 0.0);

    // Gather α_k for each phase (phase 0 = primary, α_0 = 1-Σα_k).
    std::vector<util::aligned_vector<double>*> alpha(phases_.size(), nullptr);
    for (std::size_t k = 1; k < phases_.size(); ++k) {
        const std::string n = "alpha_" + phases_[k].name;
        alpha[k] = F.find_scalar(n);
    }

    // Mixture density / viscosity.
    for (std::size_t c = 0; c < nC; ++c) {
        double rho = 0.0, mu = 0.0, asum = 0.0;
        for (std::size_t k = 1; k < phases_.size(); ++k) {
            const double a = alpha[k] ? std::clamp((*alpha[k])[c], 0.0, 1.0) : 0.0;
            rho  += a * phases_[k].rho;
            mu   += a * phases_[k].mu;
            asum += a;
        }
        const double a0 = std::clamp(1.0 - asum, 0.0, 1.0);
        rho += a0 * phases_[0].rho;
        mu  += a0 * phases_[0].mu;
        rhoM[c] = rho;
        muM [c] = mu;
    }

    // Algebraic-slip relative velocities U_kr per dispersed phase.
    // For each secondary phase k, U_kr = τ_p (ρ_k - ρ_m)/ρ_k · (g − ∇p/ρ_m approx)
    // Here we use only gravity + mixture momentum advection as the driving force.
    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(Umix->x, gUx);
    G.evaluate(Umix->y, gUy);
    G.evaluate(Umix->z, gUz);

    double maxSlip = 0.0;
    for (std::size_t k = 1; k < phases_.size(); ++k) {
        auto& Ukr = F.vector("U_kr_" + phases_[k].name, nC);
        const double rhoP = phases_[k].rho;
        const double dP   = std::max(phases_[k].diameter, 1e-12);
        for (std::size_t c = 0; c < nC; ++c) {
            // Driving force per unit mass.
            const double ax = g_[0]
                - (Umix->x[c]*gUx.x[c] + Umix->y[c]*gUx.y[c] + Umix->z[c]*gUx.z[c]);
            const double ay = g_[1]
                - (Umix->x[c]*gUy.x[c] + Umix->y[c]*gUy.y[c] + Umix->z[c]*gUy.z[c]);
            const double az = g_[2]
                - (Umix->x[c]*gUz.x[c] + Umix->y[c]*gUz.y[c] + Umix->z[c]*gUz.z[c]);

            // Stokes relaxation time with Schiller-Naumann correction.
            // Initial guess (Stokes): U_kr0 = (ρ_k-ρ_m)d²/(18μ_m) |a|
            const double mu_m = std::max(muM[c], 1e-12);
            const double tauStokes = rhoP * dP*dP / (18.0 * mu_m);
            const double drho = (rhoP - rhoM[c]) / rhoP;

            // Two-step iteration over f_d to converge on Re_p.
            double Uxr = tauStokes * drho * ax;
            double Uyr = tauStokes * drho * ay;
            double Uzr = tauStokes * drho * az;
            for (int it = 0; it < 3; ++it) {
                const double Umag = std::sqrt(Uxr*Uxr + Uyr*Uyr + Uzr*Uzr);
                const double Rep  = rhoM[c] * Umag * dP / mu_m;
                const double fd   = 1.0 + 0.15 * std::pow(std::max(Rep, 1e-12), 0.687);
                const double tau  = tauStokes / fd;
                Uxr = tau * drho * ax;
                Uyr = tau * drho * ay;
                Uzr = tau * drho * az;
            }
            Ukr.x[c] = Uxr; Ukr.y[c] = Uyr; Ukr.z[c] = Uzr;
            const double Umag = std::sqrt(Uxr*Uxr + Uyr*Uyr + Uzr*Uzr);
            if (Umag > maxSlip) maxSlip = Umag;
        }
    }
    return maxSlip;
}

}  // namespace simall::multiphase
