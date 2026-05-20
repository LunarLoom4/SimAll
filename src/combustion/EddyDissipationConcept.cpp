// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/EddyDissipationConcept.cpp
// =============================================================================
#include "combustion/EddyDissipationConcept.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::combustion
{

namespace
{
constexpr double Runi = 8.314462618; // J/(mol·K)
}

void EddyDissipationConcept::initialize(const meshing::Mesh& m,
                                        solver::FieldRegistry& F,
                                        std::vector<Species> s,
                                        std::vector<Reaction> r,
                                        EdcProps p)
{
    mesh_ = &m;
    p_ = p;
    species_ = std::move(s);
    reactions_ = std::move(r);
    const std::size_t nC = m.cells().size();
    for (const auto& sp : species_) {
        F.scalar("Y_" + sp.name, nC);
        F.scalar("wdot_" + sp.name, nC);
    }
    F.scalar("S_combustion", nC);
    F.scalar("edc_gamma", nC);
    F.scalar("edc_tau", nC);
}

void EddyDissipationConcept::laminar_rates(const std::vector<double>& Y,
                                           double rho,
                                           double T,
                                           std::vector<double>& wdot) const
{
    const std::size_t N = species_.size();
    wdot.assign(N, 0.0);
    // Mole concentrations [X_j] = ρ Y_j / W_j.
    std::vector<double> conc(N, 0.0);
    for (std::size_t j = 0; j < N; ++j) {
        const double W = std::max(species_[j].molarMass, 1e-12);
        conc[j] = rho * Y[j] / W;
    }
    for (const auto& rxn : reactions_) {
        const double kf = rxn.A * std::pow(std::max(T, 1.0), rxn.beta)
                          * std::exp(-rxn.Ea / (Runi * std::max(T, 1.0)));
        double q = kf;
        const std::size_t Nr = std::min(rxn.order.size(), N);
        for (std::size_t j = 0; j < Nr; ++j) {
            if (rxn.order[j] != 0.0)
                q *= std::pow(std::max(conc[j], 0.0), rxn.order[j]);
        }
        // ω̇_k mass rate [kg/(m³·s)] = (ν'_k - ν''_k) W_k · q
        const std::size_t Np = std::min(rxn.nuP.size(), N);
        const std::size_t Nrr = std::min(rxn.nuR.size(), N);
        for (std::size_t k = 0; k < N; ++k) {
            const double nuP = (k < Np) ? rxn.nuP[k] : 0.0;
            const double nuR = (k < Nrr) ? rxn.nuR[k] : 0.0;
            wdot[k] += (nuP - nuR) * species_[k].molarMass * q;
        }
    }
}

void EddyDissipationConcept::psr_step(const std::vector<double>& Yin,
                                      double rho,
                                      double T,
                                      double tau,
                                      std::vector<double>& Yout) const
{
    // Implicit Euler with one Newton iteration (sufficient for stiff PSR
    // when τ is one local mixing time and the linearised Jacobian is
    // approximated by ∂ω̇/∂Y ≈ ω̇/Y, the local-equilibrium relaxation).
    const std::size_t N = species_.size();
    Yout = Yin;
    if (tau <= 0.0 || N == 0)
        return;
    std::vector<double> w(N, 0.0);
    laminar_rates(Yout, rho, T, w);
    for (std::size_t k = 0; k < N; ++k) {
        // Implicit Euler: Y_new = Y_in + τ w/ρ, but stabilise w against Y_new.
        const double jac = (Yout[k] > 1e-12) ? std::abs(w[k]) / (rho * Yout[k]) : 1.0 / tau;
        const double denom = 1.0 + tau * jac;
        Yout[k] += tau * w[k] / (rho * denom);
        Yout[k] = std::clamp(Yout[k], 0.0, 1.0);
    }
    // Renormalise to preserve Σ Y = 1.
    double sum = 0.0;
    for (double y : Yout)
        sum += y;
    if (sum > 1e-12)
        for (auto& y : Yout)
            y /= sum;
}

void EddyDissipationConcept::apply(solver::FieldRegistry& F, double rho_default)
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    const std::size_t N = species_.size();
    if (N == 0)
        return;

    const auto* rhoF = F.find_scalar("rho");
    const auto* TF = F.find_scalar("T");
    const auto* kF = F.find_scalar("k");
    const auto* eF = F.find_scalar("epsilon");
    auto& Scomb = *F.find_scalar("S_combustion");
    auto& gammaOut = *F.find_scalar("edc_gamma");
    auto& tauOut = *F.find_scalar("edc_tau");

    std::vector<solver::ScalarField*> Yfield(N), Wfield(N);
    for (std::size_t k = 0; k < N; ++k) {
        Yfield[k] = F.find_scalar("Y_" + species_[k].name);
        Wfield[k] = F.find_scalar("wdot_" + species_[k].name);
    }

    std::vector<double> Yc(N), wlam(N), Yfs(N), Yc_orig(N);
    for (std::size_t c = 0; c < nC; ++c) {
        const double rho = rhoF ? (*rhoF)[c] : rho_default;
        const double T = TF ? (*TF)[c] : 300.0;
        const double kT = kF ? std::max((*kF)[c], 1e-12) : 0.1;
        const double eT = eF ? std::max((*eF)[c], 1e-12) : 1.0;

        for (std::size_t k = 0; k < N; ++k)
            Yc[k] = (*Yfield[k])[c];
        Yc_orig = Yc;

        if (p_.mode == TciMode::EDM) {
            // Mean mixing rate τ_mix^-1 ~ ε / k.
            const double mixRate = eT / kT;
            // Naïve EDM: pick the limiting reactant for each reaction.
            std::fill(wlam.begin(), wlam.end(), 0.0);
            for (const auto& rxn : reactions_) {
                double Ylim = 1e30;
                for (std::size_t j = 0; j < N; ++j) {
                    if (j >= rxn.nuR.size() || rxn.nuR[j] == 0.0)
                        continue;
                    const double s = rxn.nuR[j] * species_[j].molarMass;
                    if (s <= 0.0)
                        continue;
                    Ylim = std::min(Ylim, Yc[j] / s);
                }
                if (Ylim == 1e30)
                    continue;
                const double q = p_.A_EDM * rho * mixRate * Ylim;
                const std::size_t Np = std::min(rxn.nuP.size(), N);
                const std::size_t Nrr = std::min(rxn.nuR.size(), N);
                for (std::size_t k = 0; k < N; ++k) {
                    const double nuP = (k < Np) ? rxn.nuP[k] : 0.0;
                    const double nuR = (k < Nrr) ? rxn.nuR[k] : 0.0;
                    wlam[k] += (nuP - nuR) * species_[k].molarMass * q;
                }
            }
            gammaOut[c] = 0.0;
            tauOut[c] = 1.0 / mixRate;
        } else {
            // EDC.
            const double nu = p_.nu;
            const double gam = p_.C_gamma * std::pow(nu * eT / (kT * kT), 0.25);
            const double tau = p_.C_tau * std::sqrt(nu / eT);
            const double g3 = std::min(0.999, gam * gam * gam);
            const double prefactor = rho * gam * gam / (tau * (1.0 - g3));
            // Fine-structure PSR step.
            psr_step(Yc, rho, T, tau, Yfs);
            for (std::size_t k = 0; k < N; ++k)
                wlam[k] = prefactor * (Yfs[k] - Yc[k]);
            gammaOut[c] = gam;
            tauOut[c] = tau;
        }

        double Scomb_local = 0.0;
        for (std::size_t k = 0; k < N; ++k) {
            (*Wfield[k])[c] = wlam[k];
            Scomb_local -= species_[k].formationEnthalpy * wlam[k];
        }
        Scomb[c] = Scomb_local;
    }
}

} // namespace simall::combustion
