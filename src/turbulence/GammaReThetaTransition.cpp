// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/GammaReThetaTransition.cpp
// =============================================================================
#include "turbulence/GammaReThetaTransition.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence
{

namespace
{
// Langtry-Menter constants.
constexpr double c_a1 = 2.0;
constexpr double c_a2 = 0.06;
constexpr double c_e1 = 1.0;
constexpr double c_e2 = 50.0;
constexpr double sigma_f = 1.0;
constexpr double sigma_thetat = 2.0;
} // namespace

GammaReThetaTransition::GammaReThetaTransition() = default;

void GammaReThetaTransition::initialize(meshing::Mesh& m,
                                        solver::FieldRegistry& F,
                                        double rho,
                                        double mu)
{
    mesh_ = &m;
    rho_ = rho;
    mu_ = mu;
    const std::size_t nC = m.cells().size();
    F.scalar("gamma_int", nC);
    F.scalar("ReTheta_t", nC);
    F.scalar("gamma_eff", nC);
    F.scalar("wallDistance", nC);

    // Initialise to typical free-stream (γ = 1, Re̅θ_t = 100).
    auto& g = *F.find_scalar("gamma_int");
    auto& r = *F.find_scalar("ReTheta_t");
    auto& ge = *F.find_scalar("gamma_eff");
    for (std::size_t c = 0; c < nC; ++c) {
        g[c] = 1.0;
        r[c] = 100.0;
        ge[c] = 1.0;
    }

    S_.assign(nC, 0.0);
    W_.assign(nC, 0.0);
    Pgamma_.assign(nC, 0.0);
    ReThetaT_.assign(nC, 100.0);

    solver::LinearSolverConfig cfg;
    cfg.maxIterations = 200;
    cfg.tolerance = 1e-8;
    lin_ = solver::make_bicgstab(cfg);

    wd_ = solver::make_wall_distance_exact();
    wd_->compute(m, bcs_, *F.find_scalar("wallDistance"));

    gammaEq_ = std::make_unique<solver::ScalarTransport>(m, F, *lin_);
    ReThetaEq_ = std::make_unique<solver::ScalarTransport>(m, F, *lin_);
    gammaEq_->set_field("gamma_int");
    gammaEq_->set_diffusivity(mu_ / sigma_f); // diffusivity refined below per-cell N/A
    gammaEq_->set_density(rho_);
    gammaEq_->set_source("gamma_src");
    ReThetaEq_->set_field("ReTheta_t");
    ReThetaEq_->set_diffusivity(sigma_thetat * mu_);
    ReThetaEq_->set_density(rho_);
    ReThetaEq_->set_source("ReTheta_src");
    for (const auto& bc : bcs_) {
        solver::ScalarBC sb;
        sb.zone = bc.zone;
        sb.kind = solver::ScalarBC::Kind::Neumann;
        sb.value = 0.0;
        gammaEq_->add_bc(sb);
        ReThetaEq_->add_bc(sb);
    }
    F.scalar("gamma_src", nC);
    F.scalar("ReTheta_src", nC);
}

double GammaReThetaTransition::correlation_ReThetaT(double Tu, double /*lambda*/)
{
    // Langtry-Menter (2009), zero pressure gradient form.
    Tu = std::max(Tu, 0.027);
    if (Tu <= 1.3)
        return 1173.51 - 589.428 * Tu + 0.2196 / (Tu * Tu);
    else
        return 331.5 * std::pow(Tu - 0.5658, -0.671);
}

double GammaReThetaTransition::F_onset(
    double rho, double mu, double mut, double S, double d, double ReThetaT)
{
    if (d <= 0 || mu <= 0)
        return 0.0;
    const double ReV = rho * d * d * S / mu;
    const double RT = (mu > 0) ? mut / mu : 0.0;
    const double ReThetaC = ReThetaT - 396.035e-2; // calibration offset
    const double Fonset1 = ReV / (2.193 * std::max(1.0, ReThetaC));
    const double Fonset2 = std::min(std::max(Fonset1, std::pow(Fonset1, 4.0)), 2.0);
    const double Fonset3 = std::max(1.0 - std::pow(RT / 2.5, 3.0), 0.0);
    return std::max(Fonset2 - Fonset3, 0.0);
}

void GammaReThetaTransition::compute_strain_and_vorticity(const solver::FieldRegistry& F)
{
    const auto* U = F.find_vector("U");
    if (!U || !mesh_)
        return;
    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gU, gV, gW;
    G.evaluate(U->x, gU);
    G.evaluate(U->y, gV);
    G.evaluate(U->z, gW);
    const std::size_t nC = mesh_->cells().size();
    for (std::size_t c = 0; c < nC; ++c) {
        const double Sxx = gU.x[c], Syy = gV.y[c], Szz = gW.z[c];
        const double Sxy = 0.5 * (gU.y[c] + gV.x[c]);
        const double Sxz = 0.5 * (gU.z[c] + gW.x[c]);
        const double Syz = 0.5 * (gV.z[c] + gW.y[c]);
        S_[c] = std::sqrt(
            2.0 * (Sxx * Sxx + Syy * Syy + Szz * Szz + 2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz)));
        const double Wxy = 0.5 * (gU.y[c] - gV.x[c]);
        const double Wxz = 0.5 * (gU.z[c] - gW.x[c]);
        const double Wyz = 0.5 * (gV.z[c] - gW.y[c]);
        W_[c] = std::sqrt(2.0 * 2.0 * (Wxy * Wxy + Wxz * Wxz + Wyz * Wyz));
    }
}

void GammaReThetaTransition::compute_sources(const solver::FieldRegistry& F)
{
    const std::size_t nC = mesh_->cells().size();
    const auto& g = *F.find_scalar("gamma_int");
    const auto& rt = *F.find_scalar("ReTheta_t");
    const auto* k = F.find_scalar("k");
    const auto* mt = F.find_scalar("mu_t");
    const auto* d = F.find_scalar("wallDistance");
    const auto* U = F.find_vector("U");
    auto& gSrc = *F.find_scalar("gamma_src");
    auto& rSrc = *F.find_scalar("ReTheta_src");

    for (std::size_t c = 0; c < nC; ++c) {
        const double mut = (mt ? (*mt)[c] : 0.0);
        const double dwall = (d ? (*d)[c] : 0.0);
        const double S = S_[c];
        const double Umag =
            U ? std::sqrt(U->x[c] * U->x[c] + U->y[c] * U->y[c] + U->z[c] * U->z[c]) : 0.0;
        const double Tu = k ? std::min(100.0,
                                       100.0 * std::sqrt(2.0 * std::max(0.0, (*k)[c]) / 3.0)
                                           / std::max(1e-6, Umag))
                            : 1.0;
        const double ReThetaT_local = correlation_ReThetaT(Tu, 0.0);
        ReThetaT_[c] = ReThetaT_local;

        // γ source (Langtry-Menter): P_γ = F_length c_a1 ρ S (γ F_onset)^0.5 (1 - c_e1 γ)
        const double Flength = 100.0; // conservative default for high-Re
        const double Fonset = F_onset(rho_, mu_, mut, S, dwall, rt[c]);
        const double Pg = Flength * c_a1 * rho_ * S * std::sqrt(std::max(0.0, g[c] * Fonset))
                          * (1.0 - c_e1 * g[c]);
        // E_γ (destruction) ~ c_a2 ρ Ω γ F_turb (c_e2 γ - 1)
        const double Eg = c_a2 * rho_ * W_[c] * g[c] * (c_e2 * g[c] - 1.0);
        gSrc[c] = Pg - Eg;
        Pgamma_[c] = Pg;

        // Reθ_t source: P_θt = c_θt (ρ/t) (Reθ_T - Reθ_t)(1 - F_θt),
        // with time scale t = 500 μ / (ρ U²)
        const double timeScale = (Umag > 1e-6) ? 500.0 * mu_ / (rho_ * Umag * Umag) : 0.0;
        const double Pthetat =
            (timeScale > 1e-30) ? 0.03 * (rho_ / timeScale) * (ReThetaT_local - rt[c]) : 0.0;
        rSrc[c] = Pthetat;
    }
}

void GammaReThetaTransition::update_gamma_eff(const solver::FieldRegistry& F)
{
    const std::size_t nC = mesh_->cells().size();
    const auto& g = *F.find_scalar("gamma_int");
    auto& ge = *const_cast<solver::FieldRegistry&>(F).find_scalar("gamma_eff");
    for (std::size_t c = 0; c < nC; ++c) {
        // Effective intermittency uses separation-induced correction γ_sep
        // (LM2009 eq. 15). Simplified: γ_eff = max(γ, γ_sep) with γ_sep≤2.
        const double gsep = std::min(2.0, std::max(0.0, S_[c] / (3.235e0))); // calibration scale
        ge[c] = std::max(g[c], 0.5 * gsep);
    }
}

double GammaReThetaTransition::solve_iteration(double dt, solver::FieldRegistry& F)
{
    if (!mesh_)
        return 0.0;
    compute_strain_and_vorticity(F);
    compute_sources(F);
    const double rG = gammaEq_->solve_iteration();
    const double rR = ReThetaEq_->solve_iteration();
    update_gamma_eff(F);
    // Clamp γ to [0, 1] (physical).
    auto& g = *F.find_scalar("gamma_int");
    for (auto& gi : g)
        gi = std::clamp(gi, 0.0, 1.0);
    (void) dt;
    return std::max(rG, rR);
}

} // namespace simall::turbulence
