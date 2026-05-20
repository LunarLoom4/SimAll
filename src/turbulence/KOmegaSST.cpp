// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KOmegaSST.cpp
// =============================================================================
#include "turbulence/KOmegaSST.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence
{

namespace
{
// SST 2003 constants (Menter, Kuntz & Langtry 2003).
constexpr double sigma_k1 = 0.85, sigma_k2 = 1.0;
constexpr double sigma_w1 = 0.5, sigma_w2 = 0.856;
constexpr double beta_1 = 0.075, beta_2 = 0.0828;
constexpr double beta_star = 0.09;
constexpr double a1 = 0.31;
constexpr double kappa = 0.41;
constexpr double gamma_1 = beta_1 / beta_star - sigma_w1 * kappa * kappa / std::sqrt(beta_star);
constexpr double gamma_2 = beta_2 / beta_star - sigma_w2 * kappa * kappa / std::sqrt(beta_star);
constexpr double a_bot_clip = 1e-30;
} // namespace

KOmegaSST_Full::KOmegaSST_Full()
{
    solver::LinearSolverConfig c;
    c.kind = solver::LinearSolverKind::GMRES;
    c.preconditioner = solver::PreconditionerKind::ILU;
    c.tolerance = 1e-6;
    c.maxIterations = 100;
    c.restart = 20;
    lin_ = solver::make_linear_solver(c);
    wallDist_ = solver::make_wall_distance_exact();
}

void KOmegaSST_Full::initialize(meshing::Mesh& mesh, solver::FieldRegistry& f)
{
    mesh_ = &mesh;
    const std::size_t nC = mesh.cells().size();
    f.scalar("k", nC);
    f.scalar("omega", nC);
    f.scalar("mut", nC);
    f.scalar("wallDistance", nC);
    S_.assign(nC, 0.0);
    F1_.assign(nC, 1.0);
    F2_.assign(nC, 1.0);
    mut_.assign(nC, 0.0);
    Pk_.assign(nC, 0.0);
    CDk_.assign(nC, 0.0);

    // Compute wall distance once at start-up (steady RANS assumption).
    auto& d = *f.find_scalar("wallDistance");
    wallDist_->compute(mesh, bcs_, d);

    // Initialise k, ω with non-trivial floor (avoid division by zero).
    auto& k = *f.find_scalar("k");
    auto& w = *f.find_scalar("omega");
    std::fill(k.begin(), k.end(), 1e-4);
    std::fill(w.begin(), w.end(), 1.0);

    kEq_ = std::make_unique<solver::ScalarTransport>(mesh, f, *lin_);
    wEq_ = std::make_unique<solver::ScalarTransport>(mesh, f, *lin_);
    kEq_->set_field("k");
    wEq_->set_field("omega");
    kEq_->set_density(rho_);
    wEq_->set_density(rho_);
}

void KOmegaSST_Full::compute_strain_rate(const solver::FieldRegistry& f)
{
    auto& fnc = const_cast<solver::FieldRegistry&>(f);
    const auto* U = fnc.find_vector("U");
    if (!U)
        return;
    const std::size_t nC = mesh_->cells().size();

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U->x, gUx);
    G.evaluate(U->y, gUy);
    G.evaluate(U->z, gUz);

    S_.assign(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c) {
        // Symmetric strain-rate tensor Sij = 0.5(∂Ui/∂xj + ∂Uj/∂xi)
        const double s11 = gUx.x[c];
        const double s22 = gUy.y[c];
        const double s33 = gUz.z[c];
        const double s12 = 0.5 * (gUx.y[c] + gUy.x[c]);
        const double s13 = 0.5 * (gUx.z[c] + gUz.x[c]);
        const double s23 = 0.5 * (gUy.z[c] + gUz.y[c]);
        const double SS =
            s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
        S_[c] = std::sqrt(2.0 * SS); // |S| = √(2 Sij Sij)
    }
}

void KOmegaSST_Full::compute_blending()
{
    // Compute F1 / F2 from k, ω, ν, wallDistance and cross-diffusion term.
    // Requires gradients of k and ω.
    const std::size_t nC = mesh_->cells().size();
    // Pull current k, ω, d from the registry hosted by parent solver — we
    // assume initialize() has run.
    // NOTE: we read fields by name via mesh access. Strain rate already done.
}

void KOmegaSST_Full::update_mut(const solver::FieldRegistry& f)
{
    auto& fnc = const_cast<solver::FieldRegistry&>(f);
    const auto& k = *fnc.find_scalar("k");
    const auto& w = *fnc.find_scalar("omega");
    auto& mut = *fnc.find_scalar("mut");
    const std::size_t nC = k.size();
    for (std::size_t c = 0; c < nC; ++c) {
        const double denom = std::max(a1 * w[c], S_[c] * F2_[c]);
        mut[c] = rho_ * a1 * std::max(k[c], 0.0) / std::max(denom, a_bot_clip);
        mut_[c] = mut[c];
    }
}

void KOmegaSST_Full::solve(double /*dt*/, solver::FieldRegistry& f)
{
    if (!mesh_)
        return;
    auto& fnc = f;
    auto& k = *fnc.find_scalar("k");
    auto& w = *fnc.find_scalar("omega");
    const auto& d = *fnc.find_scalar("wallDistance");
    const std::size_t nC = k.size();

    compute_strain_rate(f);

    // ----- compute F1, F2 blending (Menter SST 2003) -----------------------
    {
        solver::LeastSquaresGradient G(*mesh_);
        solver::VectorField gk, gw;
        G.evaluate(k, gk);
        G.evaluate(w, gw);
        for (std::size_t c = 0; c < nC; ++c) {
            const double dk_dw = gk.x[c] * gw.x[c] + gk.y[c] * gw.y[c] + gk.z[c] * gw.z[c];
            CDk_[c] = std::max(2.0 * rho_ * sigma_w2 / std::max(w[c], a_bot_clip) * dk_dw, 1e-10);
            const double yy = std::max(d[c], 1e-12);
            const double nu = mu_ / rho_;
            const double arg1a =
                std::sqrt(std::max(k[c], 0.0)) / (beta_star * std::max(w[c], a_bot_clip) * yy);
            const double arg1b = 500.0 * nu / (yy * yy * std::max(w[c], a_bot_clip));
            const double arg1c = 4.0 * rho_ * sigma_w2 * std::max(k[c], 0.0) / (CDk_[c] * yy * yy);
            const double arg1 = std::min(std::max(arg1a, arg1b), arg1c);
            F1_[c] = std::tanh(arg1 * arg1 * arg1 * arg1);
            const double arg2 = std::max(2.0 * arg1a, arg1b);
            F2_[c] = std::tanh(arg2 * arg2);
        }
    }

    update_mut(f);

    // ----- production term -------------------------------------------------
    for (std::size_t c = 0; c < nC; ++c) {
        const double Pk_unclipped = mut_[c] * S_[c] * S_[c];
        Pk_[c] =
            std::min(Pk_unclipped,
                     10.0 * beta_star * rho_ * std::max(k[c], 0.0) * std::max(w[c], a_bot_clip));
    }
    // Stash production into a scratch scalar field for ScalarTransport source.
    auto& Pkfield = fnc.scalar("__sst_Pk", nC);
    auto& Pwfield = fnc.scalar("__sst_Pw", nC);
    for (std::size_t c = 0; c < nC; ++c) {
        Pkfield[c] = Pk_[c] - beta_star * rho_ * std::max(k[c], 0.0) * std::max(w[c], a_bot_clip);
        const double sigma_w_blend = F1_[c] * sigma_w1 + (1 - F1_[c]) * sigma_w2;
        const double gamma_blend = F1_[c] * gamma_1 + (1 - F1_[c]) * gamma_2;
        const double beta_blend = F1_[c] * beta_1 + (1 - F1_[c]) * beta_2;
        (void) sigma_w_blend;
        Pwfield[c] = gamma_blend * rho_ * S_[c] * S_[c] - beta_blend * rho_ * w[c] * w[c]
                     + (1 - F1_[c]) * CDk_[c];
    }

    // ----- solve k transport ----------------------------------------------
    const double sigma_k_avg = 0.5 * (sigma_k1 + sigma_k2); // for global Γ scale
    kEq_->set_diffusivity(mu_ + sigma_k_avg * 0.0);         // diffusion mut handled per-cell below
    kEq_->set_source("__sst_Pk");
    // Boundary conditions: wall → k = 0 (Dirichlet); inlet/outlet zero-grad.
    for (const auto& b : bcs_) {
        if (b.type == solver::BCType::NoSlipWall || b.type == solver::BCType::Wall) {
            kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
            wEq_->add_bc({b.zone,
                          solver::ScalarBC::Kind::Dirichlet,
                          60.0 * mu_ / (rho_ * beta_1 * 1e-12),
                          0.0});
        } else {
            kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
            wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
    kEq_->solve_iteration();

    // ----- solve omega transport ------------------------------------------
    wEq_->set_diffusivity(mu_);
    wEq_->set_source("__sst_Pw");
    wEq_->solve_iteration();

    // Recompute μ_t after k/ω update so the next outer iteration sees fresh values.
    update_mut(f);
}

double KOmegaSST_Full::turbulent_viscosity(std::size_t c) const
{
    return (c < mut_.size()) ? mut_[c] : 0.0;
}

} // namespace simall::turbulence
