// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KOmegaStandard.cpp
// =============================================================================
#include "turbulence/KOmegaStandard.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence
{

namespace
{
constexpr double beta_star = 9.0 / 100.0; // 0.09
constexpr double alpha = 5.0 / 9.0;
constexpr double beta = 3.0 / 40.0;
constexpr double sigK = 0.5;
constexpr double sigW = 0.5;
constexpr double C_lim = 7.0 / 8.0; // Wilcox 2006 stress limiter
constexpr double kFloor = 1.0e-12;
constexpr double wFloor = 1.0e-12;
} // namespace

void KOmegaStandard_Full::initialize(meshing::Mesh& m, solver::FieldRegistry& f)
{
    mesh_ = &m;
    const std::size_t nC = m.cells().size();

    solver::LinearSolverConfig c;
    c.kind = solver::LinearSolverKind::GMRES;
    c.preconditioner = solver::PreconditionerKind::ILU;
    c.tolerance = 1e-6;
    c.maxIterations = 100;
    c.restart = 20;
    lin_ = solver::make_linear_solver(c);

    f.scalar("k", nC);
    f.scalar("omega", nC);
    f.scalar("mut", nC);
    f.scalar("__kos_srcK", nC);
    f.scalar("__kos_srcW", nC);

    auto& k = *f.find_scalar("k");
    auto& w = *f.find_scalar("omega");
    std::fill(k.begin(), k.end(), 1.0e-4);
    std::fill(w.begin(), w.end(), 1.0);

    mut_.assign(nC, 0.0);

    kEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    wEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    kEq_->set_field("k");
    wEq_->set_field("omega");
    kEq_->set_density(rho_);
    wEq_->set_density(rho_);
    kEq_->set_source(std::string{"__kos_srcK"});
    wEq_->set_source(std::string{"__kos_srcW"});
    kEq_->set_urf(0.7);
    wEq_->set_urf(0.7);

    for (const auto& b : bcs_) {
        switch (b.type) {
        case solver::BCType::NoSlipWall:
        case solver::BCType::Wall:
            kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
            // ω wall:  ω = 60 ν / (β y1²)  →  large value
            wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e6, 0.0});
            break;
        case solver::BCType::VelocityInlet:
        case solver::BCType::PressureInlet:
        case solver::BCType::MassFlowInlet:
            kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e-3, 0.0});
            wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0, 0.0});
            break;
        default:
            kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
            wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
    SIMALL_LOG_INFO("Turbulence", "Wilcox k-ω initialised (cells=", nC, ")");
}

void KOmegaStandard_Full::solve(double /*dt*/, solver::FieldRegistry& f)
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U = *f.find_vector("U");
    auto& k = *f.find_scalar("k");
    auto& w = *f.find_scalar("omega");
    auto& mut = *f.find_scalar("mut");
    auto& srcK = *f.find_scalar("__kos_srcK");
    auto& srcW = *f.find_scalar("__kos_srcW");

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx);
    G.evaluate(U.y, gUy);
    G.evaluate(U.z, gUz);

    for (std::size_t c = 0; c < nC; ++c) {
        const double ki = std::max(k[c], kFloor);
        const double wi = std::max(w[c], wFloor);

        const double s11 = gUx.x[c], s22 = gUy.y[c], s33 = gUz.z[c];
        const double s12 = 0.5 * (gUx.y[c] + gUy.x[c]);
        const double s13 = 0.5 * (gUx.z[c] + gUz.x[c]);
        const double s23 = 0.5 * (gUy.z[c] + gUz.y[c]);
        const double SijSij =
            s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
        const double S2 = 2.0 * SijSij;
        const double Smag = std::sqrt(S2);

        // Wilcox 2006 stress limiter:
        const double wbar = std::max(wi, C_lim * std::sqrt(S2 / beta_star));
        mut[c] = rho_ * ki / wbar;
        mut_[c] = mut[c];

        const double Pk = mut[c] * S2;
        srcK[c] = Pk - beta_star * rho_ * ki * wi;
        srcW[c] = alpha * (wi / ki) * Pk - beta * rho_ * wi * wi;
    }

    double mutAvg = 0.0;
    for (std::size_t c = 0; c < nC; ++c)
        mutAvg += mut[c];
    mutAvg /= std::max<std::size_t>(nC, 1);
    kEq_->set_diffusivity(mu_ + sigK * mutAvg);
    wEq_->set_diffusivity(mu_ + sigW * mutAvg);
    kEq_->solve_iteration();
    wEq_->solve_iteration();

    for (std::size_t c = 0; c < nC; ++c) {
        k[c] = std::max(k[c], kFloor);
        w[c] = std::max(w[c], wFloor);
    }
}

} // namespace simall::turbulence
