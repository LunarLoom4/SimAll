// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/SpalartAllmaras.cpp
// =============================================================================
#include "turbulence/SpalartAllmaras.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence
{

namespace
{
constexpr double cb1 = 0.1355;
constexpr double cb2 = 0.622;
constexpr double sigma = 2.0 / 3.0;
constexpr double cv1 = 7.1;
constexpr double cv1_3 = cv1 * cv1 * cv1;
constexpr double kappa = 0.41;
constexpr double cw2 = 0.3;
constexpr double cw3 = 2.0;
constexpr double cw3_6 = cw3 * cw3 * cw3 * cw3 * cw3 * cw3;
constexpr double cw1 = cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
} // namespace

void SpalartAllmaras_Full::initialize(meshing::Mesh& m, solver::FieldRegistry& f)
{
    mesh_ = &m;
    fields_ = &f;
    const std::size_t nC = m.cells().size();

    solver::LinearSolverConfig c;
    c.kind = solver::LinearSolverKind::GMRES;
    c.preconditioner = solver::PreconditionerKind::ILU;
    c.tolerance = 1e-6;
    c.maxIterations = 100;
    c.restart = 20;
    lin_ = solver::make_linear_solver(c);
    wallDist_ = solver::make_wall_distance_exact();

    f.scalar("nuTilde", nC);
    f.scalar("mut", nC);
    f.scalar("wallDistance", nC);
    f.scalar("__sa_src", nC);

    auto& d = *f.find_scalar("wallDistance");
    wallDist_->compute(m, bcs_, d);

    auto& nt = *f.find_scalar("nuTilde");
    const double nu = mu_ / rho_;
    std::fill(nt.begin(), nt.end(), 3.0 * nu); // common free-stream init

    mut_.assign(nC, 0.0);
    Stilde_.assign(nC, 0.0);
    src_.assign(nC, 0.0);

    nuTildeEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    nuTildeEq_->set_field("nuTilde");
    nuTildeEq_->set_density(rho_);
    nuTildeEq_->set_source(*f.find_scalar("__sa_src"));
    nuTildeEq_->set_urf(0.7);

    for (const auto& b : bcs_) {
        // Wall: ν̃ = 0; Inlet: ν̃ = 3ν; otherwise zero gradient.
        switch (b.type) {
        case solver::BCType::NoSlipWall:
        case solver::BCType::Wall:
        case solver::BCType::MovingWall:
            nuTildeEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
            break;
        case solver::BCType::VelocityInlet:
        case solver::BCType::PressureInlet:
        case solver::BCType::MassFlowInlet:
            nuTildeEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 3.0 * nu, 0.0});
            break;
        default:
            nuTildeEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
    SIMALL_LOG_INFO("Turbulence", "Spalart-Allmaras initialized (cells=", nC, ")");
}

void SpalartAllmaras_Full::solve(double /*dt*/, solver::FieldRegistry& f)
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U = *f.find_vector("U");
    auto& nt = *f.find_scalar("nuTilde");
    auto& d = *f.find_scalar("wallDistance");
    auto& mut = *f.find_scalar("mut");
    auto& src = *f.find_scalar("__sa_src");
    const double nu = mu_ / rho_;

    // ---- vorticity magnitude Ω = |∇×U|
    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz, gNt;
    G.evaluate(U.x, gUx);
    G.evaluate(U.y, gUy);
    G.evaluate(U.z, gUz);
    G.evaluate(nt, gNt);

    for (std::size_t c = 0; c < nC; ++c) {
        const double wx = gUz.y[c] - gUy.z[c];
        const double wy = gUx.z[c] - gUz.x[c];
        const double wz = gUy.x[c] - gUx.y[c];
        const double Omega = std::sqrt(wx * wx + wy * wy + wz * wz);

        const double nT = std::max(nt[c], 0.0);
        const double chi = nT / std::max(nu, 1e-30);
        const double chi3 = chi * chi * chi;
        const double fv1 = chi3 / (chi3 + cv1_3);
        const double fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const double dist = std::max(d[c], 1e-10);
        const double Sbar = nT / (kappa * kappa * dist * dist) * fv2;
        const double Stld = Omega + Sbar;
        Stilde_[c] = Stld;

        const double r = std::min(nT / (std::max(Stld, 1e-30) * kappa * kappa * dist * dist), 10.0);
        const double g = r + cw2 * (std::pow(r, 6.0) - r);
        const double g6 = std::pow(g, 6.0);
        const double fw = g * std::pow((1.0 + cw3_6) / (g6 + cw3_6), 1.0 / 6.0);

        const double P = cb1 * Stld * nT;
        const double D = cw1 * fw * (nT * nT) / (dist * dist);
        const double crossDiff =
            (cb2 / sigma) * (gNt.x[c] * gNt.x[c] + gNt.y[c] * gNt.y[c] + gNt.z[c] * gNt.z[c]);

        src[c] = rho_ * (P - D + crossDiff);

        // μ_t update
        mut[c] = rho_ * nT * fv1;
        mut_[c] = mut[c];
    }

    // ---- transport with diffusivity (ν + ν̃)/σ (mean of cell + neighbour
    // handled internally by ScalarTransport via constant set_diffusivity;
    // we use a representative cell value approximation: σ_avg ≈ (μ + ρ·avg_nuTilde)/σ).
    double avgNt = 0.0;
    for (std::size_t c = 0; c < nC; ++c)
        avgNt += std::max(nt[c], 0.0);
    avgNt /= std::max<std::size_t>(nC, 1);
    nuTildeEq_->set_diffusivity((mu_ + rho_ * avgNt) / sigma);
    nuTildeEq_->solve_iteration();

    for (std::size_t c = 0; c < nC; ++c)
        nt[c] = std::max(nt[c], 0.0);
}

} // namespace simall::turbulence
