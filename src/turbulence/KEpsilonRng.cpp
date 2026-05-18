// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KEpsilonRng.cpp
// =============================================================================
#include "turbulence/KEpsilonRng.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence {

namespace {
// RNG k-ε constants (Yakhot, Orszag et al. 1992).
constexpr double Cmu    = 0.0845;
constexpr double sigK   = 0.7194;
constexpr double sigE   = 0.7194;
constexpr double C1eps  = 1.42;
constexpr double C2eps  = 1.68;
constexpr double eta0   = 4.38;
constexpr double beta   = 0.012;
constexpr double kFloor = 1.0e-12;
constexpr double eFloor = 1.0e-12;
}  // namespace

void KEpsilonRng_Full::initialize(meshing::Mesh& m, solver::FieldRegistry& f) {
    mesh_ = &m;
    const std::size_t nC = m.cells().size();

    solver::LinearSolverConfig c;
    c.kind = solver::LinearSolverKind::GMRES;
    c.preconditioner = solver::PreconditionerKind::ILU;
    c.tolerance = 1e-6; c.maxIterations = 100; c.restart = 20;
    lin_ = solver::make_linear_solver(c);

    f.scalar("k",            nC);
    f.scalar("epsilon",      nC);
    f.scalar("mut",          nC);
    f.scalar("__rng_srcK",   nC);
    f.scalar("__rng_srcE",   nC);

    auto& k = *f.find_scalar("k");
    auto& e = *f.find_scalar("epsilon");
    std::fill(k.begin(), k.end(), 1.0e-4);
    std::fill(e.begin(), e.end(), 1.0e-4);

    mut_.assign(nC, 0.0);
    Pk_ .assign(nC, 0.0);

    kEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    eEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    kEq_->set_field("k");          eEq_->set_field("epsilon");
    kEq_->set_density(rho_);       eEq_->set_density(rho_);
    kEq_->set_source(std::string{"__rng_srcK"});
    eEq_->set_source(std::string{"__rng_srcE"});
    kEq_->set_urf(0.7);            eEq_->set_urf(0.7);

    for (const auto& b : bcs_) {
        switch (b.type) {
            case solver::BCType::NoSlipWall:
            case solver::BCType::Wall:
                kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
                eEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann,   0.0, 0.0});
                break;
            case solver::BCType::VelocityInlet:
            case solver::BCType::PressureInlet:
            case solver::BCType::MassFlowInlet:
                kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e-3, 0.0});
                eEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e-3, 0.0});
                break;
            default:
                kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
                eEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
    SIMALL_LOG_INFO("Turbulence", "RNG k-ε initialised (cells=", nC, ")");
}

void KEpsilonRng_Full::solve(double /*dt*/, solver::FieldRegistry& f) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U = *f.find_vector("U");
    auto& k       = *f.find_scalar("k");
    auto& e       = *f.find_scalar("epsilon");
    auto& mut     = *f.find_scalar("mut");
    auto& srcK    = *f.find_scalar("__rng_srcK");
    auto& srcE    = *f.find_scalar("__rng_srcE");

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx); G.evaluate(U.y, gUy); G.evaluate(U.z, gUz);

    for (std::size_t c = 0; c < nC; ++c) {
        const double ki = std::max(k[c], kFloor);
        const double ei = std::max(e[c], eFloor);
        mut[c]  = rho_ * Cmu * ki * ki / ei;
        mut_[c] = mut[c];

        const double s11=gUx.x[c], s22=gUy.y[c], s33=gUz.z[c];
        const double s12=0.5*(gUx.y[c]+gUy.x[c]);
        const double s13=0.5*(gUx.z[c]+gUz.x[c]);
        const double s23=0.5*(gUy.z[c]+gUz.y[c]);
        const double S2 = 2.0*(s11*s11+s22*s22+s33*s33
                             + 2.0*(s12*s12+s13*s13+s23*s23));
        const double Smag = std::sqrt(S2);
        Pk_[c] = mut[c] * S2;

        // RNG η-term: η = S k / ε
        const double eta  = Smag * ki / ei;
        const double Reta = Cmu * eta*eta*eta * (1.0 - eta / eta0)
                          / (1.0 + beta * eta*eta*eta);
        const double C2star = C2eps + Reta;

        srcK[c] = Pk_[c] - rho_ * ei;
        srcE[c] = (ei / ki) * (C1eps * Pk_[c] - C2star * rho_ * ei);
    }

    double mutAvg = 0.0;
    for (std::size_t c = 0; c < nC; ++c) mutAvg += mut[c];
    mutAvg /= std::max<std::size_t>(nC, 1);
    kEq_->set_diffusivity(mu_ + mutAvg / sigK);
    eEq_->set_diffusivity(mu_ + mutAvg / sigE);
    kEq_->solve_iteration();
    eEq_->solve_iteration();

    for (std::size_t c = 0; c < nC; ++c) {
        k[c] = std::max(k[c], kFloor);
        e[c] = std::max(e[c], eFloor);
    }
}

}  // namespace simall::turbulence
