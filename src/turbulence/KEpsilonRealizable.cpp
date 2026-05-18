// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KEpsilonRealizable.cpp
// =============================================================================
#include "turbulence/KEpsilonRealizable.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence {

namespace {
constexpr double A0      = 4.04;
constexpr double sigK    = 1.0;
constexpr double sigE    = 1.2;
constexpr double C2eps   = 1.9;
constexpr double kFloor  = 1.0e-12;
constexpr double eFloor  = 1.0e-12;
}

void KEpsilonRealizable_Full::initialize(meshing::Mesh& m, solver::FieldRegistry& f) {
    mesh_ = &m;
    const std::size_t nC = m.cells().size();

    solver::LinearSolverConfig c;
    c.kind = solver::LinearSolverKind::GMRES;
    c.preconditioner = solver::PreconditionerKind::ILU;
    c.tolerance = 1e-6; c.maxIterations = 100; c.restart = 20;
    lin_ = solver::make_linear_solver(c);

    f.scalar("k",          nC);
    f.scalar("epsilon",    nC);
    f.scalar("mut",        nC);
    f.scalar("__rke_srcK", nC);
    f.scalar("__rke_srcE", nC);

    auto& k = *f.find_scalar("k");
    auto& e = *f.find_scalar("epsilon");
    std::fill(k.begin(), k.end(), 1.0e-4);
    std::fill(e.begin(), e.end(), 1.0e-4);

    mut_.assign(nC, 0.0);

    kEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    eEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    kEq_->set_field("k");          eEq_->set_field("epsilon");
    kEq_->set_density(rho_);       eEq_->set_density(rho_);
    kEq_->set_source(std::string{"__rke_srcK"});
    eEq_->set_source(std::string{"__rke_srcE"});
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
    SIMALL_LOG_INFO("Turbulence", "Realisable k-ε initialised (cells=", nC, ")");
}

void KEpsilonRealizable_Full::solve(double /*dt*/, solver::FieldRegistry& f) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U  = *f.find_vector("U");
    auto& k        = *f.find_scalar("k");
    auto& e        = *f.find_scalar("epsilon");
    auto& mut      = *f.find_scalar("mut");
    auto& srcK     = *f.find_scalar("__rke_srcK");
    auto& srcE     = *f.find_scalar("__rke_srcE");

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx); G.evaluate(U.y, gUy); G.evaluate(U.z, gUz);

    const double nu = mu_ / rho_;
    for (std::size_t c = 0; c < nC; ++c) {
        const double ki = std::max(k[c], kFloor);
        const double ei = std::max(e[c], eFloor);

        const double s11=gUx.x[c], s22=gUy.y[c], s33=gUz.z[c];
        const double s12=0.5*(gUx.y[c]+gUy.x[c]);
        const double s13=0.5*(gUx.z[c]+gUz.x[c]);
        const double s23=0.5*(gUy.z[c]+gUz.y[c]);
        const double w12=0.5*(gUx.y[c]-gUy.x[c]);
        const double w13=0.5*(gUx.z[c]-gUz.x[c]);
        const double w23=0.5*(gUy.z[c]-gUz.y[c]);

        // SijSij
        const double SijSij = s11*s11+s22*s22+s33*s33
                            + 2.0*(s12*s12+s13*s13+s23*s23);
        // Ω̃ij Ω̃ij  (no system rotation)
        const double WijWij = 2.0*(w12*w12+w13*w13+w23*w23);

        const double Smag = std::sqrt(2.0 * SijSij);
        const double Ustar = std::sqrt(SijSij + WijWij);

        // φ = (1/3) arccos( √6 W ),  W = S_ij S_jk S_ki / S̃³,  S̃=√(SijSij)
        double Stilde3 = std::pow(std::max(SijSij, 1e-30), 1.5);
        double SSS =  s11*(s11*s11 + s12*s12 + s13*s13)
                    + s22*(s12*s12 + s22*s22 + s23*s23)
                    + s33*(s13*s13 + s23*s23 + s33*s33)
                    + 2.0*( s12*(s11*s12 + s12*s22 + s13*s23)
                          + s13*(s11*s13 + s12*s23 + s13*s33)
                          + s23*(s12*s13 + s22*s23 + s23*s33));
        double Wnorm = std::sqrt(6.0) * SSS / Stilde3;
        Wnorm = std::clamp(Wnorm, -1.0, 1.0);
        const double phi = (1.0/3.0) * std::acos(Wnorm);
        const double As  = std::sqrt(6.0) * std::cos(phi);
        const double Cmu = 1.0 / (A0 + As * ki * Ustar / ei);

        mut[c]  = rho_ * Cmu * ki*ki / ei;
        mut_[c] = mut[c];

        const double Pk = mut[c] * 2.0 * SijSij;

        const double eta = Smag * ki / ei;
        const double C1  = std::max(0.43, eta / (eta + 5.0));

        srcK[c] = Pk - rho_ * ei;
        srcE[c] = rho_ * C1 * Smag * ei
                - rho_ * C2eps * ei * ei / (ki + std::sqrt(nu * ei));
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
