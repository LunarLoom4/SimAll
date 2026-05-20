// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/KKLOmegaTransition.cpp
// =============================================================================
#include "turbulence/KKLOmegaTransition.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence
{

namespace
{
// Walters & Cokljat 2008 constants.
constexpr double A0 = 4.04;
constexpr double As = 2.12;
constexpr double Av = 6.75;
constexpr double Abp = 0.6;
constexpr double Anat = 200.0;
constexpr double Ats = 200.0;
constexpr double CbpCrit = 1.2;
constexpr double Cnc = 0.1;
constexpr double CnatCrit = 1250.0;
constexpr double Cints = 0.75;
constexpr double Ctscrit = 1000.0;
constexpr double Crnat = 0.02;
constexpr double C11 = 3.4e-6;
constexpr double C12 = 1.0e-10;
constexpr double Cr = 0.12;
constexpr double Calpha_theta = 0.035;
constexpr double Css = 1.5;
constexpr double CtauL = 4360.0;
constexpr double Cw1 = 0.44;
constexpr double Cw2 = 0.92;
constexpr double Cw3 = 0.3;
constexpr double CwR = 1.5;
constexpr double Clam_x = 2.495;
constexpr double Cmu_std = 0.09;
constexpr double Prtheta = 0.85;
constexpr double SigmaK = 1.0;
constexpr double SigmaW = 1.17;
constexpr double kFloor = 1.0e-12;
constexpr double wFloor = 1.0e-12;
constexpr double kappa = 0.41;
} // namespace

void KKLOmegaTransition_Full::initialize(meshing::Mesh& m, solver::FieldRegistry& f)
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

    f.scalar("kT", nC);
    f.scalar("kL", nC);
    f.scalar("omega", nC);
    f.scalar("mut", nC);
    f.scalar("wallDistance", nC);
    f.scalar("__kkl_srcKT", nC);
    f.scalar("__kkl_srcKL", nC);
    f.scalar("__kkl_srcW", nC);

    auto& kT = *f.find_scalar("kT");
    auto& kL = *f.find_scalar("kL");
    auto& w = *f.find_scalar("omega");
    std::fill(kT.begin(), kT.end(), 1.0e-5);
    std::fill(kL.begin(), kL.end(), 1.0e-7);
    std::fill(w.begin(), w.end(), 1.0);

    mut_.assign(nC, 0.0);

    wallDist_ = solver::make_wall_distance_exact();
    wallDist_->compute(m, bcs_, *f.find_scalar("wallDistance"));

    kTEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    kLEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    wEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    kTEq_->set_field("kT");
    kLEq_->set_field("kL");
    wEq_->set_field("omega");
    kTEq_->set_density(rho_);
    kLEq_->set_density(rho_);
    wEq_->set_density(rho_);
    kTEq_->set_source(std::string{"__kkl_srcKT"});
    kLEq_->set_source(std::string{"__kkl_srcKL"});
    wEq_->set_source(std::string{"__kkl_srcW"});
    kTEq_->set_urf(0.7);
    kLEq_->set_urf(0.7);
    wEq_->set_urf(0.7);

    for (const auto& b : bcs_) {
        switch (b.type) {
        case solver::BCType::NoSlipWall:
        case solver::BCType::Wall:
            kTEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
            kLEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
            wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e6, 0.0});
            break;
        case solver::BCType::VelocityInlet:
        case solver::BCType::PressureInlet:
        case solver::BCType::MassFlowInlet:
            kTEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e-5, 0.0});
            kLEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e-7, 0.0});
            wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0, 0.0});
            break;
        default:
            kTEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
            kLEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
            wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
    SIMALL_LOG_INFO("Turbulence", "k-kL-ω transition initialised (cells=", nC, ")");
}

void KKLOmegaTransition_Full::solve(double /*dt*/, solver::FieldRegistry& f)
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U = *f.find_vector("U");
    auto& kT = *f.find_scalar("kT");
    auto& kL = *f.find_scalar("kL");
    auto& w = *f.find_scalar("omega");
    auto& mut = *f.find_scalar("mut");
    auto& d = *f.find_scalar("wallDistance");
    auto& sKT = *f.find_scalar("__kkl_srcKT");
    auto& sKL = *f.find_scalar("__kkl_srcKL");
    auto& sW = *f.find_scalar("__kkl_srcW");

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx);
    G.evaluate(U.y, gUy);
    G.evaluate(U.z, gUz);

    const double nu = mu_ / rho_;
    for (std::size_t c = 0; c < nC; ++c) {
        const double dw = std::max(d[c], 1e-12);
        const double kTv = std::max(kT[c], kFloor);
        const double kLv = std::max(kL[c], 0.0);
        const double wi = std::max(w[c], wFloor);

        // Strain & rotation.
        const double s11 = gUx.x[c], s22 = gUy.y[c], s33 = gUz.z[c];
        const double s12 = 0.5 * (gUx.y[c] + gUy.x[c]);
        const double s13 = 0.5 * (gUx.z[c] + gUz.x[c]);
        const double s23 = 0.5 * (gUy.z[c] + gUz.y[c]);
        const double w12 = 0.5 * (gUx.y[c] - gUy.x[c]);
        const double w13 = 0.5 * (gUx.z[c] - gUz.x[c]);
        const double w23 = 0.5 * (gUy.z[c] - gUz.y[c]);
        const double S = std::sqrt(
            2.0 * (s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23)));
        const double Om = std::sqrt(2.0 * (w12 * w12 + w13 * w13 + w23 * w23));

        // Effective length scales.
        const double lam_eff = std::min(Clam_x * dw, std::sqrt(kTv) / std::max(wi, 1e-12));
        const double f_w = std::pow(lam_eff / std::max(std::sqrt(kTv) / wi, 1e-12), 2.0 / 3.0);
        // Small-scale fraction
        const double Re_T = (rho_ * kTv) / (mu_ * wi);
        const double f_v = 1.0 - std::exp(-std::sqrt(Re_T) / Av);
        const double f_INT = std::min(kLv / std::max(Cints * (kLv + kTv), 1e-30), 1.0);
        const double f_SS = std::exp(-std::pow(Css * nu * Om / std::max(kTv, 1e-30), 2));
        const double kTs = f_SS * f_w * kTv;
        const double mu_ts = f_v * f_INT * Cmu_std * std::sqrt(kTs) * lam_eff * rho_;

        // Large-scale (transition production) viscosity.
        const double phi_NAT =
            std::max(rho_ * dw * std::sqrt(kLv) / std::max(mu_, 1e-30) - CnatCrit / Anat, 0.0);
        const double mu_tl_max = 0.5 * (kLv + kTv) / std::max(S, 1e-12);
        const double beta_TS =
            1.0
            - std::exp(-std::max(rho_ * dw * dw * Om / std::max(mu_, 1e-30) - Ctscrit, 0.0) / Ats);
        const double mu_tl =
            std::min(C11 * (phi_NAT * phi_NAT * phi_NAT * phi_NAT) * Om * lam_eff * lam_eff
                         + C12 * beta_TS * f_w * f_w * phi_NAT * Om * lam_eff * lam_eff,
                     mu_tl_max);

        const double mu_t = mu_ts + mu_tl;
        mut[c] = mu_t;
        mut_[c] = mu_t;

        // Production for k_T (turbulent) — based on small-scale eddy viscosity.
        const double P_kT = mu_ts * S * S;
        // Production for k_L (laminar fluctuations).
        const double P_kL = mu_tl * S * S;

        // Bypass transition production.
        const double phi_BP = std::max(rho_ * kTv / (mu_ * std::max(Om, 1e-12)) - CbpCrit, 0.0);
        const double f_NATcrit = 1.0 - std::exp(-Cnc * std::sqrt(kLv) * dw / std::max(nu, 1e-30));
        const double R_BP = Cr * phi_BP * kLv * beta_TS / std::max(f_w, 1e-12);
        const double R_NAT = Crnat * f_NATcrit * Anat * kLv * Om;

        // ω equation.
        const double D_T = 2.0 * mu_ * (std::sqrt(kTv) / dw) * (std::sqrt(kTv) / dw);
        const double D_L = 2.0 * mu_ * (std::sqrt(kLv) / dw) * (std::sqrt(kLv) / dw);

        sKT[c] = P_kT + R_BP + R_NAT - rho_ * wi * kTv - D_T;
        sKL[c] = P_kL - R_BP - R_NAT - D_L;
        sW[c] = Cw1 * (wi / std::max(kTv, kFloor)) * P_kT
                + (CwR / std::max(f_w, 1e-12) - 1.0) * (wi / std::max(kTv, kFloor)) * (R_BP + R_NAT)
                - Cw2 * rho_ * wi * wi + Cw3 * f_w * f_w * std::sqrt(kTv) / (dw * dw * dw);
        // Reference unused constants to silence linker/compile noise (the model uses several
        // auxiliary constants):
        (void) A0;
        (void) As;
        (void) Abp;
        (void) Anat;
        (void) Ats;
        (void) CbpCrit;
        (void) Calpha_theta;
        (void) CtauL;
        (void) Cmu_std;
        (void) Prtheta;
    }

    double mutAvg = 0.0;
    for (std::size_t c = 0; c < nC; ++c)
        mutAvg += mut[c];
    mutAvg /= std::max<std::size_t>(nC, 1);
    kTEq_->set_diffusivity(mu_ + mutAvg / SigmaK);
    kLEq_->set_diffusivity(mu_);
    wEq_->set_diffusivity(mu_ + mutAvg / SigmaW);
    kTEq_->solve_iteration();
    kLEq_->solve_iteration();
    wEq_->solve_iteration();

    for (std::size_t c = 0; c < nC; ++c) {
        kT[c] = std::max(kT[c], kFloor);
        kL[c] = std::max(kL[c], 0.0);
        w[c] = std::max(w[c], wFloor);
    }
}

} // namespace simall::turbulence
