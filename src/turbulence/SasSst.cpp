// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/SasSst.cpp
// =============================================================================
#include "turbulence/SasSst.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence {

namespace {
// SST inner (1) / outer (2) constants — Menter, Kuntz, Langtry 2003.
constexpr double sigK1=0.85, sigK2=1.00;
constexpr double sigW1=0.5,  sigW2=0.856;
constexpr double beta1=0.075, beta2=0.0828;
constexpr double betaStar=0.09;
constexpr double a1=0.31;
constexpr double kappa = 0.41;
const     double gamma1 = beta1/betaStar - sigW1*kappa*kappa/std::sqrt(betaStar);
const     double gamma2 = beta2/betaStar - sigW2*kappa*kappa/std::sqrt(betaStar);
// SAS constants — Menter & Egorov 2010.
constexpr double zeta2  = 3.51;
constexpr double sigPhi = 2.0/3.0;
constexpr double C_SAS  = 2.0;
constexpr double kFloor = 1.0e-12;
constexpr double wFloor = 1.0e-12;
}

void SasSst_Full::initialize(meshing::Mesh& m, solver::FieldRegistry& f) {
    mesh_ = &m;
    const std::size_t nC = m.cells().size();

    solver::LinearSolverConfig c;
    c.kind = solver::LinearSolverKind::GMRES;
    c.preconditioner = solver::PreconditionerKind::ILU;
    c.tolerance = 1e-6; c.maxIterations = 100; c.restart = 20;
    lin_ = solver::make_linear_solver(c);

    f.scalar("k", nC); f.scalar("omega", nC); f.scalar("mut", nC);
    f.scalar("wallDistance", nC);
    f.scalar("__sas_srcK", nC);
    f.scalar("__sas_srcW", nC);

    auto& kk = *f.find_scalar("k");
    auto& ww = *f.find_scalar("omega");
    std::fill(kk.begin(), kk.end(), 1.0e-4);
    std::fill(ww.begin(), ww.end(), 1.0);

    mut_.assign(nC, 0.0); F1_.assign(nC, 0.0); F2_.assign(nC, 0.0);

    wallDist_ = solver::make_wall_distance_exact();
    wallDist_->compute(m, bcs_, *f.find_scalar("wallDistance"));

    kEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    wEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    kEq_->set_field("k");          wEq_->set_field("omega");
    kEq_->set_density(rho_);       wEq_->set_density(rho_);
    kEq_->set_source(std::string{"__sas_srcK"});
    wEq_->set_source(std::string{"__sas_srcW"});
    kEq_->set_urf(0.7);            wEq_->set_urf(0.7);

    for (const auto& b : bcs_) {
        switch (b.type) {
            case solver::BCType::NoSlipWall:
            case solver::BCType::Wall:
                kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0,   0.0});
                wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e6, 0.0});
                break;
            case solver::BCType::VelocityInlet:
            case solver::BCType::PressureInlet:
            case solver::BCType::MassFlowInlet:
                kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e-3, 0.0});
                wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0,    0.0});
                break;
            default:
                kEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
                wEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
    SIMALL_LOG_INFO("Turbulence", "SAS-SST initialised (cells=", nC, ")");
}

void SasSst_Full::solve(double /*dt*/, solver::FieldRegistry& f) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U  = *f.find_vector("U");
    auto& k        = *f.find_scalar("k");
    auto& w        = *f.find_scalar("omega");
    auto& mut      = *f.find_scalar("mut");
    auto& d        = *f.find_scalar("wallDistance");
    auto& srcK     = *f.find_scalar("__sas_srcK");
    auto& srcW     = *f.find_scalar("__sas_srcW");

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx); G.evaluate(U.y, gUy); G.evaluate(U.z, gUz);
    // ∇k, ∇ω
    solver::VectorField gK, gW;
    G.evaluate(k, gK); G.evaluate(w, gW);
    // 2nd derivatives of U components (∇²U_i ≈ ∂_j ∂_j U_i).
    solver::VectorField ggUxx, ggUxy, ggUxz, ggUyx, ggUyy, ggUyz, ggUzx, ggUzy, ggUzz;
    G.evaluate(gUx.x, ggUxx); G.evaluate(gUx.y, ggUxy); G.evaluate(gUx.z, ggUxz);
    G.evaluate(gUy.x, ggUyx); G.evaluate(gUy.y, ggUyy); G.evaluate(gUy.z, ggUyz);
    G.evaluate(gUz.x, ggUzx); G.evaluate(gUz.y, ggUzy); G.evaluate(gUz.z, ggUzz);

    const double nu = mu_ / rho_;
    for (std::size_t c = 0; c < nC; ++c) {
        const double ki = std::max(k[c], kFloor);
        const double wi = std::max(w[c], wFloor);
        const double dw = std::max(d[c], 1e-12);

        // strain & rotation
        const double s11=gUx.x[c], s22=gUy.y[c], s33=gUz.z[c];
        const double s12=0.5*(gUx.y[c]+gUy.x[c]);
        const double s13=0.5*(gUx.z[c]+gUz.x[c]);
        const double s23=0.5*(gUy.z[c]+gUz.y[c]);
        const double SS = 2.0*(s11*s11+s22*s22+s33*s33
                             + 2.0*(s12*s12+s13*s13+s23*s23));
        const double Smag = std::sqrt(std::max(SS, 0.0));

        // SST blending
        const double CD_kw = std::max(2.0*rho_*sigW2*(gK.x[c]*gW.x[c]+gK.y[c]*gW.y[c]+gK.z[c]*gW.z[c])/wi, 1e-10);
        const double arg1a = std::sqrt(ki)/(betaStar*wi*dw);
        const double arg1b = 500.0*nu/(dw*dw*wi);
        const double arg1c = 4.0*rho_*sigW2*ki/(CD_kw*dw*dw);
        const double arg1  = std::min(std::max(arg1a, arg1b), arg1c);
        F1_[c] = std::tanh(arg1*arg1*arg1*arg1);
        const double arg2  = std::max(2.0*arg1a, arg1b);
        F2_[c] = std::tanh(arg2*arg2);

        mut[c]  = rho_ * a1 * ki / std::max(a1*wi, Smag*F2_[c]);
        mut_[c] = mut[c];

        const double Pk_raw = mut[c] * SS;
        const double Pk     = std::min(Pk_raw, 10.0 * betaStar * rho_ * ki * wi);

        const double sigK = F1_[c]*sigK1 + (1.0-F1_[c])*sigK2;
        (void)sigK;
        const double gamma = F1_[c]*gamma1 + (1.0-F1_[c])*gamma2;
        const double beta_ = F1_[c]*beta1  + (1.0-F1_[c])*beta2;

        srcK[c] = Pk - betaStar*rho_*ki*wi;
        const double crossDiff = 2.0*(1.0-F1_[c])*rho_*sigW2*
            (gK.x[c]*gW.x[c]+gK.y[c]*gW.y[c]+gK.z[c]*gW.z[c]) / wi;
        const double P_omega = (gamma/std::max(mut[c]/rho_, 1e-30)) * Pk
                              - beta_*rho_*wi*wi + crossDiff;

        // SAS source: Q_SAS
        // L  = √k / (Cμ^{1/4} ω)
        const double L  = std::sqrt(ki) / (std::pow(betaStar, 0.25) * wi);
        // |∇²U| = √( Σ (∂_j ∂_j U_i)² ) ≈ √( Σ (ggUxx+ggUyy+ggUzz per component)² )
        const double lUx = ggUxx.x[c] + ggUxy.y[c] + ggUxz.z[c];
        const double lUy = ggUyx.x[c] + ggUyy.y[c] + ggUyz.z[c];
        const double lUz = ggUzx.x[c] + ggUzy.y[c] + ggUzz.z[c];
        const double lapU = std::sqrt(lUx*lUx + lUy*lUy + lUz*lUz);
        const double Lvk  = kappa * Smag / std::max(lapU, 1e-12);
        const double termA = zeta2 * kappa * SS * (L/std::max(Lvk,1e-12)) * (L/std::max(Lvk,1e-12));
        const double gradKsq = gK.x[c]*gK.x[c]+gK.y[c]*gK.y[c]+gK.z[c]*gK.z[c];
        const double gradWsq = gW.x[c]*gW.x[c]+gW.y[c]*gW.y[c]+gW.z[c]*gW.z[c];
        const double termB = C_SAS * 2.0 * ki / sigPhi *
            std::max(gradKsq/(ki*ki), gradWsq/(wi*wi));
        const double Qsas = std::max(rho_*(termA - termB), 0.0);

        srcW[c] = P_omega + Qsas;
    }

    double mutAvg=0.0;
    for (std::size_t c=0;c<nC;++c) mutAvg += mut[c];
    mutAvg /= std::max<std::size_t>(nC,1);
    // Use F1-blended diffusivities at the average (per-cell anisotropic ScalarTransport not supported here).
    const double sigK_avg = 0.5*(sigK1+sigK2);
    const double sigW_avg = 0.5*(sigW1+sigW2);
    kEq_->set_diffusivity(mu_ + sigK_avg*mutAvg);
    wEq_->set_diffusivity(mu_ + sigW_avg*mutAvg);
    kEq_->solve_iteration();
    wEq_->solve_iteration();

    for (std::size_t c=0;c<nC;++c) {
        k[c] = std::max(k[c], kFloor);
        w[c] = std::max(w[c], wFloor);
    }
}

}  // namespace simall::turbulence
