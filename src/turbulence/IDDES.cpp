// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/IDDES.cpp
// =============================================================================
#include "turbulence/IDDES.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence
{

namespace
{
// Spalart-Allmaras constants.
constexpr double cb1 = 0.1355, cb2 = 0.622, sigma = 2.0 / 3.0;
constexpr double kappa = 0.41;
constexpr double cw1 = cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
constexpr double cw2 = 0.3, cw3 = 2.0;
constexpr double cv1 = 7.1;
// IDDES constants — Shur, Spalart, Strelets, Travin 2008.
constexpr double C_DES = 0.65;
constexpr double C_w = 0.15;
constexpr double C_dt1 = 20.0;
constexpr double C_dt2 = 3.0;
constexpr double C_l = 5.0;
constexpr double C_t = 1.87;
} // namespace

void IDDES_Full::compute_hmax()
{
    const auto& C = mesh_->cells();
    const auto& F = mesh_->faces();
    const auto& N = mesh_->nodes();
    const std::size_t nC = C.size();
    hmax_.assign(nC, 0.0);
    // Loop faces of each cell; node-pair maximum distance approximates h_max.
    for (std::size_t c = 0; c < nC; ++c) {
        const int fbeg = C.faceOffsets[c];
        const int fend = C.faceOffsets[c + 1];
        double hmax = 0.0;
        for (int fi = fbeg; fi < fend; ++fi) {
            const auto f = C.faceIndices[fi];
            const int nbeg = F.nodeOffsets[f];
            const int nend = F.nodeOffsets[f + 1];
            for (int i = nbeg; i < nend; ++i)
                for (int j = i + 1; j < nend; ++j) {
                    const auto a = F.nodeIndices[i], b = F.nodeIndices[j];
                    const double dx = N.x[a] - N.x[b];
                    const double dy = N.y[a] - N.y[b];
                    const double dz = N.z[a] - N.z[b];
                    const double L = std::sqrt(dx * dx + dy * dy + dz * dz);
                    if (L > hmax)
                        hmax = L;
                }
        }
        hmax_[c] = (hmax > 0.0) ? hmax : std::cbrt(std::max(C.volume[c], 1e-30));
    }
}

void IDDES_Full::initialize(meshing::Mesh& m, solver::FieldRegistry& f)
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

    f.scalar("nuTilde", nC);
    f.scalar("mut", nC);
    f.scalar("wallDistance", nC);
    f.scalar("__iddes_src", nC);

    auto& nt = *f.find_scalar("nuTilde");
    std::fill(nt.begin(), nt.end(), 3.0 * mu_ / rho_);

    mut_.assign(nC, 0.0);
    delta_.assign(nC, 0.0);
    for (std::size_t i = 0; i < nC; ++i)
        delta_[i] = std::cbrt(std::max(m.cells().volume[i], 1e-30));
    compute_hmax();

    wallDist_ = solver::make_wall_distance_exact();
    wallDist_->compute(m, bcs_, *f.find_scalar("wallDistance"));

    nuTildeEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    nuTildeEq_->set_field("nuTilde");
    nuTildeEq_->set_density(rho_);
    nuTildeEq_->set_diffusivity(mu_ / sigma);
    nuTildeEq_->set_source(std::string{"__iddes_src"});
    nuTildeEq_->set_urf(0.7);

    for (const auto& b : bcs_) {
        switch (b.type) {
        case solver::BCType::NoSlipWall:
        case solver::BCType::Wall:
            nuTildeEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
            break;
        case solver::BCType::VelocityInlet:
        case solver::BCType::PressureInlet:
        case solver::BCType::MassFlowInlet:
            nuTildeEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 3.0 * mu_ / rho_, 0.0});
            break;
        default:
            nuTildeEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
    SIMALL_LOG_INFO("Turbulence", "IDDES (SA) initialised (cells=", nC, ")");
}

void IDDES_Full::solve(double /*dt*/, solver::FieldRegistry& f)
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U = *f.find_vector("U");
    auto& nt = *f.find_scalar("nuTilde");
    auto& d = *f.find_scalar("wallDistance");
    auto& mut = *f.find_scalar("mut");
    auto& src = *f.find_scalar("__iddes_src");

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx);
    G.evaluate(U.y, gUy);
    G.evaluate(U.z, gUz);

    const double nu = mu_ / rho_;
    for (std::size_t c = 0; c < nC; ++c) {
        const double s11 = gUx.x[c], s22 = gUy.y[c], s33 = gUz.z[c];
        const double s12 = 0.5 * (gUx.y[c] + gUy.x[c]);
        const double s13 = 0.5 * (gUx.z[c] + gUz.x[c]);
        const double s23 = 0.5 * (gUy.z[c] + gUz.y[c]);
        const double w12 = 0.5 * (gUx.y[c] - gUy.x[c]);
        const double w13 = 0.5 * (gUx.z[c] - gUz.x[c]);
        const double w23 = 0.5 * (gUy.z[c] - gUz.y[c]);
        const double SijSij =
            s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
        const double WijWij = 2.0 * (w12 * w12 + w13 * w13 + w23 * w23);
        const double Smag = std::sqrt(2.0 * SijSij);
        const double Omag = std::sqrt(WijWij);
        (void) Omag;
        const double UijUij =
            std::sqrt(gUx.x[c] * gUx.x[c] + gUx.y[c] * gUx.y[c] + gUx.z[c] * gUx.z[c]
                      + gUy.x[c] * gUy.x[c] + gUy.y[c] * gUy.y[c] + gUy.z[c] * gUy.z[c]
                      + gUz.x[c] * gUz.x[c] + gUz.y[c] * gUz.y[c] + gUz.z[c] * gUz.z[c]);

        const double dw = std::max(d[c], 1e-12);
        const double nt_c = std::max(nt[c], 0.0);
        const double chi = nt_c / std::max(nu, 1e-30);
        const double chi3 = chi * chi * chi;
        const double fv1 = chi3 / (chi3 + cv1 * cv1 * cv1);
        const double fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const double mu_t = rho_ * nt_c * fv1;
        const double nu_t = mu_t / rho_;
        mut[c] = mu_t;
        mut_[c] = mu_t;

        // Shielding function f_d (DDES).
        const double rd = (nu_t + nu) / (std::max(UijUij, 1e-10) * kappa * kappa * dw * dw);
        const double fd = 1.0 - std::tanh(std::pow(8.0 * rd, 3));

        // WMLES branch — f_dt and elevation f_e.
        const double rdt = nu_t / (std::max(UijUij, 1e-10) * kappa * kappa * dw * dw);
        const double fdt = 1.0 - std::tanh(std::pow(C_dt1 * rdt, C_dt2));
        const double fB = std::min(2.0
                                       * std::exp(-9.0 * (0.25 - dw / std::max(hmax_[c], 1e-30))
                                                  * (0.25 - dw / std::max(hmax_[c], 1e-30))),
                                   1.0);
        const double fe1 = (dw / std::max(hmax_[c], 1e-30) < 0.5)
                               ? 2.0
                                     * std::exp(-11.09 * (0.45 - dw / std::max(hmax_[c], 1e-30))
                                                * (0.45 - dw / std::max(hmax_[c], 1e-30)))
                               : 2.0
                                     * std::exp(-9.0 * (0.65 - dw / std::max(hmax_[c], 1e-30))
                                                * (0.65 - dw / std::max(hmax_[c], 1e-30)));
        const double ft = std::tanh(std::pow(C_t * C_t * rdt, 3));
        const double fl = std::tanh(std::pow(C_l * C_l * rdt, 10));
        const double fe2 = 1.0 - std::max(ft, fl);
        const double psi = 1.0; // ψ low-Reynolds correction (=1 in high-Re)
        const double fe = std::max((fe1 - 1.0), 0.0) * psi * fe2;

        const double fdt_tilde = std::max((1.0 - fdt), fd);

        // IDDES length scale.
        const double Delta_IDDES =
            std::min(std::max(std::max(C_w * dw, C_w * hmax_[c]), delta_[c]), hmax_[c]);
        const double l_LES = C_DES * Delta_IDDES;
        const double l_RANS = dw;
        const double l_iddes = fdt_tilde * (1.0 + fe) * l_RANS + (1.0 - fdt_tilde) * l_LES;

        const double d2 = l_iddes * l_iddes + 1e-30;
        // SA destruction term uses l_iddes instead of d_wall.
        const double Stilde = std::max(Smag + nt_c * fv2 / (kappa * kappa * d2), 0.3 * Smag);
        const double r = std::min(nt_c / (std::max(Stilde, 1e-30) * kappa * kappa * d2), 10.0);
        const double g = r + cw2 * (std::pow(r, 6) - r);
        const double fw =
            g * std::pow((1.0 + std::pow(cw3, 6)) / (std::pow(g, 6) + std::pow(cw3, 6)), 1.0 / 6.0);
        const double P = cb1 * rho_ * Stilde * nt_c;
        const double D = cw1 * rho_ * fw * (nt_c * nt_c) / d2;
        src[c] = P - D;
    }

    nuTildeEq_->solve_iteration();
    for (auto& v : nt)
        v = std::max(v, 0.0);
}

} // namespace simall::turbulence
