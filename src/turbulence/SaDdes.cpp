// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/SaDdes.cpp
// =============================================================================
#include "turbulence/SaDdes.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence
{

namespace
{
// Spalart-Allmaras model constants.
constexpr double cb1 = 0.1355;
constexpr double cb2 = 0.622;
constexpr double sigma = 2.0 / 3.0;
constexpr double kappa = 0.41;
constexpr double cw1 = cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
constexpr double cw2 = 0.3;
constexpr double cw3 = 2.0;
constexpr double cv1 = 7.1;
} // namespace

void SaDdes_Full::initialize(meshing::Mesh& m, solver::FieldRegistry& F)
{
    mesh_ = &m;
    fields_ = &F;
    const std::size_t nC = m.cells().size();
    F.scalar("nuTilde", nC);
    F.scalar("mu_t", nC);
    F.scalar("wallDistance", nC);
    F.scalar("ddes_dHybrid", nC);
    F.scalar("nuTilde_src", nC);

    dHybrid_.assign(nC, 0.0);
    delta_.assign(nC, 0.0);
    mut_.assign(nC, 0.0);
    src_.assign(nC, 0.0);
    Smag_.assign(nC, 0.0);

    // Δ = V^{1/3}.
    const auto& C = m.cells();
    for (std::size_t c = 0; c < nC; ++c)
        delta_[c] = std::cbrt(std::max(1e-30, C.volume[c]));

    solver::LinearSolverConfig cfg;
    cfg.maxIterations = 200;
    cfg.tolerance = 1e-8;
    lin_ = solver::make_bicgstab(cfg);
    wallDist_ = solver::make_wall_distance_exact();
    wallDist_->compute(m, bcs_, *F.find_scalar("wallDistance"));

    nuTildeEq_ = std::make_unique<solver::ScalarTransport>(m, F, *lin_);
    nuTildeEq_->set_field("nuTilde");
    nuTildeEq_->set_density(rho_);
    nuTildeEq_->set_diffusivity(mu_ / sigma); // base; explicit src handles the rest
    nuTildeEq_->set_source("nuTilde_src");
    for (const auto& bc : bcs_) {
        solver::ScalarBC sb;
        sb.zone = bc.zone;
        sb.kind = solver::ScalarBC::Kind::Dirichlet;
        sb.value = 0.0;
        nuTildeEq_->add_bc(sb);
    }
}

void SaDdes_Full::compute_length_scales(const solver::FieldRegistry& F)
{
    const auto* U = F.find_vector("U");
    const auto* nuT = F.find_scalar("nuTilde");
    const auto* d = F.find_scalar("wallDistance");
    if (!mesh_ || !U || !nuT || !d)
        return;
    const std::size_t nC = mesh_->cells().size();
    const double nu = mu_ / rho_;

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gU, gV, gW;
    G.evaluate(U->x, gU);
    G.evaluate(U->y, gV);
    G.evaluate(U->z, gW);

    auto& dh = *const_cast<solver::FieldRegistry&>(F).find_scalar("ddes_dHybrid");
    for (std::size_t c = 0; c < nC; ++c) {
        const double Sxx = gU.x[c], Syy = gV.y[c], Szz = gW.z[c];
        const double Sxy = 0.5 * (gU.y[c] + gV.x[c]);
        const double Sxz = 0.5 * (gU.z[c] + gW.x[c]);
        const double Syz = 0.5 * (gV.z[c] + gW.y[c]);
        Smag_[c] = std::sqrt(2.0 * (Sxx * Sxx + Syy * Syy + Szz * Szz)
                             + 4.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz));
        // r_d = (ν_t + ν) / (sqrt(U_{i,j} U_{i,j}) κ² d²)
        const double UijUij =
            std::sqrt(gU.x[c] * gU.x[c] + gU.y[c] * gU.y[c] + gU.z[c] * gU.z[c] + gV.x[c] * gV.x[c]
                      + gV.y[c] * gV.y[c] + gV.z[c] * gV.z[c] + gW.x[c] * gW.x[c]
                      + gW.y[c] * gW.y[c] + gW.z[c] * gW.z[c]);
        const double dwall = std::max((*d)[c], 1e-12);
        const double nu_t = (mut_[c] > 0) ? mut_[c] / rho_ : 0.0;
        const double rd = (nu_t + nu) / (std::max(UijUij, 1e-10) * kappa * kappa * dwall * dwall);
        const double fd = 1.0 - std::tanh(std::pow(8.0 * rd, 3));
        const double dRANS = (*d)[c];
        const double dLES = C_DES_ * delta_[c];
        dHybrid_[c] = dRANS - fd * std::max(0.0, dRANS - dLES);
        dh[c] = dHybrid_[c];
    }
}

void SaDdes_Full::compute_sources(const solver::FieldRegistry& F)
{
    const auto* nuT = F.find_scalar("nuTilde");
    if (!nuT)
        return;
    const std::size_t nC = mesh_->cells().size();
    const double nu = mu_ / rho_;
    auto& srcF = *const_cast<solver::FieldRegistry&>(F).find_scalar("nuTilde_src");

    for (std::size_t c = 0; c < nC; ++c) {
        const double nt = std::max(0.0, (*nuT)[c]);
        const double chi = nt / std::max(nu, 1e-30);
        const double chi3 = chi * chi * chi;
        const double fv1 = chi3 / (chi3 + cv1 * cv1 * cv1);
        const double fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const double d2 = dHybrid_[c] * dHybrid_[c] + 1e-30;
        const double S = Smag_[c];
        const double Stilde = S + nt * fv2 / (kappa * kappa * d2);
        const double r = std::min(nt / (std::max(Stilde, 1e-30) * kappa * kappa * d2), 10.0);
        const double g = r + cw2 * (std::pow(r, 6) - r);
        const double fw =
            g * std::pow((1.0 + std::pow(cw3, 6)) / (std::pow(g, 6) + std::pow(cw3, 6)), 1.0 / 6.0);
        const double P = cb1 * rho_ * Stilde * nt;
        const double D = cw1 * rho_ * fw * (nt * nt) / d2;
        srcF[c] = P - D;
        mut_[c] = rho_ * nt * fv1;
        src_[c] = srcF[c];
    }
    auto& mt = *const_cast<solver::FieldRegistry&>(F).find_scalar("mu_t");
    std::copy(mut_.begin(), mut_.end(), mt.begin());
}

void SaDdes_Full::solve(double dt, solver::FieldRegistry& F)
{
    (void) dt;
    compute_length_scales(F);
    compute_sources(F);
    nuTildeEq_->solve_iteration();
    // Clamp non-negative.
    auto& nt = *F.find_scalar("nuTilde");
    for (auto& v : nt)
        v = std::max(v, 0.0);
}

} // namespace simall::turbulence
