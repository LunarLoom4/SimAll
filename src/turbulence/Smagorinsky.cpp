// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/Smagorinsky.cpp
// =============================================================================
#include "turbulence/Smagorinsky.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence {

void Smagorinsky_LES::initialize(meshing::Mesh& m, solver::FieldRegistry& f) {
    mesh_ = &m;
    const std::size_t nC = m.cells().size();
    f.scalar("mut", nC);
    if (vanDriest_) {
        f.scalar("wallDistance", nC);
        wallDist_ = solver::make_wall_distance_exact();
        wallDist_->compute(m, bcs_, *f.find_scalar("wallDistance"));
    }
    mut_.assign(nC, 0.0);
    delta_.assign(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c)
        delta_[c] = std::cbrt(std::max(m.cells().volume[c], 1e-30));
    SIMALL_LOG_INFO("Turbulence", "Smagorinsky LES initialized (Cs=", Cs_,
        ", van-Driest=", vanDriest_, ", cells=", nC, ")");
}

void Smagorinsky_LES::solve(double /*dt*/, solver::FieldRegistry& f) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U = *f.find_vector("U");
    auto& mut    = *f.find_scalar("mut");

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx);
    G.evaluate(U.y, gUy);
    G.evaluate(U.z, gUz);

    const double* d = vanDriest_ ? f.find_scalar("wallDistance")->data() : nullptr;
    const double nu = mu_ / rho_;
    const double Aplus = 26.0;

    for (std::size_t c = 0; c < nC; ++c) {
        const double S11 = gUx.x[c];
        const double S22 = gUy.y[c];
        const double S33 = gUz.z[c];
        const double S12 = 0.5 * (gUx.y[c] + gUy.x[c]);
        const double S13 = 0.5 * (gUx.z[c] + gUz.x[c]);
        const double S23 = 0.5 * (gUy.z[c] + gUz.y[c]);
        const double SS  = 2.0 * (S11*S11 + S22*S22 + S33*S33
                                + 2.0 * (S12*S12 + S13*S13 + S23*S23));
        const double Smag = std::sqrt(std::max(SS, 0.0));

        double damp = 1.0;
        if (d && d[c] > 0.0) {
            // y+ from local wall-shear estimate: u_τ ≈ sqrt(ν |S|), y+ = y u_τ / ν.
            const double utau = std::sqrt(nu * Smag);
            const double yplus = d[c] * utau / std::max(nu, 1e-30);
            damp = 1.0 - std::exp(-yplus / Aplus);
        }
        const double Ls = Cs_ * delta_[c] * damp;
        mut[c] = rho_ * Ls * Ls * Smag;
        mut_[c] = mut[c];
    }
}

}  // namespace simall::turbulence
