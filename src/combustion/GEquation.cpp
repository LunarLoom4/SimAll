// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/GEquation.cpp
// =============================================================================
#include "combustion/GEquation.hpp"

#include "core/Logger.hpp"
#include "solver/LeastSquaresGradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::combustion
{

void GEquation::initialize(const meshing::Mesh& mesh,
                           solver::FieldRegistry& fields,
                           GEquationParams params)
{
    mesh_ = &mesh;
    fields_ = &fields;
    p_ = params;
    const std::size_t nC = mesh.cells().size();
    fields.scalar("G", nC);
    fields.scalar("u_rms", nC);
    fields.scalar("progress", nC);
    fields.vector("gradG", nC);
    SIMALL_LOG_INFO("Combustion",
                    "G-equation initialised: s_L=",
                    p_.s_laminar,
                    " C_DA=",
                    p_.C_damkohler,
                    " (",
                    nC,
                    " cells)");
}

double GEquation::cell_size(meshing::CellId c) const
{
    return std::cbrt(std::max(mesh_->cells().volume[c], 1e-30));
}

double GEquation::step(double dt)
{
    if (!mesh_ || !fields_)
        return 0.0;
    const std::size_t nC = mesh_->cells().size();
    auto& G = *fields_->find_scalar("G");
    const auto& urms = *fields_->find_scalar("u_rms");
    const auto& U = *fields_->find_vector("U");
    auto& gG = *fields_->find_vector("gradG");

    // ∇G via least-squares.
    solver::LeastSquaresGradient gradOp(*mesh_);
    gradOp.evaluate(G, gG);

    double l2 = 0.0;
    for (std::size_t c = 0; c < nC; ++c) {
        const double gx = gG.x[c], gy = gG.y[c], gz = gG.z[c];
        const double magG = std::sqrt(gx * gx + gy * gy + gz * gz);
        const double udotg = U.x[c] * gx + U.y[c] * gy + U.z[c] * gz;
        // Damköhler turbulent flame speed.
        const double sT =
            p_.s_laminar * (1.0 + p_.C_damkohler * urms[c] / std::max(p_.s_laminar, 1e-30));
        const double dG = dt * (sT * magG - udotg);
        l2 += dG * dG;
        G[c] += dG;
    }
    l2 = std::sqrt(l2 / std::max<std::size_t>(nC, 1));

    if (p_.reinitInterval > 0 && (++stepCounter_ % p_.reinitInterval) == 0) {
        reinitialise();
    }
    update_progress_indicator();
    return l2;
}

void GEquation::reinitialise()
{
    // Sussman et al. 1994 pseudo-time PDE:
    //   ∂G/∂τ = sign(G_0)(1 - |∇G|)
    // Discretised explicitly with Godunov upwinding (here approximated by
    // least-squares gradient magnitude on the unstructured mesh).
    if (!mesh_ || !fields_)
        return;
    const std::size_t nC = mesh_->cells().size();
    auto& G = *fields_->find_scalar("G");
    auto& gG = *fields_->find_vector("gradG");
    solver::LeastSquaresGradient gradOp(*mesh_);

    util::aligned_vector<double> G0(G.begin(), G.end());
    for (int it = 0; it < p_.reinitSubSteps; ++it) {
        gradOp.evaluate(G, gG);
        for (std::size_t c = 0; c < nC; ++c) {
            const double h = cell_size(static_cast<meshing::CellId>(c));
            const double dtau = p_.reinitDtFactor * h;
            const double mag = std::sqrt(gG.x[c] * gG.x[c] + gG.y[c] * gG.y[c] + gG.z[c] * gG.z[c]);
            const double sg = G0[c] / std::sqrt(G0[c] * G0[c] + h * h);
            G[c] += dtau * sg * (1.0 - mag);
        }
    }
}

void GEquation::update_progress_indicator()
{
    if (!mesh_ || !fields_)
        return;
    const std::size_t nC = mesh_->cells().size();
    const auto& G = *fields_->find_scalar("G");
    auto& Pr = *fields_->find_scalar("progress");
    for (std::size_t c = 0; c < nC; ++c) {
        const double eps = cell_size(static_cast<meshing::CellId>(c));
        const double d = (p_.G_iso - G[c]) / std::max(eps, 1e-30);
        // Smoothed Heaviside (Tornberg-Engquist 2004 5-point):
        Pr[c] = 0.5 * (1.0 + std::tanh(d));
    }
}

} // namespace simall::combustion
