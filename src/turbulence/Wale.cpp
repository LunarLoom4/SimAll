// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/Wale.cpp
// =============================================================================
#include "turbulence/Wale.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence
{

void WALE_LES::initialize(meshing::Mesh& m, solver::FieldRegistry& f)
{
    mesh_ = &m;
    const std::size_t nC = m.cells().size();
    f.scalar("mut", nC);
    mut_.assign(nC, 0.0);
    delta_.assign(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c)
        delta_[c] = std::cbrt(std::max(m.cells().volume[c], 1e-30));
    SIMALL_LOG_INFO("Turbulence", "WALE LES initialized (Cw=", Cw_, ", cells=", nC, ")");
}

void WALE_LES::solve(double /*dt*/, solver::FieldRegistry& f)
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U = *f.find_vector("U");
    auto& mut = *f.find_scalar("mut");

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx);
    G.evaluate(U.y, gUy);
    G.evaluate(U.z, gUz);

    for (std::size_t c = 0; c < nC; ++c) {
        // g[i][j] = ∂u_i/∂x_j
        const double g[3][3] = {{gUx.x[c], gUx.y[c], gUx.z[c]},
                                {gUy.x[c], gUy.y[c], gUy.z[c]},
                                {gUz.x[c], gUz.y[c], gUz.z[c]}};
        // S_ij = ½ (g_ij + g_ji)
        double S[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                S[i][j] = 0.5 * (g[i][j] + g[j][i]);
        double SS = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                SS += S[i][j] * S[i][j];

        // g² = g·g  → (g²)_ij = g_ik g_kj
        double g2[3][3] = {{0}};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    g2[i][j] += g[i][k] * g[k][j];
        // S^d_ij = ½ ((g²)_ij + (g²)_ji) - ⅓ δ_ij (g²)_kk
        const double trg2 = g2[0][0] + g2[1][1] + g2[2][2];
        double Sd[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Sd[i][j] = 0.5 * (g2[i][j] + g2[j][i]) - (i == j ? trg2 / 3.0 : 0.0);
        double SdSd = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                SdSd += Sd[i][j] * Sd[i][j];

        const double num = std::pow(std::max(SdSd, 0.0), 1.5);
        const double den =
            std::pow(std::max(SS, 0.0), 2.5) + std::pow(std::max(SdSd, 0.0), 1.25) + 1e-30;
        const double Ls = Cw_ * delta_[c];
        mut[c] = rho_ * Ls * Ls * num / den;
        mut_[c] = mut[c];
    }
}

} // namespace simall::turbulence
