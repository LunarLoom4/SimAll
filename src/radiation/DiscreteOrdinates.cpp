// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/DiscreteOrdinates.cpp
// =============================================================================
#include "radiation/DiscreteOrdinates.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::radiation
{

namespace
{
constexpr double SIGMA_SB = 5.670374419e-8; // Stefan-Boltzmann [W/(m²·K⁴)]
constexpr double FOUR_PI = 12.566370614359172;
}

DiscreteOrdinates::DiscreteOrdinates(meshing::Mesh& m, solver::FieldRegistry& f, QuadratureOrder o)
    : mesh_(m), F_(f), order_(o)
{
    build_quadrature();
    const std::size_t nC = m.cells().size();
    I_.assign(omega_.size(), util::aligned_vector<double>(nC, 0.0));
    F_.scalar("G_rad", nC);
}

void DiscreteOrdinates::build_quadrature()
{
    // Carlson Level-Symmetric S_N (Lewis & Miller 1984 Table 4-1) — direction
    // cosines lie on the surface of a unit sphere with octahedral symmetry.
    // Per-octant the ordinates use the cosines below; weights w_m sum to 1 per
    // octant and total to 4π (we normalise by 4π at consumption time).
    omega_.clear();
    std::vector<std::array<double, 4>> oct; // (μ_x, μ_y, μ_z, w_octant)
    if (order_ == QuadratureOrder::S4) {
        // S4: 3 directions per octant, weights all equal to 1/3 (octant).
        const double a = 0.2958759, b = 0.9082483;
        const double w = 1.0 / 3.0;
        oct = {{a, a, b, w}, {a, b, a, w}, {b, a, a, w}};
    } else { // S6
        const double a = 0.1838670, b = 0.6950514, c = 0.9656013;
        const double w1 = 0.1761263, w2 = 0.1572071;
        oct = {{a, a, c, w1},
               {a, c, a, w1},
               {c, a, a, w1},
               {a, b, b, w2},
               {b, a, b, w2},
               {b, b, a, w2}};
    }
    for (int sx : {+1, -1})
        for (int sy : {+1, -1})
            for (int sz : {+1, -1})
                for (const auto& q : oct)
                    omega_.push_back({sx * q[0], sy * q[1], sz * q[2], q[3] * (FOUR_PI / 8.0)});
}

const WallEmission* DiscreteOrdinates::find_wall(meshing::ZoneId z) const
{
    for (const auto& w : walls_)
        if (w.zone == z)
            return &w;
    return nullptr;
}

// Single ordinate FV sweep. For arbitrary polyhedra we use first-order
// upwind cell-centred discretisation: outgoing fluxes through faces with
// (s·n_f) > 0 are at the upwind cell value (I_c itself); incoming faces
// with (s·n_f) < 0 use the neighbour cell value (or BC). The implicit
// equation per cell is solved by Jacobi sweep within the source iteration.
//
//   sum_f (s·A_f)^+ I_c - sum_f (s·A_f)^- I_nb + (κ+σ_s) V_c I_c
//       = κ I_b V_c + (σ_s/4π) Σ_m' w_m' I_m' V_c
void DiscreteOrdinates::sweep_ordinate(int m, double& maxDelta)
{
    const auto& C = mesh_.cells();
    const auto& F = mesh_.faces();
    const std::size_t nC = C.size();
    const auto* T = F_.find_scalar("T");
    const double sx = omega_[m].sx, sy = omega_[m].sy, sz = omega_[m].sz;

    // Pre-compute scattering source from current intensities.
    util::aligned_vector<double> scat(nC, 0.0);
    if (sigmaS_ > 0.0) {
        for (std::size_t mp = 0; mp < omega_.size(); ++mp) {
            const double w = omega_[mp].w / FOUR_PI;
            for (std::size_t c = 0; c < nC; ++c)
                scat[c] += w * I_[mp][c];
        }
    }

    auto& Im = I_[m];
    for (std::size_t c = 0; c < nC; ++c) {
        const double Tc = T ? (*T)[c] : 300.0;
        const double Ib = SIGMA_SB * Tc * Tc * Tc * Tc / M_PI;
        double aP = (kappa_ + sigmaS_) * C.volume[c];
        double rhs = (kappa_ * Ib + sigmaS_ * scat[c]) * C.volume[c];
        const int beg = C.faceOffsets[c], end = C.faceOffsets[c + 1];
        for (int k = beg; k < end; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            double Ax = F.areaX[fid], Ay = F.areaY[fid], Az = F.areaZ[fid];
            // Ensure area vector points outward from cell c.
            if (F.owner[fid] != c) {
                Ax = -Ax;
                Ay = -Ay;
                Az = -Az;
            }
            const double sdotA = sx * Ax + sy * Ay + sz * Az;
            if (sdotA >= 0.0) {
                aP += sdotA; // outgoing — upwind is cell c
            } else {
                const meshing::CellId nb = (F.owner[fid] == c) ? F.neighbor[fid] : F.owner[fid];
                if (nb != meshing::kBoundaryCell) {
                    rhs -= sdotA * Im[nb]; // incoming neighbour
                } else {
                    // Wall BC: diffusely emitting.
                    const auto* w = find_wall(F.boundaryZone[fid]);
                    double Iw = 0.0;
                    if (w)
                        Iw = w->emissivity * SIGMA_SB * std::pow(w->temperature, 4.0) / M_PI;
                    rhs -= sdotA * Iw;
                }
            }
        }
        const double newI = (aP > 1e-30) ? rhs / aP : Im[c];
        maxDelta = std::max(maxDelta, std::abs(newI - Im[c]));
        Im[c] = newI;
    }
}

double DiscreteOrdinates::sweep()
{
    double maxDelta = 0.0;
    for (std::size_t m = 0; m < omega_.size(); ++m)
        sweep_ordinate(static_cast<int>(m), maxDelta);
    // Refresh mean intensity G = Σ w I.
    auto& G = *F_.find_scalar("G_rad");
    const std::size_t nC = G.size();
    std::fill(G.begin(), G.end(), 0.0);
    for (std::size_t m = 0; m < omega_.size(); ++m)
        for (std::size_t c = 0; c < nC; ++c)
            G[c] += omega_[m].w * I_[m][c];
    return maxDelta;
}

void DiscreteOrdinates::compute_source(util::aligned_vector<double>& Srad) const
{
    const auto& G = *const_cast<solver::FieldRegistry&>(F_).find_scalar("G_rad");
    const auto* T = const_cast<solver::FieldRegistry&>(F_).find_scalar("T");
    const std::size_t nC = G.size();
    Srad.assign(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c) {
        const double Tc = T ? (*T)[c] : 300.0;
        const double Ib = SIGMA_SB * Tc * Tc * Tc * Tc / M_PI;
        Srad[c] = kappa_ * (G[c] - 4.0 * M_PI * Ib);
    }
}

} // namespace simall::radiation
