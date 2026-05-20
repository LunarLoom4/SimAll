// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/SurfaceTensionCsf.cpp
//
// Implementation: two least-squares gradients on the same field — the
// first delivers ∇α from which we form the unit normal n̂ = ∇α / |∇α|,
// the second computes the divergence ∇·n̂ to recover the curvature κ.
// Wall-adhesion tilt is applied by overwriting n̂ at boundary cells from
// Young's relation n̂ = n_wall cos θ + t_wall sin θ where t_wall is the
// projection of the fluid normal on the wall plane.
// =============================================================================
#include "multiphase/SurfaceTensionCsf.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase
{

void SurfaceTensionCsf::initialize(const meshing::Mesh& m, solver::FieldRegistry& F, CsfProps props)
{
    mesh_ = &m;
    p_ = std::move(props);
    const std::size_t nC = m.cells().size();
    F.vector("S_SurfaceT", nC);
    F.scalar("kappa", nC);
}

void SurfaceTensionCsf::compute(solver::FieldRegistry& F)
{
    if (!mesh_)
        return;
    const auto* alpha = F.find_scalar(p_.alphaField);
    if (!alpha)
        return;
    const std::size_t nC = mesh_->cells().size();

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gAlpha;
    G.evaluate(*alpha, gAlpha);

    // Unit normal field nHat = ∇α / |∇α|.
    util::aligned_vector<double> nx(nC), ny(nC), nz(nC);
    const double eps = 1.0e-8;
    for (std::size_t c = 0; c < nC; ++c) {
        const double mag = std::sqrt(gAlpha.x[c] * gAlpha.x[c] + gAlpha.y[c] * gAlpha.y[c]
                                     + gAlpha.z[c] * gAlpha.z[c]);
        if (mag > eps) {
            nx[c] = gAlpha.x[c] / mag;
            ny[c] = gAlpha.y[c] / mag;
            nz[c] = gAlpha.z[c] / mag;
        }
    }

    // Apply wall-adhesion BC for cells adjacent to boundary faces with a
    // specified contact angle. We tilt the cell's nHat in-place.
    if (!p_.contactAngles.empty()) {
        const auto& C = mesh_->cells();
        const auto& Ff = mesh_->faces();
        for (std::size_t f = 0; f < Ff.size(); ++f) {
            if (Ff.neighbor[f] != meshing::kBoundaryCell)
                continue;
            auto it = p_.contactAngles.find(Ff.boundaryZone[f]);
            if (it == p_.contactAngles.end())
                continue;
            const double theta = it->second;
            const meshing::CellId c = Ff.owner[f];
            double wx = Ff.areaX[f], wy = Ff.areaY[f], wz = Ff.areaZ[f];
            const double wmag = std::sqrt(wx * wx + wy * wy + wz * wz);
            if (wmag < eps)
                continue;
            wx /= wmag;
            wy /= wmag;
            wz /= wmag;
            // Wall outward from cell.
            const double sx = Ff.centroidX[f] - C.centroidX[c];
            const double sy = Ff.centroidY[f] - C.centroidY[c];
            const double sz = Ff.centroidZ[f] - C.centroidZ[c];
            if (wx * sx + wy * sy + wz * sz < 0) {
                wx = -wx;
                wy = -wy;
                wz = -wz;
            }
            // Tangent component of current normal.
            const double ndotw = nx[c] * wx + ny[c] * wy + nz[c] * wz;
            double tx = nx[c] - ndotw * wx;
            double ty = ny[c] - ndotw * wy;
            double tz = nz[c] - ndotw * wz;
            const double tmag = std::sqrt(tx * tx + ty * ty + tz * tz);
            if (tmag > eps) {
                tx /= tmag;
                ty /= tmag;
                tz /= tmag;
            }
            const double cs = std::cos(theta), sn = std::sin(theta);
            nx[c] = wx * cs + tx * sn;
            ny[c] = wy * cs + ty * sn;
            nz[c] = wz * cs + tz * sn;
        }
    }

    // Curvature κ = -∇·n̂. Compute the three component divergences via
    // independent gradients.
    solver::ScalarField scx(nx.begin(), nx.end());
    solver::ScalarField scy(ny.begin(), ny.end());
    solver::ScalarField scz(nz.begin(), nz.end());
    solver::VectorField gnx, gny, gnz;
    G.evaluate(scx, gnx);
    G.evaluate(scy, gny);
    G.evaluate(scz, gnz);

    auto& Sst = *F.find_vector("S_SurfaceT");
    auto& kappa = *F.find_scalar("kappa");
    for (std::size_t c = 0; c < nC; ++c) {
        const double divN = gnx.x[c] + gny.y[c] + gnz.z[c];
        const double k = -divN;
        kappa[c] = k;
        Sst.x[c] = p_.sigma * k * gAlpha.x[c];
        Sst.y[c] = p_.sigma * k * gAlpha.y[c];
        Sst.z[c] = p_.sigma * k * gAlpha.z[c];
    }
}

} // namespace simall::multiphase
