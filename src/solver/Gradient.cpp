// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/Gradient.cpp
// =============================================================================
#include "solver/Gradient.hpp"

#include "core/Logger.hpp"

#include <cmath>

namespace simall::solver
{

namespace
{

bool invert3x3_sym(double xx, double xy, double xz, double yy, double yz, double zz, double out[6])
{
    const double det =
        xx * (yy * zz - yz * yz) - xy * (xy * zz - yz * xz) + xz * (xy * yz - yy * xz);
    if (std::abs(det) < 1e-30) {
        out[0] = out[1] = out[2] = out[3] = out[4] = out[5] = 0.0;
        return false;
    }
    const double inv = 1.0 / det;
    out[0] = (yy * zz - yz * yz) * inv;  // xx
    out[1] = -(xy * zz - yz * xz) * inv; // xy
    out[2] = (xy * yz - yy * xz) * inv;  // xz
    out[3] = (xx * zz - xz * xz) * inv;  // yy
    out[4] = -(xx * yz - xy * xz) * inv; // yz
    out[5] = (xx * yy - xy * xy) * inv;  // zz
    return true;
}

} // namespace

LeastSquaresGradient::LeastSquaresGradient(const meshing::Mesh& m) : mesh_(m)
{
    precompute();
}

void LeastSquaresGradient::precompute()
{
    const auto& C = mesh_.cells();
    const auto& F = mesh_.faces();
    const std::size_t nC = C.size();
    mInv_.assign(nC * 6, 0.0);

    util::aligned_vector<double> Axx(nC, 0), Axy(nC, 0), Axz(nC, 0), Ayy(nC, 0), Ayz(nC, 0),
        Azz(nC, 0);

    auto accumulate = [&](meshing::CellId c, double dx, double dy, double dz) {
        const double w2 = 1.0 / (dx * dx + dy * dy + dz * dz + 1e-30);
        Axx[c] += w2 * dx * dx;
        Axy[c] += w2 * dx * dy;
        Axz[c] += w2 * dx * dz;
        Ayy[c] += w2 * dy * dy;
        Ayz[c] += w2 * dy * dz;
        Azz[c] += w2 * dz * dz;
    };

    for (std::size_t f = 0; f < F.size(); ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        if (o == meshing::kBoundaryCell)
            continue;
        if (n == meshing::kBoundaryCell) {
            // Boundary face: treat face centroid as virtual neighbour.
            const double dx = F.centroidX[f] - C.centroidX[o];
            const double dy = F.centroidY[f] - C.centroidY[o];
            const double dz = F.centroidZ[f] - C.centroidZ[o];
            accumulate(o, dx, dy, dz);
        } else {
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            accumulate(o, dx, dy, dz);
            accumulate(n, -dx, -dy, -dz);
        }
    }

    for (std::size_t c = 0; c < nC; ++c) {
        double inv[6];
        invert3x3_sym(Axx[c], Axy[c], Axz[c], Ayy[c], Ayz[c], Azz[c], inv);
        for (int k = 0; k < 6; ++k)
            mInv_[c * 6 + k] = inv[k];
    }
}

void LeastSquaresGradient::evaluate(const ScalarField& phi, VectorField& grad) const
{
    const auto& C = mesh_.cells();
    const auto& F = mesh_.faces();
    const std::size_t nC = C.size();
    grad.resize(nC);
    util::aligned_vector<double> bx(nC, 0), by(nC, 0), bz(nC, 0);

    for (std::size_t f = 0; f < F.size(); ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        if (o == meshing::kBoundaryCell)
            continue;
        if (n == meshing::kBoundaryCell) {
            // Use cell value at the boundary (zero-gradient default; Dirichlet
            // closures apply the BC value during assembly).
            const double dx = F.centroidX[f] - C.centroidX[o];
            const double dy = F.centroidY[f] - C.centroidY[o];
            const double dz = F.centroidZ[f] - C.centroidZ[o];
            const double w2 = 1.0 / (dx * dx + dy * dy + dz * dz + 1e-30);
            const double dphi = 0.0;
            bx[o] += w2 * dx * dphi;
            by[o] += w2 * dy * dphi;
            bz[o] += w2 * dz * dphi;
        } else {
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            const double w2 = 1.0 / (dx * dx + dy * dy + dz * dz + 1e-30);
            const double dphi = phi[n] - phi[o];
            bx[o] += w2 * dx * dphi;
            by[o] += w2 * dy * dphi;
            bz[o] += w2 * dz * dphi;
            bx[n] += w2 * dx * -dphi;
            by[n] += w2 * dy * -dphi;
            bz[n] += w2 * dz * -dphi;
        }
    }

    for (std::size_t c = 0; c < nC; ++c) {
        const double m0 = mInv_[c * 6 + 0], m1 = mInv_[c * 6 + 1], m2 = mInv_[c * 6 + 2];
        const double m3 = mInv_[c * 6 + 3], m4 = mInv_[c * 6 + 4], m5 = mInv_[c * 6 + 5];
        grad.x[c] = m0 * bx[c] + m1 * by[c] + m2 * bz[c];
        grad.y[c] = m1 * bx[c] + m3 * by[c] + m4 * bz[c];
        grad.z[c] = m2 * bx[c] + m4 * by[c] + m5 * bz[c];
    }
}

} // namespace simall::solver
