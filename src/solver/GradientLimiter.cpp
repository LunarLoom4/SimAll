// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/GradientLimiter.cpp
// =============================================================================
#include "solver/GradientLimiter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace simall::solver {

GradientLimiter::GradientLimiter(const meshing::Mesh& m) : mesh_(m) {}

namespace {
/// Barth-Jespersen scalar limit (min over faces).
inline double bj_face(double phi_c, double phi_min, double phi_max, double dPhi) {
    if (dPhi > 1e-30) return std::min(1.0, (phi_max - phi_c) / dPhi);
    if (dPhi < -1e-30) return std::min(1.0, (phi_min - phi_c) / dPhi);
    return 1.0;
}
/// Venkatakrishnan smooth limiter, ε² = (K h)³ with h ~ cell length scale.
inline double venkat_face(double dM, double dPhi, double eps2) {
    // dM = either (phi_max - phi_c) or (phi_min - phi_c) depending on sign of dPhi
    const double a = dM * dM + eps2;
    const double b = 2.0 * dPhi * dPhi + dPhi * dM + eps2;
    return (a + 2.0 * dPhi * dM) / std::max(b, 1e-30);
}
}  // namespace

void GradientLimiter::evaluate(const ScalarField& phi, const VectorField& g,
                               ScalarField& psi, LimiterKind kind,
                               double K) const {
    const auto& C = mesh_.cells();
    const auto& F = mesh_.faces();
    const std::size_t nC = C.size();
    if (psi.size() != nC) psi.assign(nC, 1.0);
    if (kind == LimiterKind::None) {
        std::fill(psi.begin(), psi.end(), 1.0);
        return;
    }

    // Per-cell neighbourhood extrema.
    util::aligned_vector<double> phiMin(nC, 0.0), phiMax(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c) {
        phiMin[c] = phi[c]; phiMax[c] = phi[c];
        const int fs = C.faceOffsets[c], fe = C.faceOffsets[c+1];
        for (int k = fs; k < fe; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            const meshing::CellId nb = (F.owner[fid] == c) ? F.neighbor[fid]
                                                           : F.owner[fid];
            if (nb == meshing::kBoundaryCell) continue;
            phiMin[c] = std::min(phiMin[c], phi[nb]);
            phiMax[c] = std::max(phiMax[c], phi[nb]);
        }
    }

    for (std::size_t c = 0; c < nC; ++c) {
        double psi_c = 1.0;
        const int fs = C.faceOffsets[c], fe = C.faceOffsets[c+1];
        // Cell characteristic length h = V^{1/3}.
        const double h = std::cbrt(std::max(C.volume[c], 1e-30));
        const double eps2 = std::pow(K * h, 3.0);
        for (int k = fs; k < fe; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            const double dx = F.centroidX[fid] - C.centroidX[c];
            const double dy = F.centroidY[fid] - C.centroidY[c];
            const double dz = F.centroidZ[fid] - C.centroidZ[c];
            const double dPhi = g.x[c]*dx + g.y[c]*dy + g.z[c]*dz;
            double psi_f = 1.0;
            if (kind == LimiterKind::BarthJespersen) {
                psi_f = bj_face(phi[c], phiMin[c], phiMax[c], dPhi);
            } else {  // Venkatakrishnan
                const double dM = (dPhi >= 0) ? (phiMax[c] - phi[c])
                                              : (phiMin[c] - phi[c]);
                psi_f = venkat_face(dM, dPhi, eps2);
            }
            psi_c = std::min(psi_c, std::max(0.0, psi_f));
        }
        psi[c] = psi_c;
    }
}

}  // namespace simall::solver
