// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/PeriodicPairing.cpp
// =============================================================================
#include "solver/PeriodicPairing.hpp"
#include "core/Logger.hpp"

#include <cmath>
#include <limits>

namespace simall::solver {

namespace {

util::Vec3d apply_transform(const util::Vec3d& x, const PeriodicTransform& T) {
    // Normalise axis.
    const double an = std::sqrt(T.rotationAxis.x*T.rotationAxis.x
                              + T.rotationAxis.y*T.rotationAxis.y
                              + T.rotationAxis.z*T.rotationAxis.z);
    const double inv = an > 0 ? 1.0 / an : 0.0;
    const double ux = T.rotationAxis.x * inv;
    const double uy = T.rotationAxis.y * inv;
    const double uz = T.rotationAxis.z * inv;
    const double ca = std::cos(T.rotationAngleRad);
    const double sa = std::sin(T.rotationAngleRad);
    const double oc = 1.0 - ca;

    // Rodrigues' rotation about axis through origin then translate centre back.
    const double rx = x.x - T.rotationCentre.x;
    const double ry = x.y - T.rotationCentre.y;
    const double rz = x.z - T.rotationCentre.z;
    const double R11 = ca + ux*ux*oc;
    const double R12 = ux*uy*oc - uz*sa;
    const double R13 = ux*uz*oc + uy*sa;
    const double R21 = uy*ux*oc + uz*sa;
    const double R22 = ca + uy*uy*oc;
    const double R23 = uy*uz*oc - ux*sa;
    const double R31 = uz*ux*oc - uy*sa;
    const double R32 = uz*uy*oc + ux*sa;
    const double R33 = ca + uz*uz*oc;
    return {
        R11*rx + R12*ry + R13*rz + T.rotationCentre.x + T.translation.x,
        R21*rx + R22*ry + R23*rz + T.rotationCentre.y + T.translation.y,
        R31*rx + R32*ry + R33*rz + T.rotationCentre.z + T.translation.z
    };
}

}  // namespace

std::vector<int> build_periodic_pairs(const meshing::Mesh& mesh,
                                      const std::vector<PeriodicTransform>& Ts) {
    const auto& F = mesh.faces();
    std::vector<int> twin(F.size(), -1);

    for (const auto& T : Ts) {
        std::vector<int> A, B;
        for (std::size_t f = 0; f < F.size(); ++f) {
            if (F.neighbor[f] != meshing::kBoundaryCell) continue;
            if (F.boundaryZone[f] == T.zoneA) A.push_back(static_cast<int>(f));
            else if (F.boundaryZone[f] == T.zoneB) B.push_back(static_cast<int>(f));
        }
        if (A.empty() || B.empty()) {
            SIMALL_LOG_WARN("PeriodicPairing", "zone pair (",
                T.zoneA, ",", T.zoneB, ") empty on at least one side");
            continue;
        }
        int matched = 0;
        for (int fa : A) {
            const util::Vec3d xa{ F.centroidX[fa], F.centroidY[fa], F.centroidZ[fa] };
            const util::Vec3d xt = apply_transform(xa, T);
            int best = -1; double bestd2 = std::numeric_limits<double>::max();
            for (int fb : B) {
                const double dx = F.centroidX[fb] - xt.x;
                const double dy = F.centroidY[fb] - xt.y;
                const double dz = F.centroidZ[fb] - xt.z;
                const double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 < bestd2) { bestd2 = d2; best = fb; }
            }
            if (best >= 0 && std::sqrt(bestd2) <= T.tolerance) {
                twin[fa] = best;
                twin[best] = fa;
                ++matched;
            }
        }
        SIMALL_LOG_INFO("PeriodicPairing", "zone pair (", T.zoneA, ",", T.zoneB,
            "): matched ", matched, " / ", A.size(), " faces");
    }
    return twin;
}

}  // namespace simall::solver
