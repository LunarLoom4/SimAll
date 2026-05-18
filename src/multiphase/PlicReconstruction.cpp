// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/PlicReconstruction.cpp
// =============================================================================
#include "multiphase/PlicReconstruction.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace simall::multiphase {

namespace {

inline double dot3(const util::Vec3d& a, const util::Vec3d& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline util::Vec3d cross3(const util::Vec3d& a, const util::Vec3d& b) {
    return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}

/// Signed volume of tetrahedron (a, b, c, d).
inline double tet_signed_volume(const util::Vec3d& a, const util::Vec3d& b,
                                const util::Vec3d& c, const util::Vec3d& d) {
    const util::Vec3d ab{b.x-a.x, b.y-a.y, b.z-a.z};
    const util::Vec3d ac{c.x-a.x, c.y-a.y, c.z-a.z};
    const util::Vec3d ad{d.x-a.x, d.y-a.y, d.z-a.z};
    return (1.0/6.0) * dot3(ab, cross3(ac, ad));
}

/// Cut a tetrahedron by half-space n·x < δ and return the volume of the
/// "inside" part. Handles all 16 cases by sign of the 4 vertices.
double tet_cut_volume(const std::array<util::Vec3d, 4>& v,
                      const util::Vec3d& n, double delta) {
    double s[4];
    int below[4], above[4]; int nb = 0, na = 0;
    for (int i = 0; i < 4; ++i) {
        s[i] = dot3(n, v[i]) - delta;
        if (s[i] <= 0) below[nb++] = i;
        else           above[na++] = i;
    }
    const double Vtot = std::abs(tet_signed_volume(v[0], v[1], v[2], v[3]));
    if (na == 0) return Vtot;
    if (nb == 0) return 0.0;

    auto cut = [&](int b, int a) {
        const double t = s[b] / (s[b] - s[a]);
        return util::Vec3d{ v[b].x + t*(v[a].x - v[b].x),
                            v[b].y + t*(v[a].y - v[b].y),
                            v[b].z + t*(v[a].z - v[b].z) };
    };

    if (nb == 1) {
        // Inside region = 1 small tetrahedron at the below vertex.
        const util::Vec3d& b = v[below[0]];
        const util::Vec3d p0 = cut(below[0], above[0]);
        const util::Vec3d p1 = cut(below[0], above[1]);
        const util::Vec3d p2 = cut(below[0], above[2]);
        return std::abs(tet_signed_volume(b, p0, p1, p2));
    }
    if (na == 1) {
        // Inside region = original tet minus 1 small tetrahedron at the above vertex.
        const util::Vec3d& a = v[above[0]];
        const util::Vec3d p0 = cut(below[0], above[0]);
        const util::Vec3d p1 = cut(below[1], above[0]);
        const util::Vec3d p2 = cut(below[2], above[0]);
        return Vtot - std::abs(tet_signed_volume(a, p0, p1, p2));
    }
    // nb == 2, na == 2 : inside region = prism. Decompose into 3 tetrahedra.
    const util::Vec3d& b0 = v[below[0]];
    const util::Vec3d& b1 = v[below[1]];
    const util::Vec3d p00 = cut(below[0], above[0]);
    const util::Vec3d p01 = cut(below[0], above[1]);
    const util::Vec3d p10 = cut(below[1], above[0]);
    const util::Vec3d p11 = cut(below[1], above[1]);
    const double V1 = std::abs(tet_signed_volume(b0, b1, p00, p10));
    const double V2 = std::abs(tet_signed_volume(b1, p00, p10, p11));
    const double V3 = std::abs(tet_signed_volume(b0, b1, p11, p01));
    return V1 + V2 + V3;
}

}  // namespace

double PlicReconstruction::cell_cut_volume(const meshing::Mesh& m,
                                           meshing::CellId c,
                                           const util::Vec3d& n,
                                           double delta) const {
    // Decompose cell into tetrahedra: for each face, fan from face's first
    // node through edges, then attach to cell centroid.
    const auto& C = m.cells();
    const auto& F = m.faces();
    const auto& N = m.nodes();
    const util::Vec3d xc{ C.centroidX[c], C.centroidY[c], C.centroidZ[c] };
    const int fs = C.faceOffsets[c], fe = C.faceOffsets[c+1];
    double V = 0.0;
    for (int k = fs; k < fe; ++k) {
        const meshing::FaceId fid = C.faceIndices[k];
        const int ns = F.nodeOffsets[fid], ne = F.nodeOffsets[fid+1];
        if (ne - ns < 3) continue;
        const auto firstId = F.nodeIndices[ns];
        const util::Vec3d a{ N.x[firstId], N.y[firstId], N.z[firstId] };
        for (int p = ns + 1; p + 1 < ne; ++p) {
            const auto bId = F.nodeIndices[p];
            const auto cId = F.nodeIndices[p + 1];
            const util::Vec3d b{ N.x[bId], N.y[bId], N.z[bId] };
            const util::Vec3d cc{ N.x[cId], N.y[cId], N.z[cId] };
            V += tet_cut_volume({xc, a, b, cc}, n, delta);
        }
    }
    return V;
}

void PlicReconstruction::reconstruct(const meshing::Mesh& m,
                                     const solver::ScalarField& alpha) {
    const std::size_t nC = m.cells().size();
    patches_.assign(nC, PlicPatch{});

    solver::LeastSquaresGradient G(m);
    solver::VectorField gA;
    G.evaluate(alpha, gA);

    const auto& C = m.cells();
    for (std::size_t c = 0; c < nC; ++c) {
        const double a = alpha[c];
        if (a <= 1e-6 || a >= 1.0 - 1e-6) continue;
        util::Vec3d n{ -gA.x[c], -gA.y[c], -gA.z[c] };
        const double mag = std::sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
        if (mag < 1e-30) continue;
        n.x /= mag; n.y /= mag; n.z /= mag;

        // Bisect δ on [δmin, δmax] derived from cell bbox in n direction.
        const util::Vec3d xc{ C.centroidX[c], C.centroidY[c], C.centroidZ[c] };
        const auto& F = m.faces(); const auto& N = m.nodes();
        const int fs = C.faceOffsets[c], fe = C.faceOffsets[c+1];
        double dmin =  std::numeric_limits<double>::max();
        double dmax = -std::numeric_limits<double>::max();
        for (int k = fs; k < fe; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            const int ns = F.nodeOffsets[fid], ne = F.nodeOffsets[fid+1];
            for (int p = ns; p < ne; ++p) {
                const auto v = F.nodeIndices[p];
                const double d = n.x*(N.x[v]-xc.x) + n.y*(N.y[v]-xc.y) + n.z*(N.z[v]-xc.z);
                dmin = std::min(dmin, d);
                dmax = std::max(dmax, d);
            }
        }
        const double Vtarget = a * C.volume[c];
        double lo = dmin, hi = dmax;
        for (int it = 0; it < 40; ++it) {
            const double mid = 0.5 * (lo + hi);
            const double V = cell_cut_volume(m, c, n, mid);
            if (V < Vtarget) lo = mid; else hi = mid;
            if (hi - lo < 1e-12) break;
        }
        const double delta = 0.5 * (lo + hi);
        PlicPatch P;
        P.normal = n; P.offset = delta; P.valid = true;
        // Approximate patch centroid as cell-centroid projected onto plane,
        // patch area from finite differences of cut-volume (∂V/∂δ ≈ A).
        const double dEps = 1e-4 * std::max(1e-12, dmax - dmin);
        const double Vp = cell_cut_volume(m, c, n, delta + dEps);
        const double Vm = cell_cut_volume(m, c, n, delta - dEps);
        P.area = std::max(0.0, (Vp - Vm) / (2.0 * dEps));
        // Centroid on plane closest to cell centroid:
        const double t = delta - (n.x*xc.x + n.y*xc.y + n.z*xc.z) + dot3(n, xc);
        (void)t;
        P.centroid = { xc.x + (delta - 0.0) * n.x,
                       xc.y + (delta - 0.0) * n.y,
                       xc.z + (delta - 0.0) * n.z };
        patches_[c] = P;
    }
    SIMALL_LOG_INFO("PLIC", "reconstructed (cells=", nC, ")");
}

}  // namespace simall::multiphase
