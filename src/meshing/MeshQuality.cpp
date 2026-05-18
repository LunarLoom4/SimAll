// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MeshQuality.cpp
// =============================================================================
#include "meshing/MeshQuality.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace simall::meshing {

namespace {

inline double dot3(double ax, double ay, double az,
                   double bx, double by, double bz) {
    return ax*bx + ay*by + az*bz;
}

QualityHistogram make_hist(const util::aligned_vector<double>& v) {
    QualityHistogram h;
    if (v.empty()) return h;
    h.vmin =  std::numeric_limits<double>::infinity();
    h.vmax = -std::numeric_limits<double>::infinity();
    for (double x : v) {
        if (!std::isfinite(x)) continue;
        h.vmin = std::min(h.vmin, x);
        h.vmax = std::max(h.vmax, x);
    }
    if (!std::isfinite(h.vmin)) { h.vmin = 0; h.vmax = 0; return h; }
    const double w = std::max(h.vmax - h.vmin, 1e-30);
    for (double x : v) {
        if (!std::isfinite(x)) continue;
        int idx = static_cast<int>(std::floor((x - h.vmin) / w * 10.0));
        idx = std::clamp(idx, 0, 9);
        h.bins[idx]++;
    }
    return h;
}

}  // namespace

MeshQualityReport MeshQuality::evaluate(const Mesh& m) {
    MeshQualityReport R;
    const auto& F = m.faces();
    const auto& C = m.cells();
    const std::size_t nF = F.size();
    const std::size_t nC = C.size();
    R.nFaces = nF; R.nCells = nC;
    R.skewness.assign(nF, 0.0);
    R.nonOrthoDeg.assign(nF, 0.0);
    R.aspectRatio.assign(nC, 1.0);
    R.negativeVolume.assign(nC, 0);

    for (std::size_t f = 0; f < nF; ++f) {
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
        const double Amag = std::sqrt(Ax*Ax + Ay*Ay + Az*Az);
        if (Amag < 1e-30) continue;
        const double nx = Ax/Amag, ny = Ay/Amag, nz = Az/Amag;
        const meshing::CellId o = F.owner[f], nb = F.neighbor[f];
        if (nb == meshing::kBoundaryCell) {
            R.nonOrthoDeg[f] = 0.0;
            // Boundary skewness: distance from face centroid to projection of
            // owner-centroid onto face plane, divided by √A (length scale).
            const double dx = C.centroidX[o] - F.centroidX[f];
            const double dy = C.centroidY[o] - F.centroidY[f];
            const double dz = C.centroidZ[o] - F.centroidZ[f];
            const double tn = dot3(dx, dy, dz, nx, ny, nz);
            const double tx = dx - tn*nx, ty = dy - tn*ny, tz = dz - tn*nz;
            R.skewness[f] = std::sqrt(tx*tx + ty*ty + tz*tz) / std::sqrt(Amag);
            continue;
        }
        const double dx = C.centroidX[nb] - C.centroidX[o];
        const double dy = C.centroidY[nb] - C.centroidY[o];
        const double dz = C.centroidZ[nb] - C.centroidZ[o];
        const double dmag = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (dmag < 1e-30) continue;
        const double cosa = (dx*nx + dy*ny + dz*nz) / dmag;
        const double ang = std::acos(std::clamp(std::abs(cosa), 0.0, 1.0))
                         * 180.0 / M_PI;
        R.nonOrthoDeg[f] = ang;
        // Skewness: distance from face centroid to the line connecting cell
        // centres, normalised by half the connection length.
        // Line param t s.t. closest point: pf - (po + t·d)
        const double rx = F.centroidX[f] - C.centroidX[o];
        const double ry = F.centroidY[f] - C.centroidY[o];
        const double rz = F.centroidZ[f] - C.centroidZ[o];
        const double t = (rx*dx + ry*dy + rz*dz) / (dmag*dmag);
        const double px = C.centroidX[o] + t*dx;
        const double py = C.centroidY[o] + t*dy;
        const double pz = C.centroidZ[o] + t*dz;
        const double sx = F.centroidX[f] - px;
        const double sy = F.centroidY[f] - py;
        const double sz = F.centroidZ[f] - pz;
        R.skewness[f] = std::sqrt(sx*sx + sy*sy + sz*sz) / (0.5 * dmag);
        R.maxNonOrtho = std::max(R.maxNonOrtho, ang);
        R.maxSkewness = std::max(R.maxSkewness, R.skewness[f]);
    }

    for (std::size_t c = 0; c < nC; ++c) {
        if (C.volume[c] <= 0.0) {
            R.negativeVolume[c] = 1;
            R.negativeCount++;
        }
        const int fs = C.faceOffsets[c], fe = C.faceOffsets[c+1];
        double amin = std::numeric_limits<double>::infinity();
        double amax = 0.0;
        for (int k = fs; k < fe; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            const double Amag = std::sqrt(F.areaX[fid]*F.areaX[fid]
                                        + F.areaY[fid]*F.areaY[fid]
                                        + F.areaZ[fid]*F.areaZ[fid]);
            amin = std::min(amin, Amag);
            amax = std::max(amax, Amag);
        }
        const double ar = (amin > 1e-30) ? amax / amin : 0.0;
        R.aspectRatio[c] = ar;
        R.maxAspect = std::max(R.maxAspect, ar);
    }

    R.histSkewness = make_hist(R.skewness);
    R.histNonOrtho = make_hist(R.nonOrthoDeg);
    R.histAspect   = make_hist(R.aspectRatio);
    return R;
}

std::string MeshQuality::format(const MeshQualityReport& r) {
    std::ostringstream s;
    s << "Mesh quality report\n";
    s << "  cells=" << r.nCells << "  faces=" << r.nFaces
      << "  negativeVolume=" << r.negativeCount << "\n";
    s << "  max skewness        = " << r.maxSkewness   << "\n";
    s << "  max non-orthogonality (deg) = " << r.maxNonOrtho << "\n";
    s << "  max face-area aspect ratio  = " << r.maxAspect   << "\n";
    return s.str();
}

}  // namespace simall::meshing
