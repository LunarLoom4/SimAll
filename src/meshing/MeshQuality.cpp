// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MeshQuality.cpp
// =============================================================================
#include "meshing/MeshQuality.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>

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
    R.histVolume   = make_hist(C.volume);

    if (nC > 0) {
        R.minVolume = std::numeric_limits<double>::infinity();
        R.maxVolume = -std::numeric_limits<double>::infinity();
        for (std::size_t c = 0; c < nC; ++c) {
            const double v = C.volume[c];
            R.minVolume   = std::min(R.minVolume, v);
            R.maxVolume   = std::max(R.maxVolume, v);
            R.totalVolume += v;
        }
        R.meanVolume = R.totalVolume / static_cast<double>(nC);
    }
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
    s << "  cell volume: min=" << r.minVolume
      <<              "  max=" << r.maxVolume
      <<              "  mean=" << r.meanVolume
      <<              "  total=" << r.totalVolume << "\n";
    s << format_histogram(r.histSkewness, "skewness");
    s << format_histogram(r.histNonOrtho, "non-orthogonality (deg)");
    s << format_histogram(r.histAspect,   "aspect ratio");
    s << format_histogram(r.histVolume,   "cell volume");
    return s.str();
}

std::string MeshQuality::format_histogram(const QualityHistogram& h,
                                          std::string_view        label,
                                          std::size_t             barWidth) {
    std::ostringstream s;
    s << "  histogram: " << label
      << "  [" << h.vmin << " .. " << h.vmax << "]\n";
    std::size_t peak = 0;
    for (auto n : h.bins) peak = std::max(peak, n);
    if (peak == 0) {
        s << "    (no samples)\n";
        return s.str();
    }
    const double span    = h.vmax - h.vmin;
    const double binStep = (span > 0.0) ? span / static_cast<double>(h.bins.size())
                                        : 0.0;
    for (std::size_t i = 0; i < h.bins.size(); ++i) {
        const double lo  = h.vmin + binStep * static_cast<double>(i);
        const double hi  = h.vmin + binStep * static_cast<double>(i + 1);
        const std::size_t bars =
            (barWidth * h.bins[i] + peak / 2) / peak;
        s << "    [" << std::setw(10) << lo << " .. "
                     << std::setw(10) << hi << "] "
          << std::setw(8) << h.bins[i] << " | "
          << std::string(bars, '#') << "\n";
    }
    return s.str();
}

std::string MeshQuality::to_csv_faces(const MeshQualityReport& r) {
    std::ostringstream s;
    s << "faceId,skewness,nonOrthoDeg\n";
    s.precision(17);
    for (std::size_t f = 0; f < r.skewness.size(); ++f) {
        s << f << ',' << r.skewness[f] << ',' << r.nonOrthoDeg[f] << '\n';
    }
    return s.str();
}

std::string MeshQuality::to_csv_cells(const MeshQualityReport& r) {
    std::ostringstream s;
    s << "cellId,aspectRatio,negativeVolume\n";
    s.precision(17);
    for (std::size_t c = 0; c < r.aspectRatio.size(); ++c) {
        s << c << ',' << r.aspectRatio[c] << ','
          << static_cast<int>(r.negativeVolume[c]) << '\n';
    }
    return s.str();
}

std::vector<ZoneQualityStats>
MeshQuality::per_zone_stats(const Mesh& mesh, const MeshQualityReport& r) {
    const auto& F = mesh.faces();
    const std::size_t nF = F.size();

    struct Accum {
        std::size_t n = 0;
        double sumSkew = 0.0, maxSkew = 0.0;
        double sumNonO = 0.0, maxNonO = 0.0;
    };
    std::map<ZoneId, Accum> agg;

    for (std::size_t f = 0; f < nF; ++f) {
        if (F.neighbor[f] != kBoundaryCell) continue;
        const ZoneId z = F.boundaryZone[f];
        if (z == 0) continue;                  // unassigned boundary face
        auto& a = agg[z];
        const double sk = (f < r.skewness.size())    ? r.skewness[f]    : 0.0;
        const double no = (f < r.nonOrthoDeg.size()) ? r.nonOrthoDeg[f] : 0.0;
        a.n++;
        a.sumSkew += sk;
        a.sumNonO += no;
        a.maxSkew = std::max(a.maxSkew, sk);
        a.maxNonO = std::max(a.maxNonO, no);
    }

    std::vector<ZoneQualityStats> out;
    out.reserve(agg.size());
    Mesh& m = const_cast<Mesh&>(mesh);    // find_zone is non-const-only
    for (const auto& [z, a] : agg) {
        ZoneQualityStats s;
        s.id           = z;
        if (auto* zi = m.find_zone(z)) s.name = zi->name;
        s.nFaces       = a.n;
        s.maxSkewness  = a.maxSkew;
        s.maxNonOrtho  = a.maxNonO;
        s.meanSkewness = (a.n > 0) ? a.sumSkew / static_cast<double>(a.n) : 0.0;
        s.meanNonOrtho = (a.n > 0) ? a.sumNonO / static_cast<double>(a.n) : 0.0;
        out.push_back(std::move(s));
    }
    return out;
}

std::string MeshQuality::format_per_zone(const std::vector<ZoneQualityStats>& s) {
    std::ostringstream o;
    o << "Per-zone quality\n";
    o << "  " << std::left << std::setw(20) << "zone"
      <<        std::right << std::setw(8)  << "faces"
      <<        std::setw(12) << "maxSkew"
      <<        std::setw(12) << "meanSkew"
      <<        std::setw(12) << "maxNonO"
      <<        std::setw(12) << "meanNonO" << '\n';
    if (s.empty()) {
        o << "    (no boundary zones)\n";
        return o.str();
    }
    for (const auto& z : s) {
        o << "  " << std::left << std::setw(20)
          << (z.name.empty() ? std::string("zone") + std::to_string(z.id) : z.name)
          << std::right << std::setw(8)  << z.nFaces
          << std::setw(12) << z.maxSkewness
          << std::setw(12) << z.meanSkewness
          << std::setw(12) << z.maxNonOrtho
          << std::setw(12) << z.meanNonOrtho << '\n';
    }
    return o.str();
}

}  // namespace simall::meshing
