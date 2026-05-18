// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MedialAxisDecomposition.cpp
// =============================================================================
#include "meshing/MedialAxisDecomposition.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>

namespace simall::meshing {

namespace {
inline double dist(const util::Vec3d& a, const util::Vec3d& b) {
    const double dx = a.x-b.x, dy = a.y-b.y, dz = a.z-b.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}
}  // namespace

void MedialAxisDecomposition::initialize(MedialProps props) {
    p_ = props;
    pts_.clear();
    SIMALL_LOG_INFO("Meshing",
        "MedialAxis init: clusterR=", p_.clusterRadius,
        " minBranch=", p_.branchMinLength);
}

std::size_t MedialAxisDecomposition::compute_medial_points(const Mesh& mesh) {
    const auto& C  = mesh.cells();
    const auto& Ff = mesh.faces();
    pts_.clear();
    for (std::size_t c = 0; c < C.size(); ++c) {
        // Inscribed-sphere radius ≈ min distance from cell centroid to any
        // boundary face attached to cell.
        const std::size_t off = C.faceOffsets[c];
        const std::size_t end = C.faceOffsets[c+1];
        double rMin = std::numeric_limits<double>::infinity();
        for (std::size_t k = off; k < end; ++k) {
            const std::size_t fi = C.faceIndices[k];
            if (Ff.neighbor[fi] != kBoundaryCell) continue;
            const double dx = Ff.centroidX[fi] - C.centroidX[c];
            const double dy = Ff.centroidY[fi] - C.centroidY[c];
            const double dz = Ff.centroidZ[fi] - C.centroidZ[c];
            rMin = std::min(rMin, std::sqrt(dx*dx + dy*dy + dz*dz));
        }
        if (!std::isfinite(rMin)) continue;
        pts_.push_back({ { C.centroidX[c], C.centroidY[c], C.centroidZ[c] },
                         rMin, -1 });
    }
    SIMALL_LOG_INFO("Meshing",
        "Medial points computed: ", pts_.size());
    return pts_.size();
}

std::size_t MedialAxisDecomposition::decompose(std::vector<MedialRegion>& out) {
    out.clear();
    if (pts_.empty()) return 0;
    // Greedy spatial clustering: nearest-neighbour chaining.
    const std::size_t N = pts_.size();
    std::vector<bool>  visited(N, false);
    std::vector<int>   branch(N, -1);
    int nextBranch = 0;
    for (std::size_t s = 0; s < N && nextBranch < static_cast<int>(p_.maxBranches); ++s) {
        if (visited[s]) continue;
        std::queue<std::size_t> q;
        q.push(s); visited[s] = true; branch[s] = nextBranch;
        while (!q.empty()) {
            const auto u = q.front(); q.pop();
            for (std::size_t v = 0; v < N; ++v) {
                if (visited[v]) continue;
                if (dist(pts_[u].position, pts_[v].position) <= p_.clusterRadius) {
                    visited[v] = true; branch[v] = nextBranch; q.push(v);
                }
            }
        }
        ++nextBranch;
    }
    for (std::size_t i = 0; i < N; ++i) pts_[i].branchId = branch[i];

    // Aggregate into MedialRegions (axis = AABB diagonal).
    out.assign(nextBranch, {});
    for (int b = 0; b < nextBranch; ++b) {
        out[b].branchId   = b;
        out[b].minCorner  = { 1e300,  1e300,  1e300};
        out[b].maxCorner  = {-1e300, -1e300, -1e300};
    }
    for (const auto& mp : pts_) {
        auto& r = out[mp.branchId];
        r.minCorner.x = std::min(r.minCorner.x, mp.position.x - mp.radius);
        r.minCorner.y = std::min(r.minCorner.y, mp.position.y - mp.radius);
        r.minCorner.z = std::min(r.minCorner.z, mp.position.z - mp.radius);
        r.maxCorner.x = std::max(r.maxCorner.x, mp.position.x + mp.radius);
        r.maxCorner.y = std::max(r.maxCorner.y, mp.position.y + mp.radius);
        r.maxCorner.z = std::max(r.maxCorner.z, mp.position.z + mp.radius);
    }
    // Axis = principal direction (longest box edge).
    for (auto& r : out) {
        const double dx = r.maxCorner.x - r.minCorner.x;
        const double dy = r.maxCorner.y - r.minCorner.y;
        const double dz = r.maxCorner.z - r.minCorner.z;
        const double cx = 0.5 * (r.maxCorner.x + r.minCorner.x);
        const double cy = 0.5 * (r.maxCorner.y + r.minCorner.y);
        const double cz = 0.5 * (r.maxCorner.z + r.minCorner.z);
        if (dx >= dy && dx >= dz) {
            r.axisStart = {r.minCorner.x, cy, cz};
            r.axisEnd   = {r.maxCorner.x, cy, cz};
        } else if (dy >= dz) {
            r.axisStart = {cx, r.minCorner.y, cz};
            r.axisEnd   = {cx, r.maxCorner.y, cz};
        } else {
            r.axisStart = {cx, cy, r.minCorner.z};
            r.axisEnd   = {cx, cy, r.maxCorner.z};
        }
    }
    // Drop branches shorter than min-length.
    out.erase(std::remove_if(out.begin(), out.end(), [&](const MedialRegion& r){
        return dist(r.axisStart, r.axisEnd) < p_.branchMinLength;
    }), out.end());

    SIMALL_LOG_INFO("Meshing",
        "Medial decomposition: branches=", nextBranch,
        " regions retained=", out.size());
    return out.size();
}

}  // namespace simall::meshing
