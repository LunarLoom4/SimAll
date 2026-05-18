// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/AdvancingFrontSurface.cpp
// =============================================================================
#include "meshing/AdvancingFrontSurface.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::meshing {

namespace {
inline util::Vec3d sub(const util::Vec3d& a, const util::Vec3d& b)
{ return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline double dot(const util::Vec3d& a, const util::Vec3d& b)
{ return a.x*b.x + a.y*b.y + a.z*b.z; }
inline util::Vec3d cross(const util::Vec3d& a, const util::Vec3d& b)
{ return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
inline double mag(const util::Vec3d& a) { return std::sqrt(dot(a,a)); }

bool segments_intersect_2d(const util::Vec3d& p1, const util::Vec3d& p2,
                           const util::Vec3d& p3, const util::Vec3d& p4) {
    // 2-D segment intersection in xy plane (z ignored).
    auto cross2 = [](double ax,double ay,double bx,double by){return ax*by - ay*bx;};
    const double d1 = cross2(p4.x-p3.x, p4.y-p3.y, p1.x-p3.x, p1.y-p3.y);
    const double d2 = cross2(p4.x-p3.x, p4.y-p3.y, p2.x-p3.x, p2.y-p3.y);
    const double d3 = cross2(p2.x-p1.x, p2.y-p1.y, p3.x-p1.x, p3.y-p1.y);
    const double d4 = cross2(p2.x-p1.x, p2.y-p1.y, p4.x-p1.x, p4.y-p1.y);
    if (((d1>0 && d2<0) || (d1<0 && d2>0)) &&
        ((d3>0 && d4<0) || (d3<0 && d4>0))) return true;
    return false;
}
}  // namespace

void AdvancingFrontSurface::initialize(AfsProps props, SizingFn sizing) {
    p_      = props;
    sizing_ = std::move(sizing);
    nodes_.clear(); front_.clear(); tris_.clear();
    SIMALL_LOG_INFO("Meshing",
        "AdvancingFrontSurface init: h_base=", p_.baseSize,
        " zone=", p_.outputZone);
}

NodeId AdvancingFrontSurface::find_or_create_node(const util::Vec3d& q, double tol) {
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (mag(sub(nodes_[i], q)) < tol) return static_cast<NodeId>(i);
    }
    nodes_.push_back(q);
    return static_cast<NodeId>(nodes_.size() - 1);
}

void AdvancingFrontSurface::seed_boundary(const std::vector<util::Vec3d>& loop) {
    if (loop.size() < 3) return;
    const double tol = 1e-6 * p_.baseSize;
    std::vector<NodeId> ids;
    ids.reserve(loop.size());
    for (const auto& v : loop) ids.push_back(find_or_create_node(v, tol));
    for (std::size_t i = 0; i < ids.size(); ++i) {
        const NodeId a = ids[i];
        const NodeId b = ids[(i+1) % ids.size()];
        if (a != b) front_.push_back({a, b});
    }
}

bool AdvancingFrontSurface::triangle_valid(const util::Vec3d& A,
                                           const util::Vec3d& B,
                                           const util::Vec3d& P,
                                           double h_ref) const {
    const double ab = mag(sub(B,A));
    const double bp = mag(sub(P,B));
    const double pa = mag(sub(A,P));
    if (ab < 1e-12 || bp < 1e-12 || pa < 1e-12) return false;
    // Heron-area quality measure: 4√3·A / (a²+b²+c²) ∈ [0,1].
    const double s = 0.5 * (ab + bp + pa);
    const double area2 = std::max(s*(s-ab)*(s-bp)*(s-pa), 0.0);
    const double area = std::sqrt(area2);
    const double q = 4.0 * std::sqrt(3.0) * area
                   / std::max(ab*ab + bp*bp + pa*pa, 1e-30);
    if (q < p_.minQuality) return false;
    // Aspect to local h.
    if (ab > 3*h_ref || bp > 3*h_ref || pa > 3*h_ref) return false;
    return true;
}

std::size_t AdvancingFrontSurface::generate(Mesh& outMesh) {
    if (front_.empty()) return 0;
    const double R0 = p_.searchRadiusMul * p_.baseSize;
    std::size_t iter = 0;
    while (!front_.empty() && iter < p_.maxIters) {
        ++iter;
        Edge e = front_.front();
        front_.erase(front_.begin());
        const util::Vec3d A = nodes_[e.a];
        const util::Vec3d B = nodes_[e.b];
        const util::Vec3d M = { 0.5*(A.x+B.x), 0.5*(A.y+B.y), 0.5*(A.z+B.z) };
        const double h = sizing_ ? sizing_(M.x,M.y,M.z) : p_.baseSize;
        const double L = mag(sub(B,A));
        // Construct apex P in plane: rotate (B-A) by 60° about normal.
        const util::Vec3d N{ p_.planeNormalX, p_.planeNormalY, p_.planeNormalZ };
        const util::Vec3d AB = sub(B, A);
        const util::Vec3d t  = cross(N, AB);
        const double tm = mag(t);
        if (tm < 1e-30) continue;
        const double ideal = std::sqrt(std::max(h*h - 0.25*L*L, 0.25*h*h));
        const util::Vec3d P{ M.x + (ideal/tm)*t.x,
                             M.y + (ideal/tm)*t.y,
                             M.z + (ideal/tm)*t.z };

        // Candidate snap: scan existing nodes for one inside the search disc.
        NodeId apex = static_cast<NodeId>(-1);
        double bestD = R0;
        for (std::size_t i = 0; i < nodes_.size(); ++i) {
            if (i == e.a || i == e.b) continue;
            const double d = mag(sub(nodes_[i], P));
            if (d < bestD) {
                if (triangle_valid(A, B, nodes_[i], h)) {
                    bestD = d; apex = static_cast<NodeId>(i);
                }
            }
        }
        if (apex == static_cast<NodeId>(-1)) {
            if (!triangle_valid(A, B, P, h)) continue;
            apex = find_or_create_node(P, 1e-6 * h);
        }
        tris_.push_back({e.a, e.b, apex});

        // Replace edge with two new ones; cancel reverse-duplicates in front.
        auto try_add = [&](NodeId u, NodeId v){
            for (auto it = front_.begin(); it != front_.end(); ++it) {
                if (it->a == v && it->b == u) { front_.erase(it); return; }
            }
            front_.push_back({u, v});
        };
        try_add(e.a, apex);
        try_add(apex, e.b);
    }

    // Emit mesh.
    auto& N = outMesh.nodes();
    const NodeId baseId = static_cast<NodeId>(N.size());
    for (const auto& q : nodes_) {
        N.x.push_back(q.x); N.y.push_back(q.y); N.z.push_back(q.z);
    }
    auto& F = outMesh.faces();
    const std::size_t fStart = F.size();
    for (const auto& t : tris_) {
        F.owner.push_back(static_cast<CellId>(F.size()));
        F.neighbor.push_back(kBoundaryCell);
        F.areaX.push_back(0.0); F.areaY.push_back(0.0); F.areaZ.push_back(0.0);
        F.centroidX.push_back(0.0); F.centroidY.push_back(0.0); F.centroidZ.push_back(0.0);
        F.boundaryZone.push_back(p_.outputZone);
        if (F.nodeOffsets.empty()) F.nodeOffsets.push_back(0);
        F.nodeIndices.push_back(baseId + t.a);
        F.nodeIndices.push_back(baseId + t.b);
        F.nodeIndices.push_back(baseId + t.c);
        F.nodeOffsets.push_back(static_cast<std::int32_t>(F.nodeIndices.size()));
    }
    SIMALL_LOG_INFO("Meshing",
        "AdvancingFront generated ", tris_.size(),
        " triangles in ", iter, " iters (faces ", fStart, "..", F.size(), ")");
    return tris_.size();
}

}  // namespace simall::meshing
