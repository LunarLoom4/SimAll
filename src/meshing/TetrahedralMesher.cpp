// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/TetrahedralMesher.cpp
//
// Robust 3-D Bowyer-Watson tetrahedralization with:
//   - 4×4 orient3d predicate
//   - 5×5 in-sphere predicate
//   - cavity-rebuild insertion
//   - exterior flood-fill stripping
//   - circumcentre-based Delaunay refinement for quality
//
// References:
//   - Shewchuk (2002): "Delaunay refinement algorithms for triangular mesh
//     generation", Computational Geometry 22(1-3), pp. 21-74.
//   - Si (2015): "TetGen, a Delaunay-based quality tetrahedral mesh
//     generator", ACM TOMS 41(2), Article 11.
//   - Edelsbrunner & Shah (1996): "Incremental topological flipping works
//     for regular triangulations".
// =============================================================================
#include "meshing/TetrahedralMesher.hpp"
#include "meshing/Connectivity.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace simall::meshing {

namespace {

// ----- face key (sorted triple) for adjacency -----------------------------
struct FaceKey3 {
    std::array<std::uint32_t, 3> s;
    bool operator==(const FaceKey3& o) const noexcept { return s == o.s; }
};
struct FaceKey3Hash {
    std::size_t operator()(const FaceKey3& k) const noexcept {
        std::size_t h = 1469598103934665603ull;
        for (auto v : k.s) h ^= v + 0x9e3779b97f4a7c15ull + (h << 12) + (h >> 4);
        return h;
    }
};
inline FaceKey3 make_key(std::uint32_t a, std::uint32_t b, std::uint32_t c) {
    std::array<std::uint32_t, 3> s{a, b, c};
    std::sort(s.begin(), s.end());
    return {s};
}

}  // namespace

// =========================================================== predicates
double TetrahedralMesher::orient3d(std::uint32_t a, std::uint32_t b,
                                   std::uint32_t c, std::uint32_t d) const {
    const auto& A = pts_[a]; const auto& B = pts_[b];
    const auto& C = pts_[c]; const auto& D = pts_[d];
    const double adx = A.x - D.x, ady = A.y - D.y, adz = A.z - D.z;
    const double bdx = B.x - D.x, bdy = B.y - D.y, bdz = B.z - D.z;
    const double cdx = C.x - D.x, cdy = C.y - D.y, cdz = C.z - D.z;
    return adx * (bdy * cdz - bdz * cdy)
         - ady * (bdx * cdz - bdz * cdx)
         + adz * (bdx * cdy - bdy * cdx);
}

bool TetrahedralMesher::in_sphere(std::uint32_t a, std::uint32_t b,
                                  std::uint32_t c, std::uint32_t d,
                                  std::uint32_t e) const {
    // Standard 5x5 in-sphere predicate. Positive ⇒ e is strictly inside the
    // circumsphere of (a,b,c,d), assuming (a,b,c,d) has positive orient3d.
    // Flip sign if orient3d is negative.
    const double sign = orient3d(a, b, c, d) > 0 ? 1.0 : -1.0;
    const auto& A = pts_[a]; const auto& B = pts_[b];
    const auto& C = pts_[c]; const auto& D = pts_[d];
    const auto& E = pts_[e];
    const double ax = A.x - E.x, ay = A.y - E.y, az = A.z - E.z;
    const double bx = B.x - E.x, by = B.y - E.y, bz = B.z - E.z;
    const double cx = C.x - E.x, cy = C.y - E.y, cz = C.z - E.z;
    const double dx = D.x - E.x, dy = D.y - E.y, dz = D.z - E.z;
    const double an = ax*ax + ay*ay + az*az;
    const double bn = bx*bx + by*by + bz*bz;
    const double cn = cx*cx + cy*cy + cz*cz;
    const double dn = dx*dx + dy*dy + dz*dz;

    // 4x4 cofactor expansion along the last column (norms).
    auto det3 = [](double m00,double m01,double m02,
                   double m10,double m11,double m12,
                   double m20,double m21,double m22){
        return m00*(m11*m22 - m12*m21)
             - m01*(m10*m22 - m12*m20)
             + m02*(m10*m21 - m11*m20);
    };
    const double det =
          an * det3(bx,by,bz, cx,cy,cz, dx,dy,dz)
        - bn * det3(ax,ay,az, cx,cy,cz, dx,dy,dz)
        + cn * det3(ax,ay,az, bx,by,bz, dx,dy,dz)
        - dn * det3(ax,ay,az, bx,by,bz, cx,cy,cz);
    return sign * det > 0;
}

// =========================================================== super-tet
void TetrahedralMesher::initialize_super_tet(const std::vector<util::Vec3d>& nodes) {
    pts_ = nodes;
    double xmin = nodes[0].x, xmax = xmin;
    double ymin = nodes[0].y, ymax = ymin;
    double zmin = nodes[0].z, zmax = zmin;
    for (const auto& p : nodes) {
        xmin = std::min(xmin, p.x); xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y); ymax = std::max(ymax, p.y);
        zmin = std::min(zmin, p.z); zmax = std::max(zmax, p.z);
    }
    const double cx = 0.5 * (xmin + xmax);
    const double cy = 0.5 * (ymin + ymax);
    const double cz = 0.5 * (zmin + zmax);
    const double R  = std::max({xmax - xmin, ymax - ymin, zmax - zmin}) * 50.0 + 1.0;

    superStart_ = static_cast<std::uint32_t>(pts_.size());
    // Regular tetrahedron whose circumsphere encloses the bounding box.
    pts_.push_back({cx        , cy        , cz + 3*R});
    pts_.push_back({cx + 3*R  , cy        , cz -   R});
    pts_.push_back({cx -1.5*R , cy + 2.6*R, cz -   R});
    pts_.push_back({cx -1.5*R , cy - 2.6*R, cz -   R});

    Tet t{ {superStart_, superStart_+1, superStart_+2, superStart_+3} };
    if (orient3d(t.v[0], t.v[1], t.v[2], t.v[3]) < 0) std::swap(t.v[2], t.v[3]);
    tets_.clear();
    tets_.push_back(t);
}

// =========================================================== insertion
void TetrahedralMesher::insert_point(std::uint32_t p) {
    // 1. Find every tet whose circumsphere contains p.
    std::vector<std::size_t> bad; bad.reserve(32);
    for (std::size_t i = 0; i < tets_.size(); ++i) {
        const auto& t = tets_[i];
        if (in_sphere(t.v[0], t.v[1], t.v[2], t.v[3], p)) bad.push_back(i);
    }
    if (bad.empty()) return;

    // 2. Cavity boundary: faces appearing in exactly one bad tet.
    std::unordered_map<FaceKey3, std::array<std::uint32_t,3>, FaceKey3Hash> seen;
    auto consider = [&](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
        FaceKey3 k = make_key(a, b, c);
        auto it = seen.find(k);
        if (it == seen.end()) seen.emplace(k, std::array<std::uint32_t,3>{a, b, c});
        else                   seen.erase(it);
    };
    for (std::size_t i : bad) {
        const auto& t = tets_[i];
        consider(t.v[0], t.v[1], t.v[2]);
        consider(t.v[0], t.v[1], t.v[3]);
        consider(t.v[0], t.v[2], t.v[3]);
        consider(t.v[1], t.v[2], t.v[3]);
    }

    // 3. Remove bad tets (descending index).
    std::sort(bad.begin(), bad.end(), std::greater<>());
    for (std::size_t i : bad) tets_.erase(tets_.begin() + i);

    // 4. Stitch new tets, ensuring positive orientation.
    for (auto& kv : seen) {
        const auto& f = kv.second;
        Tet nt{ {f[0], f[1], f[2], p} };
        if (orient3d(nt.v[0], nt.v[1], nt.v[2], nt.v[3]) < 0)
            std::swap(nt.v[1], nt.v[2]);
        tets_.push_back(nt);
    }
}

// =========================================================== quality / refine
double TetrahedralMesher::tet_quality_ar(const Tet& t) const {
    // Aspect ratio = (longest edge) / (3 * inradius). Lower is better;
    // values near 1.0 indicate near-regular tets.
    const auto& A = pts_[t.v[0]]; const auto& B = pts_[t.v[1]];
    const auto& C = pts_[t.v[2]]; const auto& D = pts_[t.v[3]];
    auto edge = [](const util::Vec3d& a, const util::Vec3d& b){
        const double dx = a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
        return std::sqrt(dx*dx+dy*dy+dz*dz);
    };
    const double e[6] = { edge(A,B), edge(A,C), edge(A,D),
                          edge(B,C), edge(B,D), edge(C,D) };
    const double emax = *std::max_element(e, e+6);
    const double V    = std::abs(orient3d(t.v[0], t.v[1], t.v[2], t.v[3])) / 6.0;
    auto faceArea = [&](const util::Vec3d& p, const util::Vec3d& q, const util::Vec3d& r) {
        const double ax = q.x-p.x, ay = q.y-p.y, az = q.z-p.z;
        const double bx = r.x-p.x, by = r.y-p.y, bz = r.z-p.z;
        const double nx = ay*bz - az*by;
        const double ny = az*bx - ax*bz;
        const double nz = ax*by - ay*bx;
        return 0.5 * std::sqrt(nx*nx + ny*ny + nz*nz);
    };
    const double S = faceArea(A,B,C) + faceArea(A,B,D) + faceArea(A,C,D) + faceArea(B,C,D);
    const double r = (V > 0 && S > 0) ? (3.0 * V / S) : 0;
    return (r > 0) ? emax / (3.0 * r) : 1e30;
}

util::Vec3d TetrahedralMesher::circumcentre(const Tet& t) const {
    // Solve for the circumcentre of the tetrahedron using a 3×3 linear system.
    const auto& A = pts_[t.v[0]]; const auto& B = pts_[t.v[1]];
    const auto& C = pts_[t.v[2]]; const auto& D = pts_[t.v[3]];
    const double ba[3] = { B.x-A.x, B.y-A.y, B.z-A.z };
    const double ca[3] = { C.x-A.x, C.y-A.y, C.z-A.z };
    const double da[3] = { D.x-A.x, D.y-A.y, D.z-A.z };
    const double rhs[3] = {
        0.5 * (ba[0]*ba[0]+ba[1]*ba[1]+ba[2]*ba[2]),
        0.5 * (ca[0]*ca[0]+ca[1]*ca[1]+ca[2]*ca[2]),
        0.5 * (da[0]*da[0]+da[1]*da[1]+da[2]*da[2])
    };
    const double M[3][3] = {
        {ba[0], ba[1], ba[2]},
        {ca[0], ca[1], ca[2]},
        {da[0], da[1], da[2]}
    };
    const double det = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
                      - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
                      + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    if (std::abs(det) < 1e-30) return { A.x, A.y, A.z };
    const double inv = 1.0 / det;
    auto solve_for = [&](int comp) {
        double m[3][3]; std::memcpy(m, M, sizeof(M));
        for (int r = 0; r < 3; ++r) m[r][comp] = rhs[r];
        return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
              - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
              + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])) * inv;
    };
    return { A.x + solve_for(0), A.y + solve_for(1), A.z + solve_for(2) };
}

void TetrahedralMesher::refine(TetMeshOptions opt) {
    for (int pass = 0; pass < opt.maxRefinementPasses; ++pass) {
        std::vector<util::Vec3d> newPoints;
        const double h = opt.targetEdgeLength;
        for (const Tet& t : tets_) {
            if (t.v[0] >= superStart_ || t.v[1] >= superStart_ ||
                t.v[2] >= superStart_ || t.v[3] >= superStart_) continue;
            const auto& A = pts_[t.v[0]]; const auto& B = pts_[t.v[1]];
            const auto& C = pts_[t.v[2]]; const auto& D = pts_[t.v[3]];
            auto e2 = [](const util::Vec3d& a, const util::Vec3d& b){
                const double dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
                return dx*dx+dy*dy+dz*dz;
            };
            const double emax2 = std::max({e2(A,B), e2(A,C), e2(A,D),
                                           e2(B,C), e2(B,D), e2(C,D)});
            if (emax2 > h*h * 4.0 || tet_quality_ar(t) > 4.0) {
                newPoints.push_back(circumcentre(t));
            }
        }
        if (newPoints.empty()) break;
        for (const auto& p : newPoints) {
            const std::uint32_t idx = static_cast<std::uint32_t>(pts_.size());
            pts_.push_back(p);
            insert_point(idx);
        }
    }
}

// =========================================================== strip exterior
void TetrahedralMesher::strip_exterior(
        const std::vector<std::array<std::uint32_t, 3>>& boundary) {
    // Build adjacency: for each interior face, list the two adjacent tets.
    std::unordered_map<FaceKey3, std::array<int,2>, FaceKey3Hash> adj;
    for (std::size_t i = 0; i < tets_.size(); ++i) {
        const auto& t = tets_[i];
        std::array<std::array<std::uint32_t,3>, 4> faces = {{
            {t.v[0], t.v[1], t.v[2]}, {t.v[0], t.v[1], t.v[3]},
            {t.v[0], t.v[2], t.v[3]}, {t.v[1], t.v[2], t.v[3]}
        }};
        for (auto& f : faces) {
            FaceKey3 k = make_key(f[0], f[1], f[2]);
            auto& slot = adj[k];
            if (slot[0] == 0 && slot[1] == 0) { slot = {static_cast<int>(i+1), 0}; }
            else if (slot[1] == 0) slot[1] = static_cast<int>(i+1);
        }
    }
    // Mark boundary faces as blocking (flood-fill cannot cross them).
    std::unordered_set<FaceKey3, FaceKey3Hash> block;
    for (const auto& f : boundary) block.insert(make_key(f[0], f[1], f[2]));

    // Find a seed: any tet containing a super-tet vertex is exterior.
    std::vector<char> exterior(tets_.size(), 0);
    std::queue<int> q;
    for (std::size_t i = 0; i < tets_.size(); ++i) {
        const auto& t = tets_[i];
        if (t.v[0] >= superStart_ || t.v[1] >= superStart_ ||
            t.v[2] >= superStart_ || t.v[3] >= superStart_) {
            exterior[i] = 1;
            q.push(static_cast<int>(i));
        }
    }
    while (!q.empty()) {
        const int i = q.front(); q.pop();
        const auto& t = tets_[i];
        std::array<std::array<std::uint32_t,3>, 4> faces = {{
            {t.v[0], t.v[1], t.v[2]}, {t.v[0], t.v[1], t.v[3]},
            {t.v[0], t.v[2], t.v[3]}, {t.v[1], t.v[2], t.v[3]}
        }};
        for (auto& f : faces) {
            FaceKey3 k = make_key(f[0], f[1], f[2]);
            if (block.count(k)) continue;       // boundary blocks flood-fill
            auto it = adj.find(k); if (it == adj.end()) continue;
            for (int slot : {0, 1}) {
                const int j = it->second[slot] - 1;
                if (j >= 0 && !exterior[j]) { exterior[j] = 1; q.push(j); }
            }
        }
    }
    std::vector<Tet> keep; keep.reserve(tets_.size());
    for (std::size_t i = 0; i < tets_.size(); ++i)
        if (!exterior[i]) keep.push_back(tets_[i]);
    tets_ = std::move(keep);
}

// =========================================================== driver
void TetrahedralMesher::mesh(
        const std::vector<util::Vec3d>& points,
        const std::vector<std::array<std::uint32_t, 3>>& boundaryTris,
        TetMeshOptions opt,
        Mesh& out) {
    initialize_super_tet(points);
    // Insert every boundary node first.
    for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(points.size()); ++i)
        insert_point(i);
    refine(opt);
    if (opt.recoverBoundary) strip_exterior(boundaryTris);

    // ---- materialise as polyhedral Mesh -----------------------------------
    NodeStorage ns; ns.reserve(pts_.size() - 4);
    for (std::uint32_t i = 0; i < superStart_; ++i) {
        ns.x.push_back(pts_[i].x); ns.y.push_back(pts_[i].y); ns.z.push_back(pts_[i].z);
    }
    std::vector<CellDescriptor> cells;
    cells.reserve(tets_.size());
    for (const Tet& t : tets_) {
        if (t.v[0] >= superStart_ || t.v[1] >= superStart_ ||
            t.v[2] >= superStart_ || t.v[3] >= superStart_) continue;
        CellDescriptor cd;
        // 4 triangular faces with outward-normal (CCW from outside) ordering.
        // Standard tet: v0,v1,v2,v3 with positive orient3d.
        cd.faces = {
            {t.v[0], t.v[2], t.v[1]},
            {t.v[0], t.v[1], t.v[3]},
            {t.v[0], t.v[3], t.v[2]},
            {t.v[1], t.v[2], t.v[3]}
        };
        cells.push_back(std::move(cd));
    }
    ConnectivityBuilder::build(out, ns, cells);
    SIMALL_LOG_INFO("Mesh", "Tetrahedral mesh: ", cells.size(), " tets, ",
        ns.size(), " nodes");
}

}  // namespace simall::meshing
