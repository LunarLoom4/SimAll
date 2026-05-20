// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/CutCellHexMesher.cpp
// =============================================================================
#include "meshing/CutCellHexMesher.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace simall::meshing
{

namespace
{
util::BoundingBox tri_aabb(const CutTriangle& t)
{
    util::BoundingBox a;
    for (const auto& p : t.v)
        a.expand(p);
    return a;
}
util::BoundingBox join(const util::BoundingBox& a, const util::BoundingBox& b)
{
    util::BoundingBox r = a;
    r.expand(b.min);
    r.expand(b.max);
    return r;
}
} // namespace

double CutCellHexMesher::point_tri_dist_sq(const util::Vec3d& p, const CutTriangle& t)
{
    // Eberly's projection algorithm for closest point on triangle.
    const util::Vec3d e0 = t.v[1] - t.v[0];
    const util::Vec3d e1 = t.v[2] - t.v[0];
    const util::Vec3d d = t.v[0] - p;
    const double a = e0.dot(e0);
    const double b = e0.dot(e1);
    const double c = e1.dot(e1);
    const double dD = e0.dot(d);
    const double eE = e1.dot(d);
    const double det = std::max(1e-30, a * c - b * b);
    double s = b * eE - c * dD;
    double tt = b * dD - a * eE;
    if (s + tt <= det) {
        if (s < 0) {
            if (tt < 0) {
                s = 0;
                tt = 0;
            } else {
                s = 0;
                tt = std::clamp(-eE / c, 0.0, 1.0);
            }
        } else if (tt < 0) {
            tt = 0;
            s = std::clamp(-dD / a, 0.0, 1.0);
        } else {
            s /= det;
            tt /= det;
        }
    } else {
        if (s < 0) {
            s = 0;
            tt = 1;
        } else if (tt < 0) {
            tt = 0;
            s = 1;
        } else {
            const double num = (c + eE) - (b + dD);
            if (num <= 0) {
                s = 0;
            } else {
                s = std::clamp(num / (a - 2 * b + c), 0.0, 1.0);
            }
            tt = 1 - s;
        }
    }
    const util::Vec3d q{t.v[0].x + s * e0.x + tt * e1.x,
                        t.v[0].y + s * e0.y + tt * e1.y,
                        t.v[0].z + s * e0.z + tt * e1.z};
    const util::Vec3d r = p - q;
    return r.dot(r);
}

bool CutCellHexMesher::ray_tri_intersect(const util::Vec3d& o,
                                         const util::Vec3d& d,
                                         const CutTriangle& t,
                                         double& tOut)
{
    // Möller-Trumbore.
    const util::Vec3d e1 = t.v[1] - t.v[0];
    const util::Vec3d e2 = t.v[2] - t.v[0];
    const util::Vec3d h = d.cross(e2);
    const double a = e1.dot(h);
    if (std::abs(a) < 1e-30)
        return false;
    const double f = 1.0 / a;
    const util::Vec3d s = o - t.v[0];
    const double u = f * s.dot(h);
    if (u < 0 || u > 1)
        return false;
    const util::Vec3d q = s.cross(e1);
    const double v = f * d.dot(q);
    if (v < 0 || u + v > 1)
        return false;
    tOut = f * e2.dot(q);
    return tOut > 1e-12;
}

int CutCellHexMesher::build_bvh(std::vector<BvhNode>& nodes,
                                std::vector<int>& triIdx,
                                const std::vector<CutTriangle>& tris,
                                int first,
                                int count)
{
    const int nodeIdx = static_cast<int>(nodes.size());
    nodes.emplace_back();
    util::BoundingBox bb;
    for (int i = 0; i < count; ++i)
        bb = join(bb, tri_aabb(tris[triIdx[first + i]]));
    nodes[nodeIdx].aabb = bb;
    if (count <= 6) {
        nodes[nodeIdx].triFirst = first;
        nodes[nodeIdx].triCount = count;
        return nodeIdx;
    }
    const util::Vec3d ext = bb.extent();
    int axis = (ext.x > ext.y && ext.x > ext.z) ? 0 : (ext.y > ext.z ? 1 : 2);
    const double mid = 0.5 * (bb.min.x + bb.max.x);
    (void) mid;
    auto cmp = [&](int a, int b) {
        const auto ca = tri_aabb(tris[a]).center();
        const auto cb = tri_aabb(tris[b]).center();
        return (axis == 0   ? ca.x
                : axis == 1 ? ca.y
                            : ca.z)
               < (axis == 0   ? cb.x
                  : axis == 1 ? cb.y
                              : cb.z);
    };
    std::sort(triIdx.begin() + first, triIdx.begin() + first + count, cmp);
    const int half = count / 2;
    const int L = build_bvh(nodes, triIdx, tris, first, half);
    const int R = build_bvh(nodes, triIdx, tris, first + half, count - half);
    nodes[nodeIdx].left = L;
    nodes[nodeIdx].right = R;
    return nodeIdx;
}

bool CutCellHexMesher::build(const CutCellParams& p,
                             const std::vector<CutTriangle>& tris,
                             Mesh& outMesh,
                             solver::FieldRegistry& outF)
{
    if (p.nx <= 0 || p.ny <= 0 || p.nz <= 0)
        return false;
    const double dx = (p.bboxMax.x - p.bboxMin.x) / p.nx;
    const double dy = (p.bboxMax.y - p.bboxMin.y) / p.ny;
    const double dz = (p.bboxMax.z - p.bboxMin.z) / p.nz;
    if (dx <= 0 || dy <= 0 || dz <= 0)
        return false;

    std::vector<BvhNode> bvh;
    std::vector<int> triIdx(tris.size());
    std::iota(triIdx.begin(), triIdx.end(), 0);
    if (!tris.empty()) {
        bvh.reserve(2 * tris.size());
        build_bvh(bvh, triIdx, tris, 0, static_cast<int>(tris.size()));
    }

    auto nearest = [&](const util::Vec3d& q) -> double {
        // Iterative BVH traversal for nearest squared distance.
        double best = std::numeric_limits<double>::max();
        if (bvh.empty())
            return 0.0;
        std::vector<int> stack = {0};
        while (!stack.empty()) {
            const int n = stack.back();
            stack.pop_back();
            const auto& nd = bvh[n];
            // Pruning lower-bound: squared distance from q to AABB.
            const double dxq = std::max({nd.aabb.min.x - q.x, 0.0, q.x - nd.aabb.max.x});
            const double dyq = std::max({nd.aabb.min.y - q.y, 0.0, q.y - nd.aabb.max.y});
            const double dzq = std::max({nd.aabb.min.z - q.z, 0.0, q.z - nd.aabb.max.z});
            if (dxq * dxq + dyq * dyq + dzq * dzq >= best)
                continue;
            if (nd.triCount > 0) {
                for (int i = 0; i < nd.triCount; ++i) {
                    const double d2 = point_tri_dist_sq(q, tris[triIdx[nd.triFirst + i]]);
                    if (d2 < best)
                        best = d2;
                }
            } else {
                if (nd.left >= 0)
                    stack.push_back(nd.left);
                if (nd.right >= 0)
                    stack.push_back(nd.right);
            }
        }
        return std::sqrt(best);
    };

    auto rayCount = [&](const util::Vec3d& q) -> int {
        // Robust parity test: cast +x ray and count triangle intersections.
        if (tris.empty())
            return 0;
        const util::Vec3d d{1, 0.0001234, 0.0005678}; // jitter to avoid edges
        int hits = 0;
        std::vector<int> stack = {0};
        while (!stack.empty()) {
            const int n = stack.back();
            stack.pop_back();
            const auto& nd = bvh[n];
            // AABB ray test (slab).
            double tmin = 0, tmax = std::numeric_limits<double>::max();
            for (int ax = 0; ax < 3; ++ax) {
                const double o = (ax == 0) ? q.x : (ax == 1) ? q.y : q.z;
                const double dr = (ax == 0) ? d.x : (ax == 1) ? d.y : d.z;
                const double mn = (ax == 0)   ? nd.aabb.min.x
                                  : (ax == 1) ? nd.aabb.min.y
                                              : nd.aabb.min.z;
                const double mx = (ax == 0)   ? nd.aabb.max.x
                                  : (ax == 1) ? nd.aabb.max.y
                                              : nd.aabb.max.z;
                if (std::abs(dr) < 1e-30) {
                    if (o < mn || o > mx) {
                        tmax = -1;
                        break;
                    }
                } else {
                    double t0 = (mn - o) / dr, t1 = (mx - o) / dr;
                    if (t0 > t1)
                        std::swap(t0, t1);
                    tmin = std::max(tmin, t0);
                    tmax = std::min(tmax, t1);
                    if (tmax < tmin)
                        break;
                }
            }
            if (tmax < tmin)
                continue;
            if (nd.triCount > 0) {
                for (int i = 0; i < nd.triCount; ++i) {
                    double tH;
                    if (ray_tri_intersect(q, d, tris[triIdx[nd.triFirst + i]], tH))
                        ++hits;
                }
            } else {
                if (nd.left >= 0)
                    stack.push_back(nd.left);
                if (nd.right >= 0)
                    stack.push_back(nd.right);
            }
        }
        return hits;
    };

    // Generate background grid nodes and per-cell classification.
    const int Nx = p.nx + 1, Ny = p.ny + 1, Nz = p.nz + 1;
    auto& N = outMesh.nodes();
    N.x.reserve(Nx * Ny * Nz);
    N.y.reserve(Nx * Ny * Nz);
    N.z.reserve(Nx * Ny * Nz);
    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                N.x.push_back(p.bboxMin.x + i * dx);
                N.y.push_back(p.bboxMin.y + j * dy);
                N.z.push_back(p.bboxMin.z + k * dz);
            }
    auto nodeId = [&](int i, int j, int k) { return static_cast<NodeId>(i + Nx * (j + Ny * k)); };

    // First pass: classify cells, collect kept cells.
    std::vector<int> keepIdx;
    std::vector<double> phiVec;
    std::vector<int> tagVec;
    keepIdx.reserve(static_cast<std::size_t>(p.nx) * p.ny * p.nz);
    for (int k = 0; k < p.nz; ++k)
        for (int j = 0; j < p.ny; ++j)
            for (int i = 0; i < p.nx; ++i) {
                const util::Vec3d c{p.bboxMin.x + (i + 0.5) * dx,
                                    p.bboxMin.y + (j + 0.5) * dy,
                                    p.bboxMin.z + (k + 0.5) * dz};
                const double dist = nearest(c);
                const int parity = rayCount(c) & 1;
                const double signed_d = (parity ? -dist : dist);
                int tag = (parity ? 1 : 0);
                if (dist < 0.51 * std::max({dx, dy, dz}))
                    tag = 2; // cut
                if (tag == 0 && signed_d > p.skinThickness)
                    continue;
                keepIdx.push_back(i + p.nx * (j + p.ny * k));
                phiVec.push_back(signed_d);
                tagVec.push_back(tag);
            }

    // Build cells (hex topology) — produce 6 quad faces per kept cell with
    // canonical owner/neighbour links; neighbours computed via i±1/j±1/k±1.
    auto& C = outMesh.cells();
    auto& Ff = outMesh.faces();
    const std::size_t nKept = keepIdx.size();
    C.volume.resize(nKept);
    C.centroidX.resize(nKept);
    C.centroidY.resize(nKept);
    C.centroidZ.resize(nKept);
    C.faceOffsets.assign(nKept + 1, 0);

    std::vector<int> ijkToCell(static_cast<std::size_t>(p.nx) * p.ny * p.nz, -1);
    for (std::size_t cid = 0; cid < nKept; ++cid)
        ijkToCell[keepIdx[cid]] = static_cast<int>(cid);

    // Emit unique faces by iterating directions {+x, +y, +z} per cell.
    auto add_face = [&](CellId owner,
                        CellId nb,
                        std::array<NodeId, 4> n,
                        double ax,
                        double ay,
                        double az,
                        double cx,
                        double cy,
                        double cz,
                        ZoneId zone) {
        const std::int32_t base = static_cast<std::int32_t>(Ff.nodeIndices.size());
        Ff.owner.push_back(owner);
        Ff.neighbor.push_back(nb);
        Ff.areaX.push_back(ax);
        Ff.areaY.push_back(ay);
        Ff.areaZ.push_back(az);
        Ff.centroidX.push_back(cx);
        Ff.centroidY.push_back(cy);
        Ff.centroidZ.push_back(cz);
        Ff.boundaryZone.push_back(zone);
        Ff.nodeOffsets.push_back(base);
        for (auto nn : n)
            Ff.nodeIndices.push_back(nn);
    };

    const ZoneId boundaryZone = outMesh.add_zone("auto_boundary", true);

    std::vector<std::vector<FaceId>> cellFaces(nKept);
    for (int k = 0; k < p.nz; ++k)
        for (int j = 0; j < p.ny; ++j)
            for (int i = 0; i < p.nx; ++i) {
                const int lin = i + p.nx * (j + p.ny * k);
                const int c = ijkToCell[lin];
                if (c < 0)
                    continue;
                const NodeId n000 = nodeId(i, j, k);
                const NodeId n100 = nodeId(i + 1, j, k);
                const NodeId n110 = nodeId(i + 1, j + 1, k);
                const NodeId n010 = nodeId(i, j + 1, k);
                const NodeId n001 = nodeId(i, j, k + 1);
                const NodeId n101 = nodeId(i + 1, j, k + 1);
                const NodeId n111 = nodeId(i + 1, j + 1, k + 1);
                const NodeId n011 = nodeId(i, j + 1, k + 1);
                const double cx = p.bboxMin.x + (i + 0.5) * dx;
                const double cy = p.bboxMin.y + (j + 0.5) * dy;
                const double cz = p.bboxMin.z + (k + 0.5) * dz;
                C.centroidX[c] = cx;
                C.centroidY[c] = cy;
                C.centroidZ[c] = cz;
                C.volume[c] = dx * dy * dz;

                auto neigh = [&](int di, int dj, int dk) -> CellId {
                    const int ii = i + di, jj = j + dj, kk = k + dk;
                    if (ii < 0 || ii >= p.nx || jj < 0 || jj >= p.ny || kk < 0 || kk >= p.nz)
                        return kBoundaryCell;
                    const int nb = ijkToCell[ii + p.nx * (jj + p.ny * kk)];
                    return (nb < 0) ? kBoundaryCell : static_cast<CellId>(nb);
                };
                // +x face (owned by this cell).
                {
                    const CellId nb = neigh(+1, 0, 0);
                    const FaceId fid = Ff.size();
                    add_face(static_cast<CellId>(c),
                             nb,
                             {n100, n110, n111, n101},
                             dy * dz,
                             0,
                             0,
                             p.bboxMin.x + (i + 1) * dx,
                             cy,
                             cz,
                             nb == kBoundaryCell ? boundaryZone : 0);
                    cellFaces[c].push_back(fid);
                    if (nb != kBoundaryCell)
                        cellFaces[nb].push_back(fid);
                }
                // +y face.
                {
                    const CellId nb = neigh(0, +1, 0);
                    const FaceId fid = Ff.size();
                    add_face(static_cast<CellId>(c),
                             nb,
                             {n010, n110, n111, n011},
                             0,
                             dx * dz,
                             0,
                             cx,
                             p.bboxMin.y + (j + 1) * dy,
                             cz,
                             nb == kBoundaryCell ? boundaryZone : 0);
                    cellFaces[c].push_back(fid);
                    if (nb != kBoundaryCell)
                        cellFaces[nb].push_back(fid);
                }
                // +z face.
                {
                    const CellId nb = neigh(0, 0, +1);
                    const FaceId fid = Ff.size();
                    add_face(static_cast<CellId>(c),
                             nb,
                             {n001, n101, n111, n011},
                             0,
                             0,
                             dx * dy,
                             cx,
                             cy,
                             p.bboxMin.z + (k + 1) * dz,
                             nb == kBoundaryCell ? boundaryZone : 0);
                    cellFaces[c].push_back(fid);
                    if (nb != kBoundaryCell)
                        cellFaces[nb].push_back(fid);
                }
                // -x boundary face only if no neighbour at i-1.
                if (neigh(-1, 0, 0) == kBoundaryCell) {
                    const FaceId fid = Ff.size();
                    add_face(static_cast<CellId>(c),
                             kBoundaryCell,
                             {n000, n010, n011, n001},
                             dy * dz,
                             0,
                             0,
                             p.bboxMin.x + i * dx,
                             cy,
                             cz,
                             boundaryZone);
                    cellFaces[c].push_back(fid);
                }
                if (neigh(0, -1, 0) == kBoundaryCell) {
                    const FaceId fid = Ff.size();
                    add_face(static_cast<CellId>(c),
                             kBoundaryCell,
                             {n000, n100, n101, n001},
                             0,
                             dx * dz,
                             0,
                             cx,
                             p.bboxMin.y + j * dy,
                             cz,
                             boundaryZone);
                    cellFaces[c].push_back(fid);
                }
                if (neigh(0, 0, -1) == kBoundaryCell) {
                    const FaceId fid = Ff.size();
                    add_face(static_cast<CellId>(c),
                             kBoundaryCell,
                             {n000, n100, n110, n010},
                             0,
                             0,
                             dx * dy,
                             cx,
                             cy,
                             p.bboxMin.z + k * dz,
                             boundaryZone);
                    cellFaces[c].push_back(fid);
                }
            }
    // Sentinel terminator.
    Ff.nodeOffsets.push_back(static_cast<std::int32_t>(Ff.nodeIndices.size()));

    // CSR cell→face.
    std::int32_t acc = 0;
    for (std::size_t c = 0; c < nKept; ++c) {
        C.faceOffsets[c] = acc;
        acc += static_cast<std::int32_t>(cellFaces[c].size());
    }
    C.faceOffsets[nKept] = acc;
    C.faceIndices.resize(acc);
    for (std::size_t c = 0; c < nKept; ++c) {
        for (std::size_t k = 0; k < cellFaces[c].size(); ++k)
            C.faceIndices[C.faceOffsets[c] + k] = cellFaces[c][k];
    }

    // Output classification fields.
    auto& phi = outF.scalar("phi", nKept);
    auto& tag = outF.scalar("cutTag", nKept);
    for (std::size_t c = 0; c < nKept; ++c) {
        phi[c] = phiVec[c];
        tag[c] = static_cast<double>(tagVec[c]);
    }
    return true;
}

} // namespace simall::meshing
