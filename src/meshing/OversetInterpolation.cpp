// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/OversetInterpolation.cpp
// =============================================================================
#include "meshing/OversetInterpolation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace simall::meshing
{

namespace
{
util::BoundingBox join(const util::BoundingBox& a, const util::BoundingBox& b)
{
    util::BoundingBox r = a;
    r.expand(b.min);
    r.expand(b.max);
    return r;
}
} // namespace

void OversetInterpolation::cell_aabb(CellId c, util::BoundingBox& out) const
{
    const auto& N = donor_->nodes();
    const auto& Ff = donor_->faces();
    const auto& Cc = donor_->cells();
    out = util::BoundingBox{};
    for (int k = Cc.faceOffsets[c]; k < Cc.faceOffsets[c + 1]; ++k) {
        const auto fid = Cc.faceIndices[k];
        for (int n = Ff.nodeOffsets[fid]; n < Ff.nodeOffsets[fid + 1]; ++n) {
            const auto nid = Ff.nodeIndices[n];
            out.expand(util::Vec3d{N.x[nid], N.y[nid], N.z[nid]});
        }
    }
}

bool OversetInterpolation::point_in_cell(const util::Vec3d& p, CellId c) const
{
    // Convex test: for every face, check that p is on the owner-side of
    // the face plane (signed by outward normal pointing away from owner
    // centroid). Tolerant to within 1e-9 of plane.
    const auto& Ff = donor_->faces();
    const auto& Cc = donor_->cells();
    for (int k = Cc.faceOffsets[c]; k < Cc.faceOffsets[c + 1]; ++k) {
        const auto fid = Cc.faceIndices[k];
        const util::Vec3d nrm{Ff.areaX[fid], Ff.areaY[fid], Ff.areaZ[fid]};
        const util::Vec3d fc{Ff.centroidX[fid], Ff.centroidY[fid], Ff.centroidZ[fid]};
        const util::Vec3d cc{Cc.centroidX[c], Cc.centroidY[c], Cc.centroidZ[c]};
        // Outward normal sign w.r.t. owner.
        const double sign = (Ff.owner[fid] == c) ? +1.0 : -1.0;
        const util::Vec3d nOut = nrm * sign;
        const double rhs = nOut.dot(fc - cc); // > 0 expected
        const double lhs = nOut.dot(p - cc);
        if (lhs > rhs + 1e-9 * std::max(1.0, std::abs(rhs)))
            return false;
    }
    return true;
}

int OversetInterpolation::build_bvh(int first, int count)
{
    const int idx = static_cast<int>(bvh_.size());
    bvh_.emplace_back();
    util::BoundingBox bb;
    for (int i = 0; i < count; ++i)
        bb = join(bb, cellBox_[bvhIdx_[first + i]]);
    bvh_[idx].aabb = bb;
    if (count <= 4) {
        bvh_[idx].first = first;
        bvh_[idx].count = count;
        return idx;
    }
    const util::Vec3d ext = bb.extent();
    const int ax = (ext.x > ext.y && ext.x > ext.z) ? 0 : (ext.y > ext.z ? 1 : 2);
    std::sort(bvhIdx_.begin() + first, bvhIdx_.begin() + first + count, [&](int a, int b) {
        const auto ca = cellBox_[a].center();
        const auto cb = cellBox_[b].center();
        return (ax == 0 ? ca.x : ax == 1 ? ca.y : ca.z) < (ax == 0 ? cb.x : ax == 1 ? cb.y : cb.z);
    });
    const int h = count / 2;
    const int L = build_bvh(first, h);
    const int R = build_bvh(first + h, count - h);
    bvh_[idx].left = L;
    bvh_[idx].right = R;
    return idx;
}

CellId OversetInterpolation::locate(const util::Vec3d& p) const
{
    if (bvh_.empty())
        return static_cast<CellId>(-1);
    std::vector<int> stack = {0};
    while (!stack.empty()) {
        const int n = stack.back();
        stack.pop_back();
        const auto& nd = bvh_[n];
        if (p.x < nd.aabb.min.x - 1e-12 || p.x > nd.aabb.max.x + 1e-12)
            continue;
        if (p.y < nd.aabb.min.y - 1e-12 || p.y > nd.aabb.max.y + 1e-12)
            continue;
        if (p.z < nd.aabb.min.z - 1e-12 || p.z > nd.aabb.max.z + 1e-12)
            continue;
        if (nd.count > 0) {
            for (int i = 0; i < nd.count; ++i) {
                const CellId c = static_cast<CellId>(bvhIdx_[nd.first + i]);
                if (!holeMask_.empty() && holeMask_[c])
                    continue;
                if (point_in_cell(p, c))
                    return c;
            }
        } else {
            if (nd.left >= 0)
                stack.push_back(nd.left);
            if (nd.right >= 0)
                stack.push_back(nd.right);
        }
    }
    return static_cast<CellId>(-1);
}

bool OversetInterpolation::build(const Mesh& donor,
                                 const Mesh& receptor,
                                 const std::string& holeMaskField,
                                 const solver::FieldRegistry& donorF)
{
    donor_ = &donor;
    receptor_ = &receptor;
    const std::size_t nD = donor.cells().size();
    const std::size_t nR = receptor.cells().size();

    cellBox_.resize(nD);
    bvhIdx_.resize(nD);
    std::iota(bvhIdx_.begin(), bvhIdx_.end(), 0);
    for (std::size_t c = 0; c < nD; ++c)
        cell_aabb(static_cast<CellId>(c), cellBox_[c]);
    bvh_.clear();
    bvh_.reserve(2 * nD);
    if (nD > 0)
        build_bvh(0, static_cast<int>(nD));

    holeMask_.clear();
    if (!holeMaskField.empty()) {
        const auto* m = donorF.find_scalar(holeMaskField);
        if (m && m->size() == nD) {
            holeMask_.resize(nD);
            for (std::size_t c = 0; c < nD; ++c)
                holeMask_[c] = ((*m)[c] > 0.5) ? 1u : 0u;
        }
    }

    stencils_.assign(nR, {});
    const auto& Cr = receptor.cells();
    const auto& CrF = receptor.faces();
    const auto& Cd = donor.cells();
    const auto& CdF = donor.faces();

    for (std::size_t r = 0; r < nR; ++r) {
        const util::Vec3d pc{Cr.centroidX[r], Cr.centroidY[r], Cr.centroidZ[r]};
        const CellId d = locate(pc);
        if (d == static_cast<CellId>(-1))
            continue;

        auto& st = stencils_[r];
        st.receptor = static_cast<std::int64_t>(r);
        // Donor cell + face-neighbours.
        std::array<std::int64_t, 8> ids{static_cast<std::int64_t>(d), -1, -1, -1, -1, -1, -1, -1};
        std::array<double, 8> wts{0, 0, 0, 0, 0, 0, 0, 0};
        int k = 1;
        for (int e = Cd.faceOffsets[d]; e < Cd.faceOffsets[d + 1] && k < 8; ++e) {
            const auto fid = Cd.faceIndices[e];
            const CellId nb = (CdF.owner[fid] == d) ? CdF.neighbor[fid] : CdF.owner[fid];
            if (nb == kBoundaryCell)
                continue;
            if (!holeMask_.empty() && holeMask_[nb])
                continue;
            ids[k++] = static_cast<std::int64_t>(nb);
        }
        // Inverse-distance weights.
        double wsum = 0.0;
        for (int i = 0; i < k; ++i) {
            const auto cid = ids[i];
            const double dx = Cd.centroidX[cid] - pc.x;
            const double dy = Cd.centroidY[cid] - pc.y;
            const double dz = Cd.centroidZ[cid] - pc.z;
            const double d2 = dx * dx + dy * dy + dz * dz + 1e-30;
            const double w = 1.0 / d2;
            wts[i] = w;
            wsum += w;
        }
        for (int i = 0; i < k; ++i)
            wts[i] /= wsum;
        st.donors = ids;
        st.weights = wts;
        st.count = k;
        (void) CrF;
    }
    return true;
}

void OversetInterpolation::interpolate_scalar(const solver::ScalarField& src,
                                              solver::ScalarField& dst) const
{
    for (const auto& st : stencils_) {
        if (st.receptor < 0)
            continue;
        double s = 0.0;
        for (int i = 0; i < st.count; ++i)
            s += st.weights[i] * src[st.donors[i]];
        dst[st.receptor] = s;
    }
}
void OversetInterpolation::interpolate_vector(const solver::VectorField& src,
                                              solver::VectorField& dst) const
{
    for (const auto& st : stencils_) {
        if (st.receptor < 0)
            continue;
        double sx = 0, sy = 0, sz = 0;
        for (int i = 0; i < st.count; ++i) {
            const double w = st.weights[i];
            const auto j = st.donors[i];
            sx += w * src.x[j];
            sy += w * src.y[j];
            sz += w * src.z[j];
        }
        dst.x[st.receptor] = sx;
        dst.y[st.receptor] = sy;
        dst.z[st.receptor] = sz;
    }
}

} // namespace simall::meshing
