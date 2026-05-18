// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/ParametricSurfaceDelaunay.cpp
// =============================================================================
#include "meshing/ParametricSurfaceDelaunay.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace simall::meshing {

namespace {
struct EdgeHash {
    std::size_t operator()(const ParamEdge& e) const noexcept {
        std::uint64_t a = std::min(e.a, e.b), b = std::max(e.a, e.b);
        return std::hash<std::uint64_t>{}((a << 32) ^ b);
    }
};
struct EdgeEqUndirected {
    bool operator()(const ParamEdge& x, const ParamEdge& y) const noexcept {
        return (x.a == y.a && x.b == y.b) || (x.a == y.b && x.b == y.a);
    }
};
}  // namespace

double ParametricSurfaceDelaunay::orient2d(
        std::uint32_t a, std::uint32_t b, std::uint32_t c) const {
    return (pts_[b].u - pts_[a].u) * (pts_[c].v - pts_[a].v)
         - (pts_[b].v - pts_[a].v) * (pts_[c].u - pts_[a].u);
}

bool ParametricSurfaceDelaunay::in_circle(
        std::uint32_t a, std::uint32_t b, std::uint32_t c, std::uint32_t d) const {
    const double ax = pts_[a].u - pts_[d].u, ay = pts_[a].v - pts_[d].v;
    const double bx = pts_[b].u - pts_[d].u, by = pts_[b].v - pts_[d].v;
    const double cx = pts_[c].u - pts_[d].u, cy = pts_[c].v - pts_[d].v;
    const double det =
        (ax*ax + ay*ay) * (bx*cy - cx*by)
      - (bx*bx + by*by) * (ax*cy - cx*ay)
      + (cx*cx + cy*cy) * (ax*by - bx*ay);
    return det > 0;
}

void ParametricSurfaceDelaunay::insert_point(std::uint32_t p) {
    std::vector<std::size_t> bad; bad.reserve(16);
    for (std::size_t i = 0; i < tris_.size(); ++i) {
        auto& t = tris_[i];
        std::uint32_t a = t.v[0], b = t.v[1], c = t.v[2];
        if (orient2d(a, b, c) < 0) std::swap(b, c);
        if (in_circle(a, b, c, p)) bad.push_back(i);
    }
    if (bad.empty()) return;

    std::unordered_set<ParamEdge, EdgeHash, EdgeEqUndirected> seen;
    std::vector<ParamEdge> boundary;
    auto consider = [&](std::uint32_t a, std::uint32_t b) {
        ParamEdge e{a, b};
        auto it = seen.find(e);
        if (it == seen.end()) { seen.insert(e); boundary.push_back(e); }
        else                    { seen.erase(it); std::erase(boundary, e); }
    };
    for (std::size_t i : bad) {
        const auto& t = tris_[i];
        consider(t.v[0], t.v[1]);
        consider(t.v[1], t.v[2]);
        consider(t.v[2], t.v[0]);
    }
    std::sort(bad.begin(), bad.end(), std::greater<>());
    for (std::size_t i : bad) tris_.erase(tris_.begin() + i);

    for (const ParamEdge& e : boundary) {
        ParamTriangle nt{ e.a, e.b, p };
        if (orient2d(nt.v[0], nt.v[1], nt.v[2]) < 0) std::swap(nt.v[1], nt.v[2]);
        tris_.push_back(nt);
    }
}

void ParametricSurfaceDelaunay::recover_edges(
        const std::vector<ParamEdge>& constraints) {
    for (const ParamEdge& c : constraints) {
        bool present = false;
        for (const auto& t : tris_) {
            EdgeEqUndirected eq;
            if (eq(c, {t.v[0], t.v[1]}) || eq(c, {t.v[1], t.v[2]}) ||
                eq(c, {t.v[2], t.v[0]})) { present = true; break; }
        }
        if (present) continue;
        std::uint32_t mid = static_cast<std::uint32_t>(pts_.size());
        pts_.push_back({0.5 * (pts_[c.a].u + pts_[c.b].u),
                        0.5 * (pts_[c.a].v + pts_[c.b].v)});
        insert_point(mid);
    }
}

void ParametricSurfaceDelaunay::triangulate(
        const std::vector<ParamPoint>& input,
        const std::vector<ParamEdge>& constraints) {
    pts_ = input;
    tris_.clear();
    if (pts_.size() < 3) return;

    double umin = pts_[0].u, umax = umin, vmin = pts_[0].v, vmax = vmin;
    for (const auto& p : pts_) {
        umin = std::min(umin, p.u); umax = std::max(umax, p.u);
        vmin = std::min(vmin, p.v); vmax = std::max(vmax, p.v);
    }
    const double cu = 0.5 * (umin + umax), cv = 0.5 * (vmin + vmax);
    const double R  = std::max(umax - umin, vmax - vmin) * 10.0 + 1.0;
    const std::uint32_t s0 = static_cast<std::uint32_t>(pts_.size());
    pts_.push_back({cu - 2*R, cv -   R});
    pts_.push_back({cu + 2*R, cv -   R});
    pts_.push_back({cu,        cv + 2*R});
    tris_.push_back({s0, s0 + 1, s0 + 2});

    for (std::uint32_t i = 0; i < s0; ++i) insert_point(i);
    recover_edges(constraints);

    tris_.erase(std::remove_if(tris_.begin(), tris_.end(), [&](const ParamTriangle& t) {
        return t.v[0] >= s0 || t.v[1] >= s0 || t.v[2] >= s0;
    }), tris_.end());
    pts_.resize(s0);
}

}  // namespace simall::meshing
