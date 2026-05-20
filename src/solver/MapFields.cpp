// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/MapFields.cpp
// Phase  : 23 Pass 12
// =============================================================================
#include "solver/MapFields.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace simall::solver
{

namespace
{

struct GridKey
{
    std::int64_t i, j, k;
    bool operator==(const GridKey& o) const noexcept { return i == o.i && j == o.j && k == o.k; }
};

struct GridKeyHash
{
    std::size_t operator()(const GridKey& g) const noexcept
    {
        // FNV-1a 64-bit over (i, j, k).
        std::uint64_t h = 1469598103934665603ULL;
        auto mix = [&](std::int64_t v) {
            const auto uv = static_cast<std::uint64_t>(v);
            for (int b = 0; b < 8; ++b) {
                h ^= (uv >> (b * 8)) & 0xFFu;
                h *= 1099511628211ULL;
            }
        };
        mix(g.i);
        mix(g.j);
        mix(g.k);
        return static_cast<std::size_t>(h);
    }
};

inline double sqr(double x) noexcept
{
    return x * x;
}

inline double dist2(double ax, double ay, double az, double bx, double by, double bz) noexcept
{
    return sqr(ax - bx) + sqr(ay - by) + sqr(az - bz);
}

} // namespace

// =============================================================================
MapFieldsStats map_fields(const meshing::Mesh& srcMesh,
                          const FieldRegistry& srcFields,
                          const meshing::Mesh& tgtMesh,
                          FieldRegistry& dstFields,
                          const MapFieldsOptions& opts)
{
    MapFieldsStats stats{};
    const auto& sc = srcMesh.cells();
    const auto& tc = tgtMesh.cells();
    const std::size_t nSrc = sc.size();
    const std::size_t nTgt = tc.size();
    stats.sourceCells = nSrc;
    stats.targetCells = nTgt;

    if (nSrc == 0 || nTgt == 0)
        return stats;

    // -------------------------------------------------------------------------
    // 1. Bounding box + bucket size derived from average inter-cell spacing.
    // -------------------------------------------------------------------------
    double xmin = std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double zmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    for (std::size_t c = 0; c < nSrc; ++c) {
        xmin = std::min(xmin, sc.centroidX[c]);
        xmax = std::max(xmax, sc.centroidX[c]);
        ymin = std::min(ymin, sc.centroidY[c]);
        ymax = std::max(ymax, sc.centroidY[c]);
        zmin = std::min(zmin, sc.centroidZ[c]);
        zmax = std::max(zmax, sc.centroidZ[c]);
    }
    const double bbVol = std::max(1.0e-30, (xmax - xmin) * (ymax - ymin) * (zmax - zmin));
    const double h = std::max(1.0e-12, std::cbrt(bbVol / static_cast<double>(nSrc)));
    const double inv = 1.0 / h;

    auto key_of = [&](double x, double y, double z) {
        return GridKey{static_cast<std::int64_t>(std::floor(x * inv)),
                       static_cast<std::int64_t>(std::floor(y * inv)),
                       static_cast<std::int64_t>(std::floor(z * inv))};
    };

    // -------------------------------------------------------------------------
    // 2. Bucket source cell centroids.
    // -------------------------------------------------------------------------
    std::unordered_map<GridKey, std::vector<std::uint64_t>, GridKeyHash> grid;
    grid.reserve(nSrc * 2);
    for (std::size_t c = 0; c < nSrc; ++c) {
        grid[key_of(sc.centroidX[c], sc.centroidY[c], sc.centroidZ[c])].push_back(
            static_cast<std::uint64_t>(c));
    }

    // -------------------------------------------------------------------------
    // 3. Per-target weights table. For Nearest, kNeighbors is forced to 1.
    // -------------------------------------------------------------------------
    const std::size_t kReq = (opts.method == MapMethod::Nearest)
                                 ? std::size_t{1}
                                 : std::max<std::size_t>(1, opts.kNeighbors);

    struct Pick
    {
        std::uint64_t cell;
        double w;
    };
    std::vector<std::vector<Pick>> picks(nTgt);

    constexpr int kMaxRing = 32;
    std::vector<std::pair<double, std::uint64_t>> cand;
    cand.reserve(64);

    for (std::size_t t = 0; t < nTgt; ++t) {
        const double tx = tc.centroidX[t];
        const double ty = tc.centroidY[t];
        const double tz = tc.centroidZ[t];
        const GridKey g0 = key_of(tx, ty, tz);

        cand.clear();
        int ring = 0;
        bool found = false;
        for (; ring <= kMaxRing; ++ring) {
            for (int dk = -ring; dk <= ring; ++dk)
                for (int dj = -ring; dj <= ring; ++dj)
                    for (int di = -ring; di <= ring; ++di) {
                        // Only walk the shell - inner cells were covered by smaller rings.
                        if (ring > 0 && std::abs(di) != ring && std::abs(dj) != ring
                            && std::abs(dk) != ring)
                            continue;
                        auto it = grid.find(GridKey{g0.i + di, g0.j + dj, g0.k + dk});
                        if (it == grid.end())
                            continue;
                        for (auto c : it->second) {
                            const double d2 = dist2(
                                tx, ty, tz, sc.centroidX[c], sc.centroidY[c], sc.centroidZ[c]);
                            cand.emplace_back(d2, c);
                        }
                    }
            if (cand.size() >= kReq) {
                found = true;
                break;
            }
        }

        if (!found) {
            // Fallback: scan all source cells. Rare; only when target is
            // far outside the source bbox or buckets are extremely sparse.
            ++stats.fallbackQueries;
            cand.clear();
            cand.reserve(nSrc);
            for (std::size_t c = 0; c < nSrc; ++c) {
                const double d2 =
                    dist2(tx, ty, tz, sc.centroidX[c], sc.centroidY[c], sc.centroidZ[c]);
                cand.emplace_back(d2, static_cast<std::uint64_t>(c));
            }
        }

        // Take the k nearest from cand.
        const std::size_t k = std::min(kReq, cand.size());
        std::partial_sort(cand.begin(),
                          cand.begin() + k,
                          cand.end(),
                          [](const auto& a, const auto& b) { return a.first < b.first; });

        auto& p = picks[t];
        p.reserve(k);

        if (opts.method == MapMethod::Nearest) {
            p.push_back({cand[0].second, 1.0});
            continue;
        }

        // Inverse-distance weighting.
        // If a candidate is essentially coincident, snap to it (avoid 1/0).
        constexpr double kEps2 = 1.0e-30;
        if (cand[0].first <= kEps2) {
            p.push_back({cand[0].second, 1.0});
            continue;
        }
        double wsum = 0.0;
        for (std::size_t i = 0; i < k; ++i) {
            const double d = std::sqrt(cand[i].first);
            const double w = 1.0 / std::pow(d, opts.idwPower);
            p.push_back({cand[i].second, w});
            wsum += w;
        }
        const double invW = 1.0 / wsum;
        for (auto& pi : p)
            pi.w *= invW;
    }

    // -------------------------------------------------------------------------
    // 4. Apply weights to every source field.
    // -------------------------------------------------------------------------
    for (const auto& [name, sfield] : srcFields.scalars()) {
        if (sfield.size() != nSrc)
            continue; // unsized field - skip
        auto& dfield = dstFields.scalar(name, nTgt);
        for (std::size_t t = 0; t < nTgt; ++t) {
            double acc = 0.0;
            for (const auto& pi : picks[t])
                acc += pi.w * sfield[pi.cell];
            dfield[t] = acc;
        }
        ++stats.scalarFieldsMapped;
    }
    for (const auto& [name, svec] : srcFields.vectors()) {
        if (svec.size() != nSrc)
            continue;
        auto& dvec = dstFields.vector(name, nTgt);
        for (std::size_t t = 0; t < nTgt; ++t) {
            double ax = 0.0, ay = 0.0, az = 0.0;
            for (const auto& pi : picks[t]) {
                ax += pi.w * svec.x[pi.cell];
                ay += pi.w * svec.y[pi.cell];
                az += pi.w * svec.z[pi.cell];
            }
            dvec.x[t] = ax;
            dvec.y[t] = ay;
            dvec.z[t] = az;
        }
        ++stats.vectorFieldsMapped;
    }

    return stats;
}

} // namespace simall::solver
