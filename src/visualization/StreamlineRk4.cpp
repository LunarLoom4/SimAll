// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/StreamlineRk4.cpp
// =============================================================================
#include "visualization/StreamlineRk4.hpp"

#include <cmath>

namespace simall::visualization
{

namespace
{

/// Compute one RK4 step.  Returns the (point, magnitude) after the step or
/// nullopt if any of the four samples leaves the domain.
struct StepOut
{
    util::Vec3d p;
    double mag;
};

std::optional<StepOut> rk4_step(const VectorSampler& f, const util::Vec3d& p, double h)
{
    auto k1o = f(p);
    if (!k1o)
        return std::nullopt;
    util::Vec3d k1 = *k1o;
    auto k2o = f(util::Vec3d{p.x + 0.5 * h * k1.x, p.y + 0.5 * h * k1.y, p.z + 0.5 * h * k1.z});
    if (!k2o)
        return std::nullopt;
    util::Vec3d k2 = *k2o;
    auto k3o = f(util::Vec3d{p.x + 0.5 * h * k2.x, p.y + 0.5 * h * k2.y, p.z + 0.5 * h * k2.z});
    if (!k3o)
        return std::nullopt;
    util::Vec3d k3 = *k3o;
    auto k4o = f(util::Vec3d{p.x + h * k3.x, p.y + h * k3.y, p.z + h * k3.z});
    if (!k4o)
        return std::nullopt;
    util::Vec3d k4 = *k4o;

    util::Vec3d avg{(k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x) / 6.0,
                    (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y) / 6.0,
                    (k1.z + 2.0 * k2.z + 2.0 * k3.z + k4.z) / 6.0};
    util::Vec3d np{p.x + h * avg.x, p.y + h * avg.y, p.z + h * avg.z};
    return StepOut{np, avg.norm()};
}

void trace_one_direction(const VectorSampler& sampler,
                         const util::Vec3d& seed,
                         const StreamlineConfig& cfg,
                         double sign,
                         std::vector<util::Vec3d>& outPts,
                         std::vector<double>& outMag)
{
    const double h = cfg.stepSize * sign;
    util::Vec3d p = seed;
    double length = 0.0;
    for (std::size_t i = 0; i < cfg.maxSteps; ++i) {
        auto out = rk4_step(sampler, p, h);
        if (!out)
            break;
        if (out->mag < cfg.stagnationEps)
            break;
        const double dx = out->p.x - p.x;
        const double dy = out->p.y - p.y;
        const double dz = out->p.z - p.z;
        length += std::sqrt(dx * dx + dy * dy + dz * dz);
        if (length > cfg.maxLength)
            break;
        p = out->p;
        outPts.push_back(p);
        outMag.push_back(out->mag);
    }
}

} // namespace

LineSet StreamlineRk4::trace(const VectorSampler& sampler,
                             const util::Vec3d& seed,
                             const StreamlineConfig& cfg)
{
    LineSet ls;
    if (!sampler)
        return ls;
    auto seedSample = sampler(seed);
    if (!seedSample)
        return ls;

    std::vector<util::Vec3d> fwdPts, bwdPts;
    std::vector<double> fwdMag, bwdMag;
    trace_one_direction(sampler, seed, cfg, 1.0, fwdPts, fwdMag);
    if (cfg.bidirectional) {
        trace_one_direction(sampler, seed, cfg, -1.0, bwdPts, bwdMag);
    }
    // Concatenate backward (reversed) + seed + forward.
    ls.lineOffsets.push_back(0);
    for (auto it = bwdPts.rbegin(); it != bwdPts.rend(); ++it) {
        ls.points.push_back(*it);
    }
    for (auto it = bwdMag.rbegin(); it != bwdMag.rend(); ++it) {
        ls.scalars.push_back(*it);
    }
    ls.points.push_back(seed);
    ls.scalars.push_back(seedSample->norm());
    for (auto& p : fwdPts)
        ls.points.push_back(p);
    for (auto& m : fwdMag)
        ls.scalars.push_back(m);
    ls.lineOffsets.push_back(static_cast<std::int32_t>(ls.points.size()));
    if (!cfg.colorByMagnitude)
        ls.scalars.clear();
    return ls;
}

LineSet StreamlineRk4::trace_many(const VectorSampler& sampler,
                                  const std::vector<util::Vec3d>& seeds,
                                  const StreamlineConfig& cfg)
{
    LineSet acc;
    acc.lineOffsets.push_back(0);
    for (const auto& s : seeds) {
        LineSet one = trace(sampler, s, cfg);
        if (one.points.empty())
            continue;
        // Append a single polyline (one has exactly one line, between
        // indices [0, points.size()]).
        const std::int32_t base = static_cast<std::int32_t>(acc.points.size());
        acc.points.insert(acc.points.end(), one.points.begin(), one.points.end());
        if (!one.scalars.empty()) {
            acc.scalars.insert(acc.scalars.end(), one.scalars.begin(), one.scalars.end());
        }
        acc.lineOffsets.push_back(base + static_cast<std::int32_t>(one.points.size()));
    }
    return acc;
}

} // namespace simall::visualization
