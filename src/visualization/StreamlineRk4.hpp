// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/StreamlineRk4.hpp
// Phase  : 3 / Week 14 (filters: streamline integration)
//
// Classical four-stage Runge–Kutta integration of a vector sampler from a
// seed point, with optional reverse pass for the upstream half of the
// streamline.  Termination conditions: sampler returns nullopt (out of
// domain), step count exceeded, length exceeded, or local stagnation
// (|v| < eps).
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <vector>

namespace simall::visualization
{

struct StreamlineConfig
{
    double stepSize = 0.05;      // world-space units
    std::size_t maxSteps = 4000; // per direction
    double maxLength = 1e6;      // world-space cumulative
    double stagnationEps = 1e-12;
    bool bidirectional = true; // integrate forward AND backward
    bool colorByMagnitude = true;
};

class StreamlineRk4
{
public:
    /// Trace a single streamline from `seed`.  Returns a LineSet with one
    /// polyline.  Empty when the seed is out of the sampler's domain.
    static LineSet trace(const VectorSampler& sampler,
                         const util::Vec3d& seed,
                         const StreamlineConfig& cfg = {});

    /// Trace from many seeds; outputs are concatenated into one LineSet.
    static LineSet trace_many(const VectorSampler& sampler,
                              const std::vector<util::Vec3d>& seeds,
                              const StreamlineConfig& cfg = {});
};

} // namespace simall::visualization
