// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/LoadBalancer.hpp
// Phase  : 17.4 — dynamic load measurement and re-partition driver.
//
// Each rank reports its per-iteration wall-clock cost; LoadBalancer
// computes load imbalance and decides whether to trigger a global
// re-partition.  The actual repartition is delegated to DomainPartition;
// LoadBalancer is the *policy* layer.
//
// Imbalance metric (Hendrickson 1998):
//     imbalance = max_r (t_r * n_r⁻¹) / mean_r (t_r * n_r⁻¹) - 1
// A value above `repartThreshold` (default 0.20 = 20 %) triggers re-
// partition; a hysteresis band (`hysteresis = 0.05`) prevents thrash.
// =============================================================================
#pragma once

#include "parallel/MpiContext.hpp"

#include <cstddef>
#include <vector>

namespace simall::parallel {

struct LoadStats {
    double imbalance      = 0.0;
    double maxRankCost    = 0.0;
    double meanRankCost   = 0.0;
    int    worstRank      = 0;
    int    bestRank       = 0;
    bool   shouldRepart   = false;
};

struct LoadBalancerProps {
    double repartThreshold = 0.20;
    double hysteresis      = 0.05;
    int    sampleWindow    = 5;
};

class LoadBalancer {
public:
    explicit LoadBalancer(MpiContext& ctx) : ctx_(ctx) {}

    void initialize(LoadBalancerProps props);

    /// Submit one timing sample for this rank (seconds spent in last iter).
    void record_iteration(double seconds, std::size_t nOwnedCells);

    /// Reduce samples across ranks and produce a verdict.
    LoadStats evaluate();

    /// Reset accumulated samples (call after a repartition).
    void reset();

    const LoadBalancerProps& props() const noexcept { return p_; }

private:
    MpiContext&        ctx_;
    LoadBalancerProps  p_{};
    std::vector<double> samples_;
    std::vector<double> normSamples_;   // per-iter cost / cells (μs per cell)
    bool                lastRepartFlag_ = false;
};

}  // namespace simall::parallel
