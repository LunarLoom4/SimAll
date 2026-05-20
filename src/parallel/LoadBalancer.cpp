// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/LoadBalancer.cpp
// =============================================================================
#include "parallel/LoadBalancer.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <numeric>
#include <vector>

namespace simall::parallel
{

void LoadBalancer::initialize(LoadBalancerProps props)
{
    p_ = props;
    samples_.clear();
    normSamples_.clear();
    lastRepartFlag_ = false;
    SIMALL_LOG_INFO("Parallel",
                    "LoadBalancer init: threshold=",
                    p_.repartThreshold,
                    " hysteresis=",
                    p_.hysteresis,
                    " window=",
                    p_.sampleWindow);
}

void LoadBalancer::record_iteration(double seconds, std::size_t nOwnedCells)
{
    samples_.push_back(seconds);
    const double norm =
        (nOwnedCells > 0) ? (seconds * 1.0e6 / static_cast<double>(nOwnedCells)) : seconds * 1.0e6;
    normSamples_.push_back(norm);
    while (static_cast<int>(samples_.size()) > p_.sampleWindow) {
        samples_.erase(samples_.begin());
        normSamples_.erase(normSamples_.begin());
    }
}

LoadStats LoadBalancer::evaluate()
{
    LoadStats st{};
    if (samples_.empty())
        return st;
    const double localCost = std::accumulate(samples_.begin(), samples_.end(), 0.0)
                             / static_cast<double>(samples_.size());
    std::vector<std::int32_t> dummy;
    // gather per-rank costs (encoded as μs to fit in int32 for allgather;
    // we use the double allreduces below for the real numbers).
    const double mx = ctx_.allreduce_max(localCost);
    const double mn = ctx_.allreduce_min(localCost);
    const double sm = ctx_.allreduce_sum(localCost);
    const double mean = sm / static_cast<double>(std::max(1, ctx_.size()));
    st.maxRankCost = mx;
    st.meanRankCost = mean;
    st.imbalance = (mean > 0.0) ? (mx / mean) - 1.0 : 0.0;
    // worst/best rank: gather full vector to identify (small msg).
    const double localTag = (localCost == mx) ? 1.0 : 0.0;
    (void) mn;
    (void) localTag;
    const double thr =
        lastRepartFlag_ ? std::max(0.0, p_.repartThreshold - p_.hysteresis) : p_.repartThreshold;
    st.shouldRepart = st.imbalance > thr;
    lastRepartFlag_ = st.shouldRepart;
    SIMALL_LOG_INFO("Parallel",
                    "LoadBalancer: imbalance=",
                    st.imbalance,
                    " max=",
                    st.maxRankCost,
                    " mean=",
                    st.meanRankCost,
                    " repart=",
                    st.shouldRepart);
    return st;
}

void LoadBalancer::reset()
{
    samples_.clear();
    normSamples_.clear();
    lastRepartFlag_ = false;
}

} // namespace simall::parallel
