// =============================================================================
// SimAll Beta - tests/performance/PerfHarness.hpp
// Week 19 - Micro-benchmark harness.  Header-only.  Provides a
// `run_benchmark(name, ops, body)` helper that times `body()` to within
// a configurable wall-clock budget and reports MOps/s.  The benchmark
// thresholds are intentionally generous (we just want to alarm on a
// >10× regression on a clean CI worker, not enforce a hard SLA).
// =============================================================================
#pragma once

#include "../regression/RegressionFramework.hpp"

#include <functional>

namespace simall::performance
{

struct BenchmarkResult
{
    std::string name;
    double seconds = 0.0;
    std::size_t opsPerformed = 0;
    [[nodiscard]] double mops() const
    {
        return double(opsPerformed) / std::max(seconds, 1e-12) * 1e-6;
    }
};

template <typename Body>
[[nodiscard]] inline BenchmarkResult run_benchmark(const std::string& name,
                                                   std::size_t opsPerIter,
                                                   Body body,
                                                   double budgetSec = 0.25)
{
    BenchmarkResult r{name, 0.0, 0};
    simall::regression::HiResTimer t;
    std::size_t iters = 0;
    while (t.seconds() < budgetSec) {
        body();
        r.opsPerformed += opsPerIter;
        ++iters;
        // Bail out early on extremely fast bodies to avoid millions of iterations.
        if (iters >= (1ULL << 22))
            break;
    }
    r.seconds = t.seconds();
    return r;
}

} // namespace simall::performance
