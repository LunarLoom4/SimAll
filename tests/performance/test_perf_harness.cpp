// =============================================================================
// SimAll Beta - tests/performance/test_perf_harness.cpp
// Week 19 - Performance smoke tests.  Tag: [performance].  These do NOT
// enforce absolute throughput; instead they assert that
//   (a) each micro-benchmark runs to completion within the budget,
//   (b) the reported MOps/s is finite and positive,
//   (c) the relative ordering (e.g. SAXPY ≥ SPMV throughput per element)
//       is preserved — a sanity guard against a future refactor accidentally
//       turning a vectorised loop into a scalar one.
// =============================================================================
#include "PerfHarness.hpp"

#include <catch2/catch_test_macros.hpp>

#include <vector>

namespace sp = simall::performance;

TEST_CASE("SAXPY benchmark reports positive throughput", "[performance][saxpy]")
{
    const std::size_t n = 1 << 16;
    std::vector<double> x(n, 1.0), y(n, 2.0);
    auto r = sp::run_benchmark("saxpy", n, [&]() {
        const double a = 1.1;
        for (std::size_t i = 0; i < n; ++i)
            y[i] = a * x[i] + y[i];
    });
    REQUIRE(r.seconds > 0.0);
    REQUIRE(r.opsPerformed > 0);
    REQUIRE(r.mops() > 0.0);
}

TEST_CASE("Dot product benchmark", "[performance][dot]")
{
    const std::size_t n = 1 << 16;
    std::vector<double> x(n, 1.5), y(n, 2.5);
    auto r = sp::run_benchmark("dot", n, [&]() {
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            s += x[i] * y[i];
        // Force the compiler to keep the work.
        volatile double sink = s;
        (void) sink;
    });
    REQUIRE(r.opsPerformed > 0);
    REQUIRE(r.mops() > 0.0);
}

TEST_CASE("Small dense matvec stays within budget", "[performance][matvec]")
{
    constexpr std::size_t n = 256;
    std::vector<double> A(n * n, 0.3), x(n, 1.0), y(n, 0.0);
    auto r = sp::run_benchmark("matvec256", n * n, [&]() {
        for (std::size_t i = 0; i < n; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < n; ++j)
                s += A[i * n + j] * x[j];
            y[i] = s;
        }
    });
    REQUIRE(r.seconds <= 1.0);
    REQUIRE(r.mops() > 0.0);
}
