// =============================================================================
// SimAll Beta - tests/verification/test_mms.cpp
// Week 19 - MMS verification.  Confirms that the 5-point central
// difference Laplacian converges at second order in h on the smooth sin·sin
// manufactured field (ratio of errors ≈ 4 when h halves).  Also verifies
// that the analytic residual vanishes exactly on the exact solution.
// Tag: [verification].
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "MmsSolutions.hpp"

using Catch::Matchers::WithinAbs;
namespace sv = simall::verification;

TEST_CASE("MMS heat-equation residual vanishes on exact solution",
           "[verification][mms]") {
    for (double x : {0.1, 0.3, 0.7}) {
        for (double y : {0.2, 0.5, 0.9}) {
            REQUIRE_THAT(sv::heat_residual(x, y, 0.05, 0.01),
                          WithinAbs(0.0, 1e-12));
        }
    }
}

TEST_CASE("FD5 Laplacian converges at second order", "[verification][mms]") {
    const double e1 = sv::fd5_laplacian_error(33);
    const double e2 = sv::fd5_laplacian_error(65);
    const double e3 = sv::fd5_laplacian_error(129);
    REQUIRE(e2 < e1);
    REQUIRE(e3 < e2);
    // Order-of-accuracy estimate: log2(e1/e2) should be ≈ 2.
    const double order12 = std::log2(e1 / e2);
    const double order23 = std::log2(e2 / e3);
    REQUIRE(order12 > 1.8);
    REQUIRE(order12 < 2.2);
    REQUIRE(order23 > 1.8);
    REQUIRE(order23 < 2.2);
}

TEST_CASE("MMS exact solution satisfies homogeneous Dirichlet BC",
           "[verification][mms]") {
    for (double t : {0.0, 0.1, 0.5}) {
        REQUIRE_THAT(sv::heat_exact(0.0, 0.5, t, 0.01), WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(sv::heat_exact(1.0, 0.5, t, 0.01), WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(sv::heat_exact(0.5, 0.0, t, 0.01), WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(sv::heat_exact(0.5, 1.0, t, 0.01), WithinAbs(0.0, 1e-12));
    }
}
