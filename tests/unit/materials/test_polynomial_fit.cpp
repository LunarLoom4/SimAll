// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/materials/test_polynomial_fit.cpp
// =============================================================================
#include "materials/PolynomialFit.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <stdexcept>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using simall::materials::PolynomialFit;

TEST_CASE("PolynomialFit Horner evaluation", "[materials][polyfit]")
{
    PolynomialFit p({1.0, 2.0, 3.0}, 0.0, 1000.0); // 1 + 2T + 3T^2
    CHECK_THAT(p.evaluate(0.0), WithinAbs(1.0, 1e-12));
    CHECK_THAT(p.evaluate(1.0), WithinAbs(6.0, 1e-12));
    CHECK_THAT(p.evaluate(2.0), WithinAbs(1.0 + 4.0 + 12.0, 1e-12));
    CHECK_THAT(p.derivative(0.0), WithinAbs(2.0, 1e-12));
    CHECK_THAT(p.derivative(1.0), WithinAbs(2.0 + 6.0, 1e-12));
    REQUIRE(p.order() == 2);
}

TEST_CASE("PolynomialFit clamps out-of-range", "[materials][polyfit]")
{
    PolynomialFit p({0.0, 1.0}, 10.0, 20.0, true);        // y = T
    CHECK_THAT(p.evaluate(5.0), WithinAbs(10.0, 1e-12));  // clamped to Tmin
    CHECK_THAT(p.evaluate(25.0), WithinAbs(20.0, 1e-12)); // clamped to Tmax
}

TEST_CASE("PolynomialFit validates inputs", "[materials][polyfit]")
{
    REQUIRE_THROWS_AS(PolynomialFit({}, 0.0, 1.0), std::invalid_argument);
    REQUIRE_THROWS_AS(PolynomialFit({1.0}, 10.0, 5.0), std::invalid_argument);
}
