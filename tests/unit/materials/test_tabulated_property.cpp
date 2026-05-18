// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/materials/test_tabulated_property.cpp
// =============================================================================
#include "materials/TabulatedProperty.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using simall::materials::TabulatedProperty1D;
using simall::materials::TabulatedProperty2D;
using simall::materials::TableInterp;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("TabulatedProperty1D linear interpolation", "[materials][tab1d]") {
    TabulatedProperty1D t({0.0, 1.0, 2.0}, {10.0, 20.0, 40.0});
    CHECK_THAT(t.evaluate(0.0),  WithinAbs(10.0, 1e-12));
    CHECK_THAT(t.evaluate(0.5),  WithinAbs(15.0, 1e-12));
    CHECK_THAT(t.evaluate(1.5),  WithinAbs(30.0, 1e-12));
    CHECK_THAT(t.evaluate(2.0),  WithinAbs(40.0, 1e-12));
}

TEST_CASE("TabulatedProperty1D clamps without extrapolation", "[materials][tab1d]") {
    TabulatedProperty1D t({0.0, 1.0}, {5.0, 10.0});
    CHECK_THAT(t.evaluate(-100.0), WithinAbs(5.0,  1e-12));
    CHECK_THAT(t.evaluate( 100.0), WithinAbs(10.0, 1e-12));
}

TEST_CASE("TabulatedProperty1D PCHIP is monotone for monotone data",
          "[materials][tab1d][pchip]") {
    TabulatedProperty1D t({0.0, 1.0, 2.0, 3.0, 4.0},
                          {0.0, 1.0, 4.0, 9.0, 16.0},
                          TableInterp::MonotonicCubic);
    double prev = t.evaluate(0.0);
    for (double x = 0.05; x <= 4.0; x += 0.05) {
        double v = t.evaluate(x);
        REQUIRE(v >= prev - 1e-12);
        prev = v;
    }
}

TEST_CASE("TabulatedProperty2D bilinear", "[materials][tab2d]") {
    // values[i*np + j]   for axes T = {0,1}, P = {0,2}
    TabulatedProperty2D t({0.0, 1.0}, {0.0, 2.0},
                          {0.0, 2.0,
                           1.0, 5.0});
    CHECK_THAT(t.evaluate(0.0, 0.0), WithinAbs(0.0, 1e-12));
    CHECK_THAT(t.evaluate(1.0, 2.0), WithinAbs(5.0, 1e-12));
    CHECK_THAT(t.evaluate(0.5, 1.0), WithinAbs(2.0, 1e-12));
}
