// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/core/test_unit_system.cpp
// =============================================================================
#include "core/UnitSystem.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace simall::core::units;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("UnitSystem length conversions", "[core][units]")
{
    auto v = convert(Quantity::Length, 1.0, "in", "mm");
    REQUIRE(v.has_value());
    REQUIRE_THAT(*v, WithinAbs(25.4, 1e-9));

    v = convert(Quantity::Length, 1.0, "ft", "m");
    REQUIRE(v.has_value());
    REQUIRE_THAT(*v, WithinAbs(0.3048, 1e-9));
}

TEST_CASE("UnitSystem pressure (psi → Pa)", "[core][units]")
{
    auto v = toSI(Quantity::Pressure, 14.7, "psi");
    REQUIRE(v.has_value());
    REQUIRE_THAT(*v, WithinRel(101352.9322, 1e-4));
}

TEST_CASE("UnitSystem temperature affine (C ↔ F)", "[core][units]")
{
    auto si = toSI(Quantity::Temperature, 0.0, "C");
    REQUIRE(si.has_value());
    REQUIRE_THAT(*si, WithinAbs(273.15, 1e-9));

    auto f = convert(Quantity::Temperature, 100.0, "C", "F");
    REQUIRE(f.has_value());
    REQUIRE_THAT(*f, WithinAbs(212.0, 1e-9));
}

TEST_CASE("UnitSystem returns nullopt for unknown labels", "[core][units]")
{
    auto v = toSI(Quantity::Length, 1.0, "lightyear");
    REQUIRE_FALSE(v.has_value());
}

TEST_CASE("UnitSystem system defaults", "[core][units]")
{
    double siVal = 100.0; // m
    double ft = fromSI(Quantity::Length, siVal, UnitSystem::USCS);
    REQUIRE_THAT(ft, WithinAbs(100.0 / 0.3048, 1e-9));
}
