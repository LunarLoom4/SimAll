// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/materials/test_catalogs.cpp
// =============================================================================
#include "materials/SolidsCatalog.hpp"
#include "materials/FluidsCatalog.hpp"
#include "materials/MaterialAssignment.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace simall::materials;
using Catch::Matchers::WithinRel;

TEST_CASE("SolidsCatalog contains common engineering materials",
          "[materials][solids]") {
    SolidsCatalog cat;
    REQUIRE(cat.all().size() >= 20);
    const auto* al = cat.find("Aluminum-6061");
    REQUIRE(al != nullptr);
    CHECK_THAT(al->density,      WithinRel(2700.0, 0.01));
    CHECK_THAT(al->conductivity, WithinRel(167.0,  0.01));

    REQUIRE(cat.find("DoesNotExist") == nullptr);
}

TEST_CASE("FluidsCatalog covers air + water with Sutherland",
          "[materials][fluids]") {
    FluidsCatalog cat;
    const auto* air = cat.find("Air");
    REQUIRE(air != nullptr);
    REQUIRE(air->isGas);
    CHECK_THAT(air->sutherlandVisc.mu(273.15),
               WithinRel(1.716e-5, 0.005));

    const auto* h2o = cat.find("Water");
    REQUIRE(h2o != nullptr);
    REQUIRE_FALSE(h2o->isGas);
}

TEST_CASE("MaterialAssignment resolves via database + default fallback",
          "[materials][assignment]") {
    MaterialDatabase db;
    db.add("Steel");
    db.add("Aluminum");

    MaterialAssignment a(db);
    a.setCellZone(7,  "Steel");
    a.setCellZone(11, "Aluminum");
    a.setDefault("Steel");

    REQUIRE(a.materialForCellZone(7)  != nullptr);
    REQUIRE(a.materialForCellZone(11) != nullptr);
    REQUIRE(a.materialForCellZone(99) != nullptr);   // default
    REQUIRE(a.materialForFaceZone(7)  == nullptr);   // no face assignment
}
