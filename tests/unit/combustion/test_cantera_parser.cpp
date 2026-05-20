// =============================================================================
// SimAll Beta - Unit Tests
// File   : tests/unit/combustion/test_cantera_parser.cpp
// Phase  : 11 — sanity round-trip for the in-tree Cantera YAML reader.
// =============================================================================
#include "combustion/CanteraParser.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace simall::combustion;

namespace
{

constexpr const char* kTinyMech = R"YAML(
elements:
  - H
  - O

species:
  - name: H2
    composition: {H: 2}
    thermo:
      model: NASA7
  - name: O2
    composition: {O: 2}
    thermo:
      model: NASA7
  - name: H2O
    composition: {H: 2, O: 1}
    thermo:
      model: NASA7

reactions:
  - equation: "2 H2 + O2 <=> 2 H2O"
    rate-constant: {A: 1.2e+17, b: -1.0, Ea: 0.0}
  - equation: "H2 + O2 => H2O + O"
    rate-constant: {A: 3.5e+13, b: 0.0, Ea: 7.1e+04}
)YAML";

} // namespace

TEST_CASE("CanteraParser parses a minimal H2/O2 YAML mechanism", "[combustion][cantera]")
{
    CanteraParser parser;
    auto mech = parser.parseString(kTinyMech, CanteraFormat::Yaml);

    REQUIRE(mech.elements.size() == 2);
    REQUIRE(mech.species.size() == 3);
    REQUIRE(mech.reactions.size() == 2);

    // Species composition was decoded.
    bool sawH2O = false;
    for (const auto& s : mech.species) {
        if (s.name == "H2O") {
            sawH2O = true;
            REQUIRE(s.composition.at("H") == 2);
            REQUIRE(s.composition.at("O") == 1);
        }
    }
    REQUIRE(sawH2O);

    // Equation direction is preserved.
    REQUIRE(mech.reactions[0].direction == ReactionDirection::Reversible);
    REQUIRE(mech.reactions[1].direction == ReactionDirection::Forward);

    // Stoichiometry parsed from "2 H2 + O2 <=> 2 H2O".
    REQUIRE(mech.reactions[0].reactants.at("H2") == Catch::Approx(2.0));
    REQUIRE(mech.reactions[0].reactants.at("O2") == Catch::Approx(1.0));
    REQUIRE(mech.reactions[0].products.at("H2O") == Catch::Approx(2.0));

    // Arrhenius coefficients survived the JSON-style inline parser.
    REQUIRE(mech.reactions[1].fwd.A == Catch::Approx(3.5e13));
    REQUIRE(mech.reactions[1].fwd.beta == Catch::Approx(0.0));
    REQUIRE(mech.reactions[1].fwd.Ea == Catch::Approx(7.1e4));
}
