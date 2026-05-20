// =============================================================================
// SimAll Beta - Unit Tests
// File   : tests/unit/zones/test_selection_manager.cpp
// Phase  : 25 — NamedSelection persistence and set algebra.
// =============================================================================
#include "zones/SelectionManager.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace simall::zones;

namespace
{
NamedSelection make(const std::string& name, std::vector<simall::util::PersistentId> ids)
{
    NamedSelection s;
    s.name = name;
    s.entity = SelectionEntity::Face;
    s.ids = std::move(ids);
    return s;
}
} // namespace

TEST_CASE("SelectionManager define/find/remove round-trip", "[zones][selection]")
{
    SelectionManager mgr;
    REQUIRE_FALSE(mgr.define(make("inlet", {10, 11, 12})));
    REQUIRE_FALSE(mgr.define(make("outlet", {20, 21})));
    REQUIRE(mgr.size() == 2);

    REQUIRE(mgr.find("inlet") != nullptr);
    REQUIRE(mgr.find("outlet") != nullptr);
    REQUIRE(mgr.find("walls") == nullptr);

    // Re-define overwrites.
    REQUIRE(mgr.define(make("inlet", {99})));
    REQUIRE(mgr.find("inlet")->ids.size() == 1);

    REQUIRE(mgr.rename("inlet", "primary_inlet"));
    REQUIRE(mgr.find("inlet") == nullptr);
    REQUIRE(mgr.find("primary_inlet") != nullptr);

    REQUIRE(mgr.remove("outlet"));
    REQUIRE_FALSE(mgr.remove("outlet"));
}

TEST_CASE("SelectionManager set algebra is correct", "[zones][selection]")
{
    SelectionManager mgr;
    mgr.define(make("a", {1, 2, 3, 4}));
    mgr.define(make("b", {3, 4, 5}));

    auto u = mgr.set_union("a", "b", "u");
    auto i = mgr.set_intersect("a", "b", "i");
    auto d = mgr.set_subtract("a", "b", "d");

    REQUIRE(u.has_value());
    REQUIRE(i.has_value());
    REQUIRE(d.has_value());
    REQUIRE(u->ids.size() == 5);
    REQUIRE(i->ids.size() == 2);
    REQUIRE(d->ids.size() == 2);
}

TEST_CASE("SelectionManager JSON round-trip is stable", "[zones][selection]")
{
    SelectionManager a;
    a.define(make("inlet", {10, 11, 12}));
    a.define(make("outlet", {20, 21}));
    const std::string s = a.to_json();

    SelectionManager b;
    REQUIRE(b.from_json(s));
    REQUIRE(b.size() == 2);
    REQUIRE(b.find("inlet")->ids.size() == 3);
    REQUIRE(b.find("outlet")->ids.size() == 2);
    REQUIRE(b.to_json() == s);
}
