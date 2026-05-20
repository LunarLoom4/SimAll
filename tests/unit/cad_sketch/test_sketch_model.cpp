// =============================================================================
// SimAll Beta -- Sketcher unit tests
// File   : tests/unit/cad_sketch/test_sketch_model.cpp
// Phase  : 23 Pass 23.1
// =============================================================================
#include "cad_sketch/Sketch.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace simall::cad::sketch;
using Catch::Matchers::WithinAbs;

TEST_CASE("Sketch creates parameters with sequential ids and retrievable values",
          "[cad][sketch][model]") {
    Sketch s;
    auto a = s.add_parameter(1.5);
    auto b = s.add_parameter(2.5, /*fixed=*/true, "x_pin");
    REQUIRE(a != b);
    REQUIRE(s.parameter_count() == 2);
    REQUIRE_THAT(s.get_parameter(a), WithinAbs(1.5, 1e-15));
    REQUIRE(s.parameter(b).fixed);
    REQUIRE(s.parameter(b).name == "x_pin");
}

TEST_CASE("Adding a point allocates exactly two parameters",
          "[cad][sketch][model]") {
    Sketch s;
    const auto pid = s.add_point({3.0, 4.0});
    const auto& e  = s.entity(pid);
    REQUIRE(e.kind == EntityKind::Point);
    REQUIRE(e.params.size() == 2);
    const auto v = s.point_value(pid);
    REQUIRE_THAT(v.x, WithinAbs(3.0, 1e-15));
    REQUIRE_THAT(v.y, WithinAbs(4.0, 1e-15));
}

TEST_CASE("Adding a line allocates four parameters and exposes endpoints",
          "[cad][sketch][model]") {
    Sketch s;
    const auto lid = s.add_line({0.0, 0.0}, {1.0, 2.0});
    REQUIRE(s.entity(lid).kind == EntityKind::Line);
    REQUIRE(s.entity(lid).params.size() == 4);
    const auto a = s.line_endpoint(lid, 0);
    const auto b = s.line_endpoint(lid, 1);
    REQUIRE_THAT(a.x, WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(b.y, WithinAbs(2.0, 1e-15));
}

TEST_CASE("Circle and arc store center and radius",
          "[cad][sketch][model]") {
    Sketch s;
    const auto cid = s.add_circle({1.0, 1.0}, 0.5);
    REQUIRE(s.entity(cid).kind == EntityKind::Circle);
    REQUIRE(s.entity(cid).params.size() == 3);
    REQUIRE_THAT(s.circle_radius(cid), WithinAbs(0.5, 1e-15));

    const auto aid = s.add_arc({0.0, 0.0}, 2.0, 0.0, 3.14);
    REQUIRE(s.entity(aid).kind == EntityKind::Arc);
    REQUIRE(s.entity(aid).params.size() == 5);
}

TEST_CASE("Spline stores 2*N parameters",
          "[cad][sketch][model]") {
    Sketch s;
    const auto sid = s.add_spline({{0,0}, {1,1}, {2,0}, {3,1}});
    REQUIRE(s.entity(sid).kind == EntityKind::Spline);
    REQUIRE(s.entity(sid).params.size() == 8);
}

TEST_CASE("Removing an entity drops its dependent constraints",
          "[cad][sketch][model]") {
    Sketch s;
    const auto p1 = s.add_point({0, 0});
    const auto p2 = s.add_point({1, 0});
    s.add_distance(p1, p2, 1.0);
    REQUIRE(s.constraints().size() == 1);
    REQUIRE(s.remove_entity(p2));
    REQUIRE(s.constraints().empty());
}

TEST_CASE("Fix constraint pins the underlying point parameters",
          "[cad][sketch][model]") {
    Sketch s;
    const auto p = s.add_point({7.0, 8.0});
    s.add_fix(p);
    const auto& ent = s.entity(p);
    for (auto pid : ent.params) REQUIRE(s.parameter(pid).fixed);
}
