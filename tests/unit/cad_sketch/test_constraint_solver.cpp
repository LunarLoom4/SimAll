// =============================================================================
// SimAll Beta -- Sketcher unit tests
// File   : tests/unit/cad_sketch/test_constraint_solver.cpp
// Phase  : 23 Pass 23.1
// =============================================================================
#include "cad_sketch/ConstraintSolver.hpp"
#include "cad_sketch/Sketch.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace simall::cad::sketch;
using Catch::Matchers::WithinAbs;

namespace {
constexpr double kTol = 1e-7;
}

TEST_CASE("Solver returns Empty when there are no constraints",
          "[cad][sketch][solver]") {
    Sketch s;
    s.add_point({0, 0});
    const auto rep = solve(s);
    REQUIRE(rep.status == SolveStatus::Empty);
}

TEST_CASE("Coincident constraint pulls two points together",
          "[cad][sketch][solver]") {
    Sketch s;
    const auto p1 = s.add_point({0.0, 0.0});
    const auto p2 = s.add_point({3.0, 4.0});
    s.add_fix(p1);                  // anchor p1 at origin
    s.add_coincident(p1, p2);

    const auto rep = solve(s);
    REQUIRE(rep.ok());
    const auto v = s.point_value(p2);
    REQUIRE_THAT(v.x, WithinAbs(0.0, kTol));
    REQUIRE_THAT(v.y, WithinAbs(0.0, kTol));
}

TEST_CASE("Distance constraint forces specified separation",
          "[cad][sketch][solver]") {
    Sketch s;
    const auto p1 = s.add_point({0.0, 0.0});
    const auto p2 = s.add_point({0.1, 0.0});   // start nearby to test pull
    s.add_fix(p1);
    s.add_distance(p1, p2, 5.0);

    const auto rep = solve(s);
    REQUIRE(rep.ok());
    REQUIRE_THAT(distance(s.point_value(p1), s.point_value(p2)),
                 WithinAbs(5.0, kTol));
}

TEST_CASE("Horizontal constraint zeros out a line's dy",
          "[cad][sketch][solver]") {
    Sketch s;
    const auto l = s.add_line({0.0, 0.0}, {1.0, 0.7});
    s.add_horizontal(l);

    const auto rep = solve(s);
    REQUIRE(rep.ok());
    const auto a = s.line_endpoint(l, 0);
    const auto b = s.line_endpoint(l, 1);
    REQUIRE_THAT(a.y - b.y, WithinAbs(0.0, kTol));
}

TEST_CASE("Parallel constraint aligns line directions",
          "[cad][sketch][solver]") {
    Sketch s;
    const auto l1 = s.add_line({0.0, 0.0}, {1.0, 0.0});       // pinned reference
    const auto l2 = s.add_line({0.0, 1.0}, {1.0, 1.2});       // tilted
    // Pin l1 by fixing both endpoints' parameters via Fix on point-equivalents:
    for (auto pid : s.entity(l1).params) s.fix_parameter(pid);
    s.add_parallel(l1, l2);

    const auto rep = solve(s);
    REQUIRE(rep.ok());
    const auto d1 = s.line_endpoint(l1, 1) - s.line_endpoint(l1, 0);
    const auto d2 = s.line_endpoint(l2, 1) - s.line_endpoint(l2, 0);
    REQUIRE_THAT(d1.cross(d2), WithinAbs(0.0, kTol));
}

TEST_CASE("Perpendicular constraint enforces 90 degrees",
          "[cad][sketch][solver]") {
    Sketch s;
    const auto l1 = s.add_line({0.0, 0.0}, {1.0, 0.0});
    const auto l2 = s.add_line({0.0, 0.0}, {1.0, 0.1});
    for (auto pid : s.entity(l1).params) s.fix_parameter(pid);
    // Pin l2's first endpoint at origin so only its tip is free.
    s.fix_parameter(s.entity(l2).params[0]);
    s.fix_parameter(s.entity(l2).params[1]);
    s.add_perpendicular(l1, l2);

    const auto rep = solve(s);
    REQUIRE(rep.ok());
    const auto d1 = s.line_endpoint(l1, 1) - s.line_endpoint(l1, 0);
    const auto d2 = s.line_endpoint(l2, 1) - s.line_endpoint(l2, 0);
    REQUIRE_THAT(d1.dot(d2), WithinAbs(0.0, kTol));
}

TEST_CASE("Radius constraint sets a circle's radius",
          "[cad][sketch][solver]") {
    Sketch s;
    const auto c = s.add_circle({0.0, 0.0}, 1.0);
    s.add_radius(c, 3.7);
    const auto rep = solve(s);
    REQUIRE(rep.ok());
    REQUIRE_THAT(s.circle_radius(c), WithinAbs(3.7, kTol));
}

TEST_CASE("Angle constraint sets angle between two lines",
          "[cad][sketch][solver]") {
    Sketch s;
    const auto l1 = s.add_line({0.0, 0.0}, {1.0, 0.0});
    const auto l2 = s.add_line({0.0, 0.0}, {1.0, 0.1});
    for (auto pid : s.entity(l1).params) s.fix_parameter(pid);
    s.fix_parameter(s.entity(l2).params[0]);
    s.fix_parameter(s.entity(l2).params[1]);
    const double target = 3.14159265358979323846 / 3.0;   // 60 degrees
    s.add_angle(l1, l2, target);

    const auto rep = solve(s);
    REQUIRE(rep.ok());
    const auto d1 = s.line_endpoint(l1, 1) - s.line_endpoint(l1, 0);
    const auto d2 = s.line_endpoint(l2, 1) - s.line_endpoint(l2, 0);
    const double cos_actual = d1.dot(d2) / (d1.norm() * d2.norm());
    REQUIRE_THAT(cos_actual, WithinAbs(std::cos(target), 1e-6));
}

TEST_CASE("Point-on-line constraint forces collinearity",
          "[cad][sketch][solver]") {
    Sketch s;
    const auto l = s.add_line({0.0, 0.0}, {1.0, 0.0});
    for (auto pid : s.entity(l).params) s.fix_parameter(pid);
    const auto p = s.add_point({0.5, 0.4});
    s.fix_parameter(s.entity(p).params[0]);   // free only y
    s.add_point_on_line(p, l);

    const auto rep = solve(s);
    REQUIRE(rep.ok());
    REQUIRE_THAT(s.point_value(p).y, WithinAbs(0.0, kTol));
}

TEST_CASE("Solver reports under-constrained sketch via DOF > 0",
          "[cad][sketch][solver][degeneracy]") {
    // Two free points with only a distance constraint -- one constraint,
    // four unknowns -> rank deficient (multiple solutions).
    Sketch s;
    const auto p1 = s.add_point({0.0, 0.0});
    const auto p2 = s.add_point({1.0, 0.0});
    s.add_distance(p1, p2, 1.0);

    const auto rep = solve(s);
    REQUIRE(rep.unknowns == 4);
    REQUIRE(rep.dof > 0);
    REQUIRE_FALSE(rep.degenerate_parameters.empty());
}

TEST_CASE("Solver reports over-constrained sketch when contradictory",
          "[cad][sketch][solver][degeneracy]") {
    // Same pair of points with two contradictory distance constraints.
    Sketch s;
    const auto p1 = s.add_point({0.0, 0.0});
    const auto p2 = s.add_point({1.0, 0.0});
    s.add_fix(p1);
    s.add_distance(p1, p2, 1.0);
    s.add_distance(p1, p2, 2.0);   // contradicts the first

    const auto rep = solve(s);
    REQUIRE(rep.status != SolveStatus::Converged);
    REQUIRE(rep.residual_norm > 1e-3);
}

TEST_CASE("evaluate_residuals returns zero vector when sketch is already satisfied",
          "[cad][sketch][solver]") {
    Sketch s;
    const auto p1 = s.add_point({0.0, 0.0});
    const auto p2 = s.add_point({3.0, 4.0});
    s.add_distance(p1, p2, 5.0);   // already satisfied
    const auto r = evaluate_residuals(s);
    REQUIRE(r.size() == 1);
    REQUIRE_THAT(r[0], WithinAbs(0.0, 1e-15));
}
