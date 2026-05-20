// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/cad/test_topology_graph.cpp
//
// TopologyGraph is OCC-free, so it can be exercised without an OpenCASCADE
// dependency. The cpp is compiled directly into the test binary so the
// tests do not transitively link SimAll::Cad → OpenCASCADE.
// =============================================================================
#include "cad/TopologyGraph.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace simall::cad;

TEST_CASE("TopologyGraph add/find/by_type", "[cad][topology]")
{
    TopologyGraph g;
    char dummy1, dummy2, dummy3;

    auto f1 = g.add(TopologyType::Face, &dummy1);
    auto f2 = g.add(TopologyType::Face, &dummy2);
    auto e1 = g.add(TopologyType::Edge, &dummy3);

    REQUIRE(g.size() == 3);
    REQUIRE(g.find(f1) != nullptr);
    REQUIRE(g.find(f1)->type == TopologyType::Face);
    REQUIRE(g.find(util::PersistentId{99999}) == nullptr);

    auto faces = g.by_type(TopologyType::Face);
    REQUIRE(faces.size() == 2);
    auto edges = g.by_type(TopologyType::Edge);
    REQUIRE(edges.size() == 1);
    REQUIRE(edges.front() == e1);
}

TEST_CASE("TopologyGraph link establishes parent/child", "[cad][topology]")
{
    TopologyGraph g;
    int a, b;
    auto solid = g.add(TopologyType::Solid, &a);
    auto face = g.add(TopologyType::Face, &b);

    g.link(solid, face);

    auto* p = g.find(solid);
    auto* c = g.find(face);
    REQUIRE(p);
    REQUIRE(c);
    REQUIRE(p->children.size() == 1);
    REQUIRE(p->children.front() == face);
    REQUIRE(c->parents.size() == 1);
    REQUIRE(c->parents.front() == solid);
}
