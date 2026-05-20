// =============================================================================
// SimAll Beta — IO Unit Tests
// File   : tests/unit/io/test_gmsh_reader.cpp
// Phase  : 23 Pass 1 (v4.1 coverage); Pass 6 added v2.2 coverage
//
// Coverage for the Gmsh ASCII reader.  All tests parse in-memory
// strings via parse_gmsh_msh_string() so the suite has zero disk I/O and
// runs in milliseconds.  Covers v4.1 entity-block grammar and legacy
// v2.2 flat-list grammar.
// =============================================================================
#include "io/GmshReader.hpp"
#include "io/MeshFormats.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>

using simall::io::ElementType;
using simall::io::gmsh_element_type;
using simall::io::GmshReadResult;
using simall::io::parse_gmsh_msh_string;
using simall::io::vertices_per_element;

namespace
{

// Minimal valid Gmsh v4.1 mesh: one tetrahedron (4 nodes, 1 tet, 4 tri faces)
// inside a single physical volume "fluid" and three physical surfaces
// "inlet", "outlet", "wall".  We embed entities so the boundary-patch
// derivation has metadata to consume.
constexpr const char* kSimpleTet = R"GMSH(
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
4
2 1 "inlet"
2 2 "outlet"
2 3 "wall"
3 100 "fluid"
$EndPhysicalNames
$Entities
0 0 4 1
1 0 0 0 1 1 1 1 1 0
2 0 0 0 1 1 1 1 2 0
3 0 0 0 1 1 1 1 3 0
4 0 0 0 1 1 1 0 0
1 0 0 0 1 1 1 1 100 4 1 2 3 4
$EndEntities
$Nodes
1 4 1 4
3 4 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
5 5 1 5
2 1 2 1
1 1 2 3
2 2 2 1
2 1 2 4
2 3 2 1
3 1 3 4
2 4 2 1
4 2 3 4
3 1 4 1
5 1 2 3 4
$EndElements
)GMSH";

} // namespace

TEST_CASE("gmsh_element_type translates supported codes", "[io][gmsh]")
{
    CHECK(gmsh_element_type(1) == ElementType::Bar2);
    CHECK(gmsh_element_type(2) == ElementType::Tri3);
    CHECK(gmsh_element_type(3) == ElementType::Quad4);
    CHECK(gmsh_element_type(4) == ElementType::Tetra4);
    CHECK(gmsh_element_type(5) == ElementType::Hexa8);
    CHECK(gmsh_element_type(6) == ElementType::Penta6);
    CHECK(gmsh_element_type(7) == ElementType::Pyra5);
    CHECK(gmsh_element_type(15) == ElementType::Unknown); // points: caller skips
    CHECK(gmsh_element_type(99) == ElementType::Unknown);
}

TEST_CASE("Gmsh reader parses a minimal tetrahedron mesh", "[io][gmsh]")
{
    const auto r = parse_gmsh_msh_string(kSimpleTet, "tet.msh");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);
    const auto& z = r.mesh.zones[0];

    CHECK(z.x.size() == 4);
    CHECK(z.y.size() == 4);
    CHECK(z.z.size() == 4);
    CHECK(z.x[0] == 0.0);
    CHECK(z.x[1] == 1.0);
    CHECK(z.y[2] == 1.0);
    CHECK(z.z[3] == 1.0);

    // 4 surface tri sections + 1 volume tet section = 5 sections.
    REQUIRE(z.sections.size() == 5);

    std::size_t triSecs = 0, tetSecs = 0;
    std::size_t triElems = 0, tetElems = 0;
    for (const auto& s : z.sections) {
        if (s.type == ElementType::Tri3) {
            ++triSecs;
            triElems += s.element_count();
            CHECK(s.nodes.size() == s.element_count() * 3);
        } else if (s.type == ElementType::Tetra4) {
            ++tetSecs;
            tetElems += s.element_count();
            CHECK(s.nodes.size() == s.element_count() * 4);
        }
    }
    CHECK(triSecs == 4);
    CHECK(tetSecs == 1);
    CHECK(triElems == 4);
    CHECK(tetElems == 1);

    // ImportedMesh aggregates.
    CHECK(r.mesh.total_nodes() == 4);
    CHECK(r.mesh.total_elements() == 5);
}

TEST_CASE("Gmsh reader derives boundary patches from physical surfaces", "[io][gmsh]")
{
    const auto r = parse_gmsh_msh_string(kSimpleTet, "tet.msh");
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];

    // Three physical surfaces were declared ("inlet", "outlet", "wall");
    // boundaries are sorted alphabetically.
    REQUIRE(z.boundaries.size() == 3);
    CHECK(z.boundaries[0].name == "inlet");
    CHECK(z.boundaries[1].name == "outlet");
    CHECK(z.boundaries[2].name == "wall");
    for (const auto& bp : z.boundaries) {
        CHECK(bp.bcType == "wall"); // default until caller remaps
        CHECK_FALSE(bp.faceElementIndices.empty());
        // Every referenced section must be a Tri3 (surface) section.
        for (auto si : bp.faceElementIndices) {
            REQUIRE(si < z.sections.size());
            CHECK(z.sections[si].type == ElementType::Tri3);
        }
    }
}

TEST_CASE("Gmsh reader rejects binary MSH files", "[io][gmsh][error]")
{
    const std::string bin = "$MeshFormat\n4.1 1 8\n$EndMeshFormat\n";
    const auto r = parse_gmsh_msh_string(bin);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("binary") != std::string::npos);
}

TEST_CASE("Gmsh reader rejects unsupported format versions", "[io][gmsh][error]")
{
    // v3.x was never an official release; v5.x is future and unknown.
    const std::string v3 = "$MeshFormat\n3.0 0 8\n$EndMeshFormat\n";
    const auto r3 = parse_gmsh_msh_string(v3);
    CHECK_FALSE(r3.ok);
    CHECK(r3.error.find("unsupported Gmsh format version") != std::string::npos);

    const std::string v5 = "$MeshFormat\n5.0 0 8\n$EndMeshFormat\n";
    const auto r5 = parse_gmsh_msh_string(v5);
    CHECK_FALSE(r5.ok);
    CHECK(r5.error.find("unsupported Gmsh format version") != std::string::npos);
}

TEST_CASE("Gmsh reader rejects missing $MeshFormat", "[io][gmsh][error]")
{
    const std::string bad = "$Nodes\n0 0 0 0\n$EndNodes\n";
    const auto r = parse_gmsh_msh_string(bad);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("$MeshFormat") != std::string::npos);
}

TEST_CASE("Gmsh reader reports unknown node tags", "[io][gmsh][error]")
{
    // Node 5 referenced by tet but never defined.
    const std::string bad = R"GMSH(
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 5
$EndElements
)GMSH";
    const auto r = parse_gmsh_msh_string(bad);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unknown node tag") != std::string::npos);
}

TEST_CASE("Gmsh reader skips type-15 vertex elements without failing", "[io][gmsh]")
{
    // One tet + one vertex element (type 15) attached to node 1.
    const std::string txt = R"GMSH(
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
2 2 1 2
0 1 15 1
1 1
3 1 4 1
2 1 2 3 4
$EndElements
)GMSH";
    const auto r = parse_gmsh_msh_string(txt);
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);
    // Only the tet section survives -- the vertex element is skipped.
    REQUIRE(r.mesh.zones[0].sections.size() == 1);
    CHECK(r.mesh.zones[0].sections[0].type == ElementType::Tetra4);
}

TEST_CASE("Gmsh reader handles 2D quad mesh and reports dimension==2", "[io][gmsh]")
{
    const std::string txt = R"GMSH(
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
2 1 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
$EndNodes
$Elements
1 1 1 1
2 1 3 1
1 1 2 3 4
$EndElements
)GMSH";
    const auto r = parse_gmsh_msh_string(txt);
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 2);
    REQUIRE(r.mesh.zones.size() == 1);
    REQUIRE(r.mesh.zones[0].sections.size() == 1);
    CHECK(r.mesh.zones[0].sections[0].type == ElementType::Quad4);
    CHECK(r.mesh.zones[0].sections[0].element_count() == 1);
}

TEST_CASE("Gmsh reader rejects high-order element type codes", "[io][gmsh][error]")
{
    // Type 9 = 6-node second-order triangle -- not yet supported.
    const std::string txt = R"GMSH(
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 6 1 6
2 1 0 6
1
2
3
4
5
6
0 0 0
1 0 0
0 1 0
0.5 0 0
0.5 0.5 0
0 0.5 0
$EndNodes
$Elements
1 1 1 1
2 1 9 1
1 1 2 3 4 5 6
$EndElements
)GMSH";
    const auto r = parse_gmsh_msh_string(txt);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unsupported Gmsh element type") != std::string::npos);
}

// =============================================================================
// Phase 23 Pass 6: legacy v2.2 ASCII coverage.
//
// The v2.2 grammar is a flat-list format:
//     $Nodes
//         numNodes
//         tag x y z
//     $EndNodes
//     $Elements
//         numElements
//         elm-num elm-type num-tags <tags...> <nodes...>
//     $EndElements
// Each element row carries its physical-group tag inline (tags[0]) so the
// reader synthesises one entity record per (dim, physTag) for shared
// boundary-derivation downstream.
// =============================================================================

namespace
{

// Minimal valid Gmsh v2.2 mesh: one tetrahedron (4 nodes, 1 tet, 4 tri faces)
// with physical groups 1=inlet, 2=outlet, 3=wall (surfaces) and 100=fluid
// (volume).  Each element row's first tag is the physical group.
constexpr const char* kSimpleTetV2 = R"GMSH(
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
4
2 1 "inlet"
2 2 "outlet"
2 3 "wall"
3 100 "fluid"
$EndPhysicalNames
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
5
1 2 2 1 1 1 2 3
2 2 2 2 2 1 2 4
3 2 2 3 3 1 3 4
4 2 2 3 4 2 3 4
5 4 2 100 1 1 2 3 4
$EndElements
)GMSH";

} // namespace

TEST_CASE("Gmsh v2.2 reader parses a minimal tetrahedron mesh", "[io][gmsh][v22]")
{
    const auto r = parse_gmsh_msh_string(kSimpleTetV2, "tet.msh2");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);
    const auto& z = r.mesh.zones[0];

    CHECK(z.x.size() == 4);
    CHECK(z.y.size() == 4);
    CHECK(z.z.size() == 4);
    CHECK(z.x[0] == 0.0);
    CHECK(z.x[1] == 1.0);
    CHECK(z.y[2] == 1.0);
    CHECK(z.z[3] == 1.0);

    std::size_t triSecs = 0, tetSecs = 0;
    std::size_t triElems = 0, tetElems = 0;
    for (const auto& s : z.sections) {
        if (s.type == ElementType::Tri3) {
            ++triSecs;
            triElems += s.element_count();
            CHECK(s.nodes.size() == s.element_count() * 3);
        } else if (s.type == ElementType::Tetra4) {
            ++tetSecs;
            tetElems += s.element_count();
            CHECK(s.nodes.size() == s.element_count() * 4);
        }
    }
    // 3 unique surface physical tags -> 3 Tri3 sections (group "wall" has
    // two tris -- both land in the same section).
    CHECK(triSecs == 3);
    CHECK(tetSecs == 1);
    CHECK(triElems == 4);
    CHECK(tetElems == 1);

    CHECK(r.mesh.total_nodes() == 4);
    CHECK(r.mesh.total_elements() == 5);
}

TEST_CASE("Gmsh v2.2 reader derives boundary patches from physical surfaces", "[io][gmsh][v22]")
{
    const auto r = parse_gmsh_msh_string(kSimpleTetV2, "tet.msh2");
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];

    REQUIRE(z.boundaries.size() == 3);
    CHECK(z.boundaries[0].name == "inlet");
    CHECK(z.boundaries[1].name == "outlet");
    CHECK(z.boundaries[2].name == "wall");
    for (const auto& bp : z.boundaries) {
        CHECK(bp.bcType == "wall");
        CHECK_FALSE(bp.faceElementIndices.empty());
        for (auto si : bp.faceElementIndices) {
            REQUIRE(si < z.sections.size());
            CHECK(z.sections[si].type == ElementType::Tri3);
        }
    }
}

TEST_CASE("Gmsh v2.2 reader handles 2D quad mesh and reports dimension==2", "[io][gmsh][v22]")
{
    const std::string txt = R"GMSH(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
$EndNodes
$Elements
1
1 3 2 1 1 1 2 3 4
$EndElements
)GMSH";
    const auto r = parse_gmsh_msh_string(txt);
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 2);
    REQUIRE(r.mesh.zones.size() == 1);
    REQUIRE(r.mesh.zones[0].sections.size() == 1);
    CHECK(r.mesh.zones[0].sections[0].type == ElementType::Quad4);
    CHECK(r.mesh.zones[0].sections[0].element_count() == 1);
}

TEST_CASE("Gmsh v2.2 reader skips type-15 vertex elements without failing", "[io][gmsh][v22]")
{
    const std::string txt = R"GMSH(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
2
1 15 2 0 1 1
2 4 2 100 1 1 2 3 4
$EndElements
)GMSH";
    const auto r = parse_gmsh_msh_string(txt);
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);
    REQUIRE(r.mesh.zones[0].sections.size() == 1);
    CHECK(r.mesh.zones[0].sections[0].type == ElementType::Tetra4);
}

TEST_CASE("Gmsh v2.2 reader rejects unknown node tags", "[io][gmsh][v22][error]")
{
    const std::string bad = R"GMSH(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
1
1 4 2 100 1 1 2 3 5
$EndElements
)GMSH";
    const auto r = parse_gmsh_msh_string(bad);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unknown node tag") != std::string::npos);
}

TEST_CASE("Gmsh v2.2 reader rejects truncated element rows", "[io][gmsh][v22][error]")
{
    // Tet declared but only 3 node tags supplied after the physical tag.
    const std::string bad = R"GMSH(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
1
1 4 2 100 1 1 2 3
$EndElements
)GMSH";
    const auto r = parse_gmsh_msh_string(bad);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("too few nodes") != std::string::npos);
}

TEST_CASE("Gmsh v2.2 reader still rejects binary file-type", "[io][gmsh][v22][error]")
{
    const std::string bin = "$MeshFormat\n2.2 1 8\n$EndMeshFormat\n";
    const auto r = parse_gmsh_msh_string(bin);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("binary") != std::string::npos);
}
