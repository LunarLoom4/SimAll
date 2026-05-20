// =============================================================================
// SimAll Beta - IO Unit Tests
// File   : tests/unit/io/test_gambit_reader.cpp
// Phase  : 23 Pass 7
//
// Coverage for the GAMBIT neutral file (.neu) reader.  All tests parse
// in-memory strings via parse_gambit_neu_string() so the suite has zero
// disk I/O.
// =============================================================================
#include "io/GambitReader.hpp"
#include "io/MeshFormats.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>

using simall::io::ElementType;
using simall::io::gambit_element_type;
using simall::io::GambitReadResult;
using simall::io::parse_gambit_neu_string;
using simall::io::vertices_per_element;

namespace
{

// Minimal 3D single-tet mesh: 4 nodes, 1 tet (group "fluid"), 3 BCs
// inlet/outlet/wall covering all four faces (wall holds two).
//
// Tet face conventions used by the reader (see GambitReader.cpp):
//   F1: 1,2,3   F2: 1,2,4   F3: 2,3,4   F4: 1,3,4
constexpr const char* kSimpleTet = R"NEU(
        CONTROL INFO 2.4.6
** GAMBIT NEUTRAL FILE
single-tet test case
PROGRAM:                Gambit     VERSION:  2.4.6
19 May 2026    12:00:00
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         1         3         3         3
ENDOFSECTION
   NODAL COORDINATES 2.4.6
     1   0.0   0.0   0.0
     2   1.0   0.0   0.0
     3   0.0   1.0   0.0
     4   0.0   0.0   1.0
ENDOFSECTION
      ELEMENTS/CELLS 2.4.6
        1  6  4  1  2  3  4
ENDOFSECTION
       ELEMENT GROUP 2.4.6
GROUP:    1 ELEMENTS:    1 MATERIAL:    2 NFLAGS:    1
fluid
       0
       1
ENDOFSECTION
 BOUNDARY CONDITIONS 2.4.6
                  inlet       1       1       0       6
       1       6       1
ENDOFSECTION
 BOUNDARY CONDITIONS 2.4.6
                 outlet       1       1       0       6
       1       6       3
ENDOFSECTION
 BOUNDARY CONDITIONS 2.4.6
                   wall       1       2       0       6
       1       6       2
       1       6       4
ENDOFSECTION
)NEU";

} // namespace

TEST_CASE("gambit_element_type translates supported (code, ndp) pairs", "[io][gambit]")
{
    CHECK(gambit_element_type(1, 2) == ElementType::Bar2);
    CHECK(gambit_element_type(2, 4) == ElementType::Quad4);
    CHECK(gambit_element_type(3, 3) == ElementType::Tri3);
    CHECK(gambit_element_type(4, 8) == ElementType::Hexa8);
    CHECK(gambit_element_type(5, 6) == ElementType::Penta6);
    CHECK(gambit_element_type(6, 4) == ElementType::Tetra4);
    CHECK(gambit_element_type(7, 5) == ElementType::Pyra5);
    // Wrong ndp -> Unknown (high-order detection).
    CHECK(gambit_element_type(6, 10) == ElementType::Unknown);
    CHECK(gambit_element_type(4, 20) == ElementType::Unknown);
    CHECK(gambit_element_type(99, 1) == ElementType::Unknown);
}

TEST_CASE("GAMBIT reader parses a minimal tet mesh", "[io][gambit]")
{
    const auto r = parse_gambit_neu_string(kSimpleTet, "tet.neu");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);
    const auto& z = r.mesh.zones[0];

    CHECK(z.x.size() == 4);
    CHECK(z.y.size() == 4);
    CHECK(z.z.size() == 4);
    CHECK(z.x[1] == 1.0);
    CHECK(z.y[2] == 1.0);
    CHECK(z.z[3] == 1.0);

    // One Tetra4 volume section + three Tri3 BC face sections.
    REQUIRE(z.sections.size() == 4);

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
    CHECK(triSecs == 3); // inlet / outlet / wall face sections
    CHECK(tetSecs == 1);
    CHECK(triElems == 4); // 1 + 1 + 2
    CHECK(tetElems == 1);

    CHECK(r.mesh.total_nodes() == 4);
    CHECK(r.mesh.total_elements() == 5);
}

TEST_CASE("GAMBIT reader names volume section after element group", "[io][gambit]")
{
    const auto r = parse_gambit_neu_string(kSimpleTet, "tet.neu");
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];
    bool foundFluid = false;
    for (const auto& s : z.sections) {
        if (s.type == ElementType::Tetra4) {
            foundFluid = (s.name == "fluid_Tetra4");
            break;
        }
    }
    CHECK(foundFluid);
}

TEST_CASE("GAMBIT reader derives boundary patches sorted by name with "
          "correct face connectivity",
          "[io][gambit]")
{
    const auto r = parse_gambit_neu_string(kSimpleTet, "tet.neu");
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];

    REQUIRE(z.boundaries.size() == 3);
    CHECK(z.boundaries[0].name == "inlet");
    CHECK(z.boundaries[1].name == "outlet");
    CHECK(z.boundaries[2].name == "wall");
    for (const auto& bp : z.boundaries) {
        CHECK(bp.bcType == "wall");
        REQUIRE(bp.faceElementIndices.size() == 1);
        const auto si = bp.faceElementIndices[0];
        REQUIRE(si < z.sections.size());
        CHECK(z.sections[si].type == ElementType::Tri3);
    }

    // inlet patch -> face 1 of the tet -> local nodes 1,2,3 -> NodeIdx 0,1,2
    {
        const auto si = z.boundaries[0].faceElementIndices[0];
        const auto& sec = z.sections[si];
        REQUIRE(sec.element_count() == 1);
        REQUIRE(sec.nodes.size() == 3);
        CHECK(sec.nodes[0] == 0);
        CHECK(sec.nodes[1] == 1);
        CHECK(sec.nodes[2] == 2);
    }
    // wall patch -> faces 2 and 4 -> two Tri3 elements (1,2,4) + (1,3,4)
    {
        const auto si = z.boundaries[2].faceElementIndices[0];
        const auto& sec = z.sections[si];
        REQUIRE(sec.element_count() == 2);
        REQUIRE(sec.nodes.size() == 6);
        CHECK(sec.nodes[0] == 0);
        CHECK(sec.nodes[1] == 1);
        CHECK(sec.nodes[2] == 3);
        CHECK(sec.nodes[3] == 0);
        CHECK(sec.nodes[4] == 2);
        CHECK(sec.nodes[5] == 3);
    }
}

TEST_CASE("GAMBIT reader handles 2D quad mesh (NDFCD=2)", "[io][gambit]")
{
    const std::string txt = R"NEU(
        CONTROL INFO 2.4.6
2D test
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         0         0         2         2
ENDOFSECTION
   NODAL COORDINATES 2.4.6
     1   0.0   0.0
     2   1.0   0.0
     3   1.0   1.0
     4   0.0   1.0
ENDOFSECTION
      ELEMENTS/CELLS 2.4.6
        1  2  4  1  2  3  4
ENDOFSECTION
)NEU";
    const auto r = parse_gambit_neu_string(txt);
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 2);
    REQUIRE(r.mesh.zones.size() == 1);
    REQUIRE(r.mesh.zones[0].sections.size() == 1);
    CHECK(r.mesh.zones[0].sections[0].type == ElementType::Quad4);
    CHECK(r.mesh.zones[0].sections[0].element_count() == 1);
    // For 2D NDFCD=2, z is filled with zero.
    CHECK(r.mesh.zones[0].z[0] == 0.0);
    CHECK(r.mesh.zones[0].z[3] == 0.0);
}

TEST_CASE("GAMBIT reader emits per-face-type sections for mixed-face BCs "
          "on a single wedge",
          "[io][gambit]")
{
    // Single wedge (Penta6) with one BC patch listing one quad face (F1)
    // and one tri face (F4) -> patch gets 2 face sections (Quad4 + Tri3).
    const std::string txt = R"NEU(
        CONTROL INFO 2.4.6
wedge-mixed-bc
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         6         1         0         1         3         3
ENDOFSECTION
   NODAL COORDINATES 2.4.6
     1   0.0   0.0   0.0
     2   1.0   0.0   0.0
     3   0.0   1.0   0.0
     4   0.0   0.0   1.0
     5   1.0   0.0   1.0
     6   0.0   1.0   1.0
ENDOFSECTION
      ELEMENTS/CELLS 2.4.6
        1  5  6  1  2  3  4  5  6
ENDOFSECTION
 BOUNDARY CONDITIONS 2.4.6
                mixed       1       2       0       6
       1       5       1
       1       5       4
ENDOFSECTION
)NEU";
    const auto r = parse_gambit_neu_string(txt);
    INFO(r.error);
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];
    REQUIRE(z.boundaries.size() == 1);
    const auto& bp = z.boundaries[0];
    CHECK(bp.name == "mixed");
    REQUIRE(bp.faceElementIndices.size() == 2);
    // Collect face types referenced by the patch.
    bool sawTri = false, sawQuad = false;
    for (auto si : bp.faceElementIndices) {
        REQUIRE(si < z.sections.size());
        const auto& sec = z.sections[si];
        if (sec.type == ElementType::Tri3) {
            sawTri = true;
            CHECK(sec.element_count() == 1);
        }
        if (sec.type == ElementType::Quad4) {
            sawQuad = true;
            CHECK(sec.element_count() == 1);
        }
    }
    CHECK(sawTri);
    CHECK(sawQuad);
}

TEST_CASE("GAMBIT reader rejects unsupported element types", "[io][gambit][error]")
{
    // 10-node quadratic tetrahedron: type=6, ndp=10 (not supported).
    const std::string txt = R"NEU(
        CONTROL INFO 2.4.6
bad
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
        10         1         0         0         3         3
ENDOFSECTION
   NODAL COORDINATES 2.4.6
     1  0 0 0
     2  1 0 0
     3  0 1 0
     4  0 0 1
     5  0.5 0 0
     6  0.5 0.5 0
     7  0 0.5 0
     8  0 0 0.5
     9  0.5 0 0.5
    10  0 0.5 0.5
ENDOFSECTION
      ELEMENTS/CELLS 2.4.6
        1  6  10  1 2 3 4 5 6 7 8 9 10
ENDOFSECTION
)NEU";
    const auto r = parse_gambit_neu_string(txt);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unsupported element") != std::string::npos);
}

TEST_CASE("GAMBIT reader rejects unknown node tags in elements", "[io][gambit][error]")
{
    const std::string txt = R"NEU(
        CONTROL INFO 2.4.6
bad
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         0         0         3         3
ENDOFSECTION
   NODAL COORDINATES 2.4.6
     1   0.0   0.0   0.0
     2   1.0   0.0   0.0
     3   0.0   1.0   0.0
     4   0.0   0.0   1.0
ENDOFSECTION
      ELEMENTS/CELLS 2.4.6
        1  6  4  1  2  3  9
ENDOFSECTION
)NEU";
    const auto r = parse_gambit_neu_string(txt);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unknown") != std::string::npos);
}

TEST_CASE("GAMBIT reader rejects BC cell-type mismatch", "[io][gambit][error]")
{
    // Cell 1 is a Tetra4 (type 6) but the BC entry claims type 4 (Hexa8).
    const std::string txt = R"NEU(
        CONTROL INFO 2.4.6
bad
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         0         1         3         3
ENDOFSECTION
   NODAL COORDINATES 2.4.6
     1   0.0   0.0   0.0
     2   1.0   0.0   0.0
     3   0.0   1.0   0.0
     4   0.0   0.0   1.0
ENDOFSECTION
      ELEMENTS/CELLS 2.4.6
        1  6  4  1  2  3  4
ENDOFSECTION
 BOUNDARY CONDITIONS 2.4.6
                  bad       1       1       0       6
       1       4       1
ENDOFSECTION
)NEU";
    const auto r = parse_gambit_neu_string(txt);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("cell-type mismatch") != std::string::npos);
}

TEST_CASE("GAMBIT reader rejects out-of-range face id in BC", "[io][gambit][error]")
{
    const std::string txt = R"NEU(
        CONTROL INFO 2.4.6
bad
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         0         1         3         3
ENDOFSECTION
   NODAL COORDINATES 2.4.6
     1   0.0   0.0   0.0
     2   1.0   0.0   0.0
     3   0.0   1.0   0.0
     4   0.0   0.0   1.0
ENDOFSECTION
      ELEMENTS/CELLS 2.4.6
        1  6  4  1  2  3  4
ENDOFSECTION
 BOUNDARY CONDITIONS 2.4.6
                  bad       1       1       0       6
       1       6       9
ENDOFSECTION
)NEU";
    const auto r = parse_gambit_neu_string(txt);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("invalid face id") != std::string::npos);
}

TEST_CASE("GAMBIT reader requires CONTROL INFO header first", "[io][gambit][error]")
{
    const std::string txt = "   NODAL COORDINATES 2.4.6\n"
                            "     1   0.0   0.0   0.0\n"
                            "ENDOFSECTION\n";
    const auto r = parse_gambit_neu_string(txt);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("CONTROL INFO") != std::string::npos);
}

TEST_CASE("GAMBIT reader rejects file with no element section", "[io][gambit][error]")
{
    const std::string txt = R"NEU(
        CONTROL INFO 2.4.6
nodes-only
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         0         0         0         3         3
ENDOFSECTION
   NODAL COORDINATES 2.4.6
     1   0.0   0.0   0.0
     2   1.0   0.0   0.0
     3   0.0   1.0   0.0
     4   0.0   0.0   1.0
ENDOFSECTION
)NEU";
    const auto r = parse_gambit_neu_string(txt);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("ELEMENTS/CELLS") != std::string::npos);
}
