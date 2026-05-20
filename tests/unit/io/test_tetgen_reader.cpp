// =============================================================================
// SimAll Beta - IO Unit Tests
// File   : tests/unit/io/test_tetgen_reader.cpp
// Phase  : 23 Pass 8
//
// Coverage for the TetGen `.node/.ele/.face` triple-file reader.
// All tests parse in-memory strings via parse_tetgen_strings() so the
// suite has zero disk I/O.
// =============================================================================
#include <catch2/catch_test_macros.hpp>

#include "io/MeshFormats.hpp"
#include "io/TetgenReader.hpp"

#include <string>

using simall::io::ElementType;
using simall::io::parse_tetgen_strings;
using simall::io::TetgenReadResult;

namespace {

// Two-tet mesh sharing the (1,2,3) face.  Nodes are 1-based.
//   Tet 1: 1,2,3,4   Tet 2: 1,2,3,5
// .face lists the four outer faces of tet 1 and the three outer faces of
// tet 2 (the shared face is interior and not exported).
constexpr const char* kNode2tet = R"NODE(
# 2-tet test mesh -- nodes
5 3 0 0
1   0.0  0.0  0.0
2   1.0  0.0  0.0
3   0.0  1.0  0.0
4   0.0  0.0  1.0
5   0.0  0.0 -1.0
)NODE";

constexpr const char* kEle2tet = R"ELE(
2 4 0
1   1 2 3 4
2   1 2 3 5
)ELE";

// Faces with markers: 1 = wall, 2 = inlet
constexpr const char* kFace2tet = R"FACE(
7 1
1   1 2 4   1
2   2 3 4   1
3   1 3 4   1
4   1 2 5   2
5   2 3 5   1
6   1 3 5   1
7   1 2 3   1
)FACE";

}  // namespace

TEST_CASE("TetGen reader parses minimal 2-tet mesh", "[io][tetgen]") {
    const auto r = parse_tetgen_strings(kNode2tet, kEle2tet, "", "two_tet");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);
    const auto& z = r.mesh.zones[0];
    CHECK(z.x.size() == 5);
    CHECK(z.x[1] == 1.0);
    CHECK(z.z[3] == 1.0);
    CHECK(z.z[4] == -1.0);
    REQUIRE(z.sections.size() == 1);
    CHECK(z.sections[0].type == ElementType::Tetra4);
    CHECK(z.sections[0].name == "tetgen_volume");
    CHECK(z.sections[0].element_count() == 2);
    CHECK(z.sections[0].nodes.size() == 8);
    CHECK(z.boundaries.empty());

    CHECK(r.mesh.total_nodes()    == 5);
    CHECK(r.mesh.total_elements() == 2);
}

TEST_CASE("TetGen reader splits elements per attribute into separate sections",
          "[io][tetgen]") {
    const char* node = R"NODE(
4 3 0 0
1  0 0 0
2  1 0 0
3  0 1 0
4  0 0 1
)NODE";
    // Two tets each tagged with a different region attribute (10, 20).
    const char* ele = R"ELE(
2 4 1
1   1 2 3 4   10
2   2 3 4 1   20
)ELE";
    const auto r = parse_tetgen_strings(node, ele);
    INFO(r.error);
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];
    REQUIRE(z.sections.size() == 2);
    bool sawR10 = false, sawR20 = false;
    for (const auto& s : z.sections) {
        CHECK(s.type == ElementType::Tetra4);
        CHECK(s.element_count() == 1);
        if (s.name == "region_10") sawR10 = true;
        if (s.name == "region_20") sawR20 = true;
    }
    CHECK(sawR10);
    CHECK(sawR20);
}

TEST_CASE("TetGen reader groups face markers into separate boundary patches",
          "[io][tetgen]") {
    const auto r = parse_tetgen_strings(kNode2tet, kEle2tet,
                                         kFace2tet, "two_tet");
    INFO(r.error);
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];

    // 1 volume section + 2 face sections (one per marker) = 3 sections.
    REQUIRE(z.sections.size() == 3);
    // Patches: marker_1 (6 wall faces) and marker_2 (1 inlet face).
    REQUIRE(z.boundaries.size() == 2);
    CHECK(z.boundaries[0].name == "marker_1");
    CHECK(z.boundaries[1].name == "marker_2");

    for (const auto& bp : z.boundaries) {
        CHECK(bp.bcType == "wall");
        REQUIRE(bp.faceElementIndices.size() == 1);
        const auto si = bp.faceElementIndices[0];
        REQUIRE(si < z.sections.size());
        CHECK(z.sections[si].type == ElementType::Tri3);
    }

    // marker_1 should hold 6 Tri3 elements (18 node ids).
    {
        const auto si = z.boundaries[0].faceElementIndices[0];
        const auto& sec = z.sections[si];
        CHECK(sec.element_count() == 6);
        CHECK(sec.nodes.size() == 18);
    }
    // marker_2 should hold 1 Tri3 element (3 node ids).
    {
        const auto si = z.boundaries[1].faceElementIndices[0];
        const auto& sec = z.sections[si];
        CHECK(sec.element_count() == 1);
        REQUIRE(sec.nodes.size() == 3);
        // Face was "1 2 5" with 1-based tags -> 0-based NodeIdx 0,1,4.
        CHECK(sec.nodes[0] == 0);
        CHECK(sec.nodes[1] == 1);
        CHECK(sec.nodes[2] == 4);
    }
}

TEST_CASE("TetGen reader handles 0-based node numbering", "[io][tetgen]") {
    const char* node = R"NODE(
4 3 0 0
0  0 0 0
1  1 0 0
2  0 1 0
3  0 0 1
)NODE";
    const char* ele = R"ELE(
1 4 0
0   0 1 2 3
)ELE";
    const auto r = parse_tetgen_strings(node, ele);
    INFO(r.error);
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];
    CHECK(z.x.size() == 4);
    REQUIRE(z.sections.size() == 1);
    CHECK(z.sections[0].element_count() == 1);
    REQUIRE(z.sections[0].nodes.size() == 4);
    // 0-based tag 0 -> NodeIdx 0, 1 -> 1, etc.
    CHECK(z.sections[0].nodes[0] == 0);
    CHECK(z.sections[0].nodes[3] == 3);
}

TEST_CASE("TetGen reader tolerates `#` comments and blank lines",
          "[io][tetgen]") {
    const char* node = R"NODE(
# leading comment

4 3 0 0
# nodes follow
1  0 0 0    # origin
2  1 0 0
3  0 1 0
4  0 0 1
)NODE";
    const char* ele = R"ELE(
# header
1 4 0
1   1 2 3 4   # the only tet
)ELE";
    const auto r = parse_tetgen_strings(node, ele);
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.mesh.zones[0].sections[0].element_count() == 1);
}

TEST_CASE("TetGen reader handles attribute + marker columns in node file",
          "[io][tetgen]") {
    // 2 attributes and a boundary marker per node.  We don't surface them
    // in ImportedMesh but the parser must skip them without complaint.
    const char* node = R"NODE(
4 3 2 1
1  0 0 0   1.5  2.5   0
2  1 0 0   3.5  4.5   1
3  0 1 0   5.5  6.5   1
4  0 0 1   7.5  8.5   1
)NODE";
    const char* ele = R"ELE(
1 4 0
1   1 2 3 4
)ELE";
    const auto r = parse_tetgen_strings(node, ele);
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.mesh.zones[0].x.size() == 4);
}

TEST_CASE("TetGen reader rejects 10-node parabolic tets", "[io][tetgen][error]") {
    const char* node = R"NODE(
10 3 0 0
1  0 0 0
2  1 0 0
3  0 1 0
4  0 0 1
5  0.5 0 0
6  0.5 0.5 0
7  0 0.5 0
8  0 0 0.5
9  0.5 0 0.5
10 0 0.5 0.5
)NODE";
    const char* ele = R"ELE(
1 10 0
1   1 2 3 4  5 6 7 8 9 10
)ELE";
    const auto r = parse_tetgen_strings(node, ele);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("only linear") != std::string::npos);
}

TEST_CASE("TetGen reader rejects unknown node tags in elements",
          "[io][tetgen][error]") {
    const char* node = R"NODE(
4 3 0 0
1  0 0 0
2  1 0 0
3  0 1 0
4  0 0 1
)NODE";
    const char* ele = R"ELE(
1 4 0
1   1 2 3 9
)ELE";
    const auto r = parse_tetgen_strings(node, ele);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unknown node tag") != std::string::npos);
}

TEST_CASE("TetGen reader rejects duplicate node tags",
          "[io][tetgen][error]") {
    const char* node = R"NODE(
3 3 0 0
1  0 0 0
1  1 0 0
2  0 1 0
)NODE";
    const char* ele = R"ELE(
0 4 0
)ELE";
    const auto r = parse_tetgen_strings(node, ele);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("duplicate node tag") != std::string::npos);
}

TEST_CASE("TetGen reader rejects truncated rows", "[io][tetgen][error]") {
    const char* node = R"NODE(
4 3 0 0
1  0 0 0
2  1 0 0
3  0 1 0
4  0 0 1
)NODE";
    const char* ele = R"ELE(
1 4 0
1   1 2 3
)ELE";
    const auto r = parse_tetgen_strings(node, ele);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("row too short") != std::string::npos);
}

TEST_CASE("TetGen reader rejects malformed dimension in node header",
          "[io][tetgen][error]") {
    const char* node = R"NODE(
1 7 0 0
1  0 0 0 0 0 0 0
)NODE";
    const char* ele = R"ELE(
0 4 0
)ELE";
    const auto r = parse_tetgen_strings(node, ele);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("dim must be 2 or 3") != std::string::npos);
}

TEST_CASE("TetGen reader emits a single anonymous patch when faces lack markers",
          "[io][tetgen]") {
    const char* face = R"FACE(
2 0
1   1 2 4
2   2 3 4
)FACE";
    const auto r = parse_tetgen_strings(kNode2tet, kEle2tet, face);
    INFO(r.error);
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];
    REQUIRE(z.boundaries.size() == 1);
    CHECK(z.boundaries[0].name == "all_faces");
    const auto si = z.boundaries[0].faceElementIndices[0];
    CHECK(z.sections[si].element_count() == 2);
}
