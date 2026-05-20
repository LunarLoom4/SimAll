// =============================================================================
// SimAll Beta — IO Unit Tests
// File   : tests/unit/io/test_starcd_reader.cpp
// Phase  : 23 Pass 3
// =============================================================================
#include "io/MeshFormats.hpp"
#include "io/StarCdReader.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

using simall::io::ElementType;
using simall::io::parse_starcd_strings;
using simall::io::starcd_classify_cell;
using simall::io::starcd_classify_face;

namespace
{

// Minimal tet case: 4 vertices, 1 collapsed-hex tet, 4 boundary tris on
// three patches (inlet / outlet / wall, wall has 2 faces).
constexpr const char* kVrtTet = "  1  0.0  0.0  0.0\n"
                                "  2  1.0  0.0  0.0\n"
                                "  3  0.0  1.0  0.0\n"
                                "  4  0.0  0.0  1.0\n";

constexpr const char* kCelTet = "  1   1 2 3 4 4 4 4 4   3 100\n";

constexpr const char* kBndTet = "  1   1 2 3 3   1 wall inlet\n"
                                "  2   1 2 4 4   2 wall outlet\n"
                                "  3   1 3 4 4   3 wall wall\n"
                                "  4   2 3 4 4   3 wall wall\n";

} // namespace

TEST_CASE("starcd_classify_cell maps distinct-node count to element type", "[io][starcd]")
{
    std::vector<long long> uniq;
    CHECK(starcd_classify_cell({1, 2, 3, 4, 4, 4, 4, 4}, uniq) == ElementType::Tetra4);
    CHECK(uniq.size() == 4);

    CHECK(starcd_classify_cell({1, 2, 3, 4, 5, 5, 5, 5}, uniq) == ElementType::Pyra5);
    CHECK(uniq.size() == 5);

    CHECK(starcd_classify_cell({1, 2, 3, 3, 4, 5, 6, 6}, uniq) == ElementType::Penta6);
    CHECK(uniq.size() == 6);

    CHECK(starcd_classify_cell({1, 2, 3, 4, 5, 6, 7, 8}, uniq) == ElementType::Hexa8);
    CHECK(uniq.size() == 8);

    // 7 distinct nodes are not a known collapse pattern.
    CHECK(starcd_classify_cell({1, 2, 3, 4, 5, 6, 7, 7}, uniq) == ElementType::Unknown);
}

TEST_CASE("starcd_classify_face maps 3/4 unique to Tri3/Quad4", "[io][starcd]")
{
    std::vector<long long> uniq;
    CHECK(starcd_classify_face({1, 2, 3, 3}, uniq) == ElementType::Tri3);
    CHECK(uniq.size() == 3);
    CHECK(starcd_classify_face({1, 2, 3, 4}, uniq) == ElementType::Quad4);
    CHECK(uniq.size() == 4);
    // All-equal collapse is not a valid face.
    CHECK(starcd_classify_face({1, 1, 1, 1}, uniq) == ElementType::Unknown);
}

TEST_CASE("Star-CD reader parses tet triple", "[io][starcd]")
{
    const auto r = parse_starcd_strings(kVrtTet, kCelTet, kBndTet, "tet");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);
    const auto& z = r.mesh.zones[0];

    REQUIRE(z.x.size() == 4);
    CHECK(z.x[0] == 0.0);
    CHECK(z.x[1] == 1.0);

    // 1 volume section (Tetra4) + 3 patch sections (one per patch name).
    REQUIRE(z.sections.size() == 4);
    std::size_t tetCount = 0, triCount = 0;
    for (const auto& s : z.sections) {
        if (s.type == ElementType::Tetra4) {
            tetCount += s.element_count();
            CHECK(s.nodes.size() == s.element_count() * 4);
        } else if (s.type == ElementType::Tri3) {
            triCount += s.element_count();
            CHECK(s.nodes.size() == s.element_count() * 3);
        }
    }
    CHECK(tetCount == 1);
    CHECK(triCount == 4);

    CHECK(r.mesh.total_nodes() == 4);
    CHECK(r.mesh.total_elements() == 5);
    CHECK(r.mesh.sourceFormat == "starcd_prostar");
}

TEST_CASE("Star-CD reader derives patches keyed by patch name", "[io][starcd]")
{
    const auto r = parse_starcd_strings(kVrtTet, kCelTet, kBndTet, "tet");
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
        if (bp.name == "wall") {
            // wall has 2 faces.
            CHECK(z.sections[si].element_count() == 2);
        } else {
            CHECK(z.sections[si].element_count() == 1);
        }
    }
}

TEST_CASE("Star-CD reader works without a .bnd file", "[io][starcd]")
{
    const auto r = parse_starcd_strings(kVrtTet, kCelTet, "", "tet");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.mesh.zones.size() == 1);
    CHECK(r.mesh.zones[0].boundaries.empty());
    // Only the volume section.
    CHECK(r.mesh.zones[0].sections.size() == 1);
    CHECK(r.mesh.zones[0].sections[0].type == ElementType::Tetra4);
}

TEST_CASE("Star-CD reader falls back to region_<id> when patch name absent", "[io][starcd]")
{
    // Lines with only 6 tokens (id, 4 nodes, region_id) -> no patch name.
    const char* bnd = "  1   1 2 3 3   7\n"
                      "  2   1 2 4 4   7\n"
                      "  3   1 3 4 4   9\n";
    const auto r = parse_starcd_strings(kVrtTet, kCelTet, bnd, "tet");
    INFO(r.error);
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];
    REQUIRE(z.boundaries.size() == 2); // regions 7 and 9
    CHECK(z.boundaries[0].name == "region_7");
    CHECK(z.boundaries[1].name == "region_9");
}

TEST_CASE("Star-CD reader skips blank lines and comments", "[io][starcd]")
{
    const std::string vrt = "# comment line\n"
                            "\n"
                            "  1  0  0  0\n"
                            "! pro-am style\n"
                            "  2  1  0  0\n"
                            "  3  0  1  0\n"
                            "  4  0  0  1\n";
    const auto r = parse_starcd_strings(vrt, kCelTet, "", "tet");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.mesh.zones[0].x.size() == 4);
}

TEST_CASE("Star-CD reader parses a hex cell and quad faces", "[io][starcd]")
{
    // 8 vertices of a unit cube.
    const char* vrt = "  1  0 0 0\n"
                      "  2  1 0 0\n"
                      "  3  1 1 0\n"
                      "  4  0 1 0\n"
                      "  5  0 0 1\n"
                      "  6  1 0 1\n"
                      "  7  1 1 1\n"
                      "  8  0 1 1\n";
    const char* cel = "  1   1 2 3 4 5 6 7 8   1 1\n";
    const char* bnd = "  1   1 2 3 4   1 wall bottom\n";
    const auto r = parse_starcd_strings(vrt, cel, bnd, "cube");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    const auto& z = r.mesh.zones[0];
    // 1 hex volume + 1 quad patch section.
    REQUIRE(z.sections.size() == 2);
    std::size_t hex = 0, quad = 0;
    for (const auto& s : z.sections) {
        if (s.type == ElementType::Hexa8)
            hex = s.element_count();
        if (s.type == ElementType::Quad4)
            quad = s.element_count();
    }
    CHECK(hex == 1);
    CHECK(quad == 1);
    REQUIRE(z.boundaries.size() == 1);
    CHECK(z.boundaries[0].name == "bottom");
}

TEST_CASE("Star-CD reader rejects empty .vrt", "[io][starcd][error]")
{
    const auto r = parse_starcd_strings("", kCelTet, "", "x");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("no vertices") != std::string::npos);
}

TEST_CASE("Star-CD reader rejects unknown vertex id in .cel", "[io][starcd][error]")
{
    const char* cel = "  1   1 2 3 99 99 99 99 99   3 1\n";
    const auto r = parse_starcd_strings(kVrtTet, cel, "", "bad");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unknown vertex id") != std::string::npos);
}

TEST_CASE("Star-CD reader rejects mixed-type patch", "[io][starcd][error]")
{
    // 8 vertices of a cube + 1 hex + .bnd with both a tri and a quad
    // assigned to the same patch name "mixed".
    const char* vrt = "  1  0 0 0\n"
                      "  2  1 0 0\n"
                      "  3  1 1 0\n"
                      "  4  0 1 0\n"
                      "  5  0 0 1\n"
                      "  6  1 0 1\n"
                      "  7  1 1 1\n"
                      "  8  0 1 1\n";
    const char* cel = "  1   1 2 3 4 5 6 7 8   1 1\n";
    const char* bnd = "  1   1 2 3 4   1 wall mixed\n"  // quad
                      "  2   1 2 3 3   1 wall mixed\n"; // tri
    const auto r = parse_starcd_strings(vrt, cel, bnd, "cube");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("mixes") != std::string::npos);
}

TEST_CASE("Star-CD reader rejects degenerate cell with 7 distinct nodes", "[io][starcd][error]")
{
    // 7 vertices, cell uses 1..7 with 7 distinct (last duplicated).
    const char* vrt = "  1  0 0 0\n"
                      "  2  1 0 0\n"
                      "  3  1 1 0\n"
                      "  4  0 1 0\n"
                      "  5  0 0 1\n"
                      "  6  1 0 1\n"
                      "  7  1 1 1\n";
    const char* cel = "  1   1 2 3 4 5 6 7 7   1 1\n";
    const auto r = parse_starcd_strings(vrt, cel, "", "bad");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("cannot infer element type") != std::string::npos);
}
