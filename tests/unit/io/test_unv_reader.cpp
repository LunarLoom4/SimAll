// =============================================================================
// SimAll Beta — IO Unit Tests
// File   : tests/unit/io/test_unv_reader.cpp
// Phase  : 23 Pass 2
//
// Coverage for the SDRC/I-DEAS Universal ".unv" ASCII reader.
// =============================================================================
#include "io/MeshFormats.hpp"
#include "io/UnvReader.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>

using simall::io::ElementType;
using simall::io::parse_unv_string;
using simall::io::unv_element_type;
using simall::io::unv_normalize_float_token;
using simall::io::UnvReadResult;

namespace
{

// Minimal valid UNV: 4 nodes, 1 tet, 4 surface tris, with three permanent
// groups ("inlet", "outlet", "wall"). The "wall" group holds 2 element
// references in a single 8-int row, exercising the packed-pair branch of
// the group reader.
constexpr const char* kSimpleTetUnv = R"UNV(    -1
   164
1Newtons       (N)              1Meters       (m)              1Kelvin       (K)
0.0000000000000000D+00  0.0000000000000000D+00  0.0000000000000000D+00
    -1
    -1
  2411
         1         1         1         11
   0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00
         2         1         1         11
   1.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00
         3         1         1         11
   0.00000000000000D+00  1.00000000000000D+00  0.00000000000000D+00
         4         1         1         11
   0.00000000000000D+00  0.00000000000000D+00  1.00000000000000D+00
    -1
    -1
  2412
         1       111         1         1         7         4
         1         2         3         4
         2        91         1         1         7         3
         1         2         3
         3        91         1         1         7         3
         1         2         4
         4        91         1         1         7         3
         1         3         4
         5        91         1         1         7         3
         2         3         4
    -1
    -1
  2467
         1         0         0         0         0         0         0         1
inlet
         8         2         0         0
         2         0         0         0         0         0         0         1
outlet
         8         3         0         0
         3         0         0         0         0         0         0         2
wall
         8         4         0         0         8         5         0         0
    -1
)UNV";

} // namespace

TEST_CASE("unv_element_type translates supported descriptors", "[io][unv]")
{
    CHECK(unv_element_type(11) == ElementType::Bar2);
    CHECK(unv_element_type(21) == ElementType::Bar2);
    CHECK(unv_element_type(22) == ElementType::Bar2);
    CHECK(unv_element_type(41) == ElementType::Tri3);
    CHECK(unv_element_type(91) == ElementType::Tri3);
    CHECK(unv_element_type(44) == ElementType::Quad4);
    CHECK(unv_element_type(94) == ElementType::Quad4);
    CHECK(unv_element_type(111) == ElementType::Tetra4);
    CHECK(unv_element_type(112) == ElementType::Penta6);
    CHECK(unv_element_type(115) == ElementType::Hexa8);
    CHECK(unv_element_type(42) == ElementType::Unknown);  // parabolic
    CHECK(unv_element_type(118) == ElementType::Unknown); // parabolic tet
}

TEST_CASE("unv_normalize_float_token rewrites Fortran D-exponents", "[io][unv]")
{
    CHECK(unv_normalize_float_token("1.5D+02") == "1.5e+02");
    CHECK(unv_normalize_float_token("-3.0d-01") == "-3.0e-01");
    CHECK(unv_normalize_float_token("1.0e+00") == "1.0e+00"); // pass-through
    CHECK(unv_normalize_float_token("0") == "0");
}

TEST_CASE("UNV reader parses a minimal tetrahedron mesh", "[io][unv]")
{
    const auto r = parse_unv_string(kSimpleTetUnv, "tet.unv");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);
    const auto& z = r.mesh.zones[0];

    REQUIRE(z.x.size() == 4);
    CHECK(z.x[0] == 0.0);
    CHECK(z.x[1] == 1.0);
    CHECK(z.y[2] == 1.0);
    CHECK(z.z[3] == 1.0);

    // One section per FE descriptor: descriptor 111 (tet) and 91 (shell tri).
    REQUIRE(z.sections.size() == 2);
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
    CHECK(r.mesh.sourceFormat == "ideas_unv");
}

TEST_CASE("UNV reader derives boundary patches from element groups", "[io][unv]")
{
    const auto r = parse_unv_string(kSimpleTetUnv, "tet.unv");
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

TEST_CASE("UNV reader skips unknown datasets (164 units)", "[io][unv]")
{
    // The fixture above already includes a 164 units block; the reader
    // must process the file end-to-end without diagnostics. This test
    // additionally exercises a fictional dataset number 9999.
    std::string with_extra = std::string(kSimpleTetUnv);
    // Inject an unknown dataset before the trailing newline.
    with_extra += "    -1\n  9999\nignored payload line\n42 17\n    -1\n";
    const auto r = parse_unv_string(with_extra, "extra.unv");
    INFO(r.error);
    CHECK(r.ok);
    CHECK(r.mesh.zones.size() == 1);
}

TEST_CASE("UNV reader rejects parabolic FE descriptors", "[io][unv][error]")
{
    const std::string body = "    -1\n  2411\n"
                             "         1         1         1         11\n"
                             "   0.0   0.0   0.0\n"
                             "         2         1         1         11\n"
                             "   1.0   0.0   0.0\n"
                             "         3         1         1         11\n"
                             "   0.0   1.0   0.0\n"
                             "         4         1         1         11\n"
                             "   0.5   0.0   0.0\n"
                             "         5         1         1         11\n"
                             "   0.5   0.5   0.0\n"
                             "         6         1         1         11\n"
                             "   0.0   0.5   0.0\n"
                             "    -1\n"
                             "    -1\n  2412\n"
                             "         1        42         1         1         7         6\n"
                             "         1         2         3         4         5         6\n"
                             "    -1\n";
    const auto r = parse_unv_string(body, "parabolic.unv");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unsupported FE descriptor") != std::string::npos);
}

TEST_CASE("UNV reader rejects unknown node labels in connectivity", "[io][unv][error]")
{
    // 4 nodes defined but element references node 99.
    const std::string body = "    -1\n  2411\n"
                             "         1         1         1         11\n"
                             "   0.0   0.0   0.0\n"
                             "         2         1         1         11\n"
                             "   1.0   0.0   0.0\n"
                             "         3         1         1         11\n"
                             "   0.0   1.0   0.0\n"
                             "         4         1         1         11\n"
                             "   0.0   0.0   1.0\n"
                             "    -1\n"
                             "    -1\n  2412\n"
                             "         1       111         1         1         7         4\n"
                             "         1         2         3        99\n"
                             "    -1\n";
    const auto r = parse_unv_string(body, "badnode.unv");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unknown node label") != std::string::npos);
}

TEST_CASE("UNV reader rejects unterminated dataset", "[io][unv][error]")
{
    const std::string body = "    -1\n  2411\n"
                             "         1         1         1         11\n"
                             "   0.0   0.0   0.0\n"; // missing trailing -1
    const auto r = parse_unv_string(body, "trunc.unv");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unterminated dataset") != std::string::npos);
}

TEST_CASE("UNV reader normalises Fortran D-format float coords", "[io][unv]")
{
    // Verify D-format coords are accepted by the inline parser (the fixture
    // already uses them, so this is an explicit positive assertion).
    const auto r = parse_unv_string(kSimpleTetUnv, "tet.unv");
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];
    CHECK(z.x[1] == 1.0);
    CHECK(z.y[2] == 1.0);
    CHECK(z.z[3] == 1.0);
}

TEST_CASE("UNV reader handles a 2D quad-only mesh, reports dim==2", "[io][unv]")
{
    const std::string body = "    -1\n  2411\n"
                             "         1         1         1         11\n"
                             "   0.0   0.0   0.0\n"
                             "         2         1         1         11\n"
                             "   1.0   0.0   0.0\n"
                             "         3         1         1         11\n"
                             "   1.0   1.0   0.0\n"
                             "         4         1         1         11\n"
                             "   0.0   1.0   0.0\n"
                             "    -1\n"
                             "    -1\n  2412\n"
                             "         1        44         1         1         7         4\n"
                             "         1         2         3         4\n"
                             "    -1\n";
    const auto r = parse_unv_string(body, "quad.unv");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 2);
    REQUIRE(r.mesh.zones[0].sections.size() == 1);
    CHECK(r.mesh.zones[0].sections[0].type == ElementType::Quad4);
}

TEST_CASE("UNV reader handles a beam descriptor with orientation triple", "[io][unv]")
{
    // Descriptor 21 (linear beam) inserts a 3-int beam-orientation record
    // between the header and the connectivity.
    const std::string body = "    -1\n  2411\n"
                             "         1         1         1         11\n"
                             "   0.0   0.0   0.0\n"
                             "         2         1         1         11\n"
                             "   1.0   0.0   0.0\n"
                             "    -1\n"
                             "    -1\n  2412\n"
                             "         1        21         1         1         7         2\n"
                             "         0         0         0\n"
                             "         1         2\n"
                             "    -1\n";
    const auto r = parse_unv_string(body, "beam.unv");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 1);
    REQUIRE(r.mesh.zones[0].sections.size() == 1);
    const auto& s = r.mesh.zones[0].sections[0];
    CHECK(s.type == ElementType::Bar2);
    CHECK(s.element_count() == 1);
    CHECK(s.nodes.size() == 2);
}
