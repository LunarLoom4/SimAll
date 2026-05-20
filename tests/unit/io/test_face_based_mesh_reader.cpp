// =============================================================================
// SimAll Beta — IO Unit Tests
// File   : tests/unit/io/test_face_based_mesh_reader.cpp
// Phase  : 23 Pass 4 (renamed in Pass 5 follow-up)
// =============================================================================
#include "io/FaceBasedMeshReader.hpp"
#include "io/MeshFormats.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>

using simall::io::classify_face_based_cell;
using simall::io::ElementType;
using simall::io::parse_face_based_mesh_strings;

namespace
{

// -- single-tetrahedron mesh -------------------------------------------------
// 4 points, 4 triangular boundary faces, 1 tet cell, 3 patches.

constexpr const char* kTetPoints = R"FBM(
PointHeader
{
    note "optional dictionary header is skipped";
}
4
(
(0 0 0)
(1 0 0)
(0 1 0)
(0 0 1)
)
)FBM";

constexpr const char* kTetFaces = R"FBM(
FaceHeader { }
4
(
3(0 2 1)
3(0 1 3)
3(0 3 2)
3(1 2 3)
)
)FBM";

constexpr const char* kTetOwner = R"FBM(
OwnerHeader { }
4
(
0
0
0
0
)
)FBM";

constexpr const char* kTetNeighbour = R"FBM(
NeighbourHeader { }
0
(
)
)FBM";

constexpr const char* kTetBoundary = R"FBM(
BoundaryHeader { }
3
(
    inlet
    {
        type            patch;
        nFaces          1;
        startFace       0;
    }
    outlet
    {
        type            patch;
        nFaces          1;
        startFace       1;
    }
    wall
    {
        type            wall;
        nFaces          2;
        startFace       2;
    }
)
)FBM";

// -- single-hexahedron mesh (no headers at all) ------------------------------
constexpr const char* kHexPoints = R"FBM(
8
(
(0 0 0)
(1 0 0)
(1 1 0)
(0 1 0)
(0 0 1)
(1 0 1)
(1 1 1)
(0 1 1)
)
)FBM";

constexpr const char* kHexFaces = R"FBM(
6
(
4(0 3 2 1)
4(4 5 6 7)
4(0 1 5 4)
4(1 2 6 5)
4(2 3 7 6)
4(3 0 4 7)
)
)FBM";

constexpr const char* kHexOwner = R"FBM(
6
(
0
0
0
0
0
0
)
)FBM";

constexpr const char* kHexNeighbour = R"FBM(
0
(
)
)FBM";

constexpr const char* kHexBoundary = R"FBM(
1
(
    walls
    {
        type            wall;
        nFaces          6;
        startFace       0;
    }
)
)FBM";

// -- two-hex mesh (cells stacked along z, one internal face) -----------------
constexpr const char* kTwoHexPoints = R"FBM(
12
(
(0 0 0)
(1 0 0)
(1 1 0)
(0 1 0)
(0 0 1)
(1 0 1)
(1 1 1)
(0 1 1)
(0 0 2)
(1 0 2)
(1 1 2)
(0 1 2)
)
)FBM";

constexpr const char* kTwoHexFaces = R"FBM(
11
(
4(4 5 6 7)
4(0 3 2 1)
4(0 1 5 4)
4(1 2 6 5)
4(2 3 7 6)
4(3 0 4 7)
4(8 9 10 11)
4(4 5 9 8)
4(5 6 10 9)
4(6 7 11 10)
4(7 4 8 11)
)
)FBM";

constexpr const char* kTwoHexOwner = R"FBM(
11
(
0
0
0
0
0
0
1
1
1
1
1
)
)FBM";

constexpr const char* kTwoHexNeighbour = R"FBM(
1
(
1
)
)FBM";

constexpr const char* kTwoHexBoundary = R"FBM(
1
(
    walls
    {
        type            wall;
        nFaces          10;
        startFace       1;
    }
)
)FBM";

} // namespace

TEST_CASE("classify_face_based_cell maps face signatures to element types", "[io][facebased]")
{
    CHECK(classify_face_based_cell(4, 0) == ElementType::Tetra4);
    CHECK(classify_face_based_cell(0, 6) == ElementType::Hexa8);
    CHECK(classify_face_based_cell(4, 1) == ElementType::Pyra5);
    CHECK(classify_face_based_cell(2, 3) == ElementType::Penta6);
    CHECK(classify_face_based_cell(3, 0) == ElementType::Unknown);
    CHECK(classify_face_based_cell(0, 5) == ElementType::Unknown);
}

TEST_CASE("Face-based reader parses a single-tet mesh", "[io][facebased]")
{
    const auto r = parse_face_based_mesh_strings(
        kTetPoints, kTetFaces, kTetOwner, kTetNeighbour, kTetBoundary, "tet");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);
    const auto& z = r.mesh.zones[0];
    REQUIRE(z.x.size() == 4);
    CHECK(z.x[1] == 1.0);

    REQUIRE(z.sections.size() == 4); // 1 vol + 3 patch
    std::size_t tet = 0, tri = 0;
    for (const auto& s : z.sections) {
        if (s.type == ElementType::Tetra4)
            tet += s.element_count();
        if (s.type == ElementType::Tri3)
            tri += s.element_count();
    }
    CHECK(tet == 1);
    CHECK(tri == 4);

    CHECK(r.mesh.total_nodes() == 4);
    CHECK(r.mesh.total_elements() == 5);
    CHECK(r.mesh.sourceFormat == "face_based_polymesh_dir");
}

TEST_CASE("Face-based reader preserves patch metadata", "[io][facebased]")
{
    const auto r = parse_face_based_mesh_strings(
        kTetPoints, kTetFaces, kTetOwner, kTetNeighbour, kTetBoundary, "tet");
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];
    REQUIRE(z.boundaries.size() == 3);
    CHECK(z.boundaries[0].name == "inlet");
    CHECK(z.boundaries[0].bcType == "patch");
    CHECK(z.boundaries[1].name == "outlet");
    CHECK(z.boundaries[1].bcType == "patch");
    CHECK(z.boundaries[2].name == "wall");
    CHECK(z.boundaries[2].bcType == "wall");
    for (const auto& bp : z.boundaries) {
        REQUIRE(bp.faceElementIndices.size() == 1);
        const auto si = bp.faceElementIndices[0];
        REQUIRE(si < z.sections.size());
        CHECK(z.sections[si].type == ElementType::Tri3);
        if (bp.name == "wall")
            CHECK(z.sections[si].element_count() == 2);
        else
            CHECK(z.sections[si].element_count() == 1);
    }
}

TEST_CASE("Face-based reader parses a single-hex mesh", "[io][facebased]")
{
    const auto r = parse_face_based_mesh_strings(
        kHexPoints, kHexFaces, kHexOwner, kHexNeighbour, kHexBoundary, "hex");
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.dimension == 3);
    const auto& z = r.mesh.zones[0];
    REQUIRE(z.sections.size() == 2);
    std::size_t hex = 0, quad = 0;
    for (const auto& s : z.sections) {
        if (s.type == ElementType::Hexa8)
            hex = s.element_count();
        if (s.type == ElementType::Quad4)
            quad = s.element_count();
    }
    CHECK(hex == 1);
    CHECK(quad == 6);
    REQUIRE(z.boundaries.size() == 1);
    CHECK(z.boundaries[0].name == "walls");
    CHECK(z.boundaries[0].bcType == "wall");
}

TEST_CASE("Face-based reader handles internal faces (two stacked hex cells)", "[io][facebased]")
{
    const auto r = parse_face_based_mesh_strings(
        kTwoHexPoints, kTwoHexFaces, kTwoHexOwner, kTwoHexNeighbour, kTwoHexBoundary, "twohex");
    INFO(r.error);
    REQUIRE(r.ok);
    const auto& z = r.mesh.zones[0];
    REQUIRE(z.sections.size() == 2);
    std::size_t hex = 0, quad = 0;
    for (const auto& s : z.sections) {
        if (s.type == ElementType::Hexa8)
            hex = s.element_count();
        if (s.type == ElementType::Quad4)
            quad = s.element_count();
    }
    CHECK(hex == 2);
    CHECK(quad == 10);
    REQUIRE(z.boundaries.size() == 1);
    CHECK(z.boundaries[0].name == "walls");
}

TEST_CASE("Face-based reader accepts headerless files", "[io][facebased]")
{
    const auto r = parse_face_based_mesh_strings(
        kHexPoints, kHexFaces, kHexOwner, kHexNeighbour, kHexBoundary, "hex");
    CHECK(r.ok);
}

TEST_CASE("Face-based reader rejects empty points file", "[io][facebased][error]")
{
    const auto r = parse_face_based_mesh_strings(
        "0\n(\n)\n", kTetFaces, kTetOwner, kTetNeighbour, kTetBoundary, "x");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("points file is empty") != std::string::npos);
}

TEST_CASE("Face-based reader rejects owner/faces size mismatch", "[io][facebased][error]")
{
    const char* badOwner = "3\n(\n0\n0\n0\n)\n";
    const auto r = parse_face_based_mesh_strings(
        kTetPoints, kTetFaces, badOwner, kTetNeighbour, kTetBoundary, "x");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("!= faces.size()") != std::string::npos);
}

TEST_CASE("Face-based reader rejects unsupported cell face signature", "[io][facebased][error]")
{
    const char* pts = "3\n(\n(0 0 0)\n(1 0 0)\n(0 1 0)\n)\n";
    const char* fcs = "3\n(\n3(0 1 2)\n3(0 1 2)\n3(0 1 2)\n)\n";
    const char* own = "3\n(\n0\n0\n0\n)\n";
    const auto r = parse_face_based_mesh_strings(pts, fcs, own, "", "", "bad");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("unsupported face signature") != std::string::npos);
}

TEST_CASE("Face-based reader rejects out-of-range patch face slice", "[io][facebased][error]")
{
    const char* badBoundary = "1\n(\n"
                              "    bad { type patch; nFaces 99; startFace 0; }\n"
                              ")\n";
    const auto r = parse_face_based_mesh_strings(
        kTetPoints, kTetFaces, kTetOwner, kTetNeighbour, badBoundary, "x");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("is out of bounds") != std::string::npos);
}

TEST_CASE("Face-based reader rejects non-homogeneous patch", "[io][facebased][error]")
{
    const char* pts = kHexPoints;
    const char* fcs = "7\n(\n"
                      "4(0 3 2 1)\n"
                      "4(4 5 6 7)\n"
                      "4(0 1 5 4)\n"
                      "4(1 2 6 5)\n"
                      "4(2 3 7 6)\n"
                      "4(3 0 4 7)\n"
                      "3(0 1 2)\n"
                      ")\n";
    const char* own = "7\n(\n0\n0\n0\n0\n0\n0\n0\n)\n";
    const char* bnd = "1\n(\n"
                      "    mixed { type patch; nFaces 7; startFace 0; }\n"
                      ")\n";
    const auto r = parse_face_based_mesh_strings(pts, fcs, own, "", bnd, "x");
    CHECK_FALSE(r.ok);
    CHECK((r.error.find("unsupported face signature") != std::string::npos
           || r.error.find("non-homogeneous") != std::string::npos));
}
