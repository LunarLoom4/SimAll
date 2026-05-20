// =============================================================================
// SimAll Beta - Meshing Unit Tests
// File   : tests/unit/meshing/test_blockmesh_dict.cpp
// Phase  : 23 Pass 14
//
// Exercises the OpenFOAM-style blockMeshDict parser + MultiblockHex builder:
//   - tokeniser handles // and /* */ comments
//   - convertToMeters scales vertices
//   - blocks gain divisions / optional name / simpleGrading
//   - boundary patches assigned 1-based face zones via 4-vertex set match
//   - assemble() yields a valid Mesh with the right cell count
// =============================================================================
#include "meshing/BlockMeshDict.hpp"
#include "meshing/MeshStorage.hpp"
#include "meshing/MultiblockHex.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace simall;
using Catch::Approx;

namespace
{

constexpr const char* kUnitCubeDict = R"(
/*--------------------------------*- C++ -*----------------------------------*\
| Test fixture - unit cube split into one hex block                          |
\*---------------------------------------------------------------------------*/

convertToMeters 1.0;

vertices
(
    (0 0 0)    // 0
    (1 0 0)    // 1
    (1 1 0)    // 2
    (0 1 0)    // 3
    (0 0 1)    // 4
    (1 0 1)    // 5
    (1 1 1)    // 6
    (0 1 1)    // 7
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (2 2 2) simpleGrading (1 1 1)
);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }
    walls
    {
        type wall;
        faces
        (
            (0 1 5 4)
            (3 7 6 2)
            (0 3 2 1)
            (4 5 6 7)
        );
    }
);
)";

} // namespace

// =============================================================================
TEST_CASE("blockMeshDict parser handles vertices, blocks and patches", "[meshing][blockmesh]")
{
    const auto d = meshing::parse_block_mesh_dict(kUnitCubeDict);
    REQUIRE(d.convertToMeters == Approx(1.0));
    REQUIRE(d.vertices.size() == 8);
    REQUIRE(d.vertices[6].x == Approx(1.0));
    REQUIRE(d.vertices[6].z == Approx(1.0));
    REQUIRE(d.blocks.size() == 1);
    REQUIRE(d.blocks[0].divisions[0] == 2);
    REQUIRE(d.blocks[0].divisions[1] == 2);
    REQUIRE(d.blocks[0].divisions[2] == 2);
    REQUIRE(d.blocks[0].grading[0] == Approx(1.0));
    REQUIRE(d.patches.size() == 3);
    REQUIRE(d.patches[0].name == "inlet");
    REQUIRE(d.patches[0].type == "patch");
    REQUIRE(d.patches[0].faces.size() == 1);
    REQUIRE(d.patches[2].name == "walls");
    REQUIRE(d.patches[2].type == "wall");
    REQUIRE(d.patches[2].faces.size() == 4);
}

// =============================================================================
TEST_CASE("blockMeshDict convertToMeters scales vertices", "[meshing][blockmesh]")
{
    const std::string dict = R"(
        convertToMeters 0.001;
        vertices ( (0 0 0) (1000 0 0) (1000 1000 0) (0 1000 0)
                   (0 0 1000) (1000 0 1000) (1000 1000 1000) (0 1000 1000) );
        blocks ( hex (0 1 2 3 4 5 6 7) (1 1 1) simpleGrading (1 1 1) );
    )";
    const auto d = meshing::parse_block_mesh_dict(dict);
    meshing::MultiblockHex mb;
    mb.initialize(meshing::MultiblockProps{});
    meshing::build_multiblock(d, mb);
    meshing::Mesh m;
    mb.assemble(m);
    m.compute_geometry();
    REQUIRE(m.cells().size() == 1);
    // After scaling, the cube is 1 m on a side => volume 1.0.
    REQUIRE(m.cells().volume[0] == Approx(1.0).margin(1e-9));
}

// =============================================================================
TEST_CASE("blockMeshDict builds a 2x2x2 cell mesh and tags patches", "[meshing][blockmesh]")
{
    const auto d = meshing::parse_block_mesh_dict(kUnitCubeDict);
    meshing::MultiblockHex mb;
    mb.initialize(meshing::MultiblockProps{});
    const auto stats = meshing::build_multiblock(d, mb);
    REQUIRE(stats.blocks == 1);
    REQUIRE(stats.patches == 3);
    REQUIRE(stats.patchFaces == 6);
    REQUIRE(stats.orphanFaces == 0);

    meshing::Mesh m;
    const std::size_t nCells = mb.assemble(m);
    m.compute_geometry();
    REQUIRE(nCells == 8); // 2*2*2
    REQUIRE(m.cells().size() == 8);
    for (auto v : m.cells().volume)
        REQUIRE(v > 0.0);
}

// =============================================================================
TEST_CASE("blockMeshDict parser tolerates // and /* */ comments", "[meshing][blockmesh]")
{
    const std::string dict = R"(
        // header
        convertToMeters 1.0;  // unit
        /* block of vertices */
        vertices (
            (0 0 0) (1 0 0) (1 1 0) (0 1 0)
            (0 0 1) (1 0 1) (1 1 1) (0 1 1)
        );
        blocks (
            hex (0 1 2 3 4 5 6 7) (1 1 1) /* uniform */ simpleGrading (1 1 1)
        );
    )";
    const auto d = meshing::parse_block_mesh_dict(dict);
    REQUIRE(d.vertices.size() == 8);
    REQUIRE(d.blocks.size() == 1);
    REQUIRE(d.blocks[0].divisions[0] == 1);
}

// =============================================================================
TEST_CASE("blockMeshDict reports orphan patch faces that match no block face",
          "[meshing][blockmesh]")
{
    const std::string dict = R"(
        vertices (
            (0 0 0) (1 0 0) (1 1 0) (0 1 0)
            (0 0 1) (1 0 1) (1 1 1) (0 1 1)
        );
        blocks ( hex (0 1 2 3 4 5 6 7) (1 1 1) simpleGrading (1 1 1) );
        boundary (
            ghost { type patch; faces ( (0 1 2 7) ); }
        );
    )";
    const auto d = meshing::parse_block_mesh_dict(dict);
    meshing::MultiblockHex mb;
    mb.initialize(meshing::MultiblockProps{});
    const auto stats = meshing::build_multiblock(d, mb);
    REQUIRE(stats.patchFaces == 0);
    REQUIRE(stats.orphanFaces == 1);
}

// =============================================================================
TEST_CASE("blockMeshDict parser raises on unknown top-level key", "[meshing][blockmesh]")
{
    const std::string dict = "fooBar 1;";
    REQUIRE_THROWS_AS(meshing::parse_block_mesh_dict(dict), std::runtime_error);
}

// =============================================================================
// Pass 15 - simpleGrading is now applied (not just parsed).
// =============================================================================
TEST_CASE("blockMeshDict simpleGrading produces geometric cell sizes in x",
          "[meshing][blockmesh][grading]")
{
    // 4 cells along x with G = 4 (last/first cell ratio).  Expected per-cell
    // ratio r = 4^(1/3) ~= 1.5874.  First-cell length = (r-1)/(r^4 - 1).
    const std::string dict = R"(
        vertices (
            (0 0 0) (1 0 0) (1 1 0) (0 1 0)
            (0 0 1) (1 0 1) (1 1 1) (0 1 1)
        );
        blocks ( hex (0 1 2 3 4 5 6 7) (4 1 1) simpleGrading (4 1 1) );
    )";
    const auto d = meshing::parse_block_mesh_dict(dict);
    REQUIRE(d.blocks[0].grading[0] == Approx(4.0));

    meshing::MultiblockHex mb;
    mb.initialize(meshing::MultiblockProps{});
    meshing::build_multiblock(d, mb);
    meshing::Mesh m;
    const auto nCells = mb.assemble(m);
    m.compute_geometry();
    REQUIRE(nCells == 4);

    // Locate the 4 cells by x-centroid and confirm geometric progression.
    std::vector<double> cx(m.cells().centroidX.begin(), m.cells().centroidX.end());
    std::sort(cx.begin(), cx.end());
    REQUIRE(cx.size() == 4);

    // Cell widths are differences of consecutive cell boundaries.  With u[i]
    // = (r^i-1)/(r^4-1) and r = 4^(1/3), widths w_i = u[i]-u[i-1] obey
    //   w_i = r^(i-1) * w_1
    // hence w_4 / w_1 == G == 4.
    const double r = std::pow(4.0, 1.0 / 3.0);
    const double denom = std::pow(r, 4.0) - 1.0;
    const double w1 = (r - 1.0) / denom;
    const double w4 = (std::pow(r, 4.0) - std::pow(r, 3.0)) / denom;
    REQUIRE(w4 / w1 == Approx(4.0).margin(1e-9));

    // The width of cell-i = 2 * (cx[i] - left_face_position[i]); easier check:
    // sum of widths = 1.0 (unit cube).
    const double sumVol = std::accumulate(m.cells().volume.begin(), m.cells().volume.end(), 0.0);
    REQUIRE(sumVol == Approx(1.0).margin(1e-9));

    // First-cell volume / last-cell volume = 1/4 (only x varies).
    std::vector<double> vols(m.cells().volume.begin(), m.cells().volume.end());
    std::sort(vols.begin(), vols.end());
    REQUIRE(vols.back() / vols.front() == Approx(4.0).margin(1e-9));
}

// =============================================================================
TEST_CASE("blockMeshDict simpleGrading (1 1 1) matches uniform spacing",
          "[meshing][blockmesh][grading]")
{
    const std::string dict = R"(
        vertices (
            (0 0 0) (1 0 0) (1 1 0) (0 1 0)
            (0 0 1) (1 0 1) (1 1 1) (0 1 1)
        );
        blocks ( hex (0 1 2 3 4 5 6 7) (3 1 1) simpleGrading (1 1 1) );
    )";
    meshing::MultiblockHex mb;
    mb.initialize(meshing::MultiblockProps{});
    meshing::build_multiblock(meshing::parse_block_mesh_dict(dict), mb);
    meshing::Mesh m;
    mb.assemble(m);
    m.compute_geometry();
    REQUIRE(m.cells().size() == 3);
    for (auto v : m.cells().volume)
        REQUIRE(v == Approx(1.0 / 3.0).margin(1e-12));
}

// =============================================================================
// Pass 16 - writer round-trip.
// =============================================================================
TEST_CASE("blockMeshDict writer round-trips parse -> write -> parse",
          "[meshing][blockmesh][writer]")
{
    const auto d1 = meshing::parse_block_mesh_dict(kUnitCubeDict);
    const std::string text = meshing::write_block_mesh_dict(d1);
    const auto d2 = meshing::parse_block_mesh_dict(text);

    REQUIRE(d2.convertToMeters == Approx(d1.convertToMeters));
    REQUIRE(d2.vertices.size() == d1.vertices.size());
    for (std::size_t i = 0; i < d1.vertices.size(); ++i) {
        REQUIRE(d2.vertices[i].x == Approx(d1.vertices[i].x));
        REQUIRE(d2.vertices[i].y == Approx(d1.vertices[i].y));
        REQUIRE(d2.vertices[i].z == Approx(d1.vertices[i].z));
    }
    REQUIRE(d2.blocks.size() == d1.blocks.size());
    for (std::size_t i = 0; i < d1.blocks.size(); ++i) {
        REQUIRE(d2.blocks[i].vertices == d1.blocks[i].vertices);
        REQUIRE(d2.blocks[i].divisions == d1.blocks[i].divisions);
        for (int k = 0; k < 3; ++k)
            REQUIRE(d2.blocks[i].grading[k] == Approx(d1.blocks[i].grading[k]));
    }
    REQUIRE(d2.patches.size() == d1.patches.size());
    for (std::size_t i = 0; i < d1.patches.size(); ++i) {
        REQUIRE(d2.patches[i].name == d1.patches[i].name);
        REQUIRE(d2.patches[i].type == d1.patches[i].type);
        REQUIRE(d2.patches[i].faces.size() == d1.patches[i].faces.size());
        for (std::size_t j = 0; j < d1.patches[i].faces.size(); ++j)
            REQUIRE(d2.patches[i].faces[j] == d1.patches[i].faces[j]);
    }
}

// =============================================================================
TEST_CASE("blockMeshDict writer preserves non-trivial grading and convertToMeters",
          "[meshing][blockmesh][writer]")
{
    meshing::BlockMeshDict d;
    d.convertToMeters = 0.001;
    d.vertices = {{0, 0, 0},
                  {1000, 0, 0},
                  {1000, 1000, 0},
                  {0, 1000, 0},
                  {0, 0, 1000},
                  {1000, 0, 1000},
                  {1000, 1000, 1000},
                  {0, 1000, 1000}};
    meshing::BlockMeshBlock b;
    b.vertices = {0, 1, 2, 3, 4, 5, 6, 7};
    b.divisions = {5, 3, 2};
    b.grading = {2.5, 0.5, 1.0};
    d.blocks.push_back(b);

    const auto d2 = meshing::parse_block_mesh_dict(meshing::write_block_mesh_dict(d));
    REQUIRE(d2.convertToMeters == Approx(0.001));
    REQUIRE(d2.blocks.size() == 1);
    REQUIRE(d2.blocks[0].divisions == b.divisions);
    REQUIRE(d2.blocks[0].grading[0] == Approx(2.5));
    REQUIRE(d2.blocks[0].grading[1] == Approx(0.5));
    REQUIRE(d2.blocks[0].grading[2] == Approx(1.0));
}

// =============================================================================
TEST_CASE("blockMeshDict writer round-trip yields identical assembled mesh",
          "[meshing][blockmesh][writer]")
{
    const auto d1 = meshing::parse_block_mesh_dict(kUnitCubeDict);
    const auto d2 = meshing::parse_block_mesh_dict(meshing::write_block_mesh_dict(d1));

    auto build = [](const meshing::BlockMeshDict& d) {
        meshing::MultiblockHex mb;
        mb.initialize(meshing::MultiblockProps{});
        meshing::build_multiblock(d, mb);
        meshing::Mesh m;
        mb.assemble(m);
        m.compute_geometry();
        return m;
    };
    const auto m1 = build(d1);
    const auto m2 = build(d2);
    REQUIRE(m1.cells().size() == m2.cells().size());
    REQUIRE(m1.faces().size() == m2.faces().size());
    REQUIRE(m1.nodes().size() == m2.nodes().size());

    double v1 = 0, v2 = 0;
    for (auto v : m1.cells().volume)
        v1 += v;
    for (auto v : m2.cells().volume)
        v2 += v;
    REQUIRE(v1 == Approx(v2).margin(1e-12));
}
