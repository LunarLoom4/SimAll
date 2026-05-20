// =============================================================================
// SimAll Beta - Meshing Unit Tests
// File   : tests/unit/meshing/test_mesh_ops.cpp
// Phase  : 23 Pass 11
//
// Coverage for the MeshOps utility module:
//   - transform_points / translate / rotate / mirror
//   - merge_meshes (raw concat)
//   - stitch_meshes (merge + node weld + duplicate-face dedupe)
//   - renumber_cells_cuthill_mckee (RCM)
//   - check_mesh (consolidated topology + quality report)
// =============================================================================
#include "meshing/CartesianMesher.hpp"
#include "meshing/MeshOps.hpp"
#include "meshing/MeshStorage.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace simall;
using Catch::Approx;

namespace
{

meshing::Mesh make_brick(util::Vec3d origin, util::Vec3d extent, int nx = 1, int ny = 1, int nz = 1)
{
    meshing::CartesianGridSpec s;
    s.origin = origin;
    s.extent = extent;
    s.Nx = nx;
    s.Ny = ny;
    s.Nz = nz;
    meshing::Mesh m;
    meshing::CartesianMesher{s}.generate(m);
    return m;
}

std::size_t count_boundary_faces(const meshing::Mesh& m)
{
    std::size_t k = 0;
    for (auto v : m.faces().neighbor)
        if (v == meshing::kBoundaryCell)
            ++k;
    return k;
}

} // namespace

// =============================================================================
TEST_CASE("translate shifts every node and preserves volume", "[meshing][ops][xform]")
{
    auto m = make_brick({0, 0, 0}, {1, 1, 1});
    const double v0 = std::accumulate(m.cells().volume.begin(), m.cells().volume.end(), 0.0);

    meshing::ops::translate(m, 10.0, -3.0, 2.5);

    REQUIRE(m.nodes().x.front() == Approx(10.0));
    REQUIRE(m.nodes().y.front() == Approx(-3.0));
    REQUIRE(m.nodes().z.front() == Approx(2.5));

    const double v1 = std::accumulate(m.cells().volume.begin(), m.cells().volume.end(), 0.0);
    REQUIRE(v1 == Approx(v0));
}

// =============================================================================
TEST_CASE("rotate 90deg about z maps (1,0,0) -> (0,1,0)", "[meshing][ops][xform]")
{
    auto m = make_brick({0, 0, 0}, {1, 1, 1});
    meshing::ops::rotate(m, {0, 0, 1}, M_PI / 2.0, {0, 0, 0});

    // The node that started at (1, 0, 0) should now be at (0, 1, 0).
    bool found = false;
    for (std::size_t i = 0; i < m.nodes().size(); ++i) {
        if (std::abs(m.nodes().x[i] - 0.0) < 1e-9 && std::abs(m.nodes().y[i] - 1.0) < 1e-9
            && std::abs(m.nodes().z[i] - 0.0) < 1e-9) {
            found = true;
            break;
        }
    }
    REQUIRE(found);
    // All volumes remain positive.
    for (auto v : m.cells().volume)
        REQUIRE(v > 0.0);
}

// =============================================================================
TEST_CASE("mirror across x=0 negates x and keeps volumes positive", "[meshing][ops][xform]")
{
    auto m = make_brick({0, 0, 0}, {1, 1, 1});
    meshing::ops::mirror(m, {1, 0, 0}, {0, 0, 0});

    for (std::size_t i = 0; i < m.nodes().size(); ++i) {
        REQUIRE(m.nodes().x[i] <= 1.0e-12); // mirrored into x <= 0
    }
    // Because reflection flips orientation, transform_points reverses the
    // face vertex order; volumes must remain positive.
    for (auto v : m.cells().volume)
        REQUIRE(v > 0.0);
}

// =============================================================================
TEST_CASE("merge_meshes concatenates without deduplication", "[meshing][ops][merge]")
{
    auto a = make_brick({0, 0, 0}, {1, 1, 1});
    auto b = make_brick({10, 0, 0}, {1, 1, 1}); // far away, no shared nodes

    const auto nA_nodes = a.nodes().size();
    const auto nB_nodes = b.nodes().size();
    const auto nA_faces = a.faces().size();
    const auto nB_faces = b.faces().size();
    const auto nA_cells = a.cells().size();
    const auto nB_cells = b.cells().size();

    meshing::Mesh out;
    meshing::ops::merge_meshes(a, b, out);

    REQUIRE(out.nodes().size() == nA_nodes + nB_nodes);
    REQUIRE(out.faces().size() == nA_faces + nB_faces);
    REQUIRE(out.cells().size() == nA_cells + nB_cells);

    // All faces remain boundary - meshes don't touch.
    REQUIRE(count_boundary_faces(out) == nA_faces + nB_faces);

    // All volumes remain positive.
    for (auto v : out.cells().volume)
        REQUIRE(v > 0.0);
}

// =============================================================================
TEST_CASE("stitch_meshes welds coincident nodes and dedupes shared faces", "[meshing][ops][stitch]")
{
    // Two unit cubes sharing the x=1 plane (4 coincident nodes, 1 shared face).
    auto a = make_brick({0, 0, 0}, {1, 1, 1});
    auto b = make_brick({1, 0, 0}, {1, 1, 1});

    const auto bndBefore = count_boundary_faces(a) + count_boundary_faces(b);
    const auto cellsBefore = a.cells().size() + b.cells().size();

    meshing::Mesh out;
    auto s = meshing::ops::stitch_meshes(a, b, out, {1e-9});

    REQUIRE(s.weldedNodes == 4);       // four shared nodes on the x=1 plane
    REQUIRE(s.deduplicatedFaces == 1); // the x=1 face from each cube becomes one internal face

    REQUIRE(out.cells().size() == cellsBefore); // cells preserved
    REQUIRE(out.nodes().size() == 16 - 4);      // 12 unique nodes
    // Two cubes contribute 12 boundary faces total; one pair becomes internal:
    //    surviving face count = bndBefore - 1
    REQUIRE(out.faces().size() == bndBefore - 1);
    REQUIRE(count_boundary_faces(out) == bndBefore - 2);

    // Every cell still has 6 faces (1 of which is now shared).
    for (std::size_t c = 0; c + 1 < out.cells().faceOffsets.size(); ++c) {
        const auto nFaces = out.cells().faceOffsets[c + 1] - out.cells().faceOffsets[c];
        REQUIRE(nFaces == 6);
    }
    // Volumes remain positive.
    for (auto v : out.cells().volume)
        REQUIRE(v > 0.0);
}

// =============================================================================
TEST_CASE("renumber_cells_cuthill_mckee never increases bandwidth", "[meshing][ops][rcm]")
{
    // A 4x1x1 line of cells - the natural numbering is already optimal,
    // but RCM must not make it worse.
    auto m = make_brick({0, 0, 0}, {4, 1, 1}, 4, 1, 1);

    auto stats = meshing::ops::renumber_cells_cuthill_mckee(m);
    REQUIRE(stats.cellsPermuted == m.cells().size());
    REQUIRE(stats.bandwidthAfter <= stats.bandwidthBefore);

    // Face owner/neighbor are still valid indices.
    for (std::size_t i = 0; i < m.faces().size(); ++i) {
        REQUIRE(m.faces().owner[i] < m.cells().size());
        if (m.faces().neighbor[i] != meshing::kBoundaryCell) {
            REQUIRE(m.faces().neighbor[i] < m.cells().size());
        }
    }
}

// =============================================================================
TEST_CASE("renumber on a 3x3x3 cube cuts bandwidth or matches natural order", "[meshing][ops][rcm]")
{
    auto m = make_brick({0, 0, 0}, {3, 3, 3}, 3, 3, 3);

    auto stats = meshing::ops::renumber_cells_cuthill_mckee(m);
    REQUIRE(stats.cellsPermuted == 27);
    // Reverse Cuthill-McKee on a 3x3x3 grid yields the same diagonal-front
    // ordering as the natural row-major numbering; bandwidth should not
    // increase relative to the natural numbering.
    REQUIRE(stats.bandwidthAfter <= stats.bandwidthBefore);
}

// =============================================================================
TEST_CASE("check_mesh reports OK for a healthy brick", "[meshing][ops][check]")
{
    auto m = make_brick({0, 0, 0}, {2, 2, 2}, 2, 2, 2);
    auto r = meshing::ops::check_mesh(m);
    REQUIRE(r.ok);
    REQUIRE(r.errors.empty());
    REQUIRE(r.nCells == 8);
    REQUIRE(r.nNodes > 0);
    REQUIRE(r.nFaces > 0);
    REQUIRE(r.nBoundaryFaces + r.nInteriorFaces == r.nFaces);
    REQUIRE(r.minCellVolume > 0.0);
    REQUIRE(r.nNegativeVolume == 0);
    REQUIRE(r.nZeroVolumeCells == 0);
    REQUIRE(r.nDegenerateFaces == 0);
}

// =============================================================================
TEST_CASE("check_mesh flags zero-volume cells", "[meshing][ops][check]")
{
    auto m = make_brick({0, 0, 0}, {1, 1, 1});
    // Manually corrupt the volume of cell 0 to simulate a degenerate cell.
    m.cells().volume[0] = 0.0;

    auto r = meshing::ops::check_mesh(m);
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.nZeroVolumeCells == 1);
    const auto& fmt = r.format();
    REQUIRE(fmt.find("FAIL") != std::string::npos);
    REQUIRE(fmt.find("zero-vol cells   : 1") != std::string::npos);
}

// =============================================================================
TEST_CASE("check_mesh flags negative-volume cells", "[meshing][ops][check]")
{
    auto m = make_brick({0, 0, 0}, {1, 1, 1});
    m.cells().volume[0] = -1.0;

    auto r = meshing::ops::check_mesh(m);
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.nNegativeVolume == 1);
}

// =============================================================================
// Pass 11b coverage
// =============================================================================
TEST_CASE("split_mesh extracts a contiguous half of a brick", "[meshing][ops][split]")
{
    auto src = make_brick({0, 0, 0}, {2, 1, 1}, 2, 1, 1); // 2 cells
    std::vector<std::uint8_t> keep = {1, 0};              // keep only cell 0

    meshing::Mesh out;
    auto s = meshing::ops::split_mesh(src, keep, out);

    REQUIRE(s.selectedCells == 1);
    REQUIRE(s.outputCells == 1);
    REQUIRE(out.cells().size() == 1);
    // The interface face that used to be internal is now a boundary face.
    REQUIRE(out.faces().size() == 6);
    REQUIRE(count_boundary_faces(out) == 6);
    REQUIRE(out.cells().volume[0] > 0.0);
}

// =============================================================================
TEST_CASE("split_mesh with all-zero mask emits an empty mesh", "[meshing][ops][split]")
{
    auto src = make_brick({0, 0, 0}, {1, 1, 1}, 2, 2, 2);
    std::vector<std::uint8_t> keep(src.cells().size(), 0);

    meshing::Mesh out;
    auto s = meshing::ops::split_mesh(src, keep, out);

    REQUIRE(s.selectedCells == 0);
    REQUIRE(out.cells().size() == 0);
    REQUIRE(out.nodes().size() == 0);
    REQUIRE(out.faces().size() == 0);
}

// =============================================================================
TEST_CASE("split_mesh with all-one mask round-trips cell counts", "[meshing][ops][split]")
{
    auto src = make_brick({0, 0, 0}, {2, 2, 2}, 2, 2, 2);
    std::vector<std::uint8_t> keep(src.cells().size(), 1);

    meshing::Mesh out;
    auto s = meshing::ops::split_mesh(src, keep, out);

    REQUIRE(s.selectedCells == src.cells().size());
    REQUIRE(out.cells().size() == src.cells().size());
    REQUIRE(out.nodes().size() == src.nodes().size());
    REQUIRE(out.faces().size() == src.faces().size());
    for (auto v : out.cells().volume)
        REQUIRE(v > 0.0);
}

// =============================================================================
TEST_CASE("refine_hex turns a single cube into 8 sub-cubes", "[meshing][ops][refine]")
{
    auto m = make_brick({0, 0, 0}, {1, 1, 1}, 1, 1, 1);
    const double volBefore = std::accumulate(m.cells().volume.begin(), m.cells().volume.end(), 0.0);

    auto s = meshing::ops::refine_hex(m);
    REQUIRE(s.cellsBefore == 1);
    REQUIRE(s.cellsRefined == 1);
    REQUIRE(s.rejectedCells == 0);
    REQUIRE(s.cellsAfter == 8);
    REQUIRE(m.cells().size() == 8);

    // Total volume preserved by refinement.
    const double volAfter = std::accumulate(m.cells().volume.begin(), m.cells().volume.end(), 0.0);
    REQUIRE(volAfter == Approx(volBefore));
    // Every sub-cell is 1/8 of the parent.
    for (auto v : m.cells().volume)
        REQUIRE(v == Approx(0.125));
}

// =============================================================================
TEST_CASE("refine_hex on a 2x2x2 brick produces 64 sub-cubes with shared nodes",
          "[meshing][ops][refine]")
{
    auto m = make_brick({0, 0, 0}, {1, 1, 1}, 2, 2, 2);
    const std::size_t nodesBefore = m.nodes().size();
    const double volBefore = std::accumulate(m.cells().volume.begin(), m.cells().volume.end(), 0.0);

    auto s = meshing::ops::refine_hex(m);
    REQUIRE(s.cellsBefore == 8);
    REQUIRE(s.cellsRefined == 8);
    REQUIRE(s.rejectedCells == 0);
    REQUIRE(s.cellsAfter == 64);
    REQUIRE(m.cells().size() == 64);

    // Refined brick is a 4x4x4 grid -> 5x5x5 = 125 unique nodes.
    // Edge / face / cell midpoints are deduplicated across cell boundaries
    // by the refine routine, so we must end up with exactly 125 nodes.
    REQUIRE(m.nodes().size() == 125);
    REQUIRE(m.nodes().size()
            <= nodesBefore + 8 * 19); // upper bound (12 edges + 6 faces + 1 centroid)

    const double volAfter = std::accumulate(m.cells().volume.begin(), m.cells().volume.end(), 0.0);
    REQUIRE(volAfter == Approx(volBefore));
    for (auto v : m.cells().volume)
        REQUIRE(v > 0.0);
}
