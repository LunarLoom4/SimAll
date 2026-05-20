// =============================================================================
// SimAll Beta - Meshing Unit Tests
// File   : tests/unit/meshing/test_snappy_hex_mesher.cpp
// Phase  : 23 Pass 10
//
// Drives the SnappyHexMesher orchestrator against a procedurally-built
// unit-cube STL surface and checks each pipeline stage in isolation:
//   - castellation produces a non-empty hex background mesh
//   - snapping moves boundary nodes onto the input surface (within tol)
//   - layer addition (when requested) emits a separate prism mesh
//   - mesh_background_only() honours the "no-snap, no-layer" contract
// =============================================================================
#include "meshing/MeshStorage.hpp"
#include "meshing/SnappyHexMesher.hpp"
#include "meshing/StlImporter.hpp"
#include "utilities/MathTypes.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>

using namespace simall;
using Catch::Approx;

namespace
{

// Build an in-memory unit-cube StlSurface (8 verts, 12 outward-CCW tris).
meshing::StlSurface make_unit_cube()
{
    meshing::StlSurface s;
    s.vertices = {
        {0.0, 0.0, 0.0}, // v0
        {1.0, 0.0, 0.0}, // v1
        {1.0, 1.0, 0.0}, // v2
        {0.0, 1.0, 0.0}, // v3
        {0.0, 0.0, 1.0}, // v4
        {1.0, 0.0, 1.0}, // v5
        {1.0, 1.0, 1.0}, // v6
        {0.0, 1.0, 1.0}, // v7
    };
    s.triangles = {
        // z- (outward = -z)
        {0, 3, 2},
        {0, 2, 1},
        // z+ (outward = +z)
        {4, 5, 6},
        {4, 6, 7},
        // y- (outward = -y)
        {0, 1, 5},
        {0, 5, 4},
        // y+ (outward = +y)
        {3, 7, 6},
        {3, 6, 2},
        // x- (outward = -x)
        {0, 4, 7},
        {0, 7, 3},
        // x+ (outward = +x)
        {1, 2, 6},
        {1, 6, 5},
    };
    // Per-triangle normals (not consumed by OctreeMesher or SnappyHexMesher,
    // but kept consistent with the StlSurface contract).
    s.normals.reserve(s.triangles.size());
    for (const auto& t : s.triangles) {
        const auto& a = s.vertices[t[0]];
        const auto& b = s.vertices[t[1]];
        const auto& c = s.vertices[t[2]];
        s.normals.push_back((b - a).cross(c - a).normalized());
    }
    return s;
}

// Distance from point to the surface of the unit cube [0,1]^3 (signed-distance
// magnitude collapsed to absolute residual: 0 exactly on the surface).
double cube_surface_residual(double x, double y, double z) noexcept
{
    double rx = std::min(std::abs(x), std::abs(x - 1.0));
    double ry = std::min(std::abs(y), std::abs(y - 1.0));
    double rz = std::min(std::abs(z), std::abs(z - 1.0));
    return std::min({rx, ry, rz});
}

} // namespace

// =============================================================================
TEST_CASE("SnappyHexMesher castellates a unit cube into a non-empty mesh", "[meshing][snappy]")
{
    auto cube = make_unit_cube();

    meshing::SnappyHexOptions opt{};
    opt.maxDepth = 3;
    opt.minDepth = 2;
    opt.enableSnapping = false;
    opt.nLayers = 0;

    meshing::Mesh bg, prism;
    auto stats = meshing::SnappyHexMesher{}.mesh(cube, opt, bg, prism);

    REQUIRE(stats.backgroundCells > 0);
    REQUIRE(stats.backgroundNodes > 0);
    REQUIRE(stats.snappedNodes == 0);    // snapping disabled
    REQUIRE(stats.prismLayerCells == 0); // no layers requested
    REQUIRE(prism.cells().size() == 0);

    // Every node lies in the padded bounding box [-0.05, 1.05]^3.
    const auto& nx = bg.nodes().x;
    const auto& ny = bg.nodes().y;
    const auto& nz = bg.nodes().z;
    for (std::size_t i = 0; i < nx.size(); ++i) {
        REQUIRE(nx[i] >= -0.05);
        REQUIRE(nx[i] <= 1.05);
        REQUIRE(ny[i] >= -0.05);
        REQUIRE(ny[i] <= 1.05);
        REQUIRE(nz[i] >= -0.05);
        REQUIRE(nz[i] <= 1.05);
    }
}

// =============================================================================
TEST_CASE("SnappyHexMesher snapping moves boundary nodes onto the surface", "[meshing][snappy]")
{
    auto cube = make_unit_cube();

    meshing::SnappyHexOptions noSnap{};
    noSnap.maxDepth = 3;
    noSnap.minDepth = 2;
    noSnap.enableSnapping = false;
    noSnap.nLayers = 0;

    meshing::SnappyHexOptions snap = noSnap;
    snap.enableSnapping = true;
    snap.nSnapIters = 4;
    snap.snapMaxDistFrac = 4.0; // generous so all boundary nodes participate

    // Measure boundary-residual reduction by running the same castellation
    // twice - once without snapping, once with - and comparing the average
    // distance of boundary nodes to the cube surface.
    meshing::Mesh bg0, prism0;
    auto s0 = meshing::SnappyHexMesher{}.mesh(cube, noSnap, bg0, prism0);

    meshing::Mesh bg1, prism1;
    auto s1 = meshing::SnappyHexMesher{}.mesh(cube, snap, bg1, prism1);

    REQUIRE(s1.boundaryNodes > 0);
    REQUIRE(s1.snappedNodes > 0);
    REQUIRE(s1.snapItersRun == 4);

    // Pick a reference node count from the no-snap run (same castellation
    // parameters yield the same mesh topology, modulo identical node order).
    REQUIRE(bg0.nodes().size() == bg1.nodes().size());

    double sumBefore = 0.0;
    double sumAfter = 0.0;
    std::size_t boundaryCount = 0;
    const auto& f0 = bg0.faces();
    const auto& f1 = bg1.faces();
    std::vector<bool> isBoundaryNode(bg0.nodes().size(), false);
    for (std::size_t fi = 0; fi < f0.size(); ++fi) {
        if (f0.neighbor[fi] != meshing::kBoundaryCell)
            continue;
        const auto b = f0.nodeOffsets[fi];
        const auto e = f0.nodeOffsets[fi + 1];
        for (auto k = b; k < e; ++k)
            isBoundaryNode[f0.nodeIndices[k]] = true;
    }
    for (std::size_t i = 0; i < bg0.nodes().size(); ++i) {
        if (!isBoundaryNode[i])
            continue;
        ++boundaryCount;
        sumBefore += cube_surface_residual(bg0.nodes().x[i], bg0.nodes().y[i], bg0.nodes().z[i]);
        sumAfter += cube_surface_residual(bg1.nodes().x[i], bg1.nodes().y[i], bg1.nodes().z[i]);
    }
    REQUIRE(boundaryCount > 0);
    INFO("boundaryCount=" << boundaryCount << " sumBefore=" << sumBefore
                          << " sumAfter=" << sumAfter);
    // Snapping must strictly reduce the aggregate boundary residual.
    REQUIRE(sumAfter < sumBefore);
    // And the average residual after snapping should be near machine zero
    // relative to the cube extent.
    REQUIRE((sumAfter / static_cast<double>(boundaryCount)) < 1.0e-9);
    // Faces array remains consistent after geometry refresh.
    REQUIRE(f1.size() == f0.size());
}

// =============================================================================
TEST_CASE("SnappyHexMesher layer addition produces a separate prism mesh",
          "[meshing][snappy][layers]")
{
    auto cube = make_unit_cube();

    meshing::SnappyHexOptions opt{};
    opt.maxDepth = 3;
    opt.minDepth = 2;
    opt.enableSnapping = true;
    opt.nSnapIters = 3;
    opt.snapMaxDistFrac = 4.0;
    opt.nLayers = 2;
    opt.firstLayerHeight = 0.005;
    opt.layerGrowthRatio = 1.2;

    meshing::Mesh bg, prism;
    auto stats = meshing::SnappyHexMesher{}.mesh(cube, opt, bg, prism);

    REQUIRE(stats.backgroundCells > 0);
    REQUIRE(stats.snappedNodes > 0);
    REQUIRE(stats.prismLayerCells > 0);
    REQUIRE(stats.prismLayerNodes > 0);
    // 2 layers * tris-per-wall-node.  Wall surface for a unit cube castellated
    // at depth 3 yields >= 6 quads worth of triangles; each becomes one prism
    // per layer, so prismCells must be a positive multiple of nLayers.
    REQUIRE(prism.cells().size() % static_cast<std::size_t>(opt.nLayers) == 0);
    // Background mesh is intentionally NOT modified by the layer phase.
    REQUIRE(bg.cells().size() == stats.backgroundCells);
}

// =============================================================================
TEST_CASE("SnappyHexMesher::mesh_background_only is snap-/layer-free", "[meshing][snappy]")
{
    auto cube = make_unit_cube();

    meshing::SnappyHexOptions opt{};
    opt.maxDepth = 3;
    opt.minDepth = 2;
    // These flags should be ignored by the background-only entry point.
    opt.enableSnapping = true;
    opt.nSnapIters = 5;
    opt.nLayers = 4;

    meshing::Mesh bg;
    auto stats = meshing::SnappyHexMesher{}.mesh_background_only(cube, opt, bg);

    REQUIRE(stats.backgroundCells > 0);
    REQUIRE(stats.snappedNodes == 0);
    REQUIRE(stats.snapItersRun == 0);
    REQUIRE(stats.prismLayerCells == 0);
}
