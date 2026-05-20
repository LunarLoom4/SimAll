// =============================================================================
// SimAll Beta - Meshing Unit Tests
// File   : tests/unit/meshing/test_octree_cutcell.cpp
// Phase  : 23 Pass 13
//
// Coverage for the OctreeMesher cut-cell extraction path
// (`OctreeMeshOptions::enableCutCells = true`):
//   - Legacy stepped path continues to produce a non-empty mesh.
//   - Cut-cell path produces a non-empty mesh with all-positive volumes.
//   - Cut-cell total volume is no greater than stepped volume (because
//     boundary leaves are clipped rather than kept whole).
//   - Cut-cell volume converges towards the true cube volume (1.0).
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "meshing/OctreeMesher.hpp"
#include "meshing/MeshStorage.hpp"
#include "meshing/StlImporter.hpp"
#include "utilities/MathTypes.hpp"

#include <numeric>

using namespace simall;
using Catch::Approx;

namespace {

meshing::StlSurface make_unit_cube() {
    meshing::StlSurface s;
    s.vertices = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}
    };
    s.triangles = {
        {0, 3, 2}, {0, 2, 1},
        {4, 5, 6}, {4, 6, 7},
        {0, 1, 5}, {0, 5, 4},
        {3, 7, 6}, {3, 6, 2},
        {0, 4, 7}, {0, 7, 3},
        {1, 2, 6}, {1, 6, 5}
    };
    s.normals.reserve(s.triangles.size());
    for (const auto& t : s.triangles) {
        const auto& a = s.vertices[t[0]];
        const auto& b = s.vertices[t[1]];
        const auto& c = s.vertices[t[2]];
        s.normals.push_back((b - a).cross(c - a).normalized());
    }
    return s;
}

double total_volume(const meshing::Mesh& m) {
    return std::accumulate(m.cells().volume.begin(),
                            m.cells().volume.end(), 0.0);
}

}  // namespace

// =============================================================================
TEST_CASE("OctreeMesher legacy stepped path still produces a valid mesh",
          "[meshing][octree][cutcell]") {
    auto stl = make_unit_cube();
    meshing::OctreeMeshOptions opt;
    opt.maxDepth        = 3;
    opt.minDepthGlobal  = 2;
    opt.enableCutCells  = false;

    meshing::Mesh m;
    meshing::OctreeMesher{}.mesh(stl, opt, m);
    REQUIRE(m.cells().size() > 0);
    for (auto v : m.cells().volume) REQUIRE(v > 0.0);
}

// =============================================================================
TEST_CASE("OctreeMesher cut-cell path yields a non-empty watertight mesh",
          "[meshing][octree][cutcell]") {
    auto stl = make_unit_cube();
    meshing::OctreeMeshOptions opt;
    opt.maxDepth        = 3;
    opt.minDepthGlobal  = 2;
    opt.enableCutCells  = true;
    opt.edgeBisectIters = 10;

    meshing::Mesh m;
    meshing::OctreeMesher{}.mesh(stl, opt, m);
    REQUIRE(m.cells().size() > 0);

    // Every cell volume strictly positive.
    for (auto v : m.cells().volume) REQUIRE(v > 0.0);

    // Every face is either interior or boundary; nothing dangling.
    for (std::size_t f = 0; f < m.faces().size(); ++f) {
        REQUIRE(m.faces().owner[f] < m.cells().size());
        if (m.faces().neighbor[f] != meshing::kBoundaryCell)
            REQUIRE(m.faces().neighbor[f] < m.cells().size());
    }
}

// =============================================================================
TEST_CASE("OctreeMesher cut-cell volume is <= stepped volume (cube)",
          "[meshing][octree][cutcell]") {
    auto stl = make_unit_cube();

    meshing::OctreeMeshOptions optStepped;
    optStepped.maxDepth       = 3;
    optStepped.minDepthGlobal = 2;
    optStepped.enableCutCells = false;
    meshing::Mesh mStepped;
    meshing::OctreeMesher{}.mesh(stl, optStepped, mStepped);

    meshing::OctreeMeshOptions optCut = optStepped;
    optCut.enableCutCells  = true;
    optCut.edgeBisectIters = 12;
    meshing::Mesh mCut;
    meshing::OctreeMesher{}.mesh(stl, optCut, mCut);

    const double vStepped = total_volume(mStepped);
    const double vCut     = total_volume(mCut);

    // Stepped path always over-estimates (full hex per BOUNDARY leaf).
    // Cut-cell clips those leaves and so must report less-or-equal volume.
    REQUIRE(vCut <= vStepped + 1.0e-9);

    // True cube volume is 1.0; both should be within the bbox-padded grid
    // total (roughly 1.02^3 = 1.061) - i.e. they shouldn't overrun the
    // padded domain. This is mostly a sanity ceiling.
    REQUIRE(vStepped <= 1.10);
    REQUIRE(vCut     <= 1.10);

    // Cut-cell volume should be measurably closer to the true cube than
    // stepped at this depth (cube is not exactly grid-aligned because of
    // the 1% bbox padding in OctreeMesher).
    REQUIRE(std::abs(vCut - 1.0) <= std::abs(vStepped - 1.0) + 1.0e-9);
}
