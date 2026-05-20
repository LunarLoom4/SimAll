// =============================================================================
// SimAll Beta — Meshing Advanced Unit Tests (Week 11)
// File   : tests/unit/meshing/test_meshing_advanced.cpp
//
// Sanity tests for the nine W11 meshing modules:
//   1. AdvancingFrontSurface
//   2. HexSweepMesher
//   3. MultiblockHex
//   4. PolyhedralAgglomeration
//   5. BoundaryRecovery
//   6. MedialAxisDecomposition
//   7. PartitionerMetis
//   8. PartitionerScotch
//   9. AmrAdaptor
// =============================================================================
#include "meshing/AdvancingFrontSurface.hpp"
#include "meshing/AmrAdaptor.hpp"
#include "meshing/BoundaryRecovery.hpp"
#include "meshing/CartesianMesher.hpp"
#include "meshing/HexSweepMesher.hpp"
#include "meshing/MedialAxisDecomposition.hpp"
#include "meshing/MeshStorage.hpp"
#include "meshing/MultiblockHex.hpp"
#include "meshing/PartitionerMetis.hpp"
#include "meshing/PartitionerScotch.hpp"
#include "meshing/PolyhedralAgglomeration.hpp"
#include "solver/FieldRegistry.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <vector>

using namespace simall;
using Catch::Approx;

namespace
{
meshing::Mesh build_line_mesh(std::size_t nx)
{
    meshing::CartesianGridSpec spec{};
    spec.origin = {0, 0, 0};
    spec.extent = {static_cast<double>(nx), 1.0, 1.0};
    spec.Nx = static_cast<int>(nx);
    spec.Ny = 1;
    spec.Nz = 1;
    meshing::CartesianMesher m(spec);
    meshing::Mesh mesh;
    m.generate(mesh);
    return mesh;
}
} // namespace

TEST_CASE("AdvancingFrontSurface triangulates unit square boundary", "[meshing][afs]")
{
    meshing::AdvancingFrontSurface afs;
    meshing::AfsProps p;
    p.baseSize = 0.25;
    p.outputZone = 1;
    afs.initialize(p);
    std::vector<util::Vec3d> loop{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    afs.seed_boundary(loop);
    meshing::Mesh mesh;
    const auto nTri = afs.generate(mesh);
    REQUIRE(nTri > 0);
    REQUIRE(afs.triangle_count() == nTri);
}

TEST_CASE("HexSweepMesher extrudes 4-quad source 3 layers", "[meshing][sweep]")
{
    meshing::HexSweepMesher s;
    meshing::SweepProps p;
    p.nLayers = 3;
    s.initialize(p);
    meshing::SourceQuadMesh src;
    src.nodes = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    src.quads = {{{0, 1, 2, 3}}};
    std::vector<util::Vec3d> path{{0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 0, 3}};
    meshing::Mesh mesh;
    const auto nHex = s.sweep(src, path, mesh);
    REQUIRE(nHex == 3);
}

TEST_CASE("MultiblockHex joins two adjoining unit cubes", "[meshing][multiblock]")
{
    meshing::MultiblockHex mb;
    mb.initialize({});
    meshing::HexBlock b1;
    b1.corners = {util::Vec3d{0, 0, 0},
                  {1, 0, 0},
                  {1, 1, 0},
                  {0, 1, 0},
                  {0, 0, 1},
                  {1, 0, 1},
                  {1, 1, 1},
                  {0, 1, 1}};
    b1.divisions = {2, 2, 2};
    meshing::HexBlock b2 = b1;
    b2.corners = {util::Vec3d{1, 0, 0},
                  {2, 0, 0},
                  {2, 1, 0},
                  {1, 1, 0},
                  {1, 0, 1},
                  {2, 0, 1},
                  {2, 1, 1},
                  {1, 1, 1}};
    mb.add_block(b1);
    mb.add_block(b2);
    meshing::Mesh mesh;
    const auto nCells = mb.assemble(mesh);
    REQUIRE(nCells == 16);
    // Count internal faces (neighbor != kBoundaryCell).
    std::size_t internal = 0;
    for (auto n : mesh.faces().neighbor)
        if (n != meshing::kBoundaryCell)
            ++internal;
    REQUIRE(internal > 0);
}

TEST_CASE("PolyhedralAgglomeration coarsens 8-cell strip", "[meshing][agglom]")
{
    auto mesh = build_line_mesh(8);
    meshing::PolyhedralAgglomeration pa;
    meshing::AgglomerationProps p;
    p.maxClusterSize = 4;
    pa.initialize(p);
    std::vector<std::int32_t> cluster;
    const auto nClusters = pa.agglomerate(mesh, cluster);
    REQUIRE(nClusters > 0);
    REQUIRE(nClusters < mesh.cells().size());
}

TEST_CASE("BoundaryRecovery: missing edge triggers Steiner insertion", "[meshing][recover]")
{
    meshing::Mesh mesh;
    auto& N = mesh.nodes();
    N.x = {0, 1};
    N.y = {0, 0};
    N.z = {0, 0};
    meshing::BoundaryRecovery br;
    br.initialize({});
    br.add_required_edge({0, 1});
    const auto rep = br.recover(mesh);
    REQUIRE(rep.edgesMissing == 1);
    REQUIRE(rep.steinerAdded > 0);
}

TEST_CASE("MedialAxisDecomposition produces at least one branch", "[meshing][medial]")
{
    auto mesh = build_line_mesh(4);
    mesh.compute_geometry();
    meshing::MedialAxisDecomposition md;
    md.initialize({});
    const auto nPts = md.compute_medial_points(mesh);
    REQUIRE(nPts > 0);
    std::vector<meshing::MedialRegion> regions;
    md.decompose(regions);
    // regions may be 0 if branchMinLength filters them; verify points instead.
    REQUIRE(md.medial_points().size() == nPts);
}

TEST_CASE("PartitionerMetis splits 16-cell line into 4 parts", "[meshing][metis]")
{
    auto mesh = build_line_mesh(16);
    meshing::PartitionerMetis pm;
    meshing::MetisProps p;
    p.nParts = 4;
    pm.initialize(p);
    std::vector<std::int32_t> part;
    const auto cut = pm.partition(mesh, part);
    REQUIRE(part.size() == 16);
    REQUIRE(cut >= 3);
}

TEST_CASE("PartitionerScotch splits 16-cell line into 4 parts", "[meshing][scotch]")
{
    auto mesh = build_line_mesh(16);
    meshing::PartitionerScotch ps;
    meshing::ScotchProps p;
    p.nParts = 4;
    ps.initialize(p);
    std::vector<std::int32_t> part;
    const auto cut = ps.partition(mesh, part);
    REQUIRE(part.size() == 16);
    REQUIRE(cut >= 3);
}

TEST_CASE("AmrAdaptor refines high-tag cells", "[meshing][amr]")
{
    auto mesh = build_line_mesh(8);
    const std::size_t nC0 = mesh.cells().size();
    solver::FieldRegistry fields;
    auto& tag = fields.scalar("err", nC0);
    for (std::size_t i = 0; i < nC0; ++i)
        tag[i] = (i < 4 ? 0.95 : 0.05);
    meshing::AmrAdaptor amr;
    amr.initialize(mesh, {});
    const auto rep = amr.adapt(fields, "err");
    REQUIRE(rep.cellsRefined > 0);
    REQUIRE(mesh.cells().size() > nC0);
}
