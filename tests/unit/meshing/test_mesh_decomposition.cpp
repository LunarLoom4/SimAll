// =============================================================================
// SimAll Beta - Meshing Unit Tests
// File   : tests/unit/meshing/test_mesh_decomposition.cpp
// Phase  : 23 Pass 20
//
// Coverage for ops::extract_subdomain:
//   - owned-only path on a 2x2x1 brick split half/half along x
//   - ghost-layer halo doubles the cell count on the same split
//   - rank with zero owned cells yields an empty subdomain
//   - cellRank size mismatch throws std::invalid_argument
//   - interfaceFaces equals the cut-plane face count
// =============================================================================
#include <catch2/catch_test_macros.hpp>

#include "meshing/MeshOps.hpp"
#include "meshing/CartesianMesher.hpp"
#include "meshing/MeshStorage.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

using namespace simall;

namespace {

meshing::Mesh build_brick(int Nx, int Ny, int Nz) {
    meshing::CartesianGridSpec spec;
    spec.origin = {0, 0, 0};
    spec.extent = {static_cast<double>(Nx),
                   static_cast<double>(Ny),
                   static_cast<double>(Nz)};
    spec.Nx = Nx; spec.Ny = Ny; spec.Nz = Nz;
    meshing::Mesh m;
    meshing::CartesianMesher{spec}.generate(m);
    return m;
}

// Assign rank 0 to cells with centroid.x < midX, rank 1 to the rest.
std::vector<std::int32_t> rank_split_x(const meshing::Mesh& m, double midX) {
    const auto& C = m.cells();
    std::vector<std::int32_t> r(C.size(), 0);
    for (std::size_t c = 0; c < C.size(); ++c) {
        r[c] = (C.centroidX[c] < midX) ? 0 : 1;
    }
    return r;
}

}  // namespace

// =============================================================================
TEST_CASE("extract_subdomain (owned only) on a 2x2x1 brick halved along x",
          "[meshing][decomposition]") {
    auto m = build_brick(2, 2, 1);            // 4 cells, 2 per side of x=1
    const auto rk = rank_split_x(m, 1.0);

    meshing::Mesh sub;
    auto stats = meshing::ops::extract_subdomain(m, rk, /*rank=*/0,
                                                 /*ghost=*/false, sub);

    REQUIRE(stats.ownedCells     == 2);
    REQUIRE(stats.ghostCells     == 0);
    REQUIRE(stats.interfaceFaces == 2);       // 2 cut-plane faces (z=1, y=2)
    REQUIRE(sub.cells().size()   == 2);
}

// =============================================================================
TEST_CASE("extract_subdomain ghost layer includes one-deep halo",
          "[meshing][decomposition]") {
    auto m = build_brick(2, 2, 1);
    const auto rk = rank_split_x(m, 1.0);

    meshing::Mesh sub;
    auto stats = meshing::ops::extract_subdomain(m, rk, /*rank=*/0,
                                                 /*ghost=*/true, sub);

    REQUIRE(stats.ownedCells     == 2);
    REQUIRE(stats.ghostCells     == 2);       // the other half becomes halo
    REQUIRE(stats.interfaceFaces == 2);
    REQUIRE(sub.cells().size()   == 4);       // owned + ghost
}

// =============================================================================
TEST_CASE("extract_subdomain on an absent rank yields an empty subdomain",
          "[meshing][decomposition]") {
    auto m = build_brick(2, 2, 1);
    std::vector<std::int32_t> rk(m.cells().size(), 0);  // every cell on rank 0

    meshing::Mesh sub;
    auto stats = meshing::ops::extract_subdomain(m, rk, /*rank=*/7,
                                                 /*ghost=*/true, sub);

    REQUIRE(stats.ownedCells     == 0);
    REQUIRE(stats.ghostCells     == 0);
    REQUIRE(stats.interfaceFaces == 0);
    REQUIRE(sub.cells().size()   == 0);
}

// =============================================================================
TEST_CASE("extract_subdomain rejects cellRank with wrong size",
          "[meshing][decomposition]") {
    auto m = build_brick(2, 2, 1);
    std::vector<std::int32_t> rk(m.cells().size() + 3, 0);   // bad size

    meshing::Mesh sub;
    REQUIRE_THROWS_AS(
        meshing::ops::extract_subdomain(m, rk, 0, true, sub),
        std::invalid_argument);
}
