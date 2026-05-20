// =============================================================================
// SimAll Beta - Solver Unit Tests
// File   : tests/unit/solver/test_map_fields.cpp
// Phase  : 23 Pass 12
//
// Coverage for the cross-mesh field mapper `simall::solver::map_fields`:
//   - Constant scalar field round-trip (refine).
//   - Constant scalar field round-trip (coarsen).
//   - Linear scalar field: IDW reproduces the linear trend within tolerance.
//   - Vector field mapping.
//   - Nearest method snaps to closest source cell.
//   - Unsized source fields are silently skipped.
//   - Fallback path triggers when target lies far outside source bbox.
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "solver/MapFields.hpp"
#include "solver/FieldRegistry.hpp"
#include "meshing/CartesianMesher.hpp"
#include "meshing/MeshStorage.hpp"

using namespace simall;
using Catch::Approx;

namespace {

meshing::Mesh brick(util::Vec3d origin, util::Vec3d extent,
                    int nx, int ny, int nz) {
    meshing::CartesianGridSpec s;
    s.origin = origin; s.extent = extent;
    s.Nx = nx; s.Ny = ny; s.Nz = nz;
    meshing::Mesh m;
    meshing::CartesianMesher{s}.generate(m);
    return m;
}

}  // namespace

// =============================================================================
TEST_CASE("Nearest mapping copies a constant scalar field across refinement",
          "[solver][mapfields][nearest]") {
    auto src = brick({0, 0, 0}, {1, 1, 1}, 1, 1, 1);  // 1 source cell
    auto tgt = brick({0, 0, 0}, {1, 1, 1}, 2, 2, 2);  // 8 target cells

    solver::FieldRegistry s, d;
    auto& T = s.scalar("T", src.cells().size());
    std::fill(T.begin(), T.end(), 42.0);

    solver::MapFieldsOptions opt;  // Nearest by default
    auto stats = solver::map_fields(src, s, tgt, d, opt);

    REQUIRE(stats.sourceCells        == 1);
    REQUIRE(stats.targetCells        == 8);
    REQUIRE(stats.scalarFieldsMapped == 1);
    REQUIRE(stats.vectorFieldsMapped == 0);
    REQUIRE(stats.fallbackQueries    == 0);

    const auto* T_out = d.find_scalar("T");
    REQUIRE(T_out != nullptr);
    REQUIRE(T_out->size() == 8);
    for (auto v : *T_out) REQUIRE(v == Approx(42.0));
}

// =============================================================================
TEST_CASE("Nearest mapping copies a constant scalar field across coarsening",
          "[solver][mapfields][nearest]") {
    auto src = brick({0, 0, 0}, {1, 1, 1}, 4, 4, 4);  // 64 source cells
    auto tgt = brick({0, 0, 0}, {1, 1, 1}, 2, 2, 2);  // 8 target cells

    solver::FieldRegistry s, d;
    auto& p = s.scalar("p", src.cells().size());
    std::fill(p.begin(), p.end(), 7.5);

    auto stats = solver::map_fields(src, s, tgt, d, {solver::MapMethod::Nearest});

    REQUIRE(stats.scalarFieldsMapped == 1);
    const auto* p_out = d.find_scalar("p");
    REQUIRE(p_out != nullptr);
    for (auto v : *p_out) REQUIRE(v == Approx(7.5));
}

// =============================================================================
TEST_CASE("Inverse-distance mapping reproduces a linear field accurately",
          "[solver][mapfields][idw]") {
    // Source is a fine 1D-like line; target is a coarser overlapping line.
    auto src = brick({0, 0, 0}, {1, 0.1, 0.1}, 10, 1, 1);
    auto tgt = brick({0, 0, 0}, {1, 0.1, 0.1}, 5,  1, 1);

    solver::FieldRegistry s, d;
    auto& T = s.scalar("T", src.cells().size());
    // f(x) = 3*x_centroid
    for (std::size_t c = 0; c < src.cells().size(); ++c) {
        T[c] = 3.0 * src.cells().centroidX[c];
    }

    solver::MapFieldsOptions opt;
    opt.method     = solver::MapMethod::InverseDistance;
    opt.kNeighbors = 4;
    opt.idwPower   = 2.0;
    auto stats = solver::map_fields(src, s, tgt, d, opt);
    REQUIRE(stats.scalarFieldsMapped == 1);

    const auto* T_out = d.find_scalar("T");
    REQUIRE(T_out != nullptr);
    REQUIRE(T_out->size() == tgt.cells().size());

    // IDW on a perfectly linear field reproduces the trend to a few % once
    // we are not on the boundary cell. Check interior cells exactly.
    for (std::size_t t = 0; t < tgt.cells().size(); ++t) {
        const double expected = 3.0 * tgt.cells().centroidX[t];
        REQUIRE((*T_out)[t] == Approx(expected).margin(0.25));
    }
}

// =============================================================================
TEST_CASE("Vector field mapping propagates all three components",
          "[solver][mapfields][vector]") {
    auto src = brick({0, 0, 0}, {1, 1, 1}, 2, 2, 2);
    auto tgt = brick({0, 0, 0}, {1, 1, 1}, 2, 2, 2);

    solver::FieldRegistry s, d;
    auto& U = s.vector("U", src.cells().size());
    for (std::size_t c = 0; c < src.cells().size(); ++c) {
        U.x[c] = 1.0;
        U.y[c] = 2.0;
        U.z[c] = 3.0;
    }

    auto stats = solver::map_fields(src, s, tgt, d, {solver::MapMethod::Nearest});
    REQUIRE(stats.vectorFieldsMapped == 1);

    const auto* U_out = d.find_vector("U");
    REQUIRE(U_out != nullptr);
    REQUIRE(U_out->size() == tgt.cells().size());
    for (std::size_t c = 0; c < U_out->size(); ++c) {
        REQUIRE(U_out->x[c] == Approx(1.0));
        REQUIRE(U_out->y[c] == Approx(2.0));
        REQUIRE(U_out->z[c] == Approx(3.0));
    }
}

// =============================================================================
TEST_CASE("Nearest mapping snaps each target cell to its closest source cell",
          "[solver][mapfields][nearest]") {
    // 2-cell source with distinctive values; identical-geometry target.
    auto src = brick({0, 0, 0}, {2, 1, 1}, 2, 1, 1);  // cells at x={0.5, 1.5}
    auto tgt = brick({0, 0, 0}, {2, 1, 1}, 2, 1, 1);

    solver::FieldRegistry s, d;
    auto& T = s.scalar("T", src.cells().size());
    T[0] = 10.0;  // owns x in [0,1]
    T[1] = 99.0;  // owns x in [1,2]

    solver::map_fields(src, s, tgt, d, {solver::MapMethod::Nearest});

    const auto* T_out = d.find_scalar("T");
    REQUIRE(T_out != nullptr);
    REQUIRE((*T_out)[0] == Approx(10.0));
    REQUIRE((*T_out)[1] == Approx(99.0));
}

// =============================================================================
TEST_CASE("Source fields with wrong size are silently skipped",
          "[solver][mapfields][robust]") {
    auto src = brick({0, 0, 0}, {1, 1, 1}, 2, 2, 2);
    auto tgt = brick({0, 0, 0}, {1, 1, 1}, 1, 1, 1);

    solver::FieldRegistry s, d;
    s.scalar("undersized", src.cells().size() - 1);   // wrong size
    auto& ok = s.scalar("ok", src.cells().size());
    std::fill(ok.begin(), ok.end(), 5.0);

    auto stats = solver::map_fields(src, s, tgt, d);
    REQUIRE(stats.scalarFieldsMapped == 1);
    REQUIRE(d.find_scalar("ok")        != nullptr);
    REQUIRE(d.find_scalar("undersized") == nullptr);
}

// =============================================================================
TEST_CASE("Fallback scan fires when target lies far outside source bbox",
          "[solver][mapfields][fallback]") {
    auto src = brick({0, 0, 0},   {1, 1, 1}, 1, 1, 1);  // tiny source
    auto tgt = brick({100, 0, 0}, {1, 1, 1}, 1, 1, 1);  // far away target

    solver::FieldRegistry s, d;
    auto& T = s.scalar("T", src.cells().size());
    T[0] = 13.0;

    auto stats = solver::map_fields(src, s, tgt, d);
    REQUIRE(stats.fallbackQueries    == 1);
    REQUIRE(stats.scalarFieldsMapped == 1);

    const auto* T_out = d.find_scalar("T");
    REQUIRE(T_out != nullptr);
    REQUIRE((*T_out)[0] == Approx(13.0));  // only one source cell available
}
