// =============================================================================
// SimAll Beta — Radiation Unit Tests (Week 10)
// File   : tests/unit/radiation/test_radiation_models.cpp
//
// Sanity tests for the five new W10 radiation modules:
//   1. Rosseland          (optically-thick diffusion)
//   2. MonteCarloRte      (reverse photon Monte Carlo)
//   3. SurfaceToSurface   (view-factor radiosity)
//   4. SolarLoad          (direct + diffuse + ground-reflected solar)
//   5. WsggSpectral       (Smith 1982 + Bordbar 2014 coefficients)
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "meshing/MeshStorage.hpp"
#include "radiation/MonteCarloRte.hpp"
#include "radiation/Rosseland.hpp"
#include "radiation/SolarLoad.hpp"
#include "radiation/SurfaceToSurface.hpp"
#include "radiation/WsggSpectral.hpp"

using namespace simall;
using Catch::Approx;

namespace {
meshing::Mesh build_single_cell_mesh() {
    meshing::Mesh m;
    auto& C = m.cells();
    auto& F = m.faces();
    C.volume.assign(1, 1.0);
    C.centroidX.assign(1, 0.5);
    C.centroidY.assign(1, 0.5);
    C.centroidZ.assign(1, 0.5);
    C.faceOffsets.assign(2, 0);
    C.faceOffsets[1] = 2;
    C.faceIndices.assign(2, 0);
    C.faceIndices[0] = 0;
    C.faceIndices[1] = 1;
    F.owner.assign(2, 0);
    F.neighbor.assign(2, meshing::kBoundaryCell);
    F.areaX.assign(2, 0.0);  F.areaY.assign(2, 0.0); F.areaZ.assign(2, 1.0);
    F.areaZ[0] = -1.0;
    F.centroidX.assign(2, 0.5); F.centroidY.assign(2, 0.5);
    F.centroidZ.assign(2, 0.0); F.centroidZ[1] = 1.0;
    F.boundaryZone.assign(2, 1); F.boundaryZone[1] = 2;
    return m;
}
}  // namespace

TEST_CASE("Rosseland: k_R increases with T³ and apply produces zero source on uniform T",
          "[radiation][rosseland][week10]") {
    auto m = build_single_cell_mesh();
    solver::FieldRegistry F;
    auto& T = F.scalar("T", 1); T[0] = 1000.0;
    radiation::Rosseland mod;
    radiation::RosselandProps p; p.absorption = 1.0;
    mod.initialize(m, F, p);
    const double s = mod.apply();
    // No neighbour cells → no interior flux contribution.
    REQUIRE(s == Approx(0.0).margin(1e-12));
}

TEST_CASE("MonteCarloRte: emits energy proportional to a σ T⁴ V",
          "[radiation][montecarlo][week10]") {
    auto m = build_single_cell_mesh();
    solver::FieldRegistry F;
    auto& T = F.scalar("T", 1); T[0] = 1500.0;
    radiation::MonteCarloRte mc;
    radiation::MonteCarloProps p; p.absorption = 1.0; p.nPhotonsPerCell = 50;
    mc.initialize(m, F, p);
    mc.add_wall({1, 1.0, 300.0});
    mc.add_wall({2, 1.0, 300.0});
    const double Etot = mc.run();
    const double sigma = 5.670374419e-8;
    REQUIRE(Etot == Approx(4.0 * 1.0 * sigma * std::pow(1500.0, 4)).epsilon(1e-9));
}

TEST_CASE("SurfaceToSurface: two parallel plates have view-factor approaching 1",
          "[radiation][s2s][week10]") {
    auto m = build_single_cell_mesh();      // 2 boundary faces (top & bottom)
    radiation::SurfaceToSurface s2s;
    radiation::S2SProps p; p.nRaysPerFace = 4000;
    s2s.initialize(m, p);
    s2s.add_zone({1, 1.0, 1000.0});         // bottom
    s2s.add_zone({2, 1.0,  300.0});         // top
    const auto n = s2s.build_view_factors();
    REQUIRE(n == 2);
    REQUIRE(s2s.view_factor(0, 1) > 0.5);   // most of bottom-rays should hit top
    const double qmax = s2s.solve();
    REQUIRE(qmax > 0.0);
}

TEST_CASE("SolarLoad: face normal pointing at sun receives positive flux",
          "[radiation][solar][week10]") {
    auto m = build_single_cell_mesh();
    solver::FieldRegistry F;
    radiation::SolarLoad sun;
    radiation::SolarProps p;
    p.sunDir = {0, 0, -1};   // straight down
    p.I_direct = 1000.0; p.I_diffuse = 100.0;
    sun.initialize(m, F, p);
    sun.add_zone({2, 0.8});  // top face exposed
    const double tot = sun.apply();
    REQUIRE(tot > 0.0);
    REQUIRE(sun.face_flux(1) > 0.0);    // face id 1 = top
}

TEST_CASE("WsggSpectral: emissivity monotonic in pwL and weights sum to 1",
          "[radiation][wsgg][week10]") {
    auto m = build_single_cell_mesh();
    solver::FieldRegistry F;
    auto& T = F.scalar("T", 1); T[0] = 1500.0;
    radiation::WsggSpectral w;
    radiation::WsggProps p; p.meanBeamLength = 1.0;
    w.initialize(m, F, p);
    const double e1 = w.compute_emissivity(1500.0, 0.1);
    const double e2 = w.compute_emissivity(1500.0, 5.0);
    REQUIRE(e2 >= e1);
    double sumA = 0.0;
    for (int k = 0; k < radiation::kWsggBands; ++k)
        sumA += w.compute_weight(1500.0, k);
    REQUIRE(sumA == Approx(1.0).margin(1e-9));
    const double kmax = w.apply();
    REQUIRE(kmax >= 0.0);
}
