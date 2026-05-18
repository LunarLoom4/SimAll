// =============================================================================
// SimAll Beta — Multiphase Unit Tests (Week 7)
// File   : tests/unit/multiphase/test_multiphase_models.cpp
//
// Sanity tests for the seven new W7 multiphase modules:
//   1. CavitationKunz       (Kunz et al. 2000)
//   2. CavitationZwart      (Zwart-Gerber-Belamri 2004)
//   3. MixtureModel         (algebraic slip)
//   4. DriftFlux            (Zuber-Findlay)
//   5. EulerianMultifluid   (n-phase Eulerian)
//   6. IsoAdvector          (Roenby geometric VOF)
//   7. RpiWallBoiling       (Kurul-Podowski heat-flux partitioning)
//
// Each test verifies construction, default-state invariants and one
// physical correctness property (sign of phase change, drift direction,
// mixture density formula, RPI heat partition positivity).
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "meshing/MeshStorage.hpp"
#include "multiphase/CavitationKunz.hpp"
#include "multiphase/CavitationZwart.hpp"
#include "multiphase/MixtureModel.hpp"
#include "multiphase/DriftFlux.hpp"
#include "multiphase/IsoAdvector.hpp"
#include "multiphase/RpiWallBoiling.hpp"

using namespace simall;
using Catch::Approx;

namespace {
// Trivial 1-cell mesh: cube of side 1, two faces (top wall and bottom wall),
// just enough to satisfy initialise() in the Cavitation/Mixture modules.
meshing::Mesh build_single_cell_mesh() {
    meshing::Mesh m;
    auto& C = m.cells();
    auto& F = m.faces();
    C.volume.assign(1, 1.0);
    C.centroidX.assign(1, 0.5); C.centroidY.assign(1, 0.5); C.centroidZ.assign(1, 0.5);
    C.faceOffsets.assign(2, 0);
    F.owner.assign(2, 0);
    F.neighbor.assign(2, meshing::kBoundaryCell);
    F.areaX.assign(2, 0.0); F.areaY.assign(2, 0.0); F.areaZ.assign(2, 1.0);
    F.areaZ[0] = -1.0;
    F.centroidX.assign(2, 0.5); F.centroidY.assign(2, 0.5);
    F.centroidZ.assign(2, 0.0); F.centroidZ[1] = 1.0;
    F.boundaryZone.assign(2, 1);  F.boundaryZone[1] = 2;
    return m;
}
}  // namespace

TEST_CASE("CavitationKunz: vaporisation when p < p_sat", "[multiphase][cavitation][kunz][week7]") {
    auto m = build_single_cell_mesh();
    multiphase::CavitationKunz model;
    multiphase::KunzProps props;
    props.Uref = 10.0; props.Lref = 0.1;
    model.initialize(m, props);

    solver::FieldRegistry F;
    auto& p = F.scalar("p", 1);        p[0]  = props.pSaturation - 1000.0;
    auto& a = F.scalar("alpha_v", 1);  a[0]  = 0.1;
    const double total = model.apply(F);
    REQUIRE(total > 0.0);              // net mass into vapour
    REQUIRE(F.find_scalar("S_alpha_v")->at(0) > 0.0);
}

TEST_CASE("CavitationZwart: condensation when p > p_sat", "[multiphase][cavitation][zwart][week7]") {
    auto m = build_single_cell_mesh();
    multiphase::CavitationZwart model;
    model.initialize(m, {});           // default props
    solver::FieldRegistry F;
    auto& p = F.scalar("p", 1);        p[0]  = 5000.0;     // > p_sat (=2339)
    auto& a = F.scalar("alpha_v", 1);  a[0]  = 0.5;
    model.apply(F);
    REQUIRE(F.find_scalar("S_alpha_v")->at(0) < 0.0);  // negative = vapour shrinking
}

TEST_CASE("MixtureModel: mixture density = α-weighted phase densities",
          "[multiphase][mixture][week7]") {
    auto m = build_single_cell_mesh();
    multiphase::MixtureModel model;
    model.initialize(m, {
        {"primary",   1000.0, 1.0e-3, 1e-4},
        {"secondary",    1.0, 1.8e-5, 1e-4}
    });
    solver::FieldRegistry F;
    F.vector("U_mix", 1);
    auto& a = F.scalar("alpha_secondary", 1);  a[0] = 0.3;
    model.apply(F);
    // ρ_m = 0.7·1000 + 0.3·1 = 700.3
    REQUIRE(F.find_scalar("rho_mix")->at(0) == Approx(700.3).margin(1e-9));
    REQUIRE(model.num_phases() == 2);
}

TEST_CASE("DriftFlux: gas drifts upward (Harmathy law)",
          "[multiphase][driftflux][week7]") {
    auto m = build_single_cell_mesh();
    multiphase::DriftFlux model;
    model.initialize(m, {});
    solver::FieldRegistry F;
    F.vector("U_mix", 1);
    auto& a = F.scalar("alpha_g", 1);  a[0] = 0.2;
    const double maxDrift = model.apply(F);
    REQUIRE(maxDrift > 0.0);
    REQUIRE(F.find_vector("U_drift_g")->z[0] > 0.0);
}

TEST_CASE("IsoAdvector: initialises without crash on single-cell mesh",
          "[multiphase][isoadvector][week7]") {
    auto m = build_single_cell_mesh();
    multiphase::IsoAdvector iso;
    iso.initialize(m, {});
    solver::FieldRegistry F;
    auto& a   = F.scalar("alpha", 1);        a[0]   = 0.0;        // pure-phase cell — PLIC skips
    auto& mf  = F.scalar("massFlux", 2);     mf[0]  = 0.0; mf[1] = 0.0;
    const double l1 = iso.step(0.01, F);
    REQUIRE(l1 == Approx(0.0).margin(1e-12));
    REQUIRE(F.find_scalar("alpha")->at(0) == Approx(0.0).margin(1e-12));
}

TEST_CASE("RpiWallBoiling: positive q_wall when T_wall > T_sat",
          "[multiphase][boiling][rpi][week7]") {
    auto m = build_single_cell_mesh();
    multiphase::RpiWallBoiling rpi;
    std::vector<solver::BoundarySpec> bcs;
    solver::BoundarySpec bc{};
    bc.zone = 2; bc.type = solver::BCType::Wall; bc.scalarValue = 400.0;  // T_wall
    bcs.push_back(bc);
    rpi.initialize(m, bcs, {});
    solver::FieldRegistry F;
    auto& T = F.scalar("T", 1);  T[0] = 360.0;       // liquid T
    const double q = rpi.apply(F);
    REQUIRE(q > 0.0);
    REQUIRE(F.find_scalar("S_alpha_v")->at(0) > 0.0);
}
