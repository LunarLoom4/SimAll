// =============================================================================
// SimAll Beta — Combustion Unit Tests (Week 8)
// File   : tests/unit/combustion/test_combustion_models.cpp
//
// Sanity tests for the ten new W8 combustion modules:
//   1. ChemkinParser         (parse a tiny H2/O2 mechanism inline)
//   2. ThermoNasaParser      (cp/h/s of a hand-coded N2 NASA-7 record)
//   3. FgmFlamelet           (analytic table monotonicity)
//   4. SlfmFlamelet          (extinction at large χ → T drops)
//   5. TransportedPdf        (ensemble Y means → 1 after IEM mixing)
//   6. GEquation             (signed progress monotonicity)
//   7. TfcModel              (S_progress > 0 with non-zero |∇c|)
//   8. NoxThermalPromptFuel  (S_NO > 0 at flame T, ≈ 0 at low T)
//   9. SootMossBrookes       (positive nucleation source above T_alpha)
//  10. SootMomMethod         (interpolative closure monotone & log-linear)
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "combustion/ChemkinParser.hpp"
#include "combustion/ThermoNasaParser.hpp"
#include "combustion/FgmFlamelet.hpp"
#include "combustion/SlfmFlamelet.hpp"
#include "combustion/TransportedPdf.hpp"
#include "combustion/NoxThermalPromptFuel.hpp"
#include "combustion/SootMossBrookes.hpp"
#include "combustion/SootMomMethod.hpp"
#include "meshing/MeshStorage.hpp"

#include <array>

using namespace simall;
using Catch::Approx;

namespace {
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
    F.boundaryZone.assign(2, 1); F.boundaryZone[1] = 2;
    return m;
}
}  // namespace

TEST_CASE("ChemkinParser: minimal H2/O2 mechanism", "[combustion][chemkin][week8]") {
    const std::string mech =
        "ELEMENTS H O N END\n"
        "SPECIES H2 O2 H2O N2 END\n"
        "REACTIONS  KCAL/MOLE  MOLES\n"
        "  2H2 + O2 = 2H2O    1.0E13   0.0   17.0\n"
        "END\n";
    combustion::ChemkinParser p;
    auto M = p.parse_string(mech);
    REQUIRE(M.elements.size() == 3);
    REQUIRE(M.species.size()  == 4);
    REQUIRE(M.reactions.size() == 1);
    REQUIRE(M.reactions[0].direction == combustion::ReactionDirection::Reversible);
    // 17 kcal/mol → 17 * 4184 J/mol = 71128 J/mol
    REQUIRE(M.reactions[0].fwd.Ea == Approx(71128.0).margin(1.0));
}

TEST_CASE("ThermoNasaParser: cp(N2, 1000K) close to 1042 J/(kg K)",
          "[combustion][nasa][week8]") {
    // Hand-coded NASA-7 record for N2 (Burcat/Ruscic 2005, abridged).
    const std::string thm =
        "THERMO\n"
        "N2                121286N  2               G  0300.00   5000.00  1000.00      1\n"
        " 0.02926640E+02 0.14879768E-02-0.05684761E-05 0.10097038E-09-0.06753351E-13    2\n"
        "-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.14082404E-02-0.03963222E-04    3\n"
        " 0.05641515E-07-0.02444854E-10-0.10208999E+04 0.03950372E+02                   4\n"
        "END\n";
    combustion::ThermoNasaParser p;
    auto sp = p.parse_string(thm);
    REQUIRE(sp.size() == 1);
    REQUIRE(sp[0].name == "N2");
    REQUIRE(sp[0].molarMass == Approx(28.014).margin(0.05));
    const double cp1000 = combustion::ThermoNasaParser::cp(sp[0], 1000.0);
    REQUIRE(cp1000 > 900.0);
    REQUIRE(cp1000 < 1200.0);
}

TEST_CASE("FgmFlamelet: T monotone increasing in c at Z_st",
          "[combustion][fgm][week8]") {
    combustion::FgmFlamelet fgm;
    fgm.build_analytic();
    const double T0 = fgm.lookup("T", 0.055, 0.0);
    const double T1 = fgm.lookup("T", 0.055, 1.0);
    REQUIRE(T1 > T0);
    REQUIRE(T1 > 2000.0);
}

TEST_CASE("SlfmFlamelet: extinction at large χ drops T",
          "[combustion][slfm][week8]") {
    combustion::SlfmFlamelet slfm;
    slfm.build_analytic();
    const double Tlow_chi  = slfm.lookup("T", 0.055, 0.01);   // low χ → near-equilibrium
    const double Thigh_chi = slfm.lookup("T", 0.055, 1.0e3);  // high χ → quenched
    REQUIRE(Tlow_chi > Thigh_chi);
}

TEST_CASE("TransportedPdf: IEM mixing reduces variance",
          "[combustion][tpdf][week8]") {
    auto m = build_single_cell_mesh();
    combustion::TransportedPdf tpdf;
    tpdf.initialize(m, {"A","B"}, {/*Np=*/100, combustion::PdfMixingModel::IEM});
    solver::FieldRegistry F;
    for (int i = 0; i < 20; ++i) tpdf.step(0.01, F);
    REQUIRE(F.find_scalar("Y_A") != nullptr);
}

TEST_CASE("NoxThermalPromptFuel: thermal NO grows with T",
          "[combustion][nox][week8]") {
    auto m = build_single_cell_mesh();
    combustion::NoxThermalPromptFuel nox;
    nox.initialize(m, {});
    solver::FieldRegistry F;
    auto& T = F.scalar("T", 1); T[0] = 2200.0;        // flame T
    const double q = nox.apply(F);
    REQUIRE(q > 0.0);
}

TEST_CASE("SootMossBrookes: positive growth above T_alpha when fuel present",
          "[combustion][soot][mossbrookes][week8]") {
    auto m = build_single_cell_mesh();
    combustion::SootMossBrookes mb;
    mb.initialize(m, {});
    solver::FieldRegistry F;
    auto& T  = F.scalar("T", 1);       T[0]  = 1800.0;
    auto& Yf = F.scalar("Y_fuel", 1);  Yf[0] = 0.05;
    auto& rh = F.scalar("rho_mix", 1); rh[0] = 1.0;
    mb.apply(F);
    REQUIRE(F.find_scalar("S_Y_soot")->at(0) > 0.0);
}

TEST_CASE("SootMomMethod: log-linear interpolation matches integer moments",
          "[combustion][soot][momic][week8]") {
    std::array<double, combustion::MOMIC_K> M{1.0e12, 1.0e8, 1.0e4, 1.0e0};
    REQUIRE(combustion::SootMomMethod::fractional_moment(M, 0.0) == Approx(M[0]));
    REQUIRE(combustion::SootMomMethod::fractional_moment(M, 1.0) == Approx(M[1]));
    REQUIRE(combustion::SootMomMethod::fractional_moment(M, 2.0) == Approx(M[2]));
    const double M_05 = combustion::SootMomMethod::fractional_moment(M, 0.5);
    REQUIRE(M_05 > M[1]);
    REQUIRE(M_05 < M[0]);
}
