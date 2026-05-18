// =============================================================================
// SimAll Beta — Turbulence Unit Tests (Week 6)
// File   : tests/unit/turbulence/test_turbulence_models.cpp
//
// Sanity tests covering the eight new Week-6 turbulence closures:
//   1. KEpsilonRNG          (Yakhot-Orszag)
//   2. KEpsilonRealizable   (Shih)
//   3. KOmegaStandard       (Wilcox 2006)
//   4. DynamicSmagorinsky   (Germano-Lilly)
//   5. SA-IDDES             (Shur-Spalart-Strelets-Travin)
//   6. SAS-SST              (Menter-Egorov)
//   7. LRR-RSM              (Launder-Reece-Rodi)
//   8. k-kL-ω transition    (Walters-Cokljat)
//
// Each test pulls the model from the TurbulenceRegistry and verifies that:
//   - the registry returns a non-null factory result,
//   - the model reports its expected name,
//   - the model's μ_t is non-negative before initialisation (default 0).
// Full functional tests (planar-channel decay, anisotropy invariants, IDDES
// shielding) require a fully-coupled flow case and are covered by the
// verification-suite (tests/verification/).
// =============================================================================
#include <catch2/catch_test_macros.hpp>

#include "turbulence/ITurbulenceModel.hpp"
#include "turbulence/KEpsilonRng.hpp"
#include "turbulence/KEpsilonRealizable.hpp"
#include "turbulence/KOmegaStandard.hpp"
#include "turbulence/DynamicSmagorinsky.hpp"
#include "turbulence/IDDES.hpp"
#include "turbulence/SasSst.hpp"
#include "turbulence/LRR_Reynolds_Stress.hpp"
#include "turbulence/KKLOmegaTransition.hpp"

using namespace simall::turbulence;

TEST_CASE("Week 6 turbulence models — registry contains all 8 new closures",
          "[turbulence][registry][week6]") {
    auto& reg = TurbulenceRegistry::instance();
    for (const char* key : {
        "kEpsilonRNG", "kEpsilonRealizable", "kOmegaStandard",
        "dynamicSmagorinsky", "SA-IDDES", "SAS-SST", "LRR-RSM",
        "kKLOmegaTransition"})
    {
        auto m = reg.create(key);
        REQUIRE(m);
        REQUIRE_FALSE(m->name().empty());
        // μ_t default before any initialisation must be 0 (safe value for momentum).
        REQUIRE(m->turbulent_viscosity(0) == 0.0);
    }
}

TEST_CASE("KEpsilonRNG and Realizable variants report distinct names",
          "[turbulence][week6][k-epsilon]") {
    KEpsilonRng_Full        rng;
    KEpsilonRealizable_Full real;
    REQUIRE(rng.name()  == "kEpsilonRNG");
    REQUIRE(real.name() == "kEpsilonRealizable");
}

TEST_CASE("KOmegaStandard reports Wilcox name and has direct construction",
          "[turbulence][week6][k-omega]") {
    KOmegaStandard_Full m;
    REQUIRE(m.name() == "kOmegaStandard");
    REQUIRE(m.turbulent_viscosity(0) == 0.0);
}

TEST_CASE("DynamicSmagorinsky exposes csMax setter",
          "[turbulence][week6][dynamic-smag]") {
    DynamicSmagorinsky_LES m;
    REQUIRE(m.name() == "dynamicSmagorinsky");
    m.set_cs_max(0.20);                         // physical upper bound
    REQUIRE(m.turbulent_viscosity(99) == 0.0);  // unset cells return 0
}

TEST_CASE("IDDES and SAS-SST and LRR distinct names",
          "[turbulence][week6][hybrid][rsm]") {
    IDDES_Full                    iddes;
    SasSst_Full                   sas;
    LRR_Reynolds_Stress_Full      lrr;
    KKLOmegaTransition_Full       kklw;
    REQUIRE(iddes.name() == "SA-IDDES");
    REQUIRE(sas  .name() == "SAS-SST");
    REQUIRE(lrr  .name() == "LRR-RSM");
    REQUIRE(kklw .name() == "kKLOmegaTransition");
}
