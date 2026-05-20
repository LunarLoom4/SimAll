// =============================================================================
// SimAll Beta - tests/regression/test_regression_suite.cpp
// Week 19 - Drives all 10 canonical CFD regression cases.  Each TEST_CASE
// verifies the reference benchmark numbers encoded in cases/*.hpp.  These
// are the gating checks the solver subsystem acceptance tests must beat
// (within their own per-case tolerance) before they can graduate to the
// nightly cross-platform CI.  Tagged [regression] for ctest -L filtering.
// =============================================================================
#include "cases/BackwardFacingStep.hpp"
#include "cases/Channel.hpp"
#include "cases/Cylinder.hpp"
#include "cases/FlameD.hpp"
#include "cases/LidDrivenCavity.hpp"
#include "cases/NacaAirfoil.hpp"
#include "cases/Pipe.hpp"
#include "cases/RayleighBenard.hpp"
#include "cases/Shocktube.hpp"
#include "cases/TaylorGreenVortex.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
namespace sr = simall::regression;

TEST_CASE("LDC Re=100 centreline back-flow magnitude", "[regression][ldc]")
{
    REQUIRE_THAT(sr::ldc::max_backflow_velocity(), WithinAbs(0.21090, 1e-4));
    // Lid u=1 is recovered at the top boundary sample.
    REQUIRE_THAT(sr::ldc::kGhiaU_Re100.back(), WithinAbs(1.0, 1e-12));
    // Centreline speed integral matches its tabulated value (golden number).
    const double s = sr::ldc::centerline_speed_integral();
    REQUIRE(s > 0.10);
    REQUIRE(s < 0.20);
}

TEST_CASE("BFS Armaly reattachment correlation", "[regression][bfs]")
{
    // Re=200 → X_r/h = 3.55, Re=300 → 5.30.
    REQUIRE_THAT(sr::bfs::armaly_reattachment(200.0), WithinAbs(3.55, 1e-6));
    REQUIRE_THAT(sr::bfs::armaly_reattachment(300.0), WithinAbs(5.30, 1e-6));
    REQUIRE(sr::bfs::error_against_experiment(100.0, 1.8) < 0.05);
}

TEST_CASE("Channel Reichardt profile matches log-law in inner region", "[regression][channel]")
{
    // At y+ = 30 the log-law and Reichardt should agree to <10 % (Reichardt
    // includes the additive viscous-buffer correction that vanishes asymptotically).
    const double up = sr::channel::reichardt_u_plus(30.0);
    const double ll = sr::log_law(30.0);
    REQUIRE(sr::relative_error(up, ll) < 0.10);
    REQUIRE(sr::channel::reichardt_u_plus(1e-3) < 1e-2);
}

TEST_CASE("TGV 2-D analytic kinetic-energy decay", "[regression][tgv]")
{
    const double E0 = 1.0, nu = 0.01;
    const double E1 = sr::tgv::analytic_kinetic_energy_2d(0.5, nu);
    const double E2 = sr::tgv::analytic_kinetic_energy_2d(1.0, nu);
    // exp(-2·0.01·2·0.5) = exp(-0.02), exp(-0.04)
    REQUIRE_THAT(E1, WithinAbs(std::exp(-0.02) * E0, 1e-9));
    REQUIRE_THAT(E2, WithinAbs(std::exp(-0.04) * E0, 1e-9));
    REQUIRE(sr::tgv::kBrachetPeakDissipation > 0.0);
}

TEST_CASE("NACA thin-airfoil lift slope", "[regression][naca]")
{
    const double cl = sr::naca::thin_airfoil_cl(5.0 * 3.14159265 / 180.0);
    REQUIRE_THAT(cl, WithinAbs(2.0 * 3.14159265 * (5.0 * 3.14159265 / 180.0), 1e-9));
    // NACA-0012 leading-edge thickness is positive and goes to 0 at the trailing edge.
    REQUIRE(sr::naca::naca4_halfthickness(0.5, 0.12) > 0.0);
    REQUIRE(std::abs(sr::naca::naca4_halfthickness(1.0, 0.12)) < 1e-2);
}

TEST_CASE("Sod shock-tube reference values", "[regression][shocktube]")
{
    auto r = sr::shocktube::sod_reference();
    REQUIRE_THAT(r.pStar, WithinAbs(0.30313, 1e-4));
    REQUIRE_THAT(r.uStar, WithinAbs(0.92745, 1e-4));
    REQUIRE_THAT(r.rhoLeftStar, WithinAbs(0.42632, 1e-4));
    // Speed of sound on the left (γ=1.4, ρ=1, p=1) is sqrt(1.4) ≈ 1.1832.
    REQUIRE_THAT(sr::shocktube::speed_of_sound(1.4, 1.0, 1.0), WithinAbs(1.1832, 1e-3));
}

TEST_CASE("Rayleigh-Bénard onset", "[regression][rb]")
{
    REQUIRE(sr::rb::linear_growth_rate(2000.0) > 0.0);
    REQUIRE(sr::rb::linear_growth_rate(1000.0) < 0.0);
    REQUIRE_THAT(sr::rb::linear_growth_rate(sr::rb::kRayleighCritical), WithinAbs(0.0, 1e-12));
}

TEST_CASE("Sandia Flame D centreline scaling", "[regression][flameD]")
{
    REQUIRE_THAT(sr::flameD::centerline_mixture_fraction(50.0), WithinAbs(5.4 / 50.0, 1e-12));
    REQUIRE(sr::flameD::centerline_mixture_fraction(0.0) == 1.0);
    REQUIRE(sr::flameD::kPeakTemperatureK > 1500.0);
}

TEST_CASE("Cylinder shedding correlations", "[regression][cylinder]")
{
    // At Re=100 the Roshko fit gives St ≈ 0.165.
    REQUIRE_THAT(sr::cylinder::roshko_strouhal(100.0), WithinAbs(0.165, 1e-2));
    // At Re=100, Henderson fit gives Cd ≈ 1 + 10·100^{-2/3} ≈ 1.465.
    const double cd = sr::cylinder::henderson_drag(100.0);
    REQUIRE(cd > 1.2);
    REQUIRE(cd < 1.7);
}

TEST_CASE("Pipe friction factor laminar+turbulent", "[regression][pipe]")
{
    REQUIRE_THAT(sr::pipe::darcy_friction_factor_laminar(2000.0), WithinAbs(64.0 / 2000.0, 1e-12));
    // Prandtl/Colebrook smooth-pipe at Re=1e5 → f ≈ 0.0178.
    const double f = sr::pipe::prandtl_friction_factor(1e5);
    REQUIRE(f > 0.010);
    REQUIRE(f < 0.030);
    // Parabolic profile recovers u_max = 2·u_mean at the centreline.
    REQUIRE_THAT(sr::pipe_profile(0.0, 1.0, 1.0), WithinAbs(2.0, 1e-12));
    // Zero at the wall.
    REQUIRE_THAT(sr::pipe_profile(1.0, 1.0, 1.0), WithinAbs(0.0, 1e-12));
}
