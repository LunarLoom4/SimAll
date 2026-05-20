// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/materials/test_anisotropic_conductivity.cpp
// =============================================================================
#include "materials/AnisotropicConductivity.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using simall::materials::AnisotropicConductivity;
using simall::materials::Tensor3x3;

TEST_CASE("AnisotropicConductivity diagonal heat flux", "[materials][aniso]")
{
    AnisotropicConductivity k(10.0, 20.0, 30.0);
    auto q = k.heatFlux(300.0, {1.0, 2.0, 3.0}); // ∇T
    CHECK_THAT(q[0], WithinAbs(-10.0 * 1.0, 1e-12));
    CHECK_THAT(q[1], WithinAbs(-20.0 * 2.0, 1e-12));
    CHECK_THAT(q[2], WithinAbs(-30.0 * 3.0, 1e-12));
}

TEST_CASE("Tensor3x3 similarity recovers diagonal under identity", "[materials][aniso][tensor]")
{
    AnisotropicConductivity k(10.0, 20.0, 30.0);
    auto R = Tensor3x3::identity();
    auto T = k.tensorInFrame(300.0, R);
    CHECK_THAT(T.at(0, 0), WithinAbs(10.0, 1e-12));
    CHECK_THAT(T.at(1, 1), WithinAbs(20.0, 1e-12));
    CHECK_THAT(T.at(2, 2), WithinAbs(30.0, 1e-12));
    CHECK_THAT(T.at(0, 1), WithinAbs(0.0, 1e-12));
}
