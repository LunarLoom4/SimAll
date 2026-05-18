// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/materials/test_janaf.cpp
//
// Round-trip test of the NASA-7 / JANAF parser using a hand-crafted record
// for H2 (typical CHEMKIN thermo format).
// =============================================================================
#include "materials/Janaf.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sstream>

using namespace simall::materials;
using Catch::Matchers::WithinRel;

TEST_CASE("Janaf parses a NASA-7 record (H2)", "[materials][janaf]") {
    // CHEMKIN thermo format -- one species record (4 lines, each 80 cols).
    // Columns:  name(1-18)  date(25-30)  comp(25-44)  phase(45)
    //           Tlow(46-55)  Thigh(56-65)  Tmid(66-73)
    const char* h2 =
        "H2                L 7/88H   2.   0.   0.   0.G   200.000  3500.000 1000.000    1\n"
        " 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2\n"
        "-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3\n"
        " 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01                   4\n";

    std::istringstream in(h2);
    JanafThermoParser p;
    JanafSpecies sp;
    REQUIRE(p.parseRecord(in, sp));
    REQUIRE(sp.name == "H2");
    REQUIRE(sp.Tmin == 200.0);
    REQUIRE(sp.Tmid == 1000.0);
    REQUIRE(sp.Tmax == 3500.0);

    // Assign a molecular weight so cp/h/s return per-kg quantities.
    sp.molecularWeight = 0.002016;   // kg/mol

    // cp/R at 1000 K must equal (a1 + a2*T + ... )_low.
    // For the low coeffs above:  a1=2.34433112, a2=7.98052075e-3,
    // a3=-1.94781510e-5, a4=2.01572094e-8, a5=-7.37611761e-12.
    const double T = 1000.0;
    const double cpOverR =
          2.34433112
        + 7.98052075e-3 * T
        - 1.94781510e-5 * T*T
        + 2.01572094e-8 * T*T*T
        - 7.37611761e-12 * T*T*T*T;
    const double expected = cpOverR * 8.314462618 / sp.molecularWeight;
    CHECK_THAT(sp.cp(T), WithinRel(expected, 1.0e-9));
}
