// =============================================================================
// SimAll Beta - Verification Application
// File   : applications/simall_naca/main.cpp
// Week   : 20 — NACA-0012 thin-airfoil baseline driver.  Prints the
// thin-airfoil-theory lift slope (2π/rad) and the NACA-4 half-thickness
// distribution at quarter-chord stations.  Acceptance criterion: numerical
// post-processor must reproduce Cl(α) within 2 % for α ≤ 5°.
// =============================================================================
#include "regression/cases/NacaAirfoil.hpp"

#include <cstdio>

namespace sr = simall::regression;

int main(int argc, char** argv) {
    double thickness = 0.12;   // NACA-0012
    if (argc >= 2) thickness = std::stod(argv[1]);

    std::printf("SimAll Beta - simall_naca verification driver\n");
    std::printf("=============================================\n");
    std::printf("NACA-00%02d thickness profile (10 stations):\n",
                 int(std::lround(thickness * 100.0)));
    for (int k = 0; k <= 10; ++k) {
        const double xc = double(k) / 10.0;
        std::printf("  x/c = %4.2f   yt/c = %+.6f\n", xc,
                     sr::naca::naca4_halfthickness(xc, thickness));
    }
    std::printf("\nThin-airfoil-theory lift slope dCl/dα = 2π = %.6f /rad\n",
                 2.0 * 3.14159265358979323846);
    std::printf("Worked example: α = 5° → Cl_TAT = %.6f\n",
                 sr::naca::thin_airfoil_cl(5.0 * 3.14159265358979323846 / 180.0));
    std::printf("\nAcceptance: solver Cl(α=5°) must be within 2%% of TAT.\n");
    return 0;
}
