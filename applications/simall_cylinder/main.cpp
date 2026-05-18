// =============================================================================
// SimAll Beta - Verification Application
// File   : applications/simall_cylinder/main.cpp
// Week   : 20 — Circular-cylinder shedding reference printer.  Tabulates
// the Roshko Strouhal-Re curve and the Henderson drag fit at the
// canonical Reynolds-number points used by the DES regression suite.
// =============================================================================
#include "regression/cases/Cylinder.hpp"

#include <cstdio>

namespace sr = simall::regression;

int main() {
    std::printf("SimAll Beta - simall_cylinder verification driver\n");
    std::printf("==================================================\n\n");
    std::printf("  %8s   %12s   %12s\n", "Re", "St (Roshko)", "Cd (Henderson)");
    for (double Re : {40.0, 60.0, 80.0, 100.0, 150.0, 200.0, 300.0}) {
        std::printf("  %8.1f   %12.5f   %12.5f\n", Re,
                     sr::cylinder::roshko_strouhal(Re),
                     sr::cylinder::henderson_drag(Re));
    }
    std::printf("\nAcceptance: DES result within 3%% on St and 5%% on Cd.\n");
    return 0;
}
