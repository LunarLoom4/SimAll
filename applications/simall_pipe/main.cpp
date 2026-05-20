// =============================================================================
// SimAll Beta - Verification Application
// File   : applications/simall_pipe/main.cpp
// Week   : 20 — Hagen-Poiseuille pipe-flow reference printer.  Emits the
// parabolic velocity profile, Darcy laminar friction factor, and the
// Prandtl smooth-pipe turbulent friction factor at a sweep of Re points.
// =============================================================================
#include "../../tests/regression/cases/Pipe.hpp"

#include <cstdio>

namespace sr = simall::regression;

// Inline parabolic profile (Hagen-Poiseuille): u(r) = 2·u_mean·(1 - r²/R²).
static double parabolic_profile(double r, double R, double uMean) {
    const double s = r / R;
    return 2.0 * uMean * (1.0 - s * s);
}

int main() {
    const double R     = 1.0;
    const double uMean = 1.0;

    std::printf("SimAll Beta - simall_pipe verification driver\n");
    std::printf("==============================================\n\n");
    std::printf("Parabolic profile u(r), R = %.2f, u_mean = %.2f:\n", R, uMean);
    std::printf("  %6s   %12s\n", "r/R", "u(r)");
    for (int k = 0; k <= 10; ++k) {
        const double r = double(k) * R / 10.0;
        std::printf("  %6.2f   %12.6f\n", r / R, parabolic_profile(r, R, uMean));
    }
    std::printf("\nLaminar Darcy f = 64/Re:\n");
    for (double Re : {500.0, 1000.0, 1500.0, 2000.0}) {
        std::printf("  Re = %6.0f   f = %.6f\n", Re,
                     sr::pipe::darcy_friction_factor_laminar(Re));
    }
    std::printf("\nTurbulent smooth-pipe (Prandtl) f:\n");
    for (double Re : {1.0e4, 5.0e4, 1.0e5, 5.0e5, 1.0e6}) {
        std::printf("  Re = %9.0f   f = %.6f\n", Re,
                     sr::pipe::prandtl_friction_factor(Re));
    }
    std::printf("\nAcceptance: solver f within 2%% (laminar) / 5%% (turbulent).\n");
    return 0;
}
