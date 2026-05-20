// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/PisoAlgorithm.cpp
// =============================================================================
#include "solver/PisoAlgorithm.hpp"

#include <algorithm>
#include <cmath>

namespace simall::solver
{

PisoAlgorithm::PisoAlgorithm(meshing::Mesh& mesh,
                             FieldRegistry& fields,
                             const std::vector<BoundarySpec>& boundaries,
                             ILinearSolver& momentumSolver,
                             ILinearSolver& pressureSolver,
                             PisoOptions opts)
    : SimpleAlgorithm(
          mesh, fields, boundaries, momentumSolver, pressureSolver, to_simple_opts(opts))
    , piso_(opts)
{
    if (piso_.nCorrectors < 1)
        piso_.nCorrectors = 1;
}

SimpleResiduals PisoAlgorithm::iterate()
{
    if (!sparsity_built_)
        build_sparsity();
    SimpleResiduals res{};

    // ---- Predictor: assemble & solve momentum once ----
    compute_face_fluxes();
    for (int k = 0; k < 3; ++k) {
        assemble_momentum(k);
        res.mom[k] = solve_momentum_component(k);
    }
    // Halo exchange after momentum predictor (W12).
    auto* sync = field_synchronizer();
    if (sync) {
        if (auto* U = F_.find_vector("U"))
            sync->sync_vector(U->x, U->y, U->z);
    }

    // ---- Pressure correctors ----
    for (int corr = 0; corr < piso_.nCorrectors; ++corr) {
        // Recompute face fluxes with the latest U
        compute_face_fluxes();
        // Pressure-correction equation. (Non-orthogonal correctors loop
        // updates the explicit non-orthogonal gradient inside
        // correct_fields, which we approximate by re-running the
        // assemble/solve a few extra times.)
        for (int nonOrth = 0; nonOrth <= piso_.nNonOrthCorr; ++nonOrth) {
            assemble_pressure_correction();
            const double contRes = solve_pressure_correction_and_correct();
            if (corr == piso_.nCorrectors - 1 && nonOrth == piso_.nNonOrthCorr) {
                res.cont = contRes;
            }
        }
        // Halo exchange after each pressure corrector.
        if (sync) {
            if (auto* p = F_.find_scalar("p"))
                sync->sync_scalar(*p);
            if (auto* U = F_.find_vector("U"))
                sync->sync_vector(U->x, U->y, U->z);
        }
    }
    return res;
}

} // namespace simall::solver
