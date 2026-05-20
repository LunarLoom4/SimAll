// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/PimpleAlgorithm.cpp
// Phase  : 22 Pass 1
// =============================================================================
#include "solver/PimpleAlgorithm.hpp"

#include <algorithm>
#include <cmath>

namespace simall::solver
{

PimpleAlgorithm::PimpleAlgorithm(meshing::Mesh& mesh,
                                 FieldRegistry& fields,
                                 const std::vector<BoundarySpec>& boundaries,
                                 ILinearSolver& momentumSolver,
                                 ILinearSolver& pressureSolver,
                                 PimpleOptions opts)
    : SimpleAlgorithm(
          mesh, fields, boundaries, momentumSolver, pressureSolver, to_simple_opts(opts))
    , pimple_(opts)
{
    if (pimple_.nOuterCorrectors < 1)
        pimple_.nOuterCorrectors = 1;
    if (pimple_.nCorrectors < 1)
        pimple_.nCorrectors = 1;
    if (pimple_.nNonOrthCorr < 0)
        pimple_.nNonOrthCorr = 0;
}

namespace
{
inline bool converged(const SimpleResiduals& r, double tolU, double tolP)
{
    if (tolU > 0.0) {
        const double maxMom =
            std::max({std::fabs(r.mom[0]), std::fabs(r.mom[1]), std::fabs(r.mom[2])});
        if (maxMom > tolU)
            return false;
    }
    if (tolP > 0.0 && std::fabs(r.cont) > tolP)
        return false;
    return (tolU > 0.0) || (tolP > 0.0); // tolerances only active when set
}
} // namespace

SimpleResiduals PimpleAlgorithm::iterate()
{
    if (!sparsity_built_)
        build_sparsity();
    SimpleResiduals res{};
    lastOuterIters_ = 0;

    auto* sync = field_synchronizer();

    for (int outer = 0; outer < pimple_.nOuterCorrectors; ++outer) {
        lastOuterIters_ = outer;

        // ---- Predictor: assemble & solve momentum once per outer ----
        compute_face_fluxes();
        if (pimple_.momentumPredictor) {
            for (int k = 0; k < 3; ++k) {
                assemble_momentum(k);
                res.mom[k] = solve_momentum_component(k);
            }
            if (sync) {
                if (auto* U = F_.find_vector("U"))
                    sync->sync_vector(U->x, U->y, U->z);
            }
        }

        // ---- PISO inner pressure correctors ----
        for (int corr = 0; corr < pimple_.nCorrectors; ++corr) {
            compute_face_fluxes();
            for (int nonOrth = 0; nonOrth <= pimple_.nNonOrthCorr; ++nonOrth) {
                assemble_pressure_correction();
                const double contRes = solve_pressure_correction_and_correct();
                if (corr == pimple_.nCorrectors - 1 && nonOrth == pimple_.nNonOrthCorr) {
                    res.cont = contRes;
                }
            }
            if (sync) {
                if (auto* p = F_.find_scalar("p"))
                    sync->sync_scalar(*p);
                if (auto* U = F_.find_vector("U"))
                    sync->sync_vector(U->x, U->y, U->z);
            }
        }

        // ---- Outer-loop early exit ----
        if (converged(res, pimple_.tolU, pimple_.tolP))
            break;
    }
    return res;
}

} // namespace simall::solver
