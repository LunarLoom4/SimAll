// =============================================================================
// SimAll Beta - tests/regression/cases/BackwardFacingStep.hpp
// BFS — Armaly et al., 1983 reattachment length.  For laminar Re < 400 the
// reattachment length grows roughly linearly with Re; we encode their
// curve fit and assert the prediction at Re=100, 200, 300 stays within
// the experimental scatter band.
// =============================================================================
#pragma once

#include "../RegressionFramework.hpp"

namespace simall::regression::bfs
{

// Armaly correlation fit:  X_r / h = 0.05 + 0.0175 · Re  for Re ≤ 400.
[[nodiscard]] inline double armaly_reattachment(double Re)
{
    return 0.05 + 0.0175 * Re;
}

[[nodiscard]] inline double error_against_experiment(double Re, double predicted)
{
    const double ref = armaly_reattachment(Re);
    return relative_error(predicted, ref);
}

} // namespace simall::regression::bfs
