// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/DualTimeStepping.cpp
// =============================================================================
#include "solver/DualTimeStepping.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::solver {

namespace {
double max_residual(const SimpleResiduals& r) {
    return std::max({std::abs(r.mom[0]), std::abs(r.mom[1]),
                     std::abs(r.mom[2]), std::abs(r.cont)});
}
}  // namespace

DualTimeReport DualTimeStepping::step() {
    DualTimeReport rep{};
    // First pseudo-iter rolls BDF2 history and does one outer iterate.
    SimpleResiduals res = alg_.advance_time_step(1);
    rep.iters    = 1;
    rep.finalMom = std::max({res.mom[0], res.mom[1], res.mom[2]});
    rep.finalCont= res.cont;
    const double firstRes = max_residual(res);
    if (firstRes < opt_.absTol) { rep.converged = true; return rep; }

    for (int k = 1; k < opt_.maxPseudoIters; ++k) {
        res = alg_.iterate();
        ++rep.iters;
        rep.finalMom = std::max({res.mom[0], res.mom[1], res.mom[2]});
        rep.finalCont= res.cont;
        const double cur = max_residual(res);
        if (opt_.verbose) {
            SIMALL_LOG_INFO("DualTime", "  pseudo-iter ", k,
                "  max-res = ", cur);
        }
        if (cur < opt_.absTol || cur < opt_.relTol * firstRes) {
            rep.converged = true; break;
        }
    }
    return rep;
}

DualTimeReport DualTimeStepping::solve_steady() {
    DualTimeReport rep{};
    SimpleResiduals res{};
    double firstRes = 0.0;
    for (int k = 0; k < opt_.maxPseudoIters; ++k) {
        res = alg_.iterate();
        ++rep.iters;
        rep.finalMom = std::max({res.mom[0], res.mom[1], res.mom[2]});
        rep.finalCont= res.cont;
        const double cur = max_residual(res);
        if (k == 0) firstRes = std::max(cur, 1.0e-30);
        if (opt_.verbose) {
            SIMALL_LOG_INFO("DualTime-Steady", "  iter ", k,
                "  max-res = ", cur);
        }
        if (cur < opt_.absTol || cur < opt_.relTol * firstRes) {
            rep.converged = true; break;
        }
    }
    return rep;
}

}  // namespace simall::solver
