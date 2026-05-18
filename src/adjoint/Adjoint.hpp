// Future: discrete-adjoint module. Interface only.
#pragma once
#include "solver/Solver.hpp"
namespace simall::adjoint {
class AdjointSolver {
public:
    explicit AdjointSolver(solver::Solver& primal) : primal_(primal) {}
    void compute_sensitivities();
private:
    solver::Solver& primal_;
};
}
