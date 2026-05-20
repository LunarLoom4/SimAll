// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/LinearSolvers.hpp
// Phase  : 7 — Krylov solvers with preconditioning.
//
// Concrete implementations:
//   - JacobiCG         (symmetric positive-definite systems, e.g. p-Poisson)
//   - GMRESm           (non-symmetric, restart parameter m)
//   - BiCGSTAB         (non-symmetric, no restart)
//
// Preconditioners (applied as M^-1):
//   - JacobiPreconditioner
//   - ILU0Preconditioner   (zero-fill incomplete LU, in-place CSR factorization)
//
// All solvers honour the LinearSolverConfig contract from Solver.hpp.
// =============================================================================
#pragma once

#include "CSRMatrix.hpp"
#include "Solver.hpp"

#include <memory>

namespace simall::solver
{

// ---------------- Preconditioners -----------------------------------------
class IPreconditioner
{
public:
    virtual ~IPreconditioner() = default;
    virtual void setup(const CSRMatrix& A) = 0;
    virtual void apply(const util::aligned_vector<double>& r,
                       util::aligned_vector<double>& z) const = 0;
};

std::unique_ptr<IPreconditioner> make_preconditioner(PreconditionerKind k);

// ---------------- Solvers (factory in Solver.cpp dispatches here) ---------
std::unique_ptr<ILinearSolver> make_gmres(LinearSolverConfig cfg);
std::unique_ptr<ILinearSolver> make_bicgstab(LinearSolverConfig cfg);
std::unique_ptr<ILinearSolver> make_cg(LinearSolverConfig cfg);
std::unique_ptr<ILinearSolver> make_tfqmr(LinearSolverConfig cfg);

} // namespace simall::solver
