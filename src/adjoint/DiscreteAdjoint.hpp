// =============================================================================
// SimAll Beta - Adjoint Subsystem
// File   : src/adjoint/DiscreteAdjoint.hpp
// Week   : 18
//
// Discrete adjoint solver for steady-state design sensitivity.  Given a
// converged primal linear system  R(Q; α) = 0  and a scalar functional
// J(Q; α) (drag, total pressure loss, exit temperature, ...), the
// discrete adjoint computes
//
//        ψᵀ = ∂J/∂Q · (∂R/∂Q)⁻¹
//        dJ/dα = ∂J/∂α − ψᵀ (∂R/∂α)
//
// so that *one* adjoint solve yields the sensitivity of J with respect
// to *any* number of design variables α.  The implementation:
//
//   * Builds the transpose Jacobian Aᵀ from the user-supplied operator;
//   * Solves  Aᵀ ψ = (∂J/∂Q)ᵀ  using restarted GMRES with Jacobi
//     preconditioning (matches the in-tree solver subsystem's style);
//   * Returns the per-DV gradient vector.
//
// Industrial-CFD parity: equivalent to Fluent Adjoint Solver / STAR-CCM+
// Adjoint Flow Solver, restricted to discrete (not continuous) formulation
// in this file.  Continuous adjoint lives in ContinuousAdjoint.hpp.
// =============================================================================
#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace simall::adjoint
{

/// Sparse matrix in CSR form, just enough to drive GMRES against Aᵀ.
struct CsrMatrix
{
    std::size_t n = 0;               // square dimension
    std::vector<std::size_t> rowPtr; // size n+1
    std::vector<std::size_t> colIdx;
    std::vector<double> values;

    [[nodiscard]] std::size_t nnz() const noexcept { return values.size(); }
};

/// Sparse-matrix-vector product against the *transpose* of `A`.
/// y = Aᵀ x.  Implemented as a free function so the same CSR pair can
/// drive forward (A x) and adjoint (Aᵀ x) without duplication.
void spmv_transpose(const CsrMatrix& A, const double* x, double* y) noexcept;

/// User-supplied Jacobians.  All callables are evaluated once per design
/// cycle.  `dRdAlpha` returns the (n × m) Jacobian flattened row-major.
struct DiscreteAdjointInputs
{
    const CsrMatrix* A = nullptr; // ∂R/∂Q
    std::vector<double> dJdQ;     // size n
    std::vector<double> dJdAlpha; // size m (∂J/∂α direct)
    std::function<void(std::vector<double>& out)> dRdAlphaTransposeTimes = nullptr;
    /**< Computes  out = (∂R/∂α)ᵀ ψ.  Caller may either store the full
         m×n matrix internally or evaluate matrix-free.  `out` is sized m
         on entry. */
};

struct DiscreteAdjointResult
{
    bool ok = false;
    std::string error;
    std::size_t iterations = 0;
    double residual = 0.0;
    std::vector<double> psi;      // adjoint state, size n
    std::vector<double> gradient; // dJ/dα, size m
};

struct DiscreteAdjointOptions
{
    std::size_t maxIter = 500;
    std::size_t restart = 30;
    double tol = 1e-10;
};

[[nodiscard]] DiscreteAdjointResult solve_discrete_adjoint(const DiscreteAdjointInputs& in,
                                                           std::size_t nDesignVars,
                                                           DiscreteAdjointOptions opt = {});

} // namespace simall::adjoint
