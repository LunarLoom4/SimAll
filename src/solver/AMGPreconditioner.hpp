// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/AMGPreconditioner.hpp
// Phase  : 7.4 — Classical Algebraic Multigrid (Ruge-Stüben).
//
// Implements:
//   - Strong-connection coarsening on |a_ij| > θ · max_{k≠i}|a_ik|
//   - Direct interpolation P (Ruge-Stüben 1987 §6)
//   - Galerkin coarse operator A_c = P^T A_f P
//   - Symmetric Gauss-Seidel smoothing (pre + post sweeps)
//   - V-cycle recursion down to nCoarse ≤ 32 (small-system direct solve)
//
// Used as an `IPreconditioner` inside Krylov solvers (GMRES, BiCGSTAB, CG).
// Works on arbitrary symmetric or weakly-non-symmetric CSR systems produced
// by the finite-volume discretisation (pressure-correction, energy,
// scalar-transport). For strongly non-symmetric systems we recommend
// preconditioner = ILU.
//
// References:
//   Ruge & Stüben, "Algebraic Multigrid", in Multigrid Methods (1987).
//   Trottenberg, Oosterlee & Schüller, "Multigrid", Academic Press, 2001.
// =============================================================================
#pragma once

#include "CSRMatrix.hpp"
#include "LinearSolvers.hpp"

#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <vector>

namespace simall::solver
{

class AMGPreconditioner final : public IPreconditioner
{
public:
    explicit AMGPreconditioner(double theta = 0.25,
                               int maxLevels = 10,
                               int preSweeps = 1,
                               int postSweeps = 1,
                               int coarseThreshold = 32);
    void setup(const CSRMatrix& A) override;
    void apply(const util::aligned_vector<double>& r,
               util::aligned_vector<double>& z) const override;

private:
    struct Level
    {
        CSRMatrix A;                          // operator at this level
        CSRMatrix P;                          // prolongation level→level-1 (coarse→fine)
        CSRMatrix R;                          // restriction = P^T
        util::aligned_vector<double> diagInv; // 1/A_ii for SGS
    };
    std::vector<Level> levels_;
    // Direct solve at coarsest level (dense L^TL Cholesky-like).
    std::vector<double> coarseLU_;
    int coarseN_ = 0;

    double theta_;
    int maxLevels_;
    int preSweeps_;
    int postSweeps_;
    int coarseThreshold_;

    // Helpers
    void build_strong(const CSRMatrix& A, std::vector<std::vector<int>>& S) const;
    void coarsen(const std::vector<std::vector<int>>& S,
                 std::vector<int>& cf) const; // 1 = C-point, 0 = F-point
    CSRMatrix build_prolongation(const CSRMatrix& A,
                                 const std::vector<std::vector<int>>& S,
                                 const std::vector<int>& cf,
                                 std::vector<int>& coarseId) const;
    CSRMatrix transpose(const CSRMatrix& M, int newRows) const;
    CSRMatrix triple_product_RAP(const CSRMatrix& R, const CSRMatrix& A, const CSRMatrix& P) const;
    void sgs_sweep(const Level& L,
                   util::aligned_vector<double>& x,
                   const util::aligned_vector<double>& b) const;
    void v_cycle(int lvl,
                 util::aligned_vector<double>& x,
                 const util::aligned_vector<double>& b) const;
    void direct_solve(util::aligned_vector<double>& x, const util::aligned_vector<double>& b) const;
};

} // namespace simall::solver
