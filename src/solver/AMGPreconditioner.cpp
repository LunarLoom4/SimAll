// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/AMGPreconditioner.cpp
// =============================================================================
#include "solver/AMGPreconditioner.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_map>

namespace simall::solver
{

AMGPreconditioner::AMGPreconditioner(double theta, int maxLevels, int pre, int post, int coarseThr)
    : theta_(theta)
    , maxLevels_(maxLevels)
    , preSweeps_(pre)
    , postSweeps_(post)
    , coarseThreshold_(coarseThr)
{
}

// ============================================================ helpers
namespace
{
inline int find_col(const CSRMatrix& A, int row, int col)
{
    for (int k = A.rowPtr[row]; k < A.rowPtr[row + 1]; ++k)
        if (A.colIdx[k] == col)
            return k;
    return -1;
}
inline double diag(const CSRMatrix& A, int row)
{
    const int k = find_col(A, row, row);
    return (k >= 0) ? A.values[k] : 0.0;
}
} // namespace

void AMGPreconditioner::build_strong(const CSRMatrix& A, std::vector<std::vector<int>>& S) const
{
    const int n = static_cast<int>(A.rows());
    S.assign(n, {});
    for (int i = 0; i < n; ++i) {
        // Negative-coupling strength: a_ii * theta criterion (assumes diagonally
        // dominant matrix with negative off-diagonals — true for diffusion-like
        // operators after sign convention).
        double maxOff = 0.0;
        for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k) {
            if (A.colIdx[k] == i)
                continue;
            maxOff = std::max(maxOff, -A.values[k]);
        }
        if (maxOff <= 0.0)
            continue;
        const double thr = theta_ * maxOff;
        for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k) {
            const int j = A.colIdx[k];
            if (j == i)
                continue;
            if (-A.values[k] >= thr)
                S[i].push_back(j);
        }
    }
}

void AMGPreconditioner::coarsen(const std::vector<std::vector<int>>& S, std::vector<int>& cf) const
{
    // Standard Ruge-Stüben (RS-1) coarsening: lambda_i = |S_i^T| + 2 |U_i^T|.
    // Greedy selection picks highest-lambda points as C, marks their S-neighbours
    // as F. We approximate with S_i^T = {j : i ∈ S_j}.
    const int n = static_cast<int>(S.size());
    cf.assign(n, -1);                    // -1 = undecided
    std::vector<std::vector<int>> St(n); // transpose of S
    for (int i = 0; i < n; ++i)
        for (int j : S[i])
            St[j].push_back(i);

    std::vector<int> lambda(n, 0);
    for (int i = 0; i < n; ++i)
        lambda[i] = static_cast<int>(St[i].size());

    // Repeatedly pick the highest-lambda undecided node.
    while (true) {
        int best = -1, bestLam = -1;
        for (int i = 0; i < n; ++i)
            if (cf[i] < 0 && lambda[i] > bestLam) {
                bestLam = lambda[i];
                best = i;
            }
        if (best < 0)
            break;
        cf[best] = 1; // C-point
        // Each strongly-influenced neighbour becomes F.
        for (int j : St[best]) {
            if (cf[j] < 0) {
                cf[j] = 0;
                // Boost lambda of points that influence j and are undecided.
                for (int k : S[j])
                    if (cf[k] < 0)
                        ++lambda[k];
            }
        }
        lambda[best] = -1;
    }
    // Mark any leftover as C (safety net).
    for (int i = 0; i < n; ++i)
        if (cf[i] < 0)
            cf[i] = 1;
}

CSRMatrix AMGPreconditioner::build_prolongation(const CSRMatrix& A,
                                                const std::vector<std::vector<int>>& /*S*/,
                                                const std::vector<int>& cf,
                                                std::vector<int>& coarseId) const
{
    // Direct interpolation (Ruge-Stüben 1987):
    //   if i is C : P_ii = 1
    //   if i is F : P_ij = -a_ij / a_ii  for each C neighbour j
    const int nf = static_cast<int>(A.rows());
    coarseId.assign(nf, -1);
    int nc = 0;
    for (int i = 0; i < nf; ++i)
        if (cf[i] == 1)
            coarseId[i] = nc++;

    CSRMatrix P;
    P.rowPtr.assign(nf + 1, 0);
    for (int i = 0; i < nf; ++i) {
        if (cf[i] == 1) {
            P.colIdx.push_back(coarseId[i]);
            P.values.push_back(1.0);
        } else {
            const double dii = diag(A, i);
            const double inv = (std::abs(dii) > 1e-30) ? 1.0 / dii : 0.0;
            for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k) {
                const int j = A.colIdx[k];
                if (j == i)
                    continue;
                if (cf[j] != 1)
                    continue; // only C-neighbours
                P.colIdx.push_back(coarseId[j]);
                P.values.push_back(-A.values[k] * inv);
            }
        }
        P.rowPtr[i + 1] = static_cast<int>(P.colIdx.size());
    }
    return P;
}

CSRMatrix AMGPreconditioner::transpose(const CSRMatrix& M, int newRows) const
{
    const int nr = static_cast<int>(M.rows());
    CSRMatrix T;
    T.rowPtr.assign(newRows + 1, 0);
    for (int k = 0; k < static_cast<int>(M.colIdx.size()); ++k)
        ++T.rowPtr[M.colIdx[k] + 1];
    for (int i = 0; i < newRows; ++i)
        T.rowPtr[i + 1] += T.rowPtr[i];
    T.colIdx.assign(M.colIdx.size(), 0);
    T.values.assign(M.values.size(), 0.0);
    std::vector<int> cursor(newRows, 0);
    for (int i = 0; i < nr; ++i) {
        for (int k = M.rowPtr[i]; k < M.rowPtr[i + 1]; ++k) {
            const int c = M.colIdx[k];
            const int dst = T.rowPtr[c] + cursor[c]++;
            T.colIdx[dst] = i;
            T.values[dst] = M.values[k];
        }
    }
    return T;
}

CSRMatrix AMGPreconditioner::triple_product_RAP(const CSRMatrix& R,
                                                const CSRMatrix& A,
                                                const CSRMatrix& P) const
{
    // First AP = A * P (fine x coarse).
    const int nf = static_cast<int>(A.rows());
    const int nc = static_cast<int>(
        P.colIdx.empty() ? 0 : *std::max_element(P.colIdx.begin(), P.colIdx.end()) + 1);

    CSRMatrix AP;
    AP.rowPtr.assign(nf + 1, 0);
    std::vector<double> rowAcc(nc, 0.0);
    std::vector<int> rowMarker(nc, -1);
    for (int i = 0; i < nf; ++i) {
        std::vector<int> cols;
        for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k) {
            const int j = A.colIdx[k];
            const double aij = A.values[k];
            for (int kp = P.rowPtr[j]; kp < P.rowPtr[j + 1]; ++kp) {
                const int c = P.colIdx[kp];
                const double v = aij * P.values[kp];
                if (rowMarker[c] != i) {
                    rowMarker[c] = i;
                    rowAcc[c] = v;
                    cols.push_back(c);
                } else {
                    rowAcc[c] += v;
                }
            }
        }
        std::sort(cols.begin(), cols.end());
        for (int c : cols) {
            AP.colIdx.push_back(c);
            AP.values.push_back(rowAcc[c]);
        }
        AP.rowPtr[i + 1] = static_cast<int>(AP.colIdx.size());
    }

    // Then RAP = R * AP (coarse x coarse).
    CSRMatrix RAP;
    RAP.rowPtr.assign(nc + 1, 0);
    std::vector<double> rowAcc2(nc, 0.0);
    std::vector<int> rowMarker2(nc, -1);
    for (int i = 0; i < nc; ++i) {
        std::vector<int> cols;
        for (int k = R.rowPtr[i]; k < R.rowPtr[i + 1]; ++k) {
            const int j = R.colIdx[k];
            const double rij = R.values[k];
            for (int kp = AP.rowPtr[j]; kp < AP.rowPtr[j + 1]; ++kp) {
                const int c = AP.colIdx[kp];
                const double v = rij * AP.values[kp];
                if (rowMarker2[c] != i) {
                    rowMarker2[c] = i;
                    rowAcc2[c] = v;
                    cols.push_back(c);
                } else {
                    rowAcc2[c] += v;
                }
            }
        }
        std::sort(cols.begin(), cols.end());
        for (int c : cols) {
            RAP.colIdx.push_back(c);
            RAP.values.push_back(rowAcc2[c]);
        }
        RAP.rowPtr[i + 1] = static_cast<int>(RAP.colIdx.size());
    }
    return RAP;
}

// ============================================================ setup
void AMGPreconditioner::setup(const CSRMatrix& A)
{
    levels_.clear();
    Level top;
    top.A = A;
    const int n0 = static_cast<int>(A.rows());
    top.diagInv.assign(n0, 0.0);
    for (int i = 0; i < n0; ++i) {
        const double d = diag(A, i);
        top.diagInv[i] = (std::abs(d) > 1e-30) ? 1.0 / d : 0.0;
    }
    levels_.push_back(std::move(top));

    for (int lvl = 0; lvl < maxLevels_ - 1; ++lvl) {
        const CSRMatrix& Af = levels_.back().A;
        const int nf = static_cast<int>(Af.rows());
        if (nf <= coarseThreshold_)
            break;
        std::vector<std::vector<int>> S;
        build_strong(Af, S);
        std::vector<int> cf;
        coarsen(S, cf);
        const int nc = static_cast<int>(std::count(cf.begin(), cf.end(), 1));
        if (nc == 0 || nc >= nf)
            break; // could not coarsen further
        std::vector<int> coarseId;
        CSRMatrix P = build_prolongation(Af, S, cf, coarseId);
        CSRMatrix R = transpose(P, nc);
        CSRMatrix Ac = triple_product_RAP(R, Af, P);

        levels_.back().P = std::move(P);
        levels_.back().R = std::move(R);

        Level next;
        next.A = std::move(Ac);
        next.diagInv.assign(nc, 0.0);
        for (int i = 0; i < nc; ++i) {
            const double d = diag(next.A, i);
            next.diagInv[i] = (std::abs(d) > 1e-30) ? 1.0 / d : 0.0;
        }
        levels_.push_back(std::move(next));
    }

    // Factorise the coarsest operator as a dense LU for the direct solve.
    coarseN_ = static_cast<int>(levels_.back().A.rows());
    coarseLU_.assign(coarseN_ * coarseN_, 0.0);
    const auto& Ac = levels_.back().A;
    for (int i = 0; i < coarseN_; ++i)
        for (int k = Ac.rowPtr[i]; k < Ac.rowPtr[i + 1]; ++k)
            coarseLU_[i * coarseN_ + Ac.colIdx[k]] = Ac.values[k];
    // In-place Gauss elimination with partial pivoting; pivoting indices are
    // implicit (we permute rows in-place). Stored as triangular factor.
    for (int k = 0; k < coarseN_; ++k) {
        // Partial pivot.
        int piv = k;
        double pivVal = std::abs(coarseLU_[k * coarseN_ + k]);
        for (int i = k + 1; i < coarseN_; ++i) {
            const double v = std::abs(coarseLU_[i * coarseN_ + k]);
            if (v > pivVal) {
                pivVal = v;
                piv = i;
            }
        }
        if (pivVal < 1e-30)
            continue; // singular row; skip (Krylov handles)
        if (piv != k) {
            for (int j = 0; j < coarseN_; ++j)
                std::swap(coarseLU_[k * coarseN_ + j], coarseLU_[piv * coarseN_ + j]);
        }
        const double dinv = 1.0 / coarseLU_[k * coarseN_ + k];
        for (int i = k + 1; i < coarseN_; ++i) {
            const double f = coarseLU_[i * coarseN_ + k] * dinv;
            coarseLU_[i * coarseN_ + k] = f;
            for (int j = k + 1; j < coarseN_; ++j)
                coarseLU_[i * coarseN_ + j] -= f * coarseLU_[k * coarseN_ + j];
        }
    }
    SIMALL_LOG_INFO("AMG", "Built ", levels_.size(), " levels (coarsest=", coarseN_, ")");
}

void AMGPreconditioner::sgs_sweep(const Level& L,
                                  util::aligned_vector<double>& x,
                                  const util::aligned_vector<double>& b) const
{
    const auto& A = L.A;
    const int n = static_cast<int>(A.rows());
    // Forward sweep
    for (int i = 0; i < n; ++i) {
        double s = b[i];
        double aii = 0.0;
        for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k) {
            const int j = A.colIdx[k];
            if (j == i)
                aii = A.values[k];
            else
                s -= A.values[k] * x[j];
        }
        if (std::abs(aii) > 1e-30)
            x[i] = s / aii;
    }
    // Backward sweep
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        double aii = 0.0;
        for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k) {
            const int j = A.colIdx[k];
            if (j == i)
                aii = A.values[k];
            else
                s -= A.values[k] * x[j];
        }
        if (std::abs(aii) > 1e-30)
            x[i] = s / aii;
    }
}

void AMGPreconditioner::direct_solve(util::aligned_vector<double>& x,
                                     const util::aligned_vector<double>& b) const
{
    const int n = coarseN_;
    if (n == 0)
        return;
    util::aligned_vector<double> y(n, 0.0);
    // Forward substitution (unit lower triangular implicit on the diagonal).
    for (int i = 0; i < n; ++i) {
        double s = b[i];
        for (int j = 0; j < i; ++j)
            s -= coarseLU_[i * n + j] * y[j];
        y[i] = s;
    }
    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        double s = y[i];
        for (int j = i + 1; j < n; ++j)
            s -= coarseLU_[i * n + j] * x[j];
        const double d = coarseLU_[i * n + i];
        x[i] = (std::abs(d) > 1e-30) ? s / d : 0.0;
    }
}

void AMGPreconditioner::v_cycle(int lvl,
                                util::aligned_vector<double>& x,
                                const util::aligned_vector<double>& b) const
{
    if (lvl == static_cast<int>(levels_.size()) - 1) {
        direct_solve(x, b);
        return;
    }
    const Level& L = levels_[lvl];
    // Pre-smoothing
    for (int s = 0; s < preSweeps_; ++s)
        sgs_sweep(L, x, b);

    // Residual r = b - A x
    util::aligned_vector<double> Ax(x.size(), 0.0);
    L.A.spmv(x, Ax);
    util::aligned_vector<double> r(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        r[i] = b[i] - Ax[i];

    // Restrict r -> rc
    const int nc = static_cast<int>(L.R.rows());
    util::aligned_vector<double> rc(nc, 0.0);
    L.R.spmv(r, rc);

    // Coarse solve recursion
    util::aligned_vector<double> ec(nc, 0.0);
    v_cycle(lvl + 1, ec, rc);

    // Prolongate ec -> ef; x += P ec
    util::aligned_vector<double> ef(x.size(), 0.0);
    L.P.spmv(ec, ef);
    for (std::size_t i = 0; i < x.size(); ++i)
        x[i] += ef[i];

    // Post-smoothing
    for (int s = 0; s < postSweeps_; ++s)
        sgs_sweep(L, x, b);
}

void AMGPreconditioner::apply(const util::aligned_vector<double>& r,
                              util::aligned_vector<double>& z) const
{
    if (levels_.empty()) {
        z = r;
        return;
    }
    z.assign(r.size(), 0.0);
    v_cycle(0, z, r);
}

} // namespace simall::solver
