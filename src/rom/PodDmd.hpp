// =============================================================================
// SimAll Beta - Reduced-Order Modelling Subsystem
// File   : src/rom/PodDmd.hpp
// Phase  : 21 — Snapshot-matrix Proper Orthogonal Decomposition (POD) and
// Dynamic Mode Decomposition (DMD).
//
// POD: given snapshots X = [x_1 … x_M] ∈ R^{N×M}, compute
//      modes Φ ∈ R^{N×r} from the eigendecomposition of X^T X (method of
//      snapshots, Sirovich 1987) and the projection coefficients
//      a_i(t) = Φ^T x_i.
//
// DMD: standard exact DMD of Tu et al. (2014):
//      Y = X' X^+,  via reduced SVD X = U Σ V^T ⇒  Ã = U^T X' V Σ^{-1}
//      eigenvalues & eigenvectors of Ã give continuous-time growth rates
//      and dynamic modes Φ_DMD = X' V Σ^{-1} W.
//
// All Eigen-free: pure dense LAPACK-style routines implemented in-house.
// =============================================================================
#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace simall::rom
{

class SnapshotMatrix
{
public:
    explicit SnapshotMatrix(std::size_t N = 0) : N_(N) {}

    void clear() { cols_.clear(); }
    /// Append a snapshot vector (length must equal N at construction).
    void append(const std::vector<double>& x);

    std::size_t rows() const noexcept { return N_; }
    std::size_t cols() const noexcept { return cols_.size(); }
    double at(std::size_t i, std::size_t j) const { return cols_[j][i]; }
    const std::vector<double>& column(std::size_t j) const { return cols_[j]; }

private:
    std::size_t N_;
    std::vector<std::vector<double>> cols_;
};

struct PodResult
{
    std::vector<std::vector<double>> modes;        // N × r
    std::vector<double> singularValues;            // length r
    std::vector<std::vector<double>> coefficients; // r × M  (a(t))
};

struct DmdResult
{
    std::vector<std::vector<std::complex<double>>> modes; // N × r
    std::vector<std::complex<double>> eigenvalues;        // discrete-time
    std::vector<std::complex<double>> amplitudes;         // initial fit
};

/// Compute POD using the method of snapshots. r is the truncation rank
/// (0 ⇒ keep all). Returns at most min(M, r) modes.
PodResult compute_pod(const SnapshotMatrix& X, std::size_t r = 0);

/// Compute exact DMD with rank truncation r on the snapshot matrix X.
/// X must have at least 2 columns.
DmdResult compute_dmd(const SnapshotMatrix& X, std::size_t r = 0);

} // namespace simall::rom
