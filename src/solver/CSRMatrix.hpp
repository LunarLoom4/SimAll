// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/CSRMatrix.hpp
// Phase  : 10 (SPARSE MATRIX SYSTEM)  — Section 10 of ultra-detailed spec.
//
// Compressed Sparse Row. Aligned arrays for SIMD SpMV. Thread-safe assembly
// is provided by SolverAssembler (atomic-free per-row chunking).
// =============================================================================
#pragma once

#include "utilities/AlignedAllocator.hpp"

namespace simall::solver
{

struct CSRMatrix
{
    util::aligned_vector<int> rowPtr;
    util::aligned_vector<int> colIdx;
    util::aligned_vector<double> values;

    std::size_t rows() const noexcept { return rowPtr.empty() ? 0 : rowPtr.size() - 1; }

    /// y = A * x  (single-threaded reference; parallel SpMV in solver.cpp)
    void spmv(const util::aligned_vector<double>& x, util::aligned_vector<double>& y) const
    {
        const std::size_t n = rows();
        y.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            double s = 0.0;
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k)
                s += values[k] * x[colIdx[k]];
            y[i] = s;
        }
    }
};

} // namespace simall::solver
