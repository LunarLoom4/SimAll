// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/FieldReducer.hpp
// Phase  : 17.3 — typed cross-rank reductions for solver fields.
//
// Used by the linear-solver residual-norm computation, convergence
// monitors, post-processing min/max queries, and conservation checks.
// All operations consider only the *owned* portion of the field
// (`[0, nOwned)`) so that ghost cells are NOT double-counted across
// ranks.
// =============================================================================
#pragma once

#include "parallel/MpiContext.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <cstddef>

namespace simall::parallel {

class FieldReducer {
public:
    explicit FieldReducer(MpiContext& ctx) : ctx_(ctx) {}

    /// Sum of owned values across ranks.
    double sum (const util::aligned_vector<double>& v, std::size_t nOwned) const;
    /// L2 norm (sqrt(Σ v²)) over the owned cells of all ranks.
    double l2  (const util::aligned_vector<double>& v, std::size_t nOwned) const;
    /// Minimum value across all owned cells of all ranks.
    double min (const util::aligned_vector<double>& v, std::size_t nOwned) const;
    /// Maximum value across all owned cells of all ranks.
    double max (const util::aligned_vector<double>& v, std::size_t nOwned) const;
    /// Arithmetic mean (sum/Ntotal).
    double mean(const util::aligned_vector<double>& v, std::size_t nOwned) const;

    /// Vector-field L2 norm: sqrt(Σ (vx² + vy² + vz²)).
    double l2_vector(const util::aligned_vector<double>& vx,
                     const util::aligned_vector<double>& vy,
                     const util::aligned_vector<double>& vz,
                     std::size_t nOwned) const;

    /// In-place all-reduce of a scalar (driver for tolerance comparisons).
    double reduce_sum_scalar(double local) const { return ctx_.allreduce_sum(local); }
    double reduce_max_scalar(double local) const { return ctx_.allreduce_max(local); }

private:
    MpiContext& ctx_;
};

}  // namespace simall::parallel
