// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/Parallel.hpp
// Phase  : 17 (MPI PARALLELIZATION) + threading facade
//
// Thin façade over MPI / OpenMP / TBB. The rest of the codebase calls into
// `simall::parallel` without ever including <mpi.h> directly. This keeps the
// build configurable (SIMALL_ENABLE_MPI / _OPENMP / _TBB).
// =============================================================================
#pragma once

#include <cstddef>
#include <functional>

namespace simall::parallel {

void initialize(int& argc, char**& argv);
void finalize() noexcept;

int  rank()       noexcept;
int  size()       noexcept;
bool is_master()  noexcept;

void barrier() noexcept;

/// Parallel for-loop façade. Backed by TBB when present, OpenMP otherwise,
/// serial fall-back when neither is enabled.
void parallel_for(std::size_t n, const std::function<void(std::size_t, std::size_t)>& body,
                  std::size_t grain = 1024);

/// All-reduce double sum (no-op outside MPI builds).
double all_reduce_sum(double local);
double all_reduce_max(double local);

}  // namespace simall::parallel
