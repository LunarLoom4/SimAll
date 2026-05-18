// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/Parallel.cpp
//
// Minimal serial implementation. MPI/OpenMP/TBB back-ends are switched on by
// compile-time defines wired from the root CMakeLists (SIMALL_ENABLE_*).
// =============================================================================
#include "parallel/Parallel.hpp"

#if defined(SIMALL_HAVE_MPI)
#  include <mpi.h>
#endif
#if defined(SIMALL_HAVE_OPENMP)
#  include <omp.h>
#endif
#if defined(SIMALL_HAVE_TBB)
#  include <tbb/parallel_for.h>
#  include <tbb/blocked_range.h>
#endif

namespace simall::parallel {

namespace { int g_rank = 0; int g_size = 1; }

void initialize(int& argc, char**& argv) {
#if defined(SIMALL_HAVE_MPI)
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_size);
#else
    (void)argc; (void)argv;
#endif
}

void finalize() noexcept {
#if defined(SIMALL_HAVE_MPI)
    MPI_Finalize();
#endif
}

int  rank()      noexcept { return g_rank; }
int  size()      noexcept { return g_size; }
bool is_master() noexcept { return g_rank == 0; }

void barrier() noexcept {
#if defined(SIMALL_HAVE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void parallel_for(std::size_t n,
                  const std::function<void(std::size_t, std::size_t)>& body,
                  std::size_t grain) {
#if defined(SIMALL_HAVE_TBB)
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, n, grain),
        [&](const tbb::blocked_range<std::size_t>& r) { body(r.begin(), r.end()); });
#elif defined(SIMALL_HAVE_OPENMP)
    const std::size_t step = grain;
    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(n); i += static_cast<std::ptrdiff_t>(step)) {
        std::size_t lo = static_cast<std::size_t>(i);
        std::size_t hi = std::min(lo + step, n);
        body(lo, hi);
    }
#else
    (void)grain;
    body(0, n);
#endif
}

double all_reduce_sum(double local) {
#if defined(SIMALL_HAVE_MPI)
    double g = 0; MPI_Allreduce(&local, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return g;
#else
    return local;
#endif
}
double all_reduce_max(double local) {
#if defined(SIMALL_HAVE_MPI)
    double g = 0; MPI_Allreduce(&local, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return g;
#else
    return local;
#endif
}

}  // namespace simall::parallel
