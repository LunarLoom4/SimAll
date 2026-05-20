// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/NumaPinning.hpp
// Phase  : 17.5 — per-thread NUMA / CPU affinity binding.
//
// On large multi-socket nodes (typical HPC: 2-4 NUMA domains × 24-48
// cores), failure to pin worker threads to consistent cores causes:
//   * cross-socket memory accesses costing 100 + ns extra,
//   * cache-line ping-pong on shared atomics,
//   * uneven OpenMP thread layout that defeats parallel_for grain
//     tuning.
//
// NumaPinning exposes a portable thin API that uses:
//   * Win32  SetThreadGroupAffinity / GetLogicalProcessorInformation
//   * POSIX  sched_setaffinity + sysconf(_SC_NPROCESSORS_ONLN)
//   * NUMA   numa_alloc_onnode (libnuma) when SIMALL_HAVE_NUMA is on
//
// Where back-end facilities are unavailable, the function returns
// `false` without aborting — pinning is always a soft optimization.
// =============================================================================
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace simall::parallel
{

struct NumaTopology
{
    int totalCores = 0;
    int numaNodes = 0;
    std::vector<int> coresPerNode; // size = numaNodes
};

class NumaPinning
{
public:
    /// Discover the host's NUMA layout.  Returns a best-effort estimate
    /// (one node + all cores) when no NUMA back-end is available.
    static NumaTopology discover();

    /// Pin the calling thread to a specific logical CPU index.
    /// Returns true on success, false on platform unsupported.
    static bool pin_thread_to_core(int cpuIndex);

    /// Distribute `nThreads` workers across `nNodes` NUMA nodes (round-
    /// robin) and return the per-thread core assignment.
    static std::vector<int> compute_thread_to_core_map(int nThreads);

    /// Convenience: pin the calling thread according to its index.  Uses
    /// `compute_thread_to_core_map` to compute the assignment lazily.
    static bool pin_self_by_index(int threadIndex, int nThreads);

    /// Human-readable topology summary (for log dumps).
    static std::string summarise(const NumaTopology& topo);
};

} // namespace simall::parallel
