// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/NumaPinning.cpp
// =============================================================================
#include "parallel/NumaPinning.hpp"
#include "core/Logger.hpp"

#include <sstream>
#include <thread>

#if defined(_WIN32)
#  define NOMINMAX
#  include <windows.h>
#elif defined(__linux__)
#  include <pthread.h>
#  include <sched.h>
#  include <unistd.h>
#  if defined(SIMALL_HAVE_NUMA)
#    include <numa.h>
#  endif
#endif

namespace simall::parallel {

NumaTopology NumaPinning::discover() {
    NumaTopology t{};
    t.totalCores = static_cast<int>(std::thread::hardware_concurrency());
    if (t.totalCores <= 0) t.totalCores = 1;
#if defined(__linux__) && defined(SIMALL_HAVE_NUMA)
    if (numa_available() >= 0) {
        t.numaNodes = numa_max_node() + 1;
        t.coresPerNode.assign(t.numaNodes, 0);
        for (int cpu = 0; cpu < t.totalCores; ++cpu) {
            const int node = numa_node_of_cpu(cpu);
            if (node >= 0 && node < t.numaNodes) ++t.coresPerNode[node];
        }
        return t;
    }
#endif
    t.numaNodes = 1;
    t.coresPerNode = {t.totalCores};
    return t;
}

bool NumaPinning::pin_thread_to_core(int cpuIndex) {
    if (cpuIndex < 0) return false;
#if defined(_WIN32)
    DWORD_PTR mask = static_cast<DWORD_PTR>(1) << (cpuIndex & 63);
    return SetThreadAffinityMask(GetCurrentThread(), mask) != 0;
#elif defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpuIndex, &set);
    return pthread_setaffinity_np(pthread_self(), sizeof(set), &set) == 0;
#else
    (void)cpuIndex;
    return false;
#endif
}

std::vector<int> NumaPinning::compute_thread_to_core_map(int nThreads) {
    if (nThreads <= 0) return {};
    auto topo = discover();
    if (topo.numaNodes <= 0 || topo.totalCores <= 0) {
        std::vector<int> v(nThreads);
        for (int i = 0; i < nThreads; ++i) v[i] = i;
        return v;
    }
    // Build CPU list grouped by NUMA node (assumes contiguous index ranges
    // per node — typical Linux layout; falls back gracefully otherwise).
    std::vector<int> cpus; cpus.reserve(topo.totalCores);
    int next = 0;
    for (int node = 0; node < topo.numaNodes; ++node) {
        const int n = topo.coresPerNode[node];
        for (int k = 0; k < n; ++k) cpus.push_back(next + k);
        next += n;
    }
    if (cpus.empty()) {
        for (int i = 0; i < topo.totalCores; ++i) cpus.push_back(i);
    }
    std::vector<int> assignment(nThreads, -1);
    const int nNodes = topo.numaNodes;
    for (int t = 0; t < nThreads; ++t) {
        const int node = t % nNodes;
        // Compute base CPU of this node.
        int base = 0;
        for (int k = 0; k < node && k < static_cast<int>(topo.coresPerNode.size()); ++k)
            base += topo.coresPerNode[k];
        const int slot = (t / nNodes) %
            std::max(1, node < static_cast<int>(topo.coresPerNode.size())
                       ? topo.coresPerNode[node] : 1);
        assignment[t] = (base + slot) % static_cast<int>(cpus.size());
    }
    return assignment;
}

bool NumaPinning::pin_self_by_index(int threadIndex, int nThreads) {
    auto map = compute_thread_to_core_map(nThreads);
    if (threadIndex < 0 || threadIndex >= static_cast<int>(map.size())) return false;
    return pin_thread_to_core(map[threadIndex]);
}

std::string NumaPinning::summarise(const NumaTopology& topo) {
    std::ostringstream os;
    os << "NumaTopology{ totalCores=" << topo.totalCores
       << " nodes="     << topo.numaNodes << " [";
    for (std::size_t i = 0; i < topo.coresPerNode.size(); ++i)
        os << (i ? "," : "") << topo.coresPerNode[i];
    os << "] }";
    return os.str();
}

}  // namespace simall::parallel
