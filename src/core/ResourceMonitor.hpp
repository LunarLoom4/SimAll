// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/ResourceMonitor.hpp
// Phase  : 1.3 (APPLICATION CORE → Resource Telemetry)
//
// Cross-platform process and GPU memory sampling. Returns a `Sample` struct
// with RSS (working set), virtual size, current process CPU time, and an
// optional list of per-GPU memory totals (free/used/total) collected via
// NVML when the runtime library is present. NVML is loaded dynamically so
// builds without CUDA / NVIDIA tooling work unmodified.
//
// Sampling is cheap (<100 µs on Linux) and may be called every step from
// the solver monitor without measurable overhead.
// =============================================================================
#pragma once

#include <chrono>
#include <cstdint>
#include <string>
#include <vector>

namespace simall::core
{

struct GpuMemorySample
{
    std::string name;
    std::uint64_t totalBytes = 0;
    std::uint64_t freeBytes = 0;
    std::uint64_t usedBytes = 0;
};

struct ResourceSample
{
    std::uint64_t rssBytes = 0;
    std::uint64_t virtualBytes = 0;
    double userCpuSeconds = 0.0;
    double kernelCpuSeconds = 0.0;
    std::vector<GpuMemorySample> gpus;
    std::chrono::system_clock::time_point capturedAt;
};

class ResourceMonitor
{
public:
    /// Collect a sample of the current process.
    static ResourceSample sample();

    /// Whether NVML is available (loaded and initialised). Cheap to call.
    static bool nvmlAvailable();
};

} // namespace simall::core
