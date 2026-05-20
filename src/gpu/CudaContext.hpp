// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/CudaContext.hpp
// Phase  : 18.1 — CUDA device context (W13).
//
// PURPOSE
// -------
// A process-wide handle that owns the CUDA driver/runtime resources every
// other GPU module (DeviceMemoryPool, StreamScheduler, GpuKernels, CgGpu,
// BicgstabGpu, RenderingBridge) routes through.  Wraps:
//
//   * cudaSetDevice / cudaGetDeviceProperties / cudaMemGetInfo
//   * cuBLAS handle  (cublasCreate / cublasSetStream)
//   * cuSPARSE handle (cusparseCreate / cusparseSetStream)
//   * Last-error inspection (cudaPeekAtLastError + cudaGetErrorString)
//
// DESIGN RULES
// ------------
//   1. NO <cuda_runtime.h> in this header.  Concrete CUDA types are kept
//      behind void* opaque handles; the .cpp gates everything on
//      SIMALL_HAVE_CUDA (defined by the build system when the CUDA toolkit
//      was located and ``SIMALL_ENABLE_CUDA`` is ON).
//   2. ``CudaContext::current()`` always returns a valid object; on a
//      non-CUDA build (or where no device is visible) ``is_available()``
//      returns ``false`` and the queries return sane defaults (compute
//      capability 0.0, zero memory).
//   3. Headers in this subsystem inherit the rule above — the GPU layer can
//      be linked from non-CUDA code without dragging in toolkit headers.
//
// CONCURRENCY
// -----------
//   * Construction / destruction is single-threaded (called from main()).
//   * ``synchronize()`` and ``select_device()`` are NOT thread-safe; call
//     from the driver thread only.
//   * Queries (``device_id``, ``properties``) are read-only and safe from
//     any thread once the context has been constructed.
// =============================================================================
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace simall::gpu
{

/// Compile-time + run-time CUDA availability.
///
///  ``SIMALL_HAVE_CUDA`` (compile-time):  set by build system when the
///                                        toolkit was located.
///  ``is_cuda_available()`` (run-time):   true when the host actually has
///                                        at least one usable device.
bool is_cuda_available() noexcept;

/// Description of a single visible device (mirrors the subset of
/// ``cudaDeviceProp`` SimAll actually consults).
struct DeviceProperties
{
    int device_id = -1;
    std::string name = "none";
    int compute_capability_major = 0;
    int compute_capability_minor = 0;
    int sm_count = 0; ///< multiProcessorCount
    int warp_size = 32;
    int max_threads_per_block = 0;
    std::size_t total_memory_bytes = 0;
    std::size_t free_memory_bytes = 0;
    std::size_t shared_memory_per_block = 0;
    bool unified_addressing = false; ///< Pascal+
    bool managed_memory = false;     ///< cudaMallocManaged supported
};

/// Process-wide context.  Owns one cuBLAS + one cuSPARSE handle bound to the
/// active device (re-bound by ``select_device``).
class CudaContext
{
public:
    /// Returns the singleton.  First call performs lazy initialisation:
    /// picks device 0 (or whichever environment variable
    /// ``CUDA_VISIBLE_DEVICES`` exposes), creates the cuBLAS / cuSPARSE
    /// handles, and queries properties.
    static CudaContext& current();

    /// Number of devices visible to the process.  0 means "no CUDA" / fall
    /// back to host kernels.
    int device_count() const noexcept;

    /// Currently active device id, or -1 on non-CUDA builds.
    int device_id() const noexcept;

    /// Switch the active device.  Re-creates cuBLAS / cuSPARSE handles
    /// (handles are device-affine).  No-op on non-CUDA builds.
    /// Returns the previously active device id.
    int select_device(int id);

    /// Query device properties (cached after first call per device).
    DeviceProperties properties() const;

    /// Refreshes ``free_memory_bytes`` via cudaMemGetInfo.  Cheap (a few µs).
    void refresh_memory_info();

    /// cudaDeviceSynchronize.  No-op on non-CUDA builds.
    void synchronize() noexcept;

    /// Pops the last sticky CUDA error (cudaPeekAtLastError +
    /// cudaGetErrorString).  Returns "ok" when no error is pending.
    std::string last_error() const;

    /// Opaque accessors for sub-modules (cast to ``cublasHandle_t`` /
    /// ``cusparseHandle_t`` inside .cu/.cpp gated by SIMALL_HAVE_CUDA).
    /// Returns nullptr on non-CUDA builds.
    void* cublas_handle() const noexcept;
    void* cusparse_handle() const noexcept;

    CudaContext(const CudaContext&) = delete;
    CudaContext& operator=(const CudaContext&) = delete;
    CudaContext(CudaContext&&) = delete;
    CudaContext& operator=(CudaContext&&) = delete;

private:
    CudaContext();
    ~CudaContext();

    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace simall::gpu
