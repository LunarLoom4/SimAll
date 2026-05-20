// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/DeviceMemoryPool.hpp
// Phase  : 18.2 — Bucketed device memory pool (W13).
//
// A power-of-two bucket allocator that amortises the cost of cudaMalloc /
// cudaFree across thousands of short-lived linear-algebra scratch buffers
// (CG residuals, BiCGSTAB s/t vectors, SpMV intermediates).
//
// ALGORITHM
// ---------
// Requested ``n`` bytes are rounded up to the smallest power-of-two
// ``size = 1 << k`` where ``k >= kMinShift``.  Each ``k`` maintains its
// own free-list.  ``allocate(n)``:
//   1. Round n → bucket k.
//   2. If bucket k has a free block, pop and return.
//   3. Otherwise:
//        a. If the pool's reserved memory + size <= maxBytes, fresh-alloc
//           via cudaMalloc (host malloc on non-CUDA) and return.
//        b. Else evict from larger buckets (LIFO across k+1, k+2, ...) and
//           re-attempt fresh-alloc.
//   4. If still failing, throw std::bad_alloc.
//
// ``release(p, n)`` pushes the slab onto its bucket's free-list, *without*
// fragmenting (we never split a slab — every allocation returns a whole
// power-of-two block).  Internal fragmentation is bounded at 2× (typical
// CFD workloads waste ~30 %).
//
// CORRECTNESS
// -----------
//   * The pool is process-singleton (``DeviceMemoryPool::instance()``).
//   * Thread-safe — bucket free-lists are guarded by a single mutex.
//     Cross-bucket scans (eviction) hold the lock as well.
//   * Pointers returned by ``allocate`` are 256-byte aligned (matches
//     cudaMalloc's guarantee; host fall-back uses std::aligned_alloc).
//
// NON-CUDA FALL-BACK
// ------------------
// When SIMALL_HAVE_CUDA is undefined, slabs are backed by aligned host
// memory.  This lets unit tests exercise the pool deterministically without
// a GPU.
// =============================================================================
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace simall::gpu
{

struct DeviceMemoryStats
{
    std::size_t total_reserved_bytes = 0;    ///< sum of all slabs ever allocated
    std::size_t total_in_use_bytes = 0;      ///< slabs handed out, not yet released
    std::size_t total_in_freelist_bytes = 0; ///< slabs cached in free-lists
    std::size_t slab_count = 0;
    std::size_t allocation_calls = 0; ///< cumulative count
    std::size_t freelist_hits = 0;    ///< served from cache (no cudaMalloc)
    std::size_t freelist_misses = 0;  ///< required fresh allocation
};

class DeviceMemoryPool
{
public:
    /// Minimum slab size = 2^kMinShift bytes (256 B = warp-friendly).
    static constexpr int kMinShift = 8;
    /// Maximum slab size = 2^kMaxShift bytes (1 GiB cap per slab).
    static constexpr int kMaxShift = 30;

    /// Default budget = 0 means "unlimited until cudaMalloc fails".
    explicit DeviceMemoryPool(std::size_t maxBytes = 0);
    ~DeviceMemoryPool();

    /// Process-wide pool (used by HostDeviceMirror, CgGpu, BicgstabGpu).
    static DeviceMemoryPool& instance();

    /// Returns a device-side pointer.  Size is rounded UP to the bucket
    /// granularity; pass the SAME ``bytes`` value to ``release``.
    void* allocate(std::size_t bytes);

    /// Returns a slab to its bucket free-list.  ``bytes`` MUST match the
    /// originally requested size.
    void release(void* ptr, std::size_t bytes) noexcept;

    /// Frees every cached slab (does NOT free in-use slabs).  Useful at
    /// solver shutdown to release VRAM back to the driver.
    void purge() noexcept;

    /// Snapshot of pool statistics — cheap, lock-protected.
    DeviceMemoryStats stats() const;

    /// Compute the bucket index a request would round up to.
    static int bucket_for(std::size_t bytes) noexcept;
    static std::size_t bucket_size(int k) noexcept { return std::size_t(1) << k; }

    DeviceMemoryPool(const DeviceMemoryPool&) = delete;
    DeviceMemoryPool& operator=(const DeviceMemoryPool&) = delete;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace simall::gpu
