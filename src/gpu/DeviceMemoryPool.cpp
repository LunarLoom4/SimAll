// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/DeviceMemoryPool.cpp
// Phase  : 18.2 — Bucketed device memory pool (W13).
// =============================================================================
#include "gpu/DeviceMemoryPool.hpp"
#include "gpu/CudaContext.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cstdlib>
#include <mutex>
#include <new>
#include <stdexcept>
#include <vector>

#ifdef SIMALL_HAVE_CUDA
    #include <cuda_runtime.h>
#endif

#if defined(_WIN32)
    #include <malloc.h>   // _aligned_malloc / _aligned_free
#endif

namespace simall::gpu {

namespace {

/// 256-byte alignment — matches the cudaMalloc contract and is friendly to
/// AVX-512 loads in the host fall-back.
constexpr std::size_t kAlign = 256;

inline void* aligned_host_alloc(std::size_t bytes) {
    if (bytes == 0) bytes = kAlign;
    // Round size up to alignment (some std::aligned_alloc impls require it).
    const std::size_t padded = (bytes + kAlign - 1) & ~(kAlign - 1);
#if defined(_WIN32)
    return _aligned_malloc(padded, kAlign);
#else
    return std::aligned_alloc(kAlign, padded);
#endif
}

inline void aligned_host_free(void* p) noexcept {
#if defined(_WIN32)
    _aligned_free(p);
#else
    std::free(p);
#endif
}

inline void* raw_alloc(std::size_t bytes) {
#ifdef SIMALL_HAVE_CUDA
    if (is_cuda_available()) {
        void* p = nullptr;
        if (cudaMalloc(&p, bytes) != cudaSuccess || p == nullptr) return nullptr;
        return p;
    }
#endif
    return aligned_host_alloc(bytes);
}

inline void raw_free(void* p) noexcept {
    if (!p) return;
#ifdef SIMALL_HAVE_CUDA
    if (is_cuda_available()) { cudaFree(p); return; }
#endif
    aligned_host_free(p);
}

}  // namespace

// -----------------------------------------------------------------------------
struct DeviceMemoryPool::Impl {
    mutable std::mutex                   mtx;
    std::size_t                          maxBytes = 0;
    DeviceMemoryStats                    stats{};
    /// One LIFO free-list per bucket k ∈ [kMinShift, kMaxShift].
    std::vector<std::vector<void*>>      freelist;

    Impl() : freelist(DeviceMemoryPool::kMaxShift + 1) {}
};

// -----------------------------------------------------------------------------
int DeviceMemoryPool::bucket_for(std::size_t bytes) noexcept {
    if (bytes <= bucket_size(kMinShift)) return kMinShift;
    // Smallest k such that (1 << k) >= bytes.
    int k = kMinShift;
    while (k < kMaxShift && (std::size_t(1) << k) < bytes) ++k;
    return k;
}

// -----------------------------------------------------------------------------
DeviceMemoryPool::DeviceMemoryPool(std::size_t maxBytes)
    : impl_(std::make_unique<Impl>()) {
    impl_->maxBytes = maxBytes;
}

DeviceMemoryPool::~DeviceMemoryPool() {
    purge();  // releases cached slabs; in-use slabs are leaked by intent
              // (their owners are still alive on a singleton teardown path,
              // matching Logger / CudaContext lifetime semantics).
}

DeviceMemoryPool& DeviceMemoryPool::instance() {
    static DeviceMemoryPool pool;
    return pool;
}

// -----------------------------------------------------------------------------
void* DeviceMemoryPool::allocate(std::size_t bytes) {
    if (bytes == 0) bytes = 1;
    const int k    = bucket_for(bytes);
    const std::size_t size = bucket_size(k);
    if (k > kMaxShift) throw std::bad_alloc{};

    std::lock_guard<std::mutex> lk(impl_->mtx);
    impl_->stats.allocation_calls++;

    auto& list = impl_->freelist[k];
    if (!list.empty()) {
        void* p = list.back();
        list.pop_back();
        impl_->stats.freelist_hits++;
        impl_->stats.total_in_freelist_bytes -= size;
        impl_->stats.total_in_use_bytes      += size;
        return p;
    }
    impl_->stats.freelist_misses++;

    // Try fresh allocation respecting budget.
    auto fits_budget = [&]() {
        return impl_->maxBytes == 0 ||
               impl_->stats.total_reserved_bytes + size <= impl_->maxBytes;
    };
    if (!fits_budget()) {
        // Evict larger buckets first (LIFO), smaller buckets last.
        for (int kk = kMaxShift; kk >= kMinShift && !fits_budget(); --kk) {
            auto& other = impl_->freelist[kk];
            while (!other.empty() && !fits_budget()) {
                void* p = other.back(); other.pop_back();
                raw_free(p);
                const std::size_t s = bucket_size(kk);
                impl_->stats.total_in_freelist_bytes -= s;
                impl_->stats.total_reserved_bytes    -= s;
                impl_->stats.slab_count              -= 1;
            }
        }
    }
    void* p = raw_alloc(size);
    if (!p) throw std::bad_alloc{};
    impl_->stats.total_reserved_bytes += size;
    impl_->stats.total_in_use_bytes   += size;
    impl_->stats.slab_count           += 1;
    return p;
}

// -----------------------------------------------------------------------------
void DeviceMemoryPool::release(void* ptr, std::size_t bytes) noexcept {
    if (!ptr) return;
    if (bytes == 0) bytes = 1;
    const int k    = bucket_for(bytes);
    const std::size_t size = bucket_size(k);

    std::lock_guard<std::mutex> lk(impl_->mtx);
    impl_->freelist[k].push_back(ptr);
    impl_->stats.total_in_use_bytes      -= size;
    impl_->stats.total_in_freelist_bytes += size;
}

// -----------------------------------------------------------------------------
void DeviceMemoryPool::purge() noexcept {
    std::lock_guard<std::mutex> lk(impl_->mtx);
    for (int k = kMinShift; k <= kMaxShift; ++k) {
        auto& list = impl_->freelist[k];
        const std::size_t s = bucket_size(k);
        for (void* p : list) {
            raw_free(p);
            impl_->stats.total_reserved_bytes    -= s;
            impl_->stats.total_in_freelist_bytes -= s;
            impl_->stats.slab_count              -= 1;
        }
        list.clear();
    }
}

// -----------------------------------------------------------------------------
DeviceMemoryStats DeviceMemoryPool::stats() const {
    std::lock_guard<std::mutex> lk(impl_->mtx);
    return impl_->stats;
}

}  // namespace simall::gpu
