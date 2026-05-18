// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/StreamScheduler.cpp
// Phase  : 18.3 — CUDA stream + event scheduler (W13).
// =============================================================================
#include "gpu/StreamScheduler.hpp"
#include "gpu/CudaContext.hpp"
#include "core/Logger.hpp"

#include <atomic>
#include <cstdint>
#include <mutex>
#include <vector>

#ifdef SIMALL_HAVE_CUDA
    #include <cuda_runtime.h>
#endif

namespace simall::gpu {

namespace {
/// Non-null sentinel returned on non-CUDA builds so that callers never see
/// a nullptr and accidentally fall through error checks.
constexpr void* kSerialSentinel = reinterpret_cast<void*>(static_cast<std::uintptr_t>(1));
}  // namespace

struct StreamScheduler::Impl {
    mutable std::mutex     mtx;
    std::vector<void*>     streams;       ///< owned, never null
    std::vector<bool>      streamInUse;
    std::vector<void*>     eventCache;    ///< free events (LIFO)
    std::atomic<std::size_t> rrIndex{0};  ///< round-robin cursor

#ifdef SIMALL_HAVE_CUDA
    bool cuda = false;
#endif
};

// -----------------------------------------------------------------------------
StreamScheduler::StreamScheduler(int n) : impl_(std::make_unique<Impl>()) {
    if (n < 1) n = 1;
    impl_->streams.resize(static_cast<std::size_t>(n), kSerialSentinel);
    impl_->streamInUse.assign(static_cast<std::size_t>(n), false);

#ifdef SIMALL_HAVE_CUDA
    impl_->cuda = is_cuda_available();
    if (impl_->cuda) {
        for (int i = 0; i < n; ++i) {
            cudaStream_t s = nullptr;
            if (cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking) == cudaSuccess) {
                impl_->streams[i] = static_cast<void*>(s);
            } else {
                impl_->streams[i] = kSerialSentinel;  // fall through
            }
        }
    }
#endif
    SIMALL_LOG_INFO("Gpu", "StreamScheduler initialised with ", n, " streams");
}

StreamScheduler::~StreamScheduler() {
#ifdef SIMALL_HAVE_CUDA
    if (impl_->cuda) {
        for (void* s : impl_->streams)
            if (s && s != kSerialSentinel)
                cudaStreamDestroy(static_cast<cudaStream_t>(s));
        for (void* e : impl_->eventCache)
            if (e && e != kSerialSentinel)
                cudaEventDestroy(static_cast<cudaEvent_t>(e));
    }
#endif
}

StreamScheduler& StreamScheduler::instance() {
    static StreamScheduler sched(4);
    return sched;
}

// -----------------------------------------------------------------------------
void* StreamScheduler::acquire() {
    std::lock_guard<std::mutex> lk(impl_->mtx);
    const std::size_t n = impl_->streams.size();
    // Single pass starting at the round-robin cursor; mark first available.
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t k = (impl_->rrIndex.fetch_add(1) ) % n;
        if (!impl_->streamInUse[k]) {
            impl_->streamInUse[k] = true;
            return impl_->streams[k];
        }
    }
    // Pool exhausted — return the next round-robin stream and leave the
    // in-use flag set; release tolerates double-releases.
    const std::size_t k = impl_->rrIndex.fetch_add(1) % n;
    return impl_->streams[k];
}

void StreamScheduler::release(void* stream) noexcept {
    if (!stream) return;
    std::lock_guard<std::mutex> lk(impl_->mtx);
    for (std::size_t i = 0; i < impl_->streams.size(); ++i) {
        if (impl_->streams[i] == stream && impl_->streamInUse[i]) {
            impl_->streamInUse[i] = false;
            return;
        }
    }
}

void StreamScheduler::synchronize(void* stream) noexcept {
#ifdef SIMALL_HAVE_CUDA
    if (impl_->cuda && stream && stream != kSerialSentinel)
        cudaStreamSynchronize(static_cast<cudaStream_t>(stream));
#else
    (void)stream;
#endif
}

void StreamScheduler::synchronize_all() noexcept {
#ifdef SIMALL_HAVE_CUDA
    if (!impl_->cuda) return;
    for (void* s : impl_->streams)
        if (s && s != kSerialSentinel)
            cudaStreamSynchronize(static_cast<cudaStream_t>(s));
#endif
}

void* StreamScheduler::record_event(void* stream) {
#ifdef SIMALL_HAVE_CUDA
    if (!impl_->cuda) return kSerialSentinel;
    std::lock_guard<std::mutex> lk(impl_->mtx);
    cudaEvent_t e = nullptr;
    if (!impl_->eventCache.empty()) {
        e = static_cast<cudaEvent_t>(impl_->eventCache.back());
        impl_->eventCache.pop_back();
    } else {
        if (cudaEventCreateWithFlags(&e, cudaEventDisableTiming) != cudaSuccess)
            return kSerialSentinel;
    }
    cudaEventRecord(e, stream && stream != kSerialSentinel
                       ? static_cast<cudaStream_t>(stream) : 0);
    return static_cast<void*>(e);
#else
    (void)stream;
    return kSerialSentinel;
#endif
}

void StreamScheduler::wait_event(void* waitingStream, void* event) noexcept {
#ifdef SIMALL_HAVE_CUDA
    if (!impl_->cuda) return;
    if (!event || event == kSerialSentinel) return;
    cudaStreamWaitEvent(
        waitingStream && waitingStream != kSerialSentinel
            ? static_cast<cudaStream_t>(waitingStream) : 0,
        static_cast<cudaEvent_t>(event), 0);
#else
    (void)waitingStream; (void)event;
#endif
}

void StreamScheduler::release_event(void* event) noexcept {
    if (!event || event == kSerialSentinel) return;
#ifdef SIMALL_HAVE_CUDA
    if (!impl_->cuda) return;
    std::lock_guard<std::mutex> lk(impl_->mtx);
    impl_->eventCache.push_back(event);
#endif
}

int StreamScheduler::stream_count() const noexcept {
    return static_cast<int>(impl_->streams.size());
}

}  // namespace simall::gpu
