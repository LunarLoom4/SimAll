// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/StreamScheduler.hpp
// Phase  : 18.3 — CUDA stream + event scheduler (W13).
//
// PURPOSE
// -------
// A pool of cudaStream_t objects round-robined across linear-algebra kernel
// launches so that independent operations (e.g. y = A*x and z = B*x with
// disjoint inputs) can overlap on the device.  Backed by an event-pool for
// inter-stream dependency synchronisation:
//
//     auto sA = scheduler.acquire();
//     spmv_csr(... sA ...);                 // launch on sA
//     auto eA = scheduler.record_event(sA);
//     auto sB = scheduler.acquire();
//     scheduler.wait_event(sB, eA);          // sB stalls until eA completes
//     axpy(..., sB);
//     scheduler.release_event(eA);
//     scheduler.release(sA);
//     scheduler.release(sB);
//
// CORRECTNESS
// -----------
//   * Streams are opaque (void*); .cpp casts to cudaStream_t.
//   * Acquire / release is thread-safe via a mutex (kernel launches from
//     multiple host threads are common in TBB-driven workloads).
//   * On non-CUDA builds every operation is a no-op and stream/event
//     handles are non-null sentinels so client code branches uniformly.
// =============================================================================
#pragma once

#include <cstddef>
#include <memory>

namespace simall::gpu
{

class StreamScheduler
{
public:
    /// ``n`` streams created up-front.  Common values: 2–4 for SIMPLE/PISO,
    /// 8+ for SpMV-heavy GMRES.
    explicit StreamScheduler(int n = 4);
    ~StreamScheduler();

    /// Process-wide singleton (4 streams by default).
    static StreamScheduler& instance();

    /// Round-robin pick.  Returns nullptr-sentinel on non-CUDA builds (still
    /// safe to pass to every kernel; they ignore it).
    void* acquire();

    /// Return a stream to the pool.  Optional — pool tolerates leakage,
    /// just degrades round-robin coverage.
    void release(void* stream) noexcept;

    /// cudaStreamSynchronize on a single stream.  No-op on non-CUDA.
    void synchronize(void* stream) noexcept;

    /// Block the host until every pool stream has drained.
    void synchronize_all() noexcept;

    /// Record an event on ``stream``.  Caller owns the event until it
    /// releases it via release_event.
    void* record_event(void* stream);

    /// Make ``waitingStream`` stall until ``event`` is signalled.  Does NOT
    /// block the host.
    void wait_event(void* waitingStream, void* event) noexcept;

    /// Return an event to the pool.
    void release_event(void* event) noexcept;

    /// Number of streams in the pool.
    int stream_count() const noexcept;

    StreamScheduler(const StreamScheduler&) = delete;
    StreamScheduler& operator=(const StreamScheduler&) = delete;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace simall::gpu
