// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/ThreadPool.hpp
// Phase  : 1.3 (APPLICATION CORE → Threading)
//
// Work-stealing thread pool. Each worker owns a deque of jobs; idle workers
// steal from the back of victims' deques. Optional CPU affinity pins each
// worker to a single logical processor so cache-warm jobs stay local — this
// approximates NUMA awareness without requiring libnuma. On Windows the
// affinity mask uses SetThreadAffinityMask; on Linux, pthread_setaffinity_np.
// On platforms without affinity APIs (macOS), the call is a no-op and the
// pool behaves as a plain work-stealing pool.
//
// Public API:
//   ThreadPool pool;                       // hardware_concurrency() workers
//   ThreadPool pool(n, {pinAffinity});
//   auto fut = pool.submit([]{ return 42; });
//   fut.wait();
//
// Submit places the job onto a randomly-chosen worker queue; idle workers
// rebalance by stealing FIFO-style from the front of victims' deques.
// =============================================================================
#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <random>
#include <thread>
#include <type_traits>
#include <vector>

namespace simall::core {

struct ThreadPoolConfig {
    std::size_t numThreads     = 0;        // 0 → hardware_concurrency()
    bool        pinAffinity    = false;    // pin worker i to logical CPU i
    int         numaNodeHint   = -1;       // -1 → no preference (informational)
    std::string threadNamePrefix = "simall-worker";
};

class ThreadPool {
public:
    explicit ThreadPool(ThreadPoolConfig cfg = {});
    ~ThreadPool();

    ThreadPool(const ThreadPool&)            = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    /// Submit any callable; returns std::future for its result.
    template <class F, class... Args>
    auto submit(F&& f, Args&&... args)
        -> std::future<std::invoke_result_t<F, Args...>>
    {
        using R = std::invoke_result_t<F, Args...>;
        auto task = std::make_shared<std::packaged_task<R()>>(
            [fn = std::forward<F>(f), tup = std::make_tuple(std::forward<Args>(args)...)]() mutable {
                return std::apply(std::move(fn), std::move(tup));
            });
        std::future<R> fut = task->get_future();
        pushJob([task]() { (*task)(); });
        return fut;
    }

    /// Block until all currently-queued jobs are complete.
    void waitIdle();

    std::size_t workerCount() const noexcept { return workers_.size(); }
    std::size_t pendingJobs() const noexcept { return pending_.load(std::memory_order_relaxed); }

    /// Process-wide default pool (lazy). Used by TaskGraph and the GUI.
    static ThreadPool& global();

private:
    using Job = std::function<void()>;

    struct Worker;

    void pushJob(Job j);
    bool tryPopLocal(std::size_t self, Job& out);
    bool trySteal(std::size_t self, Job& out);
    void workerLoop(std::size_t index);
    void setNameAndAffinity(std::size_t index);

    ThreadPoolConfig             cfg_;
    std::vector<std::thread>     workers_;
    std::vector<std::unique_ptr<Worker>> queues_;
    std::atomic<bool>            stopping_{false};
    std::atomic<std::size_t>     pending_{0};
    std::condition_variable      idle_cv_;
    std::mutex                   idle_mtx_;
};

}  // namespace simall::core
