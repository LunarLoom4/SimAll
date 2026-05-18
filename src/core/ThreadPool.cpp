// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/ThreadPool.cpp
// =============================================================================
#include "core/ThreadPool.hpp"

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <Windows.h>
#elif defined(__linux__)
#  include <pthread.h>
#  include <sched.h>
#  include <sys/prctl.h>
#elif defined(__APPLE__)
#  include <pthread.h>
#endif

namespace simall::core {

struct ThreadPool::Worker {
    std::mutex      m;
    std::deque<Job> q;          // back = LIFO local, front = victim steal
    std::condition_variable cv;
};

ThreadPool& ThreadPool::global() {
    static ThreadPool inst{ThreadPoolConfig{}};
    return inst;
}

ThreadPool::ThreadPool(ThreadPoolConfig cfg) : cfg_(std::move(cfg)) {
    std::size_t n = cfg_.numThreads;
    if (n == 0) {
        n = std::thread::hardware_concurrency();
        if (n == 0) n = 4;
        if (n > 2) n -= 1;  // leave one core for the GUI thread
    }
    queues_.reserve(n);
    workers_.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        queues_.emplace_back(std::make_unique<Worker>());
    }
    for (std::size_t i = 0; i < n; ++i) {
        workers_.emplace_back([this, i] { workerLoop(i); });
    }
}

ThreadPool::~ThreadPool() {
    stopping_.store(true, std::memory_order_release);
    for (auto& q : queues_) q->cv.notify_all();
    for (auto& w : workers_) {
        if (w.joinable()) w.join();
    }
}

void ThreadPool::pushJob(Job j) {
    // Randomized initial placement provides natural fan-out across worker
    // queues; idle workers later balance load via FIFO stealing.
    static thread_local std::mt19937 rng{std::random_device{}()};
    std::size_t target = rng() % queues_.size();
    pending_.fetch_add(1, std::memory_order_release);
    {
        auto& q = *queues_[target];
        std::lock_guard lk(q.m);
        q.q.push_back(std::move(j));
        q.cv.notify_one();
    }
}

bool ThreadPool::tryPopLocal(std::size_t self, Job& out) {
    auto& q = *queues_[self];
    std::lock_guard lk(q.m);
    if (q.q.empty()) return false;
    out = std::move(q.q.back());      // LIFO local for cache locality
    q.q.pop_back();
    return true;
}

bool ThreadPool::trySteal(std::size_t self, Job& out) {
    std::size_t n = queues_.size();
    static thread_local std::mt19937 rng{std::random_device{}() ^ static_cast<unsigned>(self)};
    std::uniform_int_distribution<std::size_t> dist(0, n - 1);
    for (std::size_t attempt = 0; attempt < n; ++attempt) {
        std::size_t v = dist(rng);
        if (v == self) continue;
        auto& q = *queues_[v];
        std::unique_lock lk(q.m, std::try_to_lock);
        if (!lk.owns_lock() || q.q.empty()) continue;
        out = std::move(q.q.front()); // FIFO steal from front
        q.q.pop_front();
        return true;
    }
    return false;
}

void ThreadPool::workerLoop(std::size_t index) {
    setNameAndAffinity(index);
    auto& my = *queues_[index];
    while (!stopping_.load(std::memory_order_acquire)) {
        Job j;
        if (tryPopLocal(index, j) || trySteal(index, j)) {
            try { j(); }
            catch (...) { /* swallow — TaskGraph captures via promises */ }
            if (pending_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                std::lock_guard lk(idle_mtx_);
                idle_cv_.notify_all();
            }
            continue;
        }
        std::unique_lock lk(my.m);
        my.cv.wait_for(lk, std::chrono::milliseconds(5),
                       [&]{ return !my.q.empty() || stopping_.load(); });
    }
}

void ThreadPool::waitIdle() {
    std::unique_lock lk(idle_mtx_);
    idle_cv_.wait(lk, [&]{ return pending_.load() == 0; });
}

void ThreadPool::setNameAndAffinity(std::size_t index) {
    std::string name = cfg_.threadNamePrefix + "-" + std::to_string(index);
#if defined(_WIN32)
    if (cfg_.pinAffinity) {
        DWORD_PTR mask = DWORD_PTR{1} << (index % (sizeof(DWORD_PTR) * 8));
        SetThreadAffinityMask(GetCurrentThread(), mask);
    }
    // Thread naming via SetThreadDescription (Win10+)
    wchar_t wname[64]{};
    for (std::size_t i = 0; i < name.size() && i < 63; ++i) {
        wname[i] = static_cast<wchar_t>(name[i]);
    }
    using SetNameFn = HRESULT(WINAPI*)(HANDLE, PCWSTR);
    static auto fn = reinterpret_cast<SetNameFn>(
        GetProcAddress(GetModuleHandleW(L"kernel32.dll"), "SetThreadDescription"));
    if (fn) fn(GetCurrentThread(), wname);
#elif defined(__linux__)
    if (cfg_.pinAffinity) {
        cpu_set_t set;
        CPU_ZERO(&set);
        unsigned cores = std::thread::hardware_concurrency();
        if (cores == 0) cores = 1;
        CPU_SET(static_cast<int>(index % cores), &set);
        pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
    }
    // PR_SET_NAME accepts max 15 chars + null.
    std::string short_name = name.substr(0, 15);
    prctl(PR_SET_NAME, short_name.c_str(), 0, 0, 0);
#elif defined(__APPLE__)
    pthread_setname_np(name.c_str());
#else
    (void)index;
#endif
}

}  // namespace simall::core
