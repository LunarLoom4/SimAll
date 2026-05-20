// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/core/test_thread_pool.cpp
// =============================================================================
#include "core/ThreadPool.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <vector>

using namespace simall::core;

TEST_CASE("ThreadPool runs many jobs and waits idle", "[core][threadpool]")
{
    ThreadPool pool(ThreadPoolConfig{4, false, -1, "test"});
    REQUIRE(pool.workerCount() == 4);

    constexpr int N = 1000;
    std::atomic<int> counter{0};
    for (int i = 0; i < N; ++i) {
        pool.submit([&] { counter.fetch_add(1, std::memory_order_relaxed); });
    }
    pool.waitIdle();
    REQUIRE(counter.load() == N);
}

TEST_CASE("ThreadPool returns values via futures", "[core][threadpool]")
{
    ThreadPool pool(ThreadPoolConfig{2});
    std::vector<std::future<int>> futs;
    futs.reserve(64);
    for (int i = 0; i < 64; ++i) {
        futs.push_back(pool.submit([i] { return i * i; }));
    }
    int sum = 0;
    for (int i = 0; i < 64; ++i)
        sum += futs[i].get();
    int expected = 0;
    for (int i = 0; i < 64; ++i)
        expected += i * i;
    REQUIRE(sum == expected);
}

TEST_CASE("ThreadPool global() is reusable across calls", "[core][threadpool]")
{
    auto& a = ThreadPool::global();
    auto& b = ThreadPool::global();
    REQUIRE(&a == &b);
    auto fut = a.submit([] { return 7; });
    REQUIRE(fut.get() == 7);
}
