// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/core/test_task_graph.cpp
// =============================================================================
#include "core/TaskGraph.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <stdexcept>
#include <vector>

using namespace simall::core;

TEST_CASE("TaskGraph respects dependencies", "[core][taskgraph]")
{
    ThreadPool pool(ThreadPoolConfig{4});
    TaskGraph g(pool);

    std::atomic<int> order{0};
    int a_seq = -1, b_seq = -1, c_seq = -1;

    auto a = g.add([&] { a_seq = order.fetch_add(1); });
    auto b = g.add([&] { b_seq = order.fetch_add(1); });
    auto c = g.add([&] { c_seq = order.fetch_add(1); });
    g.dependOn(c, a);
    g.dependOn(c, b);

    g.run();
    g.waitAll();

    REQUIRE(c_seq > a_seq);
    REQUIRE(c_seq > b_seq);
}

TEST_CASE("TaskGraph detects cycles", "[core][taskgraph]")
{
    TaskGraph g;
    auto a = g.add([] {});
    auto b = g.add([] {});
    g.dependOn(a, b);
    g.dependOn(b, a);
    REQUIRE_THROWS_AS(g.run(), std::runtime_error);
}

TEST_CASE("TaskGraph then() chains tasks", "[core][taskgraph]")
{
    TaskGraph g;
    std::vector<int> log;
    std::mutex m;
    auto push = [&](int v) {
        std::lock_guard lk(m);
        log.push_back(v);
    };

    auto t1 = g.add([&] { push(1); });
    auto t2 = t1.then([&] { push(2); });
    (void) t2.then([&] { push(3); });

    g.run();
    g.waitAll();
    REQUIRE(log == std::vector<int>{1, 2, 3});
}

TEST_CASE("TaskGraph propagates exceptions", "[core][taskgraph]")
{
    TaskGraph g;
    g.add([] { throw std::runtime_error("boom"); });
    g.run();
    REQUIRE_THROWS_AS(g.waitAll(), std::runtime_error);
}
