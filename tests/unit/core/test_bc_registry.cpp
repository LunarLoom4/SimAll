// =============================================================================
// SimAll Beta — Core Unit Tests
// File   : tests/unit/core/test_bc_registry.cpp
// Phase  : 22 Pass 2
//
// Validates the central BoundaryConditionRegistry (register / lookup /
// iterate / remove / clear / typed accessor / singleton).
// =============================================================================
#include "core/BoundaryConditionRegistry.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using namespace simall::core;

namespace {

struct FakeBc {
    int    kind  = 0;
    double value = 0.0;
};

BoundaryConditionEntry make_entry(int kindId, const std::string& kindName,
                                  BcZoneId zone, const std::string& variable,
                                  const std::string& subsystem,
                                  std::shared_ptr<FakeBc> bc = {}) {
    BoundaryConditionEntry e;
    e.kindId    = kindId;
    e.kindName  = kindName;
    e.zone      = zone;
    e.variable  = variable;
    e.subsystem = subsystem;
    e.opaque    = std::move(bc);
    return e;
}

}  // namespace

TEST_CASE("BoundaryConditionRegistry add/find returns inserted entry",
          "[core][bc-registry]") {
    BoundaryConditionRegistry reg;
    REQUIRE(reg.empty());

    auto bc = std::make_shared<FakeBc>(FakeBc{7, 3.14});
    const auto h = reg.add(make_entry(7, "Wall", 42, "U", "solver", bc));

    REQUIRE(h != BoundaryConditionRegistry::kInvalid);
    REQUIRE(reg.size() == 1);
    REQUIRE_FALSE(reg.empty());

    const auto* e = reg.find(h);
    REQUIRE(e != nullptr);
    CHECK(e->kindId    == 7);
    CHECK(e->kindName  == "Wall");
    CHECK(e->zone      == 42u);
    CHECK(e->variable  == "U");
    CHECK(e->subsystem == "solver");

    auto* typed = reg.as<FakeBc>(h);
    REQUIRE(typed != nullptr);
    CHECK(typed->kind  == 7);
    CHECK(typed->value == 3.14);
}

TEST_CASE("BoundaryConditionRegistry find on unknown handle returns null",
          "[core][bc-registry]") {
    BoundaryConditionRegistry reg;
    CHECK(reg.find(123) == nullptr);
    CHECK(reg.as<FakeBc>(99) == nullptr);
    CHECK(reg.first(0, "U") == BoundaryConditionRegistry::kInvalid);
}

TEST_CASE("BoundaryConditionRegistry lookup helpers", "[core][bc-registry]") {
    BoundaryConditionRegistry reg;
    const auto hA = reg.add(make_entry(1, "Wall",          1, "U", "solver"));
    const auto hB = reg.add(make_entry(2, "Inlet",         1, "p", "solver"));
    const auto hC = reg.add(make_entry(3, "MagneticInsul", 2, "B", "em"));
    const auto hD = reg.add(make_entry(4, "Fixed",         3, "u", "structural"));
    (void)hA; (void)hB; (void)hC; (void)hD;

    auto zone1 = reg.by_zone(1);
    CHECK(zone1.size() == 2);
    CHECK(zone1[0] == hA);     // insertion order preserved
    CHECK(zone1[1] == hB);

    auto varU = reg.by_variable("U");
    CHECK(varU.size() == 1);
    CHECK(varU[0] == hA);

    auto subEm = reg.by_subsystem("em");
    CHECK(subEm.size() == 1);
    CHECK(subEm[0] == hC);

    CHECK(reg.first(1, "U") == hA);
    CHECK(reg.first(1, "p") == hB);
    CHECK(reg.first(2, "B") == hC);
    CHECK(reg.first(9, "U") == BoundaryConditionRegistry::kInvalid);
}

TEST_CASE("BoundaryConditionRegistry for_each visits in insertion order",
          "[core][bc-registry]") {
    BoundaryConditionRegistry reg;
    const auto h1 = reg.add(make_entry(1, "A", 1, "U", "solver"));
    const auto h2 = reg.add(make_entry(2, "B", 2, "p", "solver"));
    const auto h3 = reg.add(make_entry(3, "C", 3, "T", "thermal"));

    std::vector<BoundaryConditionRegistry::Handle> visited;
    std::vector<std::string>                       names;
    reg.for_each([&](auto h, const BoundaryConditionEntry& e) {
        visited.push_back(h);
        names.push_back(e.kindName);
    });

    REQUIRE(visited.size() == 3);
    CHECK(visited[0] == h1);
    CHECK(visited[1] == h2);
    CHECK(visited[2] == h3);
    CHECK(names == std::vector<std::string>{"A", "B", "C"});
}

TEST_CASE("BoundaryConditionRegistry remove drops the entry and re-balances "
          "insertion order", "[core][bc-registry]") {
    BoundaryConditionRegistry reg;
    const auto h1 = reg.add(make_entry(1, "A", 1, "U", "solver"));
    const auto h2 = reg.add(make_entry(2, "B", 2, "p", "solver"));
    const auto h3 = reg.add(make_entry(3, "C", 3, "T", "thermal"));

    CHECK(reg.remove(h2));
    CHECK_FALSE(reg.remove(h2));         // double-remove is a no-op
    CHECK(reg.size() == 2);
    CHECK(reg.find(h2) == nullptr);

    std::vector<BoundaryConditionRegistry::Handle> visited;
    reg.for_each([&](auto h, const BoundaryConditionEntry&) {
        visited.push_back(h);
    });
    REQUIRE(visited.size() == 2);
    CHECK(visited[0] == h1);
    CHECK(visited[1] == h3);
}

TEST_CASE("BoundaryConditionRegistry clear resets state", "[core][bc-registry]") {
    BoundaryConditionRegistry reg;
    reg.add(make_entry(1, "A", 1, "U", "solver"));
    reg.add(make_entry(2, "B", 2, "p", "solver"));
    REQUIRE(reg.size() == 2);

    reg.clear();
    CHECK(reg.empty());
    CHECK(reg.size() == 0);
    CHECK(reg.by_zone(1).empty());
    CHECK(reg.by_variable("U").empty());
}

TEST_CASE("BoundaryConditionRegistry singleton is process-wide",
          "[core][bc-registry]") {
    auto& a = BoundaryConditionRegistry::instance();
    auto& b = BoundaryConditionRegistry::instance();
    CHECK(&a == &b);
    a.clear();
}

TEST_CASE("BoundaryConditionRegistry concurrent add is thread-safe",
          "[core][bc-registry][concurrency]") {
    BoundaryConditionRegistry reg;
    constexpr int kThreads = 8;
    constexpr int kPer     = 200;
    std::vector<std::thread> ts;
    std::atomic<int> okCount{0};
    for (int t = 0; t < kThreads; ++t) {
        ts.emplace_back([&, t] {
            for (int i = 0; i < kPer; ++i) {
                auto h = reg.add(make_entry(t, "Wall",
                                            static_cast<BcZoneId>(t * 100 + i),
                                            "U", "solver"));
                if (h != BoundaryConditionRegistry::kInvalid) ++okCount;
            }
        });
    }
    for (auto& th : ts) th.join();

    CHECK(okCount.load() == kThreads * kPer);
    CHECK(reg.size()      == static_cast<std::size_t>(kThreads * kPer));
}
