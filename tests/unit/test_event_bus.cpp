#include "core/EventBus.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("EventBus delivers typed events", "[core]")
{
    using namespace simall::core;
    struct Ping
    {
        int v;
    };
    int received = 0;
    auto id = EventBus::instance().subscribe<Ping>([&](const Ping& p) { received = p.v; });
    EventBus::instance().publish(Ping{42});
    REQUIRE(received == 42);
    EventBus::instance().unsubscribe<Ping>(id);
    EventBus::instance().publish(Ping{100});
    REQUIRE(received == 42);
}
