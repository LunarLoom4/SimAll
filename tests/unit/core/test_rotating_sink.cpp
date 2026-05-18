// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/core/test_rotating_sink.cpp
// =============================================================================
#include "core/RotatingFileSink.hpp"

#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

using namespace simall::core;

TEST_CASE("RotatingFileSink rotates when size exceeded", "[core][log]") {
    // Use a unique base path so concurrent or repeated test runs don't collide.
    auto stamp = std::to_string(
        std::chrono::steady_clock::now().time_since_epoch().count());
    auto base  = std::filesystem::temp_directory_path() /
                 ("simall_rot_test_" + stamp + ".log");
    std::filesystem::remove(base);
    for (unsigned i = 1; i <= 6; ++i) {
        auto p = base; p += "." + std::to_string(i);
        std::filesystem::remove(p);
    }

    {
        RotatingFileSink sink(base, /*maxBytes*/ 32, /*maxFiles*/ 3);
        std::string line(20, 'x');
        line.push_back('\n');
        for (int i = 0; i < 20; ++i) sink.write(line);
    }

    auto r1 = base; r1 += ".1";
    auto r2 = base; r2 += ".2";
    REQUIRE(std::filesystem::exists(base));
    REQUIRE(std::filesystem::exists(r1));
    REQUIRE(std::filesystem::exists(r2));

    // The oldest beyond maxFiles should not exist.
    auto r4 = base; r4 += ".4";
    REQUIRE_FALSE(std::filesystem::exists(r4));

    // Cleanup.
    std::filesystem::remove(base);
    for (unsigned i = 1; i <= 6; ++i) {
        auto p = base; p += "." + std::to_string(i);
        std::filesystem::remove(p);
    }
}
