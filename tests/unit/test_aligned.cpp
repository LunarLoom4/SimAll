#include "utilities/AlignedAllocator.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstdint>

TEST_CASE("aligned_vector storage is 64-byte aligned", "[utilities]")
{
    simall::util::aligned_vector<double> v(1024, 1.5);
    auto addr = reinterpret_cast<std::uintptr_t>(v.data());
    REQUIRE((addr % simall::util::kSimdAlignBytes) == 0);
    REQUIRE(v.size() == 1024);
    REQUIRE(v[0] == 1.5);
}
