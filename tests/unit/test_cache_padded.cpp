// =============================================================================
// SimAll Beta - Utilities Unit Tests
// File   : tests/unit/test_cache_padded.cpp
// Phase  : 24 Pass 2
//
// Coverage:
//   - sizeof(CachePadded<T>) is a positive multiple of the cache-line size
//     for several T's (int, double, std::array<int,5>, a 70-byte struct)
//   - alignof(CachePadded<T>) == kCacheLineBytes
//   - adjacent elements in std::vector<CachePadded<int>> are at least one
//     cache line apart (false-sharing layout invariant)
//   - value()/operator*()/operator->() expose the underlying T
//   - constructor-forwarding works for non-default-constructible types
// =============================================================================
#include <catch2/catch_test_macros.hpp>

#include "utilities/CachePadded.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using simall::util::CachePadded;
using simall::util::kCacheLineBytes;

namespace {
struct Big { char b[70]; };       // sizeof > one cache line
struct NoDefault {
    int x;
    explicit NoDefault(int v) : x(v) {}
};
}  // namespace

// =============================================================================
TEST_CASE("CachePadded<T> size is a positive multiple of the cache-line",
          "[utilities][memory][cache-padded]") {
    REQUIRE(sizeof(CachePadded<int>)                 % kCacheLineBytes == 0);
    REQUIRE(sizeof(CachePadded<double>)              % kCacheLineBytes == 0);
    REQUIRE(sizeof(CachePadded<std::array<int, 5>>)  % kCacheLineBytes == 0);
    REQUIRE(sizeof(CachePadded<Big>)                 % kCacheLineBytes == 0);
    REQUIRE(sizeof(CachePadded<Big>)                 >= sizeof(Big));
}

// =============================================================================
TEST_CASE("CachePadded<T> alignment matches the cache-line",
          "[utilities][memory][cache-padded]") {
    REQUIRE(alignof(CachePadded<int>)    == kCacheLineBytes);
    REQUIRE(alignof(CachePadded<double>) == kCacheLineBytes);
    REQUIRE(alignof(CachePadded<Big>)    == kCacheLineBytes);
}

// =============================================================================
TEST_CASE("Adjacent CachePadded elements are at least one cache-line apart",
          "[utilities][memory][cache-padded]") {
    std::vector<CachePadded<std::int64_t>> v(4);
    for (std::size_t i = 1; i < v.size(); ++i) {
        auto p0 = reinterpret_cast<std::uintptr_t>(&v[i - 1]);
        auto p1 = reinterpret_cast<std::uintptr_t>(&v[i]);
        REQUIRE((p1 - p0) >= kCacheLineBytes);
    }
}

// =============================================================================
TEST_CASE("CachePadded<T> accessors expose the wrapped value",
          "[utilities][memory][cache-padded]") {
    CachePadded<int> ci(42);
    REQUIRE(ci.value() == 42);
    REQUIRE(*ci        == 42);
    *ci = 7;
    REQUIRE(ci.value() == 7);

    CachePadded<std::string> cs("hello");
    REQUIRE(cs->size() == 5);
    REQUIRE(*cs        == "hello");
}

// =============================================================================
TEST_CASE("CachePadded<T> forwards constructor args to non-default types",
          "[utilities][memory][cache-padded]") {
    CachePadded<NoDefault> cn(99);
    REQUIRE(cn.value().x == 99);
}
