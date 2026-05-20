// =============================================================================
// SimAll Beta - Utilities Unit Tests
// File   : tests/unit/test_memory_pool.cpp
// Phase  : 24 Pass 1
//
// Coverage:
//   - MemoryArena returns suitably aligned blocks and reset() reuses them
//   - MemoryArena grows beyond its initial block
//   - ObjectPool<T> allocate/deallocate round-trip + freelist reuse
//   - ObjectPool<T> grows across slabs and tracks live_count/slab_count
//   - ObjectPool<T>::construct/destroy invoke ctor/dtor exactly once
// =============================================================================
#include <catch2/catch_test_macros.hpp>

#include "utilities/MemoryPool.hpp"

#include <cstdint>
#include <vector>

using simall::util::MemoryArena;
using simall::util::ObjectPool;

// =============================================================================
TEST_CASE("MemoryArena returns aligned allocations and reset reuses memory",
          "[utilities][memory][arena]") {
    MemoryArena a(4096);
    void* p1 = a.allocate(64, 64);
    void* p2 = a.allocate(64, 64);
    REQUIRE(reinterpret_cast<std::uintptr_t>(p1) % 64 == 0);
    REQUIRE(reinterpret_cast<std::uintptr_t>(p2) % 64 == 0);
    REQUIRE(p1 != p2);

    a.reset();
    void* p3 = a.allocate(64, 64);
    REQUIRE(p3 == p1);                   // first slot of first block reused
}

// =============================================================================
TEST_CASE("MemoryArena grows when a request exceeds the current block",
          "[utilities][memory][arena]") {
    MemoryArena a(256);
    void* small = a.allocate(64, 8);
    void* big   = a.allocate(4096, 16); // forces a new block
    REQUIRE(small != nullptr);
    REQUIRE(big   != nullptr);
    REQUIRE(reinterpret_cast<std::uintptr_t>(big) % 16 == 0);
}

// =============================================================================
TEST_CASE("ObjectPool<T> allocate/deallocate reuses slots from the freelist",
          "[utilities][memory][objectpool]") {
    ObjectPool<std::uint64_t> pool(/*slot_capacity=*/8);
    std::uint64_t* a = pool.allocate();
    std::uint64_t* b = pool.allocate();
    REQUIRE(a != b);
    REQUIRE(pool.live_count() == 2);

    pool.deallocate(a);
    REQUIRE(pool.live_count() == 1);
    std::uint64_t* c = pool.allocate();
    REQUIRE(c == a);                     // freelist LIFO reuse
    REQUIRE(pool.live_count() == 2);
}

// =============================================================================
TEST_CASE("ObjectPool<T> grows across slabs once initial capacity is exhausted",
          "[utilities][memory][objectpool]") {
    ObjectPool<std::uint64_t> pool(/*slot_capacity=*/4);
    std::vector<std::uint64_t*> ptrs;
    for (int i = 0; i < 10; ++i) ptrs.push_back(pool.allocate());

    REQUIRE(pool.live_count()  == 10);
    REQUIRE(pool.slab_count()  >= 3);    // 4 + 4 + at least 2 more
}

// =============================================================================
namespace {
struct Counted {
    static int  alive;
    int         tag;
    Counted(int t) : tag(t) { ++alive; }
    ~Counted()              { --alive; }
};
int Counted::alive = 0;
}  // namespace

TEST_CASE("ObjectPool<T>::construct/destroy invoke ctor and dtor exactly once",
          "[utilities][memory][objectpool]") {
    Counted::alive = 0;
    ObjectPool<Counted> pool(4);
    auto* a = pool.construct(7);
    auto* b = pool.construct(9);
    REQUIRE(Counted::alive == 2);
    REQUIRE(a->tag == 7);
    REQUIRE(b->tag == 9);

    pool.destroy(a);
    pool.destroy(b);
    REQUIRE(Counted::alive == 0);
    REQUIRE(pool.live_count() == 0);
}
