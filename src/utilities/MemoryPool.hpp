// =============================================================================
// SimAll Beta - Utilities Subsystem
// File   : src/utilities/MemoryPool.hpp
// Phase  : 24 (MEMORY ARCHITECTURE)
//
// Bump-arena memory pool for transient solver scratch allocations. Used by
// the gradient-reconstruction kernel, AMG setup, and ghost-cell packers to
// eliminate per-iteration malloc traffic on hot paths. Thread-local arenas
// are managed by parallel::ThreadContext (see src/parallel).
// =============================================================================
#pragma once

#include "AlignedAllocator.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <utility>
#include <vector>

namespace simall::util {

class MemoryArena {
public:
    explicit MemoryArena(std::size_t initial_bytes = 1 << 20)  // 1 MiB
        : block_size_(initial_bytes) { add_block(initial_bytes); }

    MemoryArena(const MemoryArena&)            = delete;
    MemoryArena& operator=(const MemoryArena&) = delete;

    [[nodiscard]] void* allocate(std::size_t bytes,
                                 std::size_t alignment = alignof(std::max_align_t)) {
        auto& blk = blocks_.back();
        std::uintptr_t base = reinterpret_cast<std::uintptr_t>(blk.data.get()) + blk.offset;
        std::uintptr_t aligned = (base + alignment - 1) & ~(alignment - 1);
        std::size_t need = (aligned - base) + bytes;

        if (blk.offset + need > blk.capacity) {
            std::size_t new_cap = std::max(block_size_ * 2, bytes + alignment);
            add_block(new_cap);
            return allocate(bytes, alignment);
        }
        blk.offset += need;
        return reinterpret_cast<void*>(aligned);
    }

    /// Drop all allocations; backing memory retained (cheap reuse).
    void reset() noexcept {
        for (auto& b : blocks_) b.offset = 0;
    }

    /// Free everything except the first block.
    void shrink() {
        if (blocks_.size() > 1) blocks_.resize(1);
        blocks_.front().offset = 0;
    }

private:
    struct Block {
        std::unique_ptr<std::byte, void(*)(void*)> data;
        std::size_t capacity;
        std::size_t offset;
    };

    void add_block(std::size_t bytes) {
        void* mem = detail::aligned_alloc_bytes(bytes, kCacheLineBytes);
        blocks_.push_back(Block{
            std::unique_ptr<std::byte, void(*)(void*)>(
                static_cast<std::byte*>(mem),
                [](void* p){ detail::aligned_free_bytes(p); }),
            bytes, 0});
    }

    std::size_t        block_size_;
    std::vector<Block> blocks_;
};

// =============================================================================
// ObjectPool<T> — slab + freelist allocator for many same-typed objects.
// Phase 24 Pass 1.
//
// Complements `MemoryArena` (good for transient batch allocations cleared
// in bulk) with O(1) individual allocate/deallocate of a fixed object
// type.  Slabs are cache-line-aligned and never released until the pool
// is destroyed, so steady-state operation is fragmentation-free.  Not
// thread-safe — wrap with an external mutex or use one pool per thread.
// =============================================================================
template <typename T>
class ObjectPool {
public:
    static_assert(sizeof(T) >= sizeof(void*),
        "ObjectPool requires sizeof(T) >= sizeof(void*) for the freelist link");

    explicit ObjectPool(std::size_t slot_capacity = 1024)
        : slot_capacity_(slot_capacity == 0 ? 1 : slot_capacity) {}

    ObjectPool(const ObjectPool&)            = delete;
    ObjectPool& operator=(const ObjectPool&) = delete;

    ~ObjectPool() { /* slabs free via unique_ptr; T destructors are caller's job */ }

    /// Allocate raw storage for one T (no constructor invoked).
    [[nodiscard]] T* allocate() {
        if (free_head_ == nullptr) grow();
        Node* n = free_head_;
        free_head_ = n->next;
        ++live_count_;
        return reinterpret_cast<T*>(n);
    }

    /// Return raw storage to the pool (no destructor invoked).
    void deallocate(T* p) noexcept {
        if (!p) return;
        auto* n  = reinterpret_cast<Node*>(p);
        n->next  = free_head_;
        free_head_ = n;
        --live_count_;
    }

    /// allocate() + placement-new T(args...).
    template <typename... Args>
    [[nodiscard]] T* construct(Args&&... args) {
        T* p = allocate();
        ::new (static_cast<void*>(p)) T(std::forward<Args>(args)...);
        return p;
    }

    /// Manual destroy + deallocate.
    void destroy(T* p) noexcept {
        if (!p) return;
        p->~T();
        deallocate(p);
    }

    std::size_t live_count()    const noexcept { return live_count_; }
    std::size_t slab_count()    const noexcept { return slabs_.size(); }
    std::size_t slot_capacity() const noexcept { return slot_capacity_; }

private:
    union Node { Node* next; alignas(T) std::byte storage[sizeof(T)]; };

    void grow() {
        const std::size_t alignment = std::max(alignof(T), kCacheLineBytes);
        const std::size_t bytes     = slot_capacity_ * sizeof(Node);
        void* mem = detail::aligned_alloc_bytes(bytes, alignment);
        slabs_.emplace_back(static_cast<std::byte*>(mem),
                            [](void* p){ detail::aligned_free_bytes(p); });
        // Thread the new slab into the freelist (in reverse so first slot is head).
        auto* nodes = static_cast<Node*>(mem);
        for (std::size_t i = slot_capacity_; i-- > 0;) {
            nodes[i].next = free_head_;
            free_head_    = &nodes[i];
        }
    }

    using SlabPtr = std::unique_ptr<std::byte, void(*)(void*)>;
    std::size_t          slot_capacity_;
    Node*                free_head_ = nullptr;
    std::size_t          live_count_ = 0;
    std::vector<SlabPtr> slabs_;
};

}  // namespace simall::util
