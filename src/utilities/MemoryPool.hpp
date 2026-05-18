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
#include <cstddef>
#include <cstdint>
#include <memory>
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

}  // namespace simall::util
