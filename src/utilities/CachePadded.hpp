// =============================================================================
// SimAll Beta - Utilities Subsystem
// File   : src/utilities/CachePadded.hpp
// Phase  : 24 Pass 2  (cache-aware structures)
//
// `CachePadded<T>` is a trivial wrapper that forces an aligned + size-padded
// layout so adjacent elements never share a CPU cache line.  Use for per-
// thread counters, residual accumulators, per-rank stats buffers, etc.,
// where independent threads write to neighbouring entries.  Without this
// wrapper, those writes invalidate the same cache line on every neighbour
// core ("false sharing") and serialise what looks like parallel work.
//
// Layout invariant: `sizeof(CachePadded<T>) % kCacheLineBytes == 0` and
// `alignof(CachePadded<T>) == kCacheLineBytes`.  Pre-condition:
// `sizeof(T) <= kCacheLineBytes` is *not* required; larger T's just round
// up to the next cache-line multiple.
// =============================================================================
#pragma once

#include "AlignedAllocator.hpp"

#include <cstddef>
#include <type_traits>
#include <utility>

namespace simall::util {

template <typename T>
class alignas(kCacheLineBytes) CachePadded {
public:
    constexpr CachePadded() noexcept(std::is_nothrow_default_constructible_v<T>)
        = default;

    template <typename... Args,
              typename = std::enable_if_t<std::is_constructible_v<T, Args&&...>>>
    constexpr explicit CachePadded(Args&&... args)
        noexcept(std::is_nothrow_constructible_v<T, Args&&...>)
        : value_(std::forward<Args>(args)...) {}

    constexpr       T& value()       noexcept { return value_; }
    constexpr const T& value() const noexcept { return value_; }

    constexpr       T* operator->()       noexcept { return &value_; }
    constexpr const T* operator->() const noexcept { return &value_; }

    constexpr       T& operator*()        noexcept { return value_; }
    constexpr const T& operator*()  const noexcept { return value_; }

private:
    T value_{};
    // Pad up to the next cache-line multiple.  When sizeof(T) is already
    // a multiple, the trailing array length evaluates to 0 (zero-length
    // arrays are illegal in standard C++, so we collapse to a 1-byte
    // tail and let the alignas rule handle the rounding).
    static constexpr std::size_t kPad =
        (kCacheLineBytes - (sizeof(T) % kCacheLineBytes)) % kCacheLineBytes;
    [[maybe_unused]] std::byte pad_[kPad == 0 ? 1 : kPad]{};
};

static_assert(alignof(CachePadded<int>) == kCacheLineBytes,
              "CachePadded<int> must be cache-line aligned");
static_assert(sizeof(CachePadded<int>) % kCacheLineBytes == 0,
              "CachePadded<int> size must be a cache-line multiple");

}  // namespace simall::util
