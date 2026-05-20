// =============================================================================
// SimAll Beta - Utilities Subsystem
// File   : src/utilities/AlignedAllocator.hpp
// Phase  : 24 (MEMORY ARCHITECTURE)
//
// 64-byte-aligned STL-compatible allocator used by all solver-critical SoA
// arrays. Alignment is mandated by Section 6.5 of the ultra-detailed spec
// (AVX-512 friendly). The allocator is stateless and noexcept-friendly so it
// composes with std::pmr and Eigen Maps without overhead.
// =============================================================================
#pragma once

#include <cstddef>
#include <cstdlib>
#include <limits>
#include <new>
#include <type_traits>
#include <vector>

namespace simall::util
{

inline constexpr std::size_t kCacheLineBytes = 64;
inline constexpr std::size_t kSimdAlignBytes = 64; // AVX-512 friendly

namespace detail
{

inline void* aligned_alloc_bytes(std::size_t bytes, std::size_t alignment)
{
    if (bytes == 0)
        return nullptr;
#if defined(_MSC_VER)
    void* p = _aligned_malloc(bytes, alignment);
    if (!p)
        throw std::bad_alloc{};
    return p;
#else
    void* p = nullptr;
    if (posix_memalign(&p, alignment, bytes) != 0)
        throw std::bad_alloc{};
    return p;
#endif
}

inline void aligned_free_bytes(void* p) noexcept
{
#if defined(_MSC_VER)
    _aligned_free(p);
#else
    std::free(p);
#endif
}

} // namespace detail

template <typename T, std::size_t Alignment = kSimdAlignBytes> class AlignedAllocator
{
    static_assert(Alignment >= alignof(T), "Alignment must satisfy the alignment of T");
    static_assert((Alignment & (Alignment - 1)) == 0, "Alignment must be a power of two");

public:
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using propagate_on_container_move_assignment = std::true_type;
    using is_always_equal = std::true_type;

    template <typename U> struct rebind
    {
        using other = AlignedAllocator<U, Alignment>;
    };

    constexpr AlignedAllocator() noexcept = default;

    template <typename U> constexpr AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept
    {
    }

    [[nodiscard]] T* allocate(size_type n)
    {
        if (n > std::numeric_limits<size_type>::max() / sizeof(T))
            throw std::bad_array_new_length{};
        return static_cast<T*>(detail::aligned_alloc_bytes(n * sizeof(T), Alignment));
    }

    void deallocate(T* p, size_type /*n*/) noexcept { detail::aligned_free_bytes(p); }

    template <typename U, std::size_t A2>
    constexpr bool operator==(const AlignedAllocator<U, A2>&) const noexcept
    {
        return Alignment == A2;
    }
    template <typename U, std::size_t A2>
    constexpr bool operator!=(const AlignedAllocator<U, A2>&) const noexcept
    {
        return Alignment != A2;
    }
};

/// SoA-friendly aligned vector. Used pervasively by mesh & solver storage.
template <typename T> using aligned_vector = std::vector<T, AlignedAllocator<T>>;

} // namespace simall::util
