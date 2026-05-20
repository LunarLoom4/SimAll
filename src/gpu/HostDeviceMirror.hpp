// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/HostDeviceMirror.hpp
// Phase  : 18.4 — RAII host/device buffer mirror (W13).
//
// A typed buffer that simultaneously exists on the host and (when available)
// on the device, with explicit one-way synchronisation primitives.  Built on
// top of DeviceMemoryPool so allocation cost is amortised.
//
// USAGE
// -----
//     HostDeviceMirror<double> b(n);
//     std::iota(b.host_begin(), b.host_end(), 0.0);
//     b.to_device();                       // H→D copy
//     gpu_axpy(2.0, x.device_ptr(), b.device_ptr(), n, stream);
//     b.to_host();                         // D→H copy
//     for (double v : b.host_view()) ...;
//
// SEMANTICS
// ---------
//   * Construction/resize allocates BOTH host and device storage in one
//     call.  Host storage is std::vector<T> for safe copy/move semantics;
//     device storage comes from DeviceMemoryPool.
//   * No automatic synchronisation — callers MUST invoke ``to_device`` /
//     ``to_host`` explicitly.  This avoids surprise PCIe traffic.
//   * Copy operations are explicit too (``operator=`` is deleted; use
//     ``clone()``).  Move is supported and cheap.
//   * On non-CUDA builds the "device" pointer aliases the host pointer and
//     to_device / to_host are no-ops, so client code stays identical.
//
// THREAD-SAFETY
// -------------
//   Single-owner.  Callers must serialise access; mirror objects are NOT
//   intended to be shared across host threads.  Multiple mirrors can however
//   be enqueued onto different streams concurrently.
// =============================================================================
#pragma once

#include "gpu/CudaContext.hpp"
#include "gpu/DeviceMemoryPool.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef SIMALL_HAVE_CUDA
#include <cuda_runtime.h>
#endif

namespace simall::gpu
{

template <typename T> class HostDeviceMirror
{
    static_assert(std::is_trivially_copyable_v<T>,
                  "HostDeviceMirror<T> requires a trivially-copyable element type "
                  "(POD numerics or aggregates).");

public:
    HostDeviceMirror() = default;

    explicit HostDeviceMirror(std::size_t n) { resize(n); }

    HostDeviceMirror(const HostDeviceMirror&) = delete;
    HostDeviceMirror& operator=(const HostDeviceMirror&) = delete;

    HostDeviceMirror(HostDeviceMirror&& o) noexcept
        : host_(std::move(o.host_))
        , device_(o.device_)
        , deviceBytes_(o.deviceBytes_)
        , ownsDevice_(o.ownsDevice_)
    {
        o.device_ = nullptr;
        o.deviceBytes_ = 0;
        o.ownsDevice_ = false;
    }

    HostDeviceMirror& operator=(HostDeviceMirror&& o) noexcept
    {
        if (this != &o) {
            free_device();
            host_ = std::move(o.host_);
            device_ = o.device_;
            deviceBytes_ = o.deviceBytes_;
            ownsDevice_ = o.ownsDevice_;
            o.device_ = nullptr;
            o.deviceBytes_ = 0;
            o.ownsDevice_ = false;
        }
        return *this;
    }

    ~HostDeviceMirror() { free_device(); }

    /// Resize host + device storage to ``n`` elements.  Existing contents
    /// are NOT preserved across resizes (mirrors std::vector::clear()).
    void resize(std::size_t n)
    {
        host_.assign(n, T{});
        const std::size_t bytes = n * sizeof(T);
#ifdef SIMALL_HAVE_CUDA
        if (is_cuda_available() && n > 0) {
            free_device();
            device_ = DeviceMemoryPool::instance().allocate(bytes);
            deviceBytes_ = bytes;
            ownsDevice_ = true;
            return;
        }
#endif
        // Serial fall-back: device pointer aliases host data.
        free_device();
        device_ = host_.data();
        deviceBytes_ = bytes;
        ownsDevice_ = false;
    }

    /// Host→device copy.  No-op on non-CUDA builds.
    void to_device()
    {
#ifdef SIMALL_HAVE_CUDA
        if (ownsDevice_ && device_ && !host_.empty())
            cudaMemcpy(device_, host_.data(), deviceBytes_, cudaMemcpyHostToDevice);
#endif
    }

    /// Asynchronous host→device on a specific stream.
    void to_device_async(void* stream)
    {
#ifdef SIMALL_HAVE_CUDA
        if (ownsDevice_ && device_ && !host_.empty())
            cudaMemcpyAsync(device_,
                            host_.data(),
                            deviceBytes_,
                            cudaMemcpyHostToDevice,
                            static_cast<cudaStream_t>(stream));
#else
        (void) stream;
#endif
    }

    /// Device→host copy.  No-op on non-CUDA builds.
    void to_host()
    {
#ifdef SIMALL_HAVE_CUDA
        if (ownsDevice_ && device_ && !host_.empty())
            cudaMemcpy(host_.data(), device_, deviceBytes_, cudaMemcpyDeviceToHost);
#endif
    }

    void to_host_async(void* stream)
    {
#ifdef SIMALL_HAVE_CUDA
        if (ownsDevice_ && device_ && !host_.empty())
            cudaMemcpyAsync(host_.data(),
                            device_,
                            deviceBytes_,
                            cudaMemcpyDeviceToHost,
                            static_cast<cudaStream_t>(stream));
#else
        (void) stream;
#endif
    }

    /// Memset device storage to zero.
    void zero_device()
    {
#ifdef SIMALL_HAVE_CUDA
        if (ownsDevice_ && device_)
            cudaMemset(device_, 0, deviceBytes_);
        else if (!host_.empty())
            std::memset(host_.data(), 0, deviceBytes_);
#else
        if (!host_.empty())
            std::memset(host_.data(), 0, deviceBytes_);
#endif
    }

    std::size_t size() const noexcept { return host_.size(); }
    std::size_t size_in_bytes() const noexcept { return deviceBytes_; }
    bool empty() const noexcept { return host_.empty(); }
    bool owns_device() const noexcept { return ownsDevice_; }

    T* host_ptr() noexcept { return host_.data(); }
    const T* host_ptr() const noexcept { return host_.data(); }
    void* device_ptr() noexcept { return device_; }
    const void* device_ptr() const noexcept { return device_; }

    T* host_begin() noexcept { return host_.data(); }
    T* host_end() noexcept { return host_.data() + host_.size(); }
    const T* host_begin() const noexcept { return host_.data(); }
    const T* host_end() const noexcept { return host_.data() + host_.size(); }

    /// std::vector reference for ergonomic algorithms.
    std::vector<T>& host_view() noexcept { return host_; }
    const std::vector<T>& host_view() const noexcept { return host_; }

    /// Deep copy (host + device).  Useful for back-up / restart vectors.
    HostDeviceMirror clone() const
    {
        HostDeviceMirror out(host_.size());
        std::copy(host_.begin(), host_.end(), out.host_.begin());
        out.to_device();
        return out;
    }

private:
    void free_device() noexcept
    {
#ifdef SIMALL_HAVE_CUDA
        if (ownsDevice_ && device_)
            DeviceMemoryPool::instance().release(device_, deviceBytes_);
#endif
        device_ = nullptr;
        deviceBytes_ = 0;
        ownsDevice_ = false;
    }

    std::vector<T> host_;
    void* device_ = nullptr;
    std::size_t deviceBytes_ = 0;
    bool ownsDevice_ = false;
};

} // namespace simall::gpu
