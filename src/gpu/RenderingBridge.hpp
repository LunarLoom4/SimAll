// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/RenderingBridge.hpp
// Phase  : 18.7 — GPU ↔ Visualization interop (W13).
//
// Provides a zero-copy pathway for fields living on the device to be drawn
// by VTK / Vulkan / OpenGL renderers in src/visualization without round-
// tripping through host memory.
//
// CONCEPTS
// --------
//   * **Buffer registration**: a HostDeviceMirror<float> (positions, scalars,
//     vectors) is registered once; the bridge tracks a stable BufferId.
//   * **Acquire / release**: the renderer calls ``acquire(id)`` before a
//     draw call to obtain a pointer it can bind to a graphics resource
//     (cudaGraphicsResourceMap on CUDA, raw host pointer on the fall-back).
//     ``release(id)`` returns control to the simulation kernel.
//   * **Synchronisation**: the bridge inserts an event on the simulation
//     stream and makes the rendering acquire wait on it (no host sync on
//     the hot loop).
//
// FALL-BACK
// ---------
// When SIMALL_HAVE_CUDA is undefined, ``acquire`` simply returns the host
// pointer of the mirror; ``release`` is a no-op.  This lets visualization
// code be written against a single API regardless of GPU presence.
//
// CUDA-GL INTEROP
// ---------------
// True zero-copy interop requires the rendering API's buffer object to be
// pre-registered via cudaGraphicsGLRegisterBuffer (or
// cudaGraphicsVulkanGetMemoryAcquire).  Because the visualization layer is
// not built yet (W14), this milestone provides:
//   * the registration plumbing (ID generation, mirror tracking),
//   * a ``request_gl_interop()`` placeholder that records intent and
//     reports a friendly diagnostic when the renderer asks for a GL
//     resource we cannot yet provide.
// The actual GL resource registration is wired in W14 when SceneGraph
// is created.
// =============================================================================
#pragma once

#include "gpu/HostDeviceMirror.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace simall::gpu
{

using BufferId = std::uint64_t;
constexpr BufferId kInvalidBufferId = 0;

enum class BufferKind
{
    Scalar,  ///< 1 float per cell/point
    Vector3, ///< 3 floats per cell/point (xyz)
    Indices, ///< 32-bit indices (e.g. triangle list)
};

struct BufferDescriptor
{
    BufferId id = kInvalidBufferId;
    BufferKind kind = BufferKind::Scalar;
    std::string name; ///< human-readable for diagnostics / pickers
    std::size_t element_count = 0;
    bool gl_interop_requested = false;
    bool gl_interop_active = false; ///< true once W14 wires GL
};

class RenderingBridge
{
public:
    RenderingBridge();
    ~RenderingBridge();

    static RenderingBridge& instance();

    /// Register a typed mirror with the bridge.  The mirror MUST outlive
    /// any pending ``acquire`` calls.  Returns a stable BufferId.
    BufferId register_scalar(const std::string& name, HostDeviceMirror<float>* m);
    BufferId register_vector3(const std::string& name, HostDeviceMirror<float>* m);
    BufferId register_indices(const std::string& name, HostDeviceMirror<int>* m);

    /// Drop a buffer (renderer will get a null pointer on next acquire).
    void unregister(BufferId id);

    /// Inform the bridge that the renderer wants GL/Vulkan interop on this
    /// buffer.  No-op on builds without a graphics interop backend.
    void request_gl_interop(BufferId id);

    /// Pointer hand-off to renderer.  Returns:
    ///   - on CUDA builds: device pointer (so a graphics shader bound via
    ///     interop can read it directly);
    ///   - on host builds: host pointer.
    void* acquire(BufferId id);

    /// Signal that the renderer is done; the simulation can mutate the
    /// buffer again.  Drops any inserted synchronisation event.
    void release(BufferId id);

    /// Make the simulation kernel running on ``simStream`` complete before
    /// the renderer acquires ``id``.  Lightweight (records an event).
    void insert_barrier(BufferId id, void* simStream);

    /// Snapshot for property panels.
    std::vector<BufferDescriptor> list_buffers() const;

    RenderingBridge(const RenderingBridge&) = delete;
    RenderingBridge& operator=(const RenderingBridge&) = delete;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace simall::gpu
