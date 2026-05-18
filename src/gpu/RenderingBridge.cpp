// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/RenderingBridge.cpp
// Phase  : 18.7 — GPU ↔ Visualization interop (W13).
// =============================================================================
#include "gpu/RenderingBridge.hpp"
#include "gpu/StreamScheduler.hpp"
#include "core/Logger.hpp"

#include <atomic>
#include <mutex>
#include <unordered_map>

namespace simall::gpu {

namespace {

struct Entry {
    BufferDescriptor desc;
    /// Exactly one of these is non-null depending on kind.
    HostDeviceMirror<float>* fmirror = nullptr;
    HostDeviceMirror<int>*   imirror = nullptr;
    /// Event the simulation last recorded (so the renderer can wait on it).
    void* pendingEvent = nullptr;
};

}  // namespace

struct RenderingBridge::Impl {
    mutable std::mutex                          mtx;
    std::unordered_map<BufferId, Entry>         buffers;
    std::atomic<BufferId>                       nextId{1};
};

// -----------------------------------------------------------------------------
RenderingBridge::RenderingBridge() : impl_(std::make_unique<Impl>()) {}
RenderingBridge::~RenderingBridge() = default;

RenderingBridge& RenderingBridge::instance() {
    static RenderingBridge bridge;
    return bridge;
}

// -----------------------------------------------------------------------------
BufferId RenderingBridge::register_scalar(const std::string& name,
                                          HostDeviceMirror<float>* m)
{
    if (!m) return kInvalidBufferId;
    std::lock_guard<std::mutex> lk(impl_->mtx);
    const BufferId id = impl_->nextId.fetch_add(1);
    Entry e{};
    e.desc.id   = id;
    e.desc.kind = BufferKind::Scalar;
    e.desc.name = name;
    e.desc.element_count = m->size();
    e.fmirror   = m;
    impl_->buffers.emplace(id, e);
    SIMALL_LOG_INFO("Gpu", "RenderingBridge: registered scalar buffer '",
        name, "' (id=", id, ", n=", m->size(), ")");
    return id;
}

BufferId RenderingBridge::register_vector3(const std::string& name,
                                           HostDeviceMirror<float>* m)
{
    if (!m) return kInvalidBufferId;
    std::lock_guard<std::mutex> lk(impl_->mtx);
    const BufferId id = impl_->nextId.fetch_add(1);
    Entry e{};
    e.desc.id   = id;
    e.desc.kind = BufferKind::Vector3;
    e.desc.name = name;
    e.desc.element_count = m->size() / 3;
    e.fmirror   = m;
    impl_->buffers.emplace(id, e);
    SIMALL_LOG_INFO("Gpu", "RenderingBridge: registered vector3 buffer '",
        name, "' (id=", id, ", n=", e.desc.element_count, ")");
    return id;
}

BufferId RenderingBridge::register_indices(const std::string& name,
                                           HostDeviceMirror<int>* m)
{
    if (!m) return kInvalidBufferId;
    std::lock_guard<std::mutex> lk(impl_->mtx);
    const BufferId id = impl_->nextId.fetch_add(1);
    Entry e{};
    e.desc.id   = id;
    e.desc.kind = BufferKind::Indices;
    e.desc.name = name;
    e.desc.element_count = m->size();
    e.imirror   = m;
    impl_->buffers.emplace(id, e);
    SIMALL_LOG_INFO("Gpu", "RenderingBridge: registered index buffer '",
        name, "' (id=", id, ", n=", m->size(), ")");
    return id;
}

// -----------------------------------------------------------------------------
void RenderingBridge::unregister(BufferId id) {
    std::lock_guard<std::mutex> lk(impl_->mtx);
    auto it = impl_->buffers.find(id);
    if (it == impl_->buffers.end()) return;
    if (it->second.pendingEvent)
        StreamScheduler::instance().release_event(it->second.pendingEvent);
    impl_->buffers.erase(it);
}

// -----------------------------------------------------------------------------
void RenderingBridge::request_gl_interop(BufferId id) {
    std::lock_guard<std::mutex> lk(impl_->mtx);
    auto it = impl_->buffers.find(id);
    if (it == impl_->buffers.end()) return;
    it->second.desc.gl_interop_requested = true;
    // W14 will instantiate the actual cudaGraphicsGLRegisterBuffer call
    // when SceneGraph asks for a renderable resource.
    SIMALL_LOG_INFO("Gpu", "RenderingBridge: GL interop requested on buffer ",
        id, " (will activate in visualization layer W14)");
}

// -----------------------------------------------------------------------------
void* RenderingBridge::acquire(BufferId id) {
    std::lock_guard<std::mutex> lk(impl_->mtx);
    auto it = impl_->buffers.find(id);
    if (it == impl_->buffers.end()) return nullptr;

    // Wait on the simulation-side event before handing the buffer over.
    if (it->second.pendingEvent) {
        // Render side typically runs on stream 0 (GL/VK presentation).
        StreamScheduler::instance().wait_event(nullptr, it->second.pendingEvent);
    }
    if (it->second.fmirror)
        return it->second.fmirror->device_ptr();
    if (it->second.imirror)
        return it->second.imirror->device_ptr();
    return nullptr;
}

void RenderingBridge::release(BufferId id) {
    std::lock_guard<std::mutex> lk(impl_->mtx);
    auto it = impl_->buffers.find(id);
    if (it == impl_->buffers.end()) return;
    if (it->second.pendingEvent) {
        StreamScheduler::instance().release_event(it->second.pendingEvent);
        it->second.pendingEvent = nullptr;
    }
}

void RenderingBridge::insert_barrier(BufferId id, void* simStream) {
    std::lock_guard<std::mutex> lk(impl_->mtx);
    auto it = impl_->buffers.find(id);
    if (it == impl_->buffers.end()) return;
    // Replace any stale event from a previous frame.
    if (it->second.pendingEvent)
        StreamScheduler::instance().release_event(it->second.pendingEvent);
    it->second.pendingEvent = StreamScheduler::instance().record_event(simStream);
}

// -----------------------------------------------------------------------------
std::vector<BufferDescriptor> RenderingBridge::list_buffers() const {
    std::lock_guard<std::mutex> lk(impl_->mtx);
    std::vector<BufferDescriptor> out;
    out.reserve(impl_->buffers.size());
    for (const auto& [id, e] : impl_->buffers) out.push_back(e.desc);
    return out;
}

}  // namespace simall::gpu
