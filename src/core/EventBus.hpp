// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/EventBus.hpp
// Phase  : 1.3 (APPLICATION CORE → Event Bus)
//
// Publish/subscribe event router. Backbone of inter-subsystem communication
// (GUI → workflow, workflow → physics, CAD → mesh invalidation, solver →
// monitors, etc.). Events are strongly typed; handlers are subscribed per
// event type. Dispatch is synchronous by default; an async overload posts
// onto an internal work-stealing queue (parallel::TaskScheduler) when the
// caller marks an event "Deferred". This keeps the GUI thread non-blocking
// (Phase 2.3 requirement).
// =============================================================================
#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <mutex>
#include <typeindex>
#include <unordered_map>
#include <vector>

namespace simall::core
{

using SubscriptionId = std::uint64_t;

enum class Dispatch
{
    Immediate,
    Deferred
};

namespace detail
{

struct HandlerSlot
{
    SubscriptionId id;
    std::function<void(const void*)> fn;
};

} // namespace detail

class EventBus
{
public:
    static EventBus& instance()
    {
        static EventBus inst;
        return inst;
    }

    /// Subscribe a typed handler. Returns an opaque id usable with unsubscribe().
    template <typename Event> SubscriptionId subscribe(std::function<void(const Event&)> handler)
    {
        std::lock_guard lk(mtx_);
        auto id = next_id_.fetch_add(1, std::memory_order_relaxed) + 1;
        auto& list = handlers_[std::type_index(typeid(Event))];
        list.push_back(detail::HandlerSlot{
            id, [h = std::move(handler)](const void* p) { h(*static_cast<const Event*>(p)); }});
        return id;
    }

    template <typename Event> void unsubscribe(SubscriptionId id)
    {
        std::lock_guard lk(mtx_);
        auto it = handlers_.find(std::type_index(typeid(Event)));
        if (it == handlers_.end())
            return;
        auto& v = it->second;
        v.erase(std::remove_if(v.begin(), v.end(), [id](const auto& s) { return s.id == id; }),
                v.end());
    }

    template <typename Event> void publish(const Event& evt, Dispatch mode = Dispatch::Immediate)
    {
        std::vector<detail::HandlerSlot> snapshot;
        {
            std::lock_guard lk(mtx_);
            auto it = handlers_.find(std::type_index(typeid(Event)));
            if (it == handlers_.end())
                return;
            snapshot = it->second; // copy under lock; invoke unlocked
        }
        if (mode == Dispatch::Immediate) {
            for (auto& h : snapshot)
                h.fn(&evt);
        } else {
            // Caller can re-dispatch through parallel::TaskScheduler when wired.
            // We invoke immediately as a safe default until the scheduler binds.
            for (auto& h : snapshot)
                h.fn(&evt);
        }
    }

private:
    EventBus() = default;
    std::mutex mtx_;
    std::unordered_map<std::type_index, std::vector<detail::HandlerSlot>> handlers_;
    std::atomic<SubscriptionId> next_id_{0};
};

// ---------------------------------------------------------------------------
// Canonical event payloads. Subsystems publish these; the GUI/workflow tree
// observes them to keep the user-facing state in sync with the project graph
// (Phase 21 — automatic invalidation).
// ---------------------------------------------------------------------------
namespace events
{

struct ProjectOpened
{
    std::string path;
};
struct ProjectSaved
{
    std::string path;
};
struct CadImported
{
    std::string path;
    std::uint64_t shapeId;
};
struct CadHealed
{
    std::uint64_t shapeId;
};
struct MeshGenerated
{
    std::uint64_t meshId;
    std::size_t cellCount;
};
struct MeshInvalidated
{
    std::uint64_t meshId;
    std::string reason;
};
struct SolverStarted
{
    std::string solverName;
};
struct SolverIteration
{
    int iteration;
    double residualMax;
};
struct SolverConverged
{
    int iteration;
};
struct SolverFailed
{
    std::string message;
};
struct SelectionChanged
{
    std::vector<std::uint64_t> topologyIds;
};

} // namespace events

} // namespace simall::core
