// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/ServiceLocator.hpp
// Phase  : 1.3 (APPLICATION CORE → Service Locator)
//
// Type-indexed global subsystem registry. Subsystems register themselves once
// during application bootstrap and are retrieved by interface type. The
// locator owns the subsystem lifetime via std::unique_ptr and tears them down
// in reverse registration order — critical for clean shutdown of GUI ↔ Solver
// ↔ CAD dependencies.
//
// THIS IS NOT A SUBSTITUTE FOR DEPENDENCY INJECTION. Subsystems that have a
// true compile-time dependency take collaborators by constructor. The locator
// is only used at the application-shell boundary (GUI, scripting, plugins).
// =============================================================================
#pragma once

#include <cassert>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <typeindex>
#include <unordered_map>
#include <vector>

namespace simall::core
{

class ServiceLocator
{
public:
    static ServiceLocator& instance()
    {
        static ServiceLocator inst;
        return inst;
    }

    template <typename Interface, typename Impl, typename... Args>
    Interface& provide(Args&&... args)
    {
        static_assert(std::is_base_of_v<Interface, Impl>, "Impl must inherit from Interface");
        std::lock_guard lk(mtx_);
        auto idx = std::type_index(typeid(Interface));
        if (services_.count(idx))
            throw std::logic_error("Service already registered");
        auto owned = std::make_unique<Impl>(std::forward<Args>(args)...);
        Interface* raw = owned.get();
        // Erase to void* via a custom deleter that invokes the real destructor.
        auto* erased = static_cast<void*>(owned.release());
        auto deleter = [](void* p) { delete static_cast<Impl*>(p); };
        services_.emplace(idx, Slot{erased, std::move(deleter)});
        order_.push_back(idx);
        return *raw;
    }

    template <typename Interface> Interface& get() const
    {
        std::lock_guard lk(mtx_);
        auto it = services_.find(std::type_index(typeid(Interface)));
        if (it == services_.end())
            throw std::logic_error("Service not registered");
        return *static_cast<Interface*>(it->second.ptr);
    }

    template <typename Interface> Interface* try_get() const noexcept
    {
        std::lock_guard lk(mtx_);
        auto it = services_.find(std::type_index(typeid(Interface)));
        return it == services_.end() ? nullptr : static_cast<Interface*>(it->second.ptr);
    }

    void shutdown()
    {
        std::lock_guard lk(mtx_);
        // Reverse-order destruction: GUI dies before solver dies before CAD.
        for (auto it = order_.rbegin(); it != order_.rend(); ++it) {
            auto& slot = services_.at(*it);
            slot.deleter(slot.ptr);
        }
        services_.clear();
        order_.clear();
    }

    ~ServiceLocator() { shutdown(); }

private:
    ServiceLocator() = default;
    ServiceLocator(const ServiceLocator&) = delete;
    ServiceLocator& operator=(const ServiceLocator&) = delete;

    struct Slot
    {
        void* ptr;
        std::function<void(void*)> deleter;
    };

    mutable std::mutex mtx_;
    std::unordered_map<std::type_index, Slot> services_;
    std::vector<std::type_index> order_;
};

} // namespace simall::core
