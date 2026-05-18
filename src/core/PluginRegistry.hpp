// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/PluginRegistry.hpp
// Phase  : 17 — host-side dispatch of loaded plugins by category.
//
// `PluginLoader` discovers shared libraries and instantiates `IPlugin`
// objects.  `PluginRegistry` is the *post-load* layer that lets host
// subsystems (post-processor toolbar, result-writer menu, BC list, …)
// retrieve the slice of currently loaded plugins relevant to them, without
// having to type-test every IPlugin pointer.
//
// A plugin participates by deriving from `plugins::CategorisedPlugin`
// (or any class that exposes a `category()` accessor); the registry uses
// runtime type information to bucket it.  Generic IPlugin instances that
// do not advertise a category are placed under `PluginCategory::Generic`.
// =============================================================================
#pragma once

#include "plugins/IPlugin.hpp"
#include "plugins/SimAllPluginSdk.hpp"

#include <mutex>
#include <unordered_map>
#include <vector>

namespace simall::core {

class PluginRegistry {
public:
    using Category = plugins::PluginCategory;

    /// Register a non-owning pointer (the loader keeps ownership) under the
    /// supplied category.  Idempotent on duplicate (plugin, category) pairs.
    void register_plugin(plugins::IPlugin* p, Category cat);

    /// Convenience: query category from a CategorisedPlugin (if available)
    /// and register; otherwise place in Generic.
    void register_auto(plugins::IPlugin* p);

    /// Forget a plugin pointer (call from PluginLoader::unload before the
    /// underlying object is destroyed).
    void unregister(plugins::IPlugin* p);

    /// All currently registered plugins in a given category, in registration
    /// order.  Returned vector is a copy — safe to iterate while plugins
    /// load/unload.
    std::vector<plugins::IPlugin*> by_category(Category cat) const;

    /// All plugins, regardless of category.
    std::vector<plugins::IPlugin*> all() const;

    /// Total count.
    std::size_t size() const;
    void        clear();

    /// Process-wide singleton.  The UI and CLI both consult this instance.
    static PluginRegistry& instance();

private:
    mutable std::mutex                                            mtx_;
    std::unordered_map<Category, std::vector<plugins::IPlugin*>> buckets_;
};

}  // namespace simall::core
