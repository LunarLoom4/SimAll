// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/PluginLoader.hpp
// Phase  : 1.3 (APPLICATION CORE → Plugin Host)
//
// Discovers and loads SimAll plugins (shared libraries) at runtime. Each
// plugin must export `CreatePlugin()` returning a heap-allocated IPlugin
// and `GetPluginAbi()` returning the ABI integer it was built against.
// Plugins are loaded in a sandboxed wrapper: load → ABI check → CreatePlugin
// → on_load() inside a try/catch. Any failure unloads the library and logs
// the failure (does not abort the host).
//
// Discovery: scanDirectory(path) walks a directory and loads every file
// matching the platform-specific shared-library extension (`*.dll`,
// `*.simallplugin`, `*.so`, `*.dylib`).
// =============================================================================
#pragma once

#include "plugins/IPlugin.hpp"

#include <filesystem>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace simall::core
{

struct LoadedPlugin
{
    std::filesystem::path path;
    void* handle = nullptr; // OS module handle
    std::unique_ptr<plugins::IPlugin> plugin;
    std::string name;
    std::string version;
    int abi = 0;
};

class PluginLoader
{
public:
    PluginLoader() = default;
    ~PluginLoader();

    PluginLoader(const PluginLoader&) = delete;
    PluginLoader& operator=(const PluginLoader&) = delete;

    /// Load a single plugin file. Returns the assigned cookie on success,
    /// or 0 on failure. Failure reason is logged via SIMALL_LOG_ERROR.
    std::size_t load(const std::filesystem::path& file);

    /// Scan a directory (non-recursive) and load every shared library
    /// matching `*.simallplugin` plus the platform's native ext as fallback.
    /// Returns the number of plugins successfully loaded.
    std::size_t scanDirectory(const std::filesystem::path& dir);

    /// Unload one plugin by cookie. Calls on_unload() before freeing.
    bool unload(std::size_t cookie);

    /// Unload all plugins.
    void unloadAll();

    /// Snapshot of currently loaded plugins.
    std::vector<LoadedPlugin> snapshot() const;

private:
    mutable std::mutex mtx_;
    std::vector<std::unique_ptr<LoadedPlugin>> plugins_;
};

} // namespace simall::core
