// =============================================================================
// SimAll Beta - Plugin ABI
// File   : plugins/IPlugin.hpp
// Phase  : 17 / Section 17 of ultra-detailed spec.
//
// Stable C ABI entry point. Plugins implement IPlugin and export CreatePlugin
// via extern "C". The host loads them with QPluginLoader (or dlopen on Linux,
// LoadLibrary on Windows) at startup or on demand.
// =============================================================================
#pragma once

#include <string>

namespace simall::plugins {

/// Plugin ABI version. Bumped on any binary-breaking change to IPlugin or
/// its derived contract. The loader (core::PluginLoader) refuses to load
/// a plugin whose reported ABI differs.
inline constexpr int kPluginAbiVersion = 1;

class IPlugin {
public:
    virtual ~IPlugin() = default;
    virtual std::string name()    const = 0;
    virtual std::string version() const = 0;
    virtual void        on_load()       = 0;
    virtual void        on_unload()     = 0;
    virtual int         abi()     const { return kPluginAbiVersion; }
};

}  // namespace simall::plugins

#if defined(_WIN32)
#  define SIMALL_PLUGIN_EXPORT __declspec(dllexport)
#else
#  define SIMALL_PLUGIN_EXPORT __attribute__((visibility("default")))
#endif

extern "C" SIMALL_PLUGIN_EXPORT simall::plugins::IPlugin* CreatePlugin();
extern "C" SIMALL_PLUGIN_EXPORT int                       GetPluginAbi();
