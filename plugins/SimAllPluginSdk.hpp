// =============================================================================
// SimAll Beta - Plugin SDK
// File   : plugins/SimAllPluginSdk.hpp
// Week   : 20 — Public-facing developer SDK header.
//
// Wraps the low-level IPlugin ABI in a small set of convenience macros and
// a default base class so out-of-tree plugin authors can write:
//
//     #include <simall/SimAllPluginSdk.hpp>
//
//     class MyFilter : public simall::plugins::PluginBase {
//     public:
//         std::string name()    const override { return "MyFilter"; }
//         std::string version() const override { return "1.0.0";    }
//         void on_load()    override { /* register */ }
//         void on_unload()  override { /* clean up */ }
//     };
//
//     SIMALL_DECLARE_PLUGIN(MyFilter)
//
// This header is **header-only**, depends only on IPlugin.hpp, and is the
// supported entry point for 3rd-party extensions.  Direct inheritance from
// IPlugin remains supported but is documented as "advanced".
// =============================================================================
#pragma once

#include "IPlugin.hpp"

#include <string>

namespace simall::plugins {

/// Convenience base class — provides storage for cached name / version, lets
/// the author override only the lifecycle hooks they care about.
class PluginBase : public IPlugin {
public:
    PluginBase() = default;
    ~PluginBase() override = default;

    PluginBase(const PluginBase&)            = delete;
    PluginBase& operator=(const PluginBase&) = delete;

    std::string name()    const override { return "unnamed"; }
    std::string version() const override { return "0.0.0";   }
    void        on_load()    override {}
    void        on_unload()  override {}
    int         abi()      const override { return kPluginAbiVersion; }
};

/// Three categories the host uses to slot the plugin into the right registry.
enum class PluginCategory { Solver, PostProcessor, Writer, Reader, Generic };

class CategorisedPlugin : public PluginBase {
public:
    explicit CategorisedPlugin(PluginCategory cat) : category_(cat) {}
    [[nodiscard]] PluginCategory category() const { return category_; }
private:
    PluginCategory category_;
};

}  // namespace simall::plugins

/// Define `CreatePlugin` and `GetPluginAbi` for a class deriving from
/// `simall::plugins::IPlugin`.  Use once per .cpp file at namespace scope.
#define SIMALL_DECLARE_PLUGIN(ClassName)                                       \
    extern "C" SIMALL_PLUGIN_EXPORT simall::plugins::IPlugin* CreatePlugin() { \
        return new ClassName();                                                \
    }                                                                          \
    extern "C" SIMALL_PLUGIN_EXPORT int GetPluginAbi() {                       \
        return simall::plugins::kPluginAbiVersion;                             \
    }
