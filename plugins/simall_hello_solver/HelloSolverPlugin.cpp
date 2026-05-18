// =============================================================================
// SimAll Beta - Reference Plugin
// File   : plugins/simall_hello_solver/HelloSolverPlugin.cpp
// Week   : 20 — Reference plugin #3.  Smallest possible "Solver" plugin —
// proves the Solver category dispatch and registers a no-op PDE step that
// downstream tutorials walk through.  Demonstrates how a 3rd-party physics
// module attaches to the kernel without modifying SimAll source.
// =============================================================================
#include "../SimAllPluginSdk.hpp"

#include <iostream>

namespace {

class HelloSolverPlugin final : public simall::plugins::CategorisedPlugin {
public:
    HelloSolverPlugin()
        : CategorisedPlugin(simall::plugins::PluginCategory::Solver) {}

    std::string name()    const override { return "simall_hello_solver"; }
    std::string version() const override { return "1.0.0"; }

    void on_load() override {
        std::clog << "[plugin] hello-solver attached.  This plugin is the\n"
                     "          minimal worked example for the SDK tutorial.\n";
    }
    void on_unload() override {
        std::clog << "[plugin] hello-solver detached\n";
    }
};

}  // namespace

SIMALL_DECLARE_PLUGIN(HelloSolverPlugin)
