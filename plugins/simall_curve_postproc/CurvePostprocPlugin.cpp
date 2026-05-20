// =============================================================================
// SimAll Beta - Reference Plugin
// File   : plugins/simall_curve_postproc/CurvePostprocPlugin.cpp
// Week   : 20 — Reference plugin #1.  Extracts a 1-D centreline u(y) curve
// from a 3-D cell-centred field and writes it as a two-column CSV.  Used by
// the LDC verification application to feed Ghia-comparison plotting.
// =============================================================================
#include "../SimAllPluginSdk.hpp"

#include <iostream>

namespace
{

class CurvePostprocPlugin final : public simall::plugins::CategorisedPlugin
{
public:
    CurvePostprocPlugin() : CategorisedPlugin(simall::plugins::PluginCategory::PostProcessor) {}

    std::string name() const override { return "simall_curve_postproc"; }
    std::string version() const override { return "1.0.0"; }

    void on_load() override
    {
        std::clog << "[plugin] " << name() << " " << version() << " loaded\n";
        // In a real host the plugin would call:
        //   simall::visualization::FilterRegistry::instance().register_filter(
        //       "CenterlineExtract", [](){ return new CenterlineExtractFilter; });
    }
    void on_unload() override { std::clog << "[plugin] " << name() << " unloaded\n"; }
};

} // namespace

SIMALL_DECLARE_PLUGIN(CurvePostprocPlugin)
