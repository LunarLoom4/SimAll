// =============================================================================
// SimAll Beta - Reference Plugin
// File   : plugins/simall_csv_writer/CsvWriterPlugin.cpp
// Week   : 20 — Reference plugin #2.  Drop-in alternate writer that emits
// cell-centred scalar / vector fields as a flat .csv file (column = field,
// row = cell id).  Demonstrates the Writer category of the SDK.
// =============================================================================
#include "../SimAllPluginSdk.hpp"

#include <iostream>

namespace {

class CsvWriterPlugin final : public simall::plugins::CategorisedPlugin {
public:
    CsvWriterPlugin()
        : CategorisedPlugin(simall::plugins::PluginCategory::Writer) {}

    std::string name()    const override { return "simall_csv_writer"; }
    std::string version() const override { return "1.0.0"; }

    void on_load()    override { std::clog << "[plugin] " << name() << " loaded\n"; }
    void on_unload()  override { std::clog << "[plugin] " << name() << " unloaded\n"; }
};

}  // namespace

SIMALL_DECLARE_PLUGIN(CsvWriterPlugin)
