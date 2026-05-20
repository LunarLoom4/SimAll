// =============================================================================
// SimAll Beta — Core Unit Tests
// File   : tests/unit/core/test_plugin_category.cpp
// Phase  : 22 Pass 4
//
// Validates the Phase 22.4 PluginCategory expansion:
//   * Legacy ordinals (Solver=0, PostProcessor=1, Writer=2, Reader=3,
//     Generic=4) remain frozen for ABI stability.
//   * Every enum value [0 .. _Count-1] maps to a non-empty, non-"unknown"
//     name and to a recognised subsystem path (or empty for cross-cutting).
//   * Names are unique across the enum.
//   * core::PluginRegistry buckets correctly under the new categories.
// =============================================================================
#include "core/PluginRegistry.hpp"
#include "plugins/IPlugin.hpp"
#include "plugins/SimAllPluginSdk.hpp"

#include <catch2/catch_test_macros.hpp>

#include <set>
#include <string>
#include <string_view>
#include <vector>

using simall::plugins::category_name;
using simall::plugins::category_subsystem_path;
using simall::plugins::PluginCategory;

namespace
{

constexpr int kCount = static_cast<int>(PluginCategory::_Count);

PluginCategory cat(int i)
{
    return static_cast<PluginCategory>(i);
}

} // namespace

TEST_CASE("PluginCategory legacy ordinals are frozen", "[core][plugin][abi]")
{
    CHECK(static_cast<int>(PluginCategory::Solver) == 0);
    CHECK(static_cast<int>(PluginCategory::PostProcessor) == 1);
    CHECK(static_cast<int>(PluginCategory::Writer) == 2);
    CHECK(static_cast<int>(PluginCategory::Reader) == 3);
    CHECK(static_cast<int>(PluginCategory::Generic) == 4);
}

TEST_CASE("PluginCategory _Count matches expanded enum size", "[core][plugin]")
{
    // Legacy 5 + Phase 22.4 expansion 30 = 35.
    CHECK(kCount == 35);
}

TEST_CASE("PluginCategory category_name() is total and non-empty", "[core][plugin]")
{
    for (int i = 0; i < kCount; ++i) {
        const char* n = category_name(cat(i));
        REQUIRE(n != nullptr);
        const std::string_view sv{n};
        INFO("ordinal " << i);
        CHECK_FALSE(sv.empty());
        CHECK(sv != "unknown");
    }
    // Out-of-range / sentinel paths return "unknown" defensively.
    CHECK(std::string(category_name(PluginCategory::_Count)) == "unknown");
}

TEST_CASE("PluginCategory category_name() values are unique", "[core][plugin]")
{
    std::set<std::string> seen;
    for (int i = 0; i < kCount; ++i) {
        const auto [it, inserted] = seen.emplace(category_name(cat(i)));
        INFO("ordinal " << i << " name '" << *it << "'");
        CHECK(inserted);
    }
    CHECK(seen.size() == static_cast<std::size_t>(kCount));
}

TEST_CASE("PluginCategory subsystem-aligned categories report src/ path", "[core][plugin]")
{
    struct Row
    {
        PluginCategory c;
        const char* path;
    };
    const Row rows[] = {
        {PluginCategory::Solver, "src/solver"},
        {PluginCategory::Cad, "src/cad"},
        {PluginCategory::Meshing, "src/meshing"},
        {PluginCategory::Materials, "src/materials"},
        {PluginCategory::Turbulence, "src/turbulence"},
        {PluginCategory::Multiphase, "src/multiphase"},
        {PluginCategory::Combustion, "src/combustion"},
        {PluginCategory::Radiation, "src/radiation"},
        {PluginCategory::HeatTransfer, "src/heat_transfer"},
        {PluginCategory::ElectroMagnetics, "src/emag"},
        {PluginCategory::Acoustics, "src/acoustics"},
        {PluginCategory::Porous, "src/porous"},
        {PluginCategory::Particles, "src/particles"},
        {PluginCategory::Adjoint, "src/adjoint"},
        {PluginCategory::Optimization, "src/optimization"},
        {PluginCategory::Rom, "src/rom"},
        {PluginCategory::Io, "src/io"},
        {PluginCategory::Scripting, "src/scripting"},
        {PluginCategory::Parallel, "src/parallel"},
        {PluginCategory::Gpu, "src/gpu"},
        {PluginCategory::Gui, "src/gui"},
        {PluginCategory::Amr, "src/amr"},
        {PluginCategory::Rotating, "src/rotating"},
        {PluginCategory::Dynamics, "src/dynamics"},
        {PluginCategory::Morphing, "src/morphing"},
        {PluginCategory::ImmersedBoundary, "src/ibm"},
        {PluginCategory::Visualization, "src/visualization"},
        {PluginCategory::Zones, "src/zones"},
        {PluginCategory::Utility, "src/utilities"},
        {PluginCategory::Core, "src/core"},
    };
    for (const auto& r : rows) {
        INFO("category " << category_name(r.c));
        CHECK(std::string(category_subsystem_path(r.c)) == r.path);
    }
}

TEST_CASE("PluginCategory cross-cutting categories report empty path", "[core][plugin]")
{
    for (auto c : {PluginCategory::PostProcessor,
                   PluginCategory::Writer,
                   PluginCategory::Reader,
                   PluginCategory::Generic,
                   PluginCategory::BoundaryCondition}) {
        INFO("category " << category_name(c));
        CHECK(std::string(category_subsystem_path(c)).empty());
    }
}

// -- Mini integration with PluginRegistry: ensures the registry buckets
//    plugins under the newly-added categories without losing entries.
namespace
{

class StubPlugin : public simall::plugins::IPlugin
{
public:
    std::string name() const override { return "stub"; }
    std::string version() const override { return "0.0.0"; }
    void on_load() override {}
    void on_unload() override {}
    int abi() const override { return simall::plugins::kPluginAbiVersion; }
};

} // namespace

TEST_CASE("PluginRegistry accepts new Phase 22.4 categories", "[core][plugin][registry]")
{
    simall::core::PluginRegistry reg;
    StubPlugin a, b, c;

    reg.register_plugin(&a, PluginCategory::Turbulence);
    reg.register_plugin(&b, PluginCategory::ElectroMagnetics);
    reg.register_plugin(&c, PluginCategory::BoundaryCondition);

    CHECK(reg.size() == 3);
    CHECK(reg.by_category(PluginCategory::Turbulence).size() == 1);
    CHECK(reg.by_category(PluginCategory::ElectroMagnetics).size() == 1);
    CHECK(reg.by_category(PluginCategory::BoundaryCondition).size() == 1);
    CHECK(reg.by_category(PluginCategory::Cad).empty());

    reg.unregister(&b);
    CHECK(reg.size() == 2);
    CHECK(reg.by_category(PluginCategory::ElectroMagnetics).empty());
}
