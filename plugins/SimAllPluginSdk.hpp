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

namespace simall::plugins
{

/// Convenience base class — provides storage for cached name / version, lets
/// the author override only the lifecycle hooks they care about.
class PluginBase : public IPlugin
{
public:
    PluginBase() = default;
    ~PluginBase() override = default;

    PluginBase(const PluginBase&) = delete;
    PluginBase& operator=(const PluginBase&) = delete;

    std::string name() const override { return "unnamed"; }
    std::string version() const override { return "0.0.0"; }
    void on_load() override {}
    void on_unload() override {}
    int abi() const override { return kPluginAbiVersion; }
};

/// Plugin category. The first five values (Solver, PostProcessor, Writer,
/// Reader, Generic) are LEGACY (pre-Phase-22.4) and MUST keep these
/// ordinals to preserve ABI for already-compiled plugins. The remaining
/// values were added in Phase 22 Pass 4 to give first-class slotting to
/// every one of the 31 SimAll subsystems plus the cross-cutting
/// BoundaryCondition role (which spans solver / em / structural / thermal /
/// acoustics / electrochemistry, see src/core/BoundaryConditionRegistry).
///
/// IMPORTANT: append new values at the END only. Never reorder or remove
/// — the integer value is part of the plugin ABI.
enum class PluginCategory
{
    // -------- Legacy (Week-20 SDK release) — ordinals 0..4 frozen --------
    Solver = 0,
    PostProcessor = 1,
    Writer = 2,
    Reader = 3,
    Generic = 4,

    // -------- Phase 22 Pass 4: subsystem-aligned slots (one per src/) ----
    Cad = 5,
    Meshing = 6,
    Materials = 7,
    Turbulence = 8,
    Multiphase = 9,
    Combustion = 10,
    Radiation = 11,
    HeatTransfer = 12,
    ElectroMagnetics = 13,
    Acoustics = 14,
    Porous = 15,
    Particles = 16,
    Adjoint = 17,
    Optimization = 18,
    Rom = 19,
    Io = 20,
    Scripting = 21,
    Parallel = 22,
    Gpu = 23,
    Gui = 24,
    Amr = 25,
    Rotating = 26,
    Dynamics = 27,
    Morphing = 28,
    ImmersedBoundary = 29,
    Visualization = 30,
    Zones = 31,
    Utility = 32,
    Core = 33,
    BoundaryCondition = 34,

    // Sentinel — always equal to (last_real_value + 1). Tests rely on this
    // to verify category_name() / category_subsystem_path() are total.
    _Count = 35,
};

/// Stable identifier ("solver", "meshing", …) suitable for log lines, .ini
/// keys, JSON manifests, and plugin-discovery filters. Returns "unknown"
/// for any value outside the declared enum (defensive — should not happen
/// after Phase 22.4 — but protects against ABI-mismatched plugins).
[[nodiscard]] inline const char* category_name(PluginCategory c) noexcept
{
    switch (c) {
    case PluginCategory::Solver:
        return "solver";
    case PluginCategory::PostProcessor:
        return "post_processor";
    case PluginCategory::Writer:
        return "writer";
    case PluginCategory::Reader:
        return "reader";
    case PluginCategory::Generic:
        return "generic";
    case PluginCategory::Cad:
        return "cad";
    case PluginCategory::Meshing:
        return "meshing";
    case PluginCategory::Materials:
        return "materials";
    case PluginCategory::Turbulence:
        return "turbulence";
    case PluginCategory::Multiphase:
        return "multiphase";
    case PluginCategory::Combustion:
        return "combustion";
    case PluginCategory::Radiation:
        return "radiation";
    case PluginCategory::HeatTransfer:
        return "heat_transfer";
    case PluginCategory::ElectroMagnetics:
        return "emag";
    case PluginCategory::Acoustics:
        return "acoustics";
    case PluginCategory::Porous:
        return "porous";
    case PluginCategory::Particles:
        return "particles";
    case PluginCategory::Adjoint:
        return "adjoint";
    case PluginCategory::Optimization:
        return "optimization";
    case PluginCategory::Rom:
        return "rom";
    case PluginCategory::Io:
        return "io";
    case PluginCategory::Scripting:
        return "scripting";
    case PluginCategory::Parallel:
        return "parallel";
    case PluginCategory::Gpu:
        return "gpu";
    case PluginCategory::Gui:
        return "gui";
    case PluginCategory::Amr:
        return "amr";
    case PluginCategory::Rotating:
        return "rotating";
    case PluginCategory::Dynamics:
        return "dynamics";
    case PluginCategory::Morphing:
        return "morphing";
    case PluginCategory::ImmersedBoundary:
        return "ibm";
    case PluginCategory::Visualization:
        return "visualization";
    case PluginCategory::Zones:
        return "zones";
    case PluginCategory::Utility:
        return "utility";
    case PluginCategory::Core:
        return "core";
    case PluginCategory::BoundaryCondition:
        return "boundary_condition";
    case PluginCategory::_Count:
        return "unknown";
    }
    return "unknown";
}

/// Workspace-relative source directory under src/ that owns the canonical
/// implementation for this category. Empty string for cross-cutting roles
/// (Generic / Writer / Reader / PostProcessor / BoundaryCondition) and for
/// roles that don't have a single owning directory (Utility = utilities/).
[[nodiscard]] inline const char* category_subsystem_path(PluginCategory c) noexcept
{
    switch (c) {
    case PluginCategory::Solver:
        return "src/solver";
    case PluginCategory::Cad:
        return "src/cad";
    case PluginCategory::Meshing:
        return "src/meshing";
    case PluginCategory::Materials:
        return "src/materials";
    case PluginCategory::Turbulence:
        return "src/turbulence";
    case PluginCategory::Multiphase:
        return "src/multiphase";
    case PluginCategory::Combustion:
        return "src/combustion";
    case PluginCategory::Radiation:
        return "src/radiation";
    case PluginCategory::HeatTransfer:
        return "src/heat_transfer";
    case PluginCategory::ElectroMagnetics:
        return "src/emag";
    case PluginCategory::Acoustics:
        return "src/acoustics";
    case PluginCategory::Porous:
        return "src/porous";
    case PluginCategory::Particles:
        return "src/particles";
    case PluginCategory::Adjoint:
        return "src/adjoint";
    case PluginCategory::Optimization:
        return "src/optimization";
    case PluginCategory::Rom:
        return "src/rom";
    case PluginCategory::Io:
        return "src/io";
    case PluginCategory::Scripting:
        return "src/scripting";
    case PluginCategory::Parallel:
        return "src/parallel";
    case PluginCategory::Gpu:
        return "src/gpu";
    case PluginCategory::Gui:
        return "src/gui";
    case PluginCategory::Amr:
        return "src/amr";
    case PluginCategory::Rotating:
        return "src/rotating";
    case PluginCategory::Dynamics:
        return "src/dynamics";
    case PluginCategory::Morphing:
        return "src/morphing";
    case PluginCategory::ImmersedBoundary:
        return "src/ibm";
    case PluginCategory::Visualization:
        return "src/visualization";
    case PluginCategory::Zones:
        return "src/zones";
    case PluginCategory::Utility:
        return "src/utilities";
    case PluginCategory::Core:
        return "src/core";
    // Cross-cutting / no single owner:
    case PluginCategory::PostProcessor:
    case PluginCategory::Writer:
    case PluginCategory::Reader:
    case PluginCategory::Generic:
    case PluginCategory::BoundaryCondition:
    case PluginCategory::_Count:
        return "";
    }
    return "";
}

class CategorisedPlugin : public PluginBase
{
public:
    explicit CategorisedPlugin(PluginCategory cat) : category_(cat) {}
    [[nodiscard]] PluginCategory category() const { return category_; }

private:
    PluginCategory category_;
};

} // namespace simall::plugins

/// Define `CreatePlugin` and `GetPluginAbi` for a class deriving from
/// `simall::plugins::IPlugin`.  Use once per .cpp file at namespace scope.
#define SIMALL_DECLARE_PLUGIN(ClassName)                                                           \
    extern "C" SIMALL_PLUGIN_EXPORT simall::plugins::IPlugin* CreatePlugin()                       \
    {                                                                                              \
        return new ClassName();                                                                    \
    }                                                                                              \
    extern "C" SIMALL_PLUGIN_EXPORT int GetPluginAbi()                                             \
    {                                                                                              \
        return simall::plugins::kPluginAbiVersion;                                                 \
    }
