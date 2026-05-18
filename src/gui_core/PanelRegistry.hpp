// =============================================================================
// SimAll Beta — GUI Core (Qt-free)
// File   : src/gui_core/PanelRegistry.hpp
//
// Decoupling layer that lets the dock framework instantiate panels by string
// id without depending on the concrete (Qt) widget classes.  GUI plugins,
// scripted UIs, and unit tests can all register factories at startup.
//
// The factory returns a `void*` (since GuiCore must not depend on QWidget),
// and the Qt-side `DockManager` casts it back to QWidget*.  Type-erased
// factories also let us swap "panel" implementations between desktop Qt and
// a future web/imgui front-end without touching the registration site.
// =============================================================================
#pragma once

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::gui_core {

struct PanelInfo {
    std::string id;                 // canonical, e.g. "Console"
    std::string displayName;        // shown in tab/menus
    std::string category;           // "Workspace", "Diagnostics", "Solver"
    std::string iconResource;       // optional QResource path
    bool        singletonPerWindow = true;
};

using PanelFactory = std::function<void*()>;  // returns a QWidget* opaquely

class PanelRegistry {
public:
    static PanelRegistry& instance();

    void register_panel(PanelInfo info, PanelFactory factory);
    [[nodiscard]] bool   has(std::string_view id) const;
    [[nodiscard]] const PanelInfo*    info(std::string_view id) const;
    [[nodiscard]] void* create(std::string_view id) const;
    [[nodiscard]] std::vector<PanelInfo> list() const;
    void clear();

private:
    PanelRegistry() = default;
    struct Entry { PanelInfo info; PanelFactory factory; };
    std::unordered_map<std::string, Entry> entries_;
};

}  // namespace simall::gui_core
