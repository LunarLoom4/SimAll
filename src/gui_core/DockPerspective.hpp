// =============================================================================
// SimAll Beta — GUI Core (Qt-free)
// File   : src/gui_core/DockPerspective.hpp
//
// Pure-data description of a window layout: where each dockable panel lives,
// whether it's floating, visible, tabbed with neighbours, and (for floating
// docks) its window rect.  Persisted to disk so users can save / restore
// named workspaces ("Mesh", "Solve", "Post-process") and so the application
// can restart with the user's last layout.
//
// Qt has its own opaque QByteArray-based saveState(); that is also stored
// (as a base64 blob inside `qtBlob`) so the precise pixel splitter ratios
// survive.  The structured fields above mirror the same data so non-Qt
// consumers — automation, headless renderers, the W16 Python REPL — can
// query and mutate the layout without spinning up a QApplication.
// =============================================================================
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace simall::gui_core {

enum class DockArea : uint8_t {
    Left = 0,
    Right,
    Top,
    Bottom,
    Floating,
    Central,    // central widget — only one panel may claim this
    Hidden
};

struct DockPlacement {
    std::string panelId;            // matches PanelRegistry key
    DockArea    area    = DockArea::Hidden;
    bool        visible = true;
    bool        floating = false;
    // Floating window rect (ignored when not floating).
    int         x = 0, y = 0, w = 0, h = 0;
    // Tab index inside the dock area; -1 if not tabbed.
    int         tabGroupId = -1;
    int         tabIndex   = 0;
};

struct DockPerspective {
    std::string                 name;        // "Default", "Meshing", "Post"
    std::vector<DockPlacement>  placements;
    std::string                 qtBlob;      // base64(QMainWindow::saveState())

    [[nodiscard]] const DockPlacement* find(std::string_view panelId) const noexcept;
    [[nodiscard]] DockPlacement*       find(std::string_view panelId)       noexcept;
    DockPlacement& upsert(const std::string& panelId);
};

// -- text serialization (INI-flavoured, deterministic) ----------------------
// Format:
//     [perspective]
//     name=Default
//     qt=BASE64
//     [panel]
//     id=Console
//     area=Bottom
//     visible=1
//     floating=0
//     rect=0,0,0,0
//     tab=2,0
//     [panel]
//     ...
[[nodiscard]] std::string  serialize(const DockPerspective& p);
[[nodiscard]] std::optional<DockPerspective> deserialize(std::string_view txt);

[[nodiscard]] std::string  area_to_string(DockArea a);
[[nodiscard]] std::optional<DockArea> area_from_string(std::string_view s);

}  // namespace simall::gui_core
