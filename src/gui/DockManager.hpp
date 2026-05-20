// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/DockManager.hpp
//
// Thin manager around QMainWindow's docking API.  Responsibilities:
//   * Add named QDockWidget instances and remember which logical panel id
//     each one carries.
//   * Save / restore a complete layout to a `gui_core::DockPerspective`,
//     including QMainWindow's opaque saveState() blob (base64) so splitter
//     ratios survive.
//   * Maintain a library of named perspectives ("Default", "Mesh", "Solve")
//     and switch between them.
//   * Persist the active perspective on shutdown via QSettings.
// =============================================================================
#pragma once

#include "gui_core/DockPerspective.hpp"

#include <QObject>
#include <QString>
#include <unordered_map>

class QMainWindow;
class QDockWidget;
class QWidget;

namespace simall::gui
{

class DockManager : public QObject
{
    Q_OBJECT
public:
    explicit DockManager(QMainWindow* main, QObject* parent = nullptr);

    // Add `panel` (any QWidget) to a dock area; panel id must match
    // gui_core::PanelRegistry registration so save/restore can round-trip.
    QDockWidget* add_panel(const QString& panelId,
                           const QString& title,
                           QWidget* panel,
                           gui_core::DockArea area = gui_core::DockArea::Right);

    void show_panel(const QString& panelId, bool visible = true);
    [[nodiscard]] bool is_visible(const QString& panelId) const;

    // -- perspectives ---------------------------------------------------------
    [[nodiscard]] gui_core::DockPerspective capture(const std::string& name) const;
    void apply_perspective(const gui_core::DockPerspective& p);

    void save_perspective(const std::string& name);
    bool load_perspective(const std::string& name);
    [[nodiscard]] std::vector<std::string> perspective_names() const;

    // QSettings persistence (organization/app set by main).
    void persist_to_settings(const QString& settingsKey) const;
    void restore_from_settings(const QString& settingsKey);

private:
    QMainWindow* main_;
    std::unordered_map<QString, QDockWidget*> docks_;
    std::unordered_map<std::string, gui_core::DockPerspective> perspectives_;
};

} // namespace simall::gui
