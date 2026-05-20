// =============================================================================
// SimAll Beta — GUI / Workbench panel
// File   : src/gui/workbench/WorkbenchSchematicView.hpp
// Phase  : 22 Pass 22.4
//
// Qt6 QGraphicsView that renders a workbench::Schematic as an ANSYS-
// Workbench-style node graph: each Cell is a rounded card with a status-
// coloured header strip, left-edge Input ports, right-edge Output ports,
// and a body label.  Links are smooth cubic Beziers between ports.
//
// The widget is *passive*: it does not mutate the Schematic directly.
// Every user gesture (drag-link, right-click action, F2 rename) is
// translated into a Qt signal that the host (MainWindow) maps onto a
// workbench::Commands factory + core::CommandHistory push.  This keeps
// undo/redo universal -- a click in the panel goes through exactly the
// same code path as a Python script that constructs commands manually.
//
// Layout: cells are placed by topological-order column + intra-column
// row.  A future pass can persist user-edited positions, but for now an
// auto-layout is computed every refresh().
// =============================================================================
#pragma once

#include "workbench/CellLink.hpp"
#include "workbench/Workbench.hpp"

#include <optional>
#include <QGraphicsView>
#include <QPointer>
#include <QString>
#include <unordered_map>

class QGraphicsScene;
class QMouseEvent;
class QContextMenuEvent;
class QGraphicsLineItem;

namespace simall::workbench
{
class Schematic;
class StateMachine;
class WorkflowEngine;
} // namespace simall::workbench

namespace simall::gui::workbench
{

class CellNodeItem;
class LinkEdgeItem;

class WorkbenchSchematicView : public QGraphicsView
{
    Q_OBJECT
public:
    explicit WorkbenchSchematicView(QWidget* parent = nullptr);
    ~WorkbenchSchematicView() override;

    // Non-owning bind.  All three pointers must outlive this widget; pass
    // `nullptr` for any of them to clear and stop rendering.
    void set_sources(simall::workbench::Schematic* s,
                     simall::workbench::StateMachine* sm,
                     simall::workbench::WorkflowEngine* eng);

    // Rebuild scene from the bound schematic.  Cheap enough to call on
    // every command (linear in cell + link count).
    void refresh();

signals:
    // Drag-link: user dropped an output port onto an input port.
    void connectRequested(simall::workbench::CellLink link);

    // Context-menu / Delete-key on a selected link.
    void disconnectRequested(simall::workbench::CellLink link);

    // Context-menu "Add cell" on empty canvas; host should mint a fresh
    // unique label.
    void cellAddRequested(simall::workbench::CellKind kind);

    // Context-menu "Remove cell" on a cell card.
    void cellRemoveRequested(simall::workbench::CellId id);

    // F2 / context-menu rename.
    void cellRenameRequested(simall::workbench::CellId id, QString newLabel);

    // Context-menu "Mark UpToDate / Failed / RefreshRequired".
    void cellStateChangeRequested(simall::workbench::CellId id,
                                  simall::workbench::CellState newState);

    // Context-menu "Run cell" -- host runs the adapter bound to the
    // cell via CellAdapterRegistry on a worker thread, then refreshes
    // the view.  (Pass 22.5.)
    void cellRunRequested(simall::workbench::CellId id);

    // Context-menu "Bind adapter..." -- host shows a picker listing
    // registry keys and pushes make_set_cell_adapter_command on confirm.
    // (Pass 22.5.)
    void cellBindAdapterRequested(simall::workbench::CellId id);

protected:
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;
    void contextMenuEvent(QContextMenuEvent*) override;
    void keyPressEvent(QKeyEvent*) override;

private:
    struct PortHit
    {
        simall::workbench::CellId cell{};
        simall::workbench::PortId port{};
        simall::workbench::PortDirection direction{};
        QPointF scene_pos;
    };

    [[nodiscard]] std::optional<PortHit> hit_port(const QPointF& scenePos) const;
    [[nodiscard]] CellNodeItem* hit_cell(const QPointF& scenePos) const;

    void rebuild_scene();
    void layout_cells();

    // Drag-link transient state.
    bool dragging_link_{false};
    PortHit drag_origin_{};
    QGraphicsLineItem* drag_preview_{nullptr};

    QGraphicsScene* scene_{nullptr};
    simall::workbench::Schematic* schematic_{nullptr};
    simall::workbench::StateMachine* state_{nullptr};
    simall::workbench::WorkflowEngine* engine_{nullptr};

    std::unordered_map<simall::workbench::CellId, CellNodeItem*> cell_items_;
};

} // namespace simall::gui::workbench
