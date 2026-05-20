// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/MainWindow.hpp
// Phase  : 2.1 (MAIN WINDOW), 2.2 (PROPERTY EDITOR), 2.3 (THREADING MODEL)
//
// QMainWindow shell: top ribbon, central viewport, left workflow tree, right
// property editor, bottom console/residual/diagnostics docks. Layout matches
// Section 3.2 of the ultra-detailed spec exactly. Dark theme is applied via
// QPalette + QSS using the colours from Section 3.6.
// =============================================================================
#pragma once

#include "core/Application.hpp"
#include <QMainWindow>
#include <memory>

class QTabWidget;
class QTreeView;
class QDockWidget;
class QPlainTextEdit;

namespace simall::visualization { class Viewport; }
namespace simall::cad           { class CadKernel; class ShapeHandle; struct TriangleMesh; }
namespace simall::meshing       { class Mesh; }
namespace simall::workbench {
    class Schematic;
    class StateMachine;
    class WorkflowEngine;
    class ChangeJournal;
    class CellAdapterRegistry;
}

namespace simall::gui {

class WorkflowTree;
class PropertyEditor;
class RibbonBar;
class ResidualPlot;
class DockManager;
class ConsolePanel;
class DiagnosticsPanel;
class MeshStatsPanel;
class SolverMonitorPanel;
class PythonConsolePanel;

namespace workbench { class WorkbenchSchematicView; }

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(core::Application& app, QWidget* parent = nullptr);
    ~MainWindow() override;

private slots:
    void on_import_cad();
    void on_open_project();
    void on_save_project();
    void on_generate_surface_mesh();
    void on_run_solver();
    void on_preferences();
    void on_about();

private:
    void build_ribbon();
    void build_docks();
    void build_workbench_dock();
    void build_status();
    void apply_dark_theme();
    void wire_events();
    void wire_ribbon_actions();
    void save_default_perspectives();

    core::Application&                  app_;
    RibbonBar*                          ribbon_      = nullptr;
    visualization::Viewport*            viewport_    = nullptr;
    WorkflowTree*                       workflow_    = nullptr;
    PropertyEditor*                     properties_  = nullptr;
    QPlainTextEdit*                     console_     = nullptr;
    ResidualPlot*                       residuals_   = nullptr;

    DockManager*                        docks_       = nullptr;
    ConsolePanel*                       consolePanel_= nullptr;
    DiagnosticsPanel*                   diagPanel_   = nullptr;
    MeshStatsPanel*                     meshPanel_   = nullptr;
    SolverMonitorPanel*                 solverPanel_ = nullptr;
    PythonConsolePanel*                 pythonPanel_ = nullptr;
    std::unique_ptr<cad::CadKernel>     cad_kernel_;

    // Workbench (Phase 22 Pass 22.4).  Lazily constructed in
    // build_workbench_dock(); destroyed in reverse order with the
    // MainWindow.  WorkbenchSchematicView is owned by its QDockWidget,
    // so we keep only a raw pointer.
    std::unique_ptr<simall::workbench::Schematic>      wbSchematic_;
    std::unique_ptr<simall::workbench::StateMachine>   wbState_;
    std::unique_ptr<simall::workbench::WorkflowEngine> wbEngine_;
    std::unique_ptr<simall::workbench::ChangeJournal>  wbJournal_;
    std::unique_ptr<simall::workbench::CellAdapterRegistry> wbAdapters_;
    workbench::WorkbenchSchematicView*                 wbView_ = nullptr;

    // Pass 22.5b -- shared pipeline state threaded between the four
    // built-in CellAdapters (cad.import -> mesh.surface -> solver.run ->
    // results.load).  Each adapter may run on a worker thread so all
    // access is mutex-guarded.  Pointers are used (rather than
    // std::optional) because cad::ShapeHandle / meshing::Mesh are move-
    // only / forward-declared here, and unique_ptr is the cleanest way
    // to keep the heavy headers out of MainWindow.hpp.
    struct WorkbenchPipeline;
    std::unique_ptr<WorkbenchPipeline>                 wbPipeline_;
};

}  // namespace simall::gui
