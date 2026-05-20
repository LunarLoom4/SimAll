// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/MainWindow.cpp
// =============================================================================
#include "gui/MainWindow.hpp"

#include "cad/CadKernel.hpp"
#include "core/Command.hpp"
#include "core/EventBus.hpp"
#include "core/Logger.hpp"
#include "gui/ConsolePanel.hpp"
#include "gui/DiagnosticsPanel.hpp"
#include "gui/Dialogs.hpp"
#include "gui/DockManager.hpp"
#include "gui/MeshStatsPanel.hpp"
#include "gui/PropertyEditor.hpp"
#include "gui/PythonConsolePanel.hpp"
#include "gui/ResidualPlot.hpp"
#include "gui/RibbonBar.hpp"
#include "gui/RibbonContent.hpp"
#include "gui/SolverMonitorPanel.hpp"
#include "gui/ThemeManager.hpp"
#include "gui/workbench/WorkbenchSchematicView.hpp"
#include "gui/WorkflowTree.hpp"
#include "io/Hdf5ResultStore.hpp"
#include "io/Project.hpp"
#include "meshing/Mesher.hpp"
#include "meshing/MeshStorage.hpp"
#include "scripting/PyBindings.hpp"
#include "scripting/PythonEngine.hpp"
#include "scripting/Repl.hpp"
#include "visualization/Viewport.hpp"
#include "workbench/CellAdapter.hpp"
#include "workbench/CellAdapterRegistry.hpp"
#include "workbench/ChangeJournal.hpp"
#include "workbench/Commands.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/Workbench.hpp"
#include "workbench/WorkflowEngine.hpp"

#include <QTabWidget>
#include <QtConcurrent>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <functional>
#include <mutex>
#include <QApplication>
#include <QDockWidget>
#include <QFileDialog>
#include <QFuture>
#include <QFutureWatcher>
#include <QInputDialog>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QPointer>
#include <QStatusBar>
#include <QVBoxLayout>
#include <QWidget>

namespace simall::gui
{

// ---------------------------------------------------------------------------
// WorkbenchPipeline -- the shared state that the four built-in
// CellAdapters (cad.import -> mesh.surface -> solver.run -> results.load)
// hand off to each other.  Defined here in the .cpp so MainWindow.hpp
// can stay free of cad / meshing / io heavy headers.  All access is
// mutex-guarded because adapters run on QtConcurrent worker threads.
// ---------------------------------------------------------------------------
struct MainWindow::WorkbenchPipeline
{
    std::mutex mu;
    std::unique_ptr<cad::ShapeHandle> shape;
    std::unique_ptr<cad::TriangleMesh> tri;
    std::unique_ptr<meshing::Mesh> mesh;
    std::vector<double> last_field;
    std::string last_field_name;
    std::atomic<std::uint64_t> next_id{1};
};

MainWindow::MainWindow(core::Application& app, QWidget* parent)
    : QMainWindow(parent)
    , app_(app)
    , cad_kernel_(std::make_unique<cad::CadKernel>())
    , wbPipeline_(std::make_unique<WorkbenchPipeline>())
{
    setWindowTitle("SimAll Beta — Industrial Multiphysics Platform");
    resize(1600, 1000);
    apply_dark_theme();
    build_ribbon();
    viewport_ = new visualization::Viewport(this);
    setCentralWidget(viewport_);
    build_docks();
    build_workbench_dock();
    build_status();
    wire_events();
}

MainWindow::~MainWindow() = default;

void MainWindow::apply_dark_theme()
{
    ThemeManager::instance().apply(Theme::Dark);
}

void MainWindow::build_ribbon()
{
    ribbon_ = new RibbonBar(this);
    auto* container = new QWidget(this);
    auto* lay = new QVBoxLayout(container);
    lay->setContentsMargins(0, 0, 0, 0);
    lay->addWidget(ribbon_);
    setMenuWidget(container);

    // Full 13-tab default population (Section 3.8).
    auto actions = populate_default_ribbon(ribbon_);

    auto bind = [&actions, this](const QString& name, void (MainWindow::*slot)()) {
        if (auto* a = actions[name])
            connect(a, &QAction::triggered, this, slot);
    };
    bind("Open Project", &MainWindow::on_open_project);
    bind("Save Project", &MainWindow::on_save_project);
    bind("Import CAD", &MainWindow::on_import_cad);
    bind("Surface Mesh", &MainWindow::on_generate_surface_mesh);
    bind("Start Solver", &MainWindow::on_run_solver);
    bind("Preferences", &MainWindow::on_preferences);
}

void MainWindow::wire_ribbon_actions()
{ /* reserved for richer bindings */
}

void MainWindow::build_docks()
{
    workflow_ = new WorkflowTree(this);
    workflow_->rebuild_from(dynamic_cast<io::Project*>(app_.project())
                                ? dynamic_cast<io::Project*>(app_.project())->root()
                                : nullptr);
    auto* dWorkflow = new QDockWidget("Project", this);
    dWorkflow->setWidget(workflow_);
    dWorkflow->setMinimumWidth(280);
    addDockWidget(Qt::LeftDockWidgetArea, dWorkflow);

    docks_ = new DockManager(this, this);

    workflow_ = new WorkflowTree(this);
    workflow_->rebuild_from(dynamic_cast<io::Project*>(app_.project())
                                ? dynamic_cast<io::Project*>(app_.project())->root()
                                : nullptr);
    docks_->add_panel("Project", "Project", workflow_, gui_core::DockArea::Left);

    properties_ = new PropertyEditor(this);
    docks_->add_panel("Properties", "Properties", properties_, gui_core::DockArea::Right);

    // Legacy console kept for backwards compatibility — primary console is
    // now `consolePanel_` (richer filter / level UI).
    console_ = new QPlainTextEdit(this);
    console_->setReadOnly(true);
    console_->setStyleSheet("QPlainTextEdit{background:#1B1B1B;color:#D4D4D4;"
                            "font-family:Consolas,monospace;font-size:10pt}");

    consolePanel_ = new ConsolePanel(this);
    diagPanel_ = new DiagnosticsPanel(this);
    meshPanel_ = new MeshStatsPanel(this);
    solverPanel_ = new SolverMonitorPanel(this);
    pythonPanel_ = new PythonConsolePanel(this);
    residuals_ = new ResidualPlot(this);

    // Wire the embedded Python interpreter (Week 16) into the panel.  When
    // SimAll was built without pybind11 the engine returns a friendly
    // "not built in" string, leaving the panel functional but read-only.
    {
        static scripting::PythonEngine s_pyEngine;
        scripting::python::install_output([this](const std::string& s) {
            QString q = QString::fromStdString(s);
            QMetaObject::invokeMethod(
                pythonPanel_, [this, q] { pythonPanel_->post_output(q); }, Qt::QueuedConnection);
        });
        pythonPanel_->set_executor([](const QString& line) -> QString {
            static scripting::Repl s_repl;
            auto r = s_repl.feed(line.toStdString());
            if (r.state == scripting::ReplState::NeedsMore)
                return QStringLiteral("... ");
            if (r.state == scripting::ReplState::SyntaxError)
                return QString::fromStdString(r.message);
            std::string buf = s_repl.take_buffer();
            s_repl.push_history(buf);
            return QString::fromStdString(scripting::python::is_built_in()
                                              ? scripting::python::exec(buf)
                                              : std::string{"Python scripting not built in."});
        });
    }

    auto* dConsole =
        docks_->add_panel("Console", "Console", consolePanel_, gui_core::DockArea::Bottom);
    auto* dDiag =
        docks_->add_panel("Diagnostics", "Diagnostics", diagPanel_, gui_core::DockArea::Bottom);
    auto* dMesh =
        docks_->add_panel("MeshStats", "Mesh stats", meshPanel_, gui_core::DockArea::Bottom);
    auto* dSolver = docks_->add_panel(
        "SolverMonitor", "Solver monitor", solverPanel_, gui_core::DockArea::Bottom);
    auto* dPython =
        docks_->add_panel("PythonConsole", "Python", pythonPanel_, gui_core::DockArea::Bottom);
    auto* dRes =
        docks_->add_panel("Residuals", "Residuals", residuals_, gui_core::DockArea::Bottom);

    tabifyDockWidget(dConsole, dDiag);
    tabifyDockWidget(dDiag, dMesh);
    tabifyDockWidget(dMesh, dSolver);
    tabifyDockWidget(dSolver, dPython);
    tabifyDockWidget(dPython, dRes);

    save_default_perspectives();
}

// ---------------------------------------------------------------------------
// wire_events -- subscribe to the core::EventBus topics that the docks
// care about.  All UI mutations are marshalled onto the GUI thread via
// QMetaObject::invokeMethod(...,Qt::QueuedConnection) so worker threads
// can publish freely.
// ---------------------------------------------------------------------------
void MainWindow::wire_events()
{
    core::EventBus::instance().subscribe<core::events::SolverIteration>(
        [this](const core::events::SolverIteration& e) {
            QMetaObject::invokeMethod(
                this,
                [this, e] {
                    if (residuals_)
                        residuals_->append(e.iteration, e.residualMax);
                    if (solverPanel_) {
                        solverPanel_->append_residual("residual_max", e.iteration, e.residualMax);
                        solverPanel_->set_iteration(e.iteration, 0.0, 0.0);
                    }
                    if (consolePanel_) {
                        consolePanel_->post(
                            2,
                            "Solver",
                            QString("iter %1  r=%2").arg(e.iteration).arg(e.residualMax));
                    }
                },
                Qt::QueuedConnection);
        });

    core::EventBus::instance().subscribe<core::events::CadImported>(
        [this](const core::events::CadImported& e) {
            QMetaObject::invokeMethod(
                this,
                [this, e] {
                    if (console_)
                        console_->appendPlainText(
                            QString("CAD imported: %1").arg(QString::fromStdString(e.path)));
                    if (consolePanel_)
                        consolePanel_->post(
                            2, "CAD", QString("Imported %1").arg(QString::fromStdString(e.path)));
                },
                Qt::QueuedConnection);
        });
}

// ---------------------------------------------------------------------------
// save_default_perspectives -- captures three named DockManager layouts
// (Default / Meshing / Solving) for the View menu's "Perspectives"
// switcher.  Called at the tail of build_docks() so all panels exist.
// ---------------------------------------------------------------------------
void MainWindow::save_default_perspectives()
{
    docks_->save_perspective("Default");
    // Mesh-focused: hide solver/python/residuals.
    docks_->show_panel("PythonConsole", false);
    docks_->show_panel("SolverMonitor", false);
    docks_->show_panel("Residuals", false);
    docks_->save_perspective("Meshing");
    // Solve-focused.
    docks_->show_panel("SolverMonitor", true);
    docks_->show_panel("Residuals", true);
    docks_->show_panel("MeshStats", false);
    docks_->save_perspective("Solving");
    // Restore Default afterwards so the window opens in the canonical layout.
    docks_->load_perspective("Default");
}

// ---------------------------------------------------------------------------
// build_status -- minimal status-bar bootstrap.  Adapter / journal
// status messages are pushed in via build_workbench_dock().
// ---------------------------------------------------------------------------
void MainWindow::build_status()
{
    statusBar()->showMessage("Ready");
}

// ---------------- slots ----------------

void MainWindow::on_open_project()
{
    QString f =
        QFileDialog::getOpenFileName(this, "Open SimAll Project", {}, "SimAll Project (*.simall)");
    if (f.isEmpty())
        return;
    try {
        app_.open_project(f.toStdString());
    } catch (const std::exception& ex) {
        QMessageBox::critical(this, "Open failed", ex.what());
        return;
    }
    if (auto* p = dynamic_cast<io::Project*>(app_.project()))
        workflow_->rebuild_from(p->root());
}

void MainWindow::on_save_project()
{
    QString f =
        QFileDialog::getSaveFileName(this, "Save SimAll Project", {}, "SimAll Project (*.simall)");
    if (f.isEmpty())
        return;
    try {
        app_.save_project(f.toStdString());
    } catch (const std::exception& ex) {
        QMessageBox::critical(this, "Save failed", ex.what());
    }
}

void MainWindow::on_import_cad()
{
    QString f = QFileDialog::getOpenFileName(
        this, "Import CAD", {}, "CAD Files (*.step *.stp *.iges *.igs *.stl)");
    if (f.isEmpty())
        return;

    // Phase 2.3 — heavy work off GUI thread.
    auto* watcher = new QFutureWatcher<cad::TriangleMesh>(this);
    connect(watcher, &QFutureWatcher<cad::TriangleMesh>::finished, this, [this, watcher] {
        auto tri = watcher->result();
        viewport_->show_cad_mesh(tri);
        watcher->deleteLater();
        statusBar()->showMessage(QString("CAD imported: %1 triangles").arg(tri.triangles.size()),
                                 5000);
    });
    const std::string path = f.toStdString();
    QFuture<cad::TriangleMesh> fut = QtConcurrent::run([this, path] {
        auto shape = cad_kernel_->import(path);
        cad_kernel_->heal(shape);
        auto tri = cad_kernel_->tessellate(shape);
        core::EventBus::instance().publish(core::events::CadImported{path, 1});
        return tri;
    });
    watcher->setFuture(fut);
}

void MainWindow::on_generate_surface_mesh()
{
    QMessageBox::information(this,
                             "Mesh",
                             "Surface meshing dispatch: wired through meshing::SurfaceMesher.\n"
                             "Concrete advancing-front kernel hand-off is Phase 5.2.");
}

void MainWindow::on_run_solver()
{
    // Demonstrates non-blocking dispatch (Phase 2.3) — runs the solver in a
    // worker; residual events propagate back to the plot.
    if (solverPanel_) {
        solverPanel_->register_variable("residual_max", QColor(0x00, 0x7A, 0xCC));
        solverPanel_->set_run_state("Running");
        solverPanel_->clear_history();
    }
    QtConcurrent::run([] {
        for (int i = 1; i <= 100; ++i) {
            core::EventBus::instance().publish(
                core::events::SolverIteration{i, std::exp(-i * 0.05)});
        }
        core::EventBus::instance().publish(core::events::SolverConverged{100});
    });
}

void MainWindow::on_preferences()
{
    PreferencesDialog dlg(this);
    if (dlg.exec() == QDialog::Accepted) {
        if (consolePanel_)
            consolePanel_->post(2, "Preferences", "Preferences updated.");
    }
}

void MainWindow::on_about()
{
    AboutDialog dlg(this);
    dlg.exec();
}

// ---------------------------------------------------------------------------
// build_workbench_dock (Phase 22 Pass 22.4)
//
// Constructs the per-window Schematic / StateMachine / WorkflowEngine /
// ChangeJournal trio that backs the schematic panel, then instantiates
// the WorkbenchSchematicView and wires its six gesture signals onto
// workbench::make_*_command factories pushed through
// `app_.history()` (core::CommandHistory).  Going through the command
// history means every panel edit is reversible and Python-scriptable
// via the same Commands.hpp entry points.
//
// The schematic starts empty -- the user populates it via right-click
// "Add cell" or by loading a `.swb.json` (Pass 22.3 io::load_schematic).
// ---------------------------------------------------------------------------
void MainWindow::build_workbench_dock()
{
    wbSchematic_ = std::make_unique<simall::workbench::Schematic>();
    wbState_ = std::make_unique<simall::workbench::StateMachine>(*wbSchematic_);
    wbEngine_ = std::make_unique<simall::workbench::WorkflowEngine>(*wbSchematic_, *wbState_);
    wbJournal_ = std::make_unique<simall::workbench::ChangeJournal>();

    wbView_ = new workbench::WorkbenchSchematicView(this);
    wbView_->set_sources(wbSchematic_.get(), wbState_.get(), wbEngine_.get());

    // Helper: push a command and refresh the view.
    auto push = [this](std::unique_ptr<simall::core::ICommand> cmd) {
        if (!cmd)
            return;
        try {
            app_.history().execute(std::move(cmd));
        } catch (const std::exception& ex) {
            if (consolePanel_) {
                consolePanel_->post(1, "Workbench", QString("command failed: %1").arg(ex.what()));
            }
        }
        wbView_->refresh();
    };

    // Auto-generate a unique-ish label for "Add cell" gestures.
    auto mint_label = [this](simall::workbench::CellKind k) {
        const QString base =
            QString::fromUtf8(std::string(simall::workbench::to_string(k)).c_str());
        const int n = static_cast<int>(wbSchematic_->size()) + 1;
        return QString("%1 %2").arg(base).arg(n).toStdString();
    };

    // Drag-link
    connect(wbView_,
            &workbench::WorkbenchSchematicView::connectRequested,
            this,
            [this, push](simall::workbench::CellLink link) {
                push(simall::workbench::make_connect_command(*wbEngine_, link, wbJournal_.get()));
            });

    // Disconnect link (context menu / Delete key)
    connect(wbView_,
            &workbench::WorkbenchSchematicView::disconnectRequested,
            this,
            [this, push](simall::workbench::CellLink link) {
                push(
                    simall::workbench::make_disconnect_command(*wbEngine_, link, wbJournal_.get()));
            });

    // Add cell (context menu)
    connect(wbView_,
            &workbench::WorkbenchSchematicView::cellAddRequested,
            this,
            [this, push, mint_label](simall::workbench::CellKind k) {
                push(simall::workbench::make_add_cell_command(
                    *wbSchematic_, k, mint_label(k), nullptr, wbJournal_.get()));
            });

    // Remove cell (context menu / Delete key)
    connect(
        wbView_,
        &workbench::WorkbenchSchematicView::cellRemoveRequested,
        this,
        [this, push](simall::workbench::CellId id) {
            push(simall::workbench::make_remove_cell_command(*wbSchematic_, id, wbJournal_.get()));
        });

    // Rename (F2 / context menu)
    connect(wbView_,
            &workbench::WorkbenchSchematicView::cellRenameRequested,
            this,
            [this, push](simall::workbench::CellId id, QString lbl) {
                push(simall::workbench::make_set_cell_label_command(
                    *wbSchematic_, id, lbl.toStdString(), wbJournal_.get()));
            });

    // State override
    connect(wbView_,
            &workbench::WorkbenchSchematicView::cellStateChangeRequested,
            this,
            [this, push](simall::workbench::CellId id, simall::workbench::CellState s) {
                push(simall::workbench::make_set_cell_state_command(
                    *wbSchematic_, id, s, wbJournal_.get()));
            });

    // Mirror journal entries into the status bar / console panel.
    wbJournal_->set_sink([this](const simall::workbench::ChangeEntry& e) {
        const QString arrow = e.direction == simall::workbench::ChangeDirection::Do ? "->" : "<-";
        const QString line = QString("[%1] %2 %3")
                                 .arg(e.sequence)
                                 .arg(arrow)
                                 .arg(QString::fromStdString(e.description));
        QMetaObject::invokeMethod(
            this,
            [this, line] {
                if (statusBar())
                    statusBar()->showMessage(line, 4000);
                if (consolePanel_)
                    consolePanel_->post(2, "Workbench", line);
            },
            Qt::QueuedConnection);
    });

    // -----------------------------------------------------------------
    // Pass 22.5 / 22.5b -- CellAdapter wiring.  We build the registry
    // and register four built-in adapters that drive the existing CAD /
    // Meshing / Solver / Results subsystems for real.  Inter-adapter
    // state (the imported ShapeHandle, the generated Mesh, the loaded
    // field) lives in wbPipeline_ which is mutex-guarded because
    // adapter bodies run on QtConcurrent worker threads.
    //
    // Convention: the per-cell free-text *label* doubles as the
    // adapter's primary input (file path for cad.import / results.load,
    // ignored for mesh.surface / solver.run which consume pipeline
    // state).  A future pass can promote per-cell properties to a
    // richer key/value bag.
    // -----------------------------------------------------------------
    wbAdapters_ = std::make_unique<simall::workbench::CellAdapterRegistry>();

    // -- "cad.import" -----------------------------------------------------
    // Reads cell.label() as a CAD file path, calls CadKernel::import +
    // heal + tessellate, stashes the result in wbPipeline_ for downstream
    // adapters, and refreshes the central viewport on the GUI thread.
    // Publishes core::events::CadImported + CadHealed for any subscriber
    // (the dock panels, the journal sink, scripted listeners).
    wbAdapters_->register_factory("cad.import", [this] {
        return std::make_unique<simall::workbench::FunctionalCellAdapter>(
            simall::workbench::CellKind::Geometry,
            "cad.import",
            [this](simall::workbench::Cell& cell,
                   simall::workbench::ExecutionContext& ctx,
                   std::string& err) {
                const std::string path = std::string(cell.label());
                if (path.empty()) {
                    err = "cad.import: cell label is empty; set it to a CAD file path";
                    return false;
                }
                if (!std::filesystem::exists(path)) {
                    err = "cad.import: file not found: " + path;
                    return false;
                }
                ctx.info(std::string("cad.import: importing ") + path);
                ctx.report_progress(0.10);
                if (ctx.is_cancelled()) {
                    err = "cancelled";
                    return false;
                }

                auto shape = std::make_unique<cad::ShapeHandle>(cad_kernel_->import(path));
                ctx.report_progress(0.40);
                if (ctx.is_cancelled()) {
                    err = "cancelled";
                    return false;
                }

                cad_kernel_->heal(*shape);
                ctx.report_progress(0.70);
                if (ctx.is_cancelled()) {
                    err = "cancelled";
                    return false;
                }

                auto tri = std::make_unique<cad::TriangleMesh>(cad_kernel_->tessellate(*shape));
                ctx.report_progress(0.95);

                const std::uint64_t sid = wbPipeline_->next_id.fetch_add(1);
                {
                    std::scoped_lock lk(wbPipeline_->mu);
                    wbPipeline_->shape = std::move(shape);
                    wbPipeline_->tri = std::move(tri);
                    wbPipeline_->mesh.reset(); // invalidate downstream
                }
                // Marshal viewport refresh onto the GUI thread.  QPointer
                // guards against MainWindow being destroyed mid-flight.
                QPointer<MainWindow> self(this);
                QMetaObject::invokeMethod(
                    this,
                    [self] {
                        if (!self)
                            return;
                        if (self->viewport_ && self->wbPipeline_->tri) {
                            std::scoped_lock lk(self->wbPipeline_->mu);
                            self->viewport_->show_cad_mesh(*self->wbPipeline_->tri);
                        }
                    },
                    Qt::QueuedConnection);

                core::EventBus::instance().publish(core::events::CadImported{path, sid});
                core::EventBus::instance().publish(core::events::CadHealed{sid});
                ctx.report_progress(1.0);
                return true;
            });
    });

    // -- "mesh.surface" ---------------------------------------------------
    // Consumes the ShapeHandle stashed by cad.import, runs the real
    // meshing::SurfaceMesher with default SurfaceMeshParams, stashes the
    // resulting meshing::Mesh in wbPipeline_, and refreshes the viewport.
    wbAdapters_->register_factory("mesh.surface", [this] {
        return std::make_unique<simall::workbench::FunctionalCellAdapter>(
            simall::workbench::CellKind::Mesh,
            "mesh.surface",
            [this](simall::workbench::Cell&,
                   simall::workbench::ExecutionContext& ctx,
                   std::string& err) {
                // Snapshot the shape pointer under the lock; release lock
                // before the (potentially long) meshing run.
                cad::ShapeHandle* shape_ptr = nullptr;
                {
                    std::scoped_lock lk(wbPipeline_->mu);
                    shape_ptr = wbPipeline_->shape.get();
                }
                if (!shape_ptr) {
                    err = "mesh.surface: no CAD shape available; "
                          "run a Geometry cell bound to cad.import first";
                    return false;
                }
                ctx.info("mesh.surface: running advancing-front surface mesher");
                ctx.report_progress(0.10);
                if (ctx.is_cancelled()) {
                    err = "cancelled";
                    return false;
                }

                meshing::SurfaceMeshParams params{};
                meshing::SurfaceMesher mesher(*shape_ptr, params);
                auto mesh = std::make_unique<meshing::Mesh>();
                try {
                    mesher.execute(*mesh);
                } catch (const std::exception& ex) {
                    err = std::string("mesh.surface: ") + ex.what();
                    return false;
                }
                ctx.report_progress(0.90);
                const std::size_t n_cells = mesh->cells().size();

                const std::uint64_t mid = wbPipeline_->next_id.fetch_add(1);
                {
                    std::scoped_lock lk(wbPipeline_->mu);
                    wbPipeline_->mesh = std::move(mesh);
                }
                QPointer<MainWindow> self(this);
                QMetaObject::invokeMethod(
                    this,
                    [self] {
                        if (!self)
                            return;
                        if (self->viewport_ && self->wbPipeline_->mesh) {
                            std::scoped_lock lk(self->wbPipeline_->mu);
                            self->viewport_->show_mesh(*self->wbPipeline_->mesh);
                        }
                    },
                    Qt::QueuedConnection);

                core::EventBus::instance().publish(core::events::MeshGenerated{mid, n_cells});
                ctx.info(std::string("mesh.surface: ") + std::to_string(n_cells) + " cells");
                ctx.report_progress(1.0);
                return true;
            });
    });

    // -- "solver.run" -----------------------------------------------------
    // NOTE -- a fully wired solver::Solver requires a
    // materials::MaterialDatabase + PhysicsConfig + boundary specs that
    // MainWindow does not yet own as members (those will be plumbed by
    // a later phase that adds the Setup-cell UI).  This adapter
    // therefore publishes the *real* SolverStarted / SolverIteration /
    // SolverConverged event stream that the dock panels and journal
    // subscribe to, with a placeholder exponential-decay residual.
    // When the Setup-cell UI lands, swap the inner loop for
    //     solver::Solver s(*mesh, materials);
    //     s.configure(cfg); s.add_boundary(bc); s.initialize();
    //     s.run(maxIter, target);
    // and the rest of the adapter (cancel, progress, event publish)
    // stays unchanged.
    wbAdapters_->register_factory("solver.run", [this] {
        return std::make_unique<simall::workbench::FunctionalCellAdapter>(
            simall::workbench::CellKind::Solution,
            "solver.run",
            [this](simall::workbench::Cell&,
                   simall::workbench::ExecutionContext& ctx,
                   std::string& err) {
                bool have_mesh = false;
                {
                    std::scoped_lock lk(wbPipeline_->mu);
                    have_mesh = static_cast<bool>(wbPipeline_->mesh);
                }
                if (!have_mesh) {
                    err = "solver.run: no mesh available; "
                          "run a Mesh cell bound to mesh.surface first";
                    return false;
                }
                ctx.info("solver.run: starting (placeholder residuals "
                         "until Setup-cell physics config lands)");
                core::EventBus::instance().publish(
                    core::events::SolverStarted{"placeholder-SIMPLE"});

                constexpr int kMaxIter = 100;
                for (int i = 1; i <= kMaxIter; ++i) {
                    if (ctx.is_cancelled()) {
                        err = "cancelled at iteration " + std::to_string(i);
                        core::EventBus::instance().publish(core::events::SolverFailed{err});
                        return false;
                    }
                    core::EventBus::instance().publish(
                        core::events::SolverIteration{i, std::exp(-i * 0.05)});
                    if ((i % 10) == 0) {
                        ctx.report_progress(static_cast<double>(i) / kMaxIter);
                    }
                }
                core::EventBus::instance().publish(core::events::SolverConverged{kMaxIter});
                ctx.report_progress(1.0);
                return true;
            });
    });

    // -- "results.load" ---------------------------------------------------
    // Reads cell.label() as a native SRS results-store path, opens it
    // via io::Hdf5ResultStore, and surfaces the first f64 cell-located
    // dataset it finds onto the viewport via show_scalar_field.
    wbAdapters_->register_factory("results.load", [this] {
        return std::make_unique<simall::workbench::FunctionalCellAdapter>(
            simall::workbench::CellKind::Results,
            "results.load",
            [this](simall::workbench::Cell& cell,
                   simall::workbench::ExecutionContext& ctx,
                   std::string& err) {
                const std::string path = std::string(cell.label());
                if (path.empty()) {
                    err = "results.load: cell label is empty; set it "
                          "to a .srs / .h5 results-store path";
                    return false;
                }
                if (!std::filesystem::exists(path)) {
                    err = "results.load: file not found: " + path;
                    return false;
                }
                ctx.info(std::string("results.load: opening ") + path);
                io::Hdf5ResultStore store;
                if (!store.load_native(path)) {
                    err = "results.load: load_native failed for " + path;
                    return false;
                }
                ctx.report_progress(0.40);
                if (ctx.is_cancelled()) {
                    err = "cancelled";
                    return false;
                }

                // Walk the group tree depth-first for the first F64
                // dataset; record its full path + values.
                std::string found_path;
                io::Dataset* found = nullptr;
                std::function<void(io::Group&, const std::string&)> walk;
                walk = [&](io::Group& g, const std::string& prefix) {
                    if (found)
                        return;
                    for (auto& [name, ds] : g.datasets) {
                        if (ds.dtype == io::StoreDtype::F64) {
                            found_path = prefix + "/" + name;
                            found = &ds;
                            return;
                        }
                    }
                    for (auto& [name, child] : g.children) {
                        if (!child)
                            continue;
                        walk(*child, prefix + "/" + name);
                        if (found)
                            return;
                    }
                };
                walk(store.root(), "");
                if (!found) {
                    err = "results.load: no f64 dataset found in " + path;
                    return false;
                }
                // Decode the F64 payload (LE-packed bytes -> doubles).
                const std::size_t n = static_cast<std::size_t>(found->element_count());
                std::vector<double> values(n);
                if (n > 0) {
                    std::memcpy(values.data(), found->bytes.data(), n * sizeof(double));
                }
                double vmin = 0.0, vmax = 1.0;
                if (!values.empty()) {
                    vmin = *std::min_element(values.begin(), values.end());
                    vmax = *std::max_element(values.begin(), values.end());
                    if (vmax <= vmin)
                        vmax = vmin + 1.0;
                }
                const std::string field_name = found_path;
                {
                    std::scoped_lock lk(wbPipeline_->mu);
                    wbPipeline_->last_field = values;
                    wbPipeline_->last_field_name = field_name;
                }
                QPointer<MainWindow> self(this);
                QMetaObject::invokeMethod(
                    this,
                    [self, field_name, values, vmin, vmax] {
                        if (!self || !self->viewport_)
                            return;
                        self->viewport_->show_scalar_field(field_name, values, vmin, vmax);
                    },
                    Qt::QueuedConnection);

                ctx.info(std::string("results.load: surfaced ") + field_name + " ("
                         + std::to_string(n) + " values)");
                ctx.report_progress(1.0);
                return true;
            });
    });

    // -- Run cell (context menu) -- dispatch on a worker thread.
    connect(
        wbView_,
        &workbench::WorkbenchSchematicView::cellRunRequested,
        this,
        [this](simall::workbench::CellId id) {
            if (consolePanel_) {
                consolePanel_->post(2, "Workbench", QString("Run requested for cell %1").arg(id));
            }
            // Capture the registry + engine pointers by raw pointer --
            // their lifetimes are bounded by MainWindow which owns the
            // QtConcurrent task indirectly through Qt's event loop.
            auto* eng = wbEngine_.get();
            auto* reg = wbAdapters_.get();
            QPointer<MainWindow> guard{this};
            QtConcurrent::run([eng, reg, id, guard] {
                simall::workbench::ExecutionContext ctx;
                ctx.set_log_sink([guard](simall::workbench::AdapterLogLevel lvl,
                                         std::string_view msg) {
                    const int lvl_i = static_cast<int>(lvl);
                    const QString qm = QString::fromUtf8(msg.data(), static_cast<int>(msg.size()));
                    if (!guard)
                        return;
                    QMetaObject::invokeMethod(
                        guard.data(),
                        [g = guard, lvl_i, qm] {
                            if (g && g->consolePanel_)
                                g->consolePanel_->post(lvl_i, "Adapter", qm);
                        },
                        Qt::QueuedConnection);
                });
                ctx.set_progress_sink([guard](double f) {
                    if (!guard)
                        return;
                    const int pct = static_cast<int>(f * 100.0);
                    QMetaObject::invokeMethod(
                        guard.data(),
                        [g = guard, pct] {
                            if (g && g->statusBar())
                                g->statusBar()->showMessage(
                                    QString("Workbench progress: %1%").arg(pct), 1000);
                        },
                        Qt::QueuedConnection);
                });
                eng->refresh_one(id, *reg, ctx);
                // Refresh the view on the GUI thread.
                if (guard) {
                    QMetaObject::invokeMethod(
                        guard.data(),
                        [g = guard] {
                            if (g && g->wbView_)
                                g->wbView_->refresh();
                        },
                        Qt::QueuedConnection);
                }
            });
        });

    // -- Bind adapter (context menu) -- show a picker against the
    // registry keys and push make_set_cell_adapter_command on confirm.
    connect(wbView_,
            &workbench::WorkbenchSchematicView::cellBindAdapterRequested,
            this,
            [this, push](simall::workbench::CellId id) {
                const auto keys = wbAdapters_->keys();
                QStringList items;
                items.reserve(static_cast<int>(keys.size()) + 1);
                items << QString(""); // empty = unbind
                for (const auto& k : keys)
                    items << QString::fromStdString(k);

                simall::workbench::Cell* c = wbSchematic_->cell(id);
                const QString current = c ? QString::fromStdString(c->adapter_id()) : QString{};
                const int curIdx = std::max(0, items.indexOf(current));
                bool ok = false;
                const QString picked = QInputDialog::getItem(this,
                                                             tr("Bind adapter"),
                                                             tr("Adapter id (empty = unbind):"),
                                                             items,
                                                             curIdx,
                                                             /*editable=*/true,
                                                             &ok);
                if (!ok)
                    return;
                push(simall::workbench::make_set_cell_adapter_command(
                    *wbSchematic_, id, picked.toStdString(), wbJournal_.get()));
            });

    // Dock it on the left, tabbed with the existing Project tree if any.
    if (docks_) {
        docks_->add_panel(
            "WorkbenchSchematic", "Workbench schematic", wbView_, gui_core::DockArea::Left);
    } else {
        auto* d = new QDockWidget("Workbench schematic", this);
        d->setWidget(wbView_);
        addDockWidget(Qt::LeftDockWidgetArea, d);
    }
}

} // namespace simall::gui
