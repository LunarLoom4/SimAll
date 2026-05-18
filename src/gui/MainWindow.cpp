// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/MainWindow.cpp
// =============================================================================
#include "gui/MainWindow.hpp"
#include "gui/RibbonBar.hpp"
#include "gui/RibbonContent.hpp"
#include "gui/WorkflowTree.hpp"
#include "gui/PropertyEditor.hpp"
#include "gui/ResidualPlot.hpp"
#include "gui/ThemeManager.hpp"
#include "gui/DockManager.hpp"
#include "gui/ConsolePanel.hpp"
#include "gui/DiagnosticsPanel.hpp"
#include "gui/MeshStatsPanel.hpp"
#include "gui/SolverMonitorPanel.hpp"
#include "gui/PythonConsolePanel.hpp"
#include "gui/Dialogs.hpp"

#include "scripting/PyBindings.hpp"
#include "scripting/PythonEngine.hpp"
#include "scripting/Repl.hpp"

#include "core/EventBus.hpp"
#include "core/Logger.hpp"
#include "io/Project.hpp"
#include "cad/CadKernel.hpp"
#include "visualization/Viewport.hpp"

#include <QApplication>
#include <QDockWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QStatusBar>
#include <QTabWidget>
#include <QVBoxLayout>
#include <QWidget>
#include <QtConcurrent>
#include <QFuture>
#include <QFutureWatcher>

namespace simall::gui {

MainWindow::MainWindow(core::Application& app, QWidget* parent)
    : QMainWindow(parent), app_(app),
      cad_kernel_(std::make_unique<cad::CadKernel>())
{
    setWindowTitle("SimAll Beta — Industrial Multiphysics Platform");
    resize(1600, 1000);
    apply_dark_theme();
    build_ribbon();
    viewport_ = new visualization::Viewport(this);
    setCentralWidget(viewport_);
    build_docks();
    build_status();
    wire_events();
}

MainWindow::~MainWindow() = default;

void MainWindow::apply_dark_theme() {
    ThemeManager::instance().apply(Theme::Dark);
}

void MainWindow::build_ribbon() {
    ribbon_ = new RibbonBar(this);
    auto* container = new QWidget(this);
    auto* lay = new QVBoxLayout(container);
    lay->setContentsMargins(0,0,0,0);
    lay->addWidget(ribbon_);
    setMenuWidget(container);

    // Full 13-tab default population (Section 3.8).
    auto actions = populate_default_ribbon(ribbon_);

    auto bind = [&actions, this](const QString& name, void (MainWindow::*slot)()){
        if (auto* a = actions[name]) connect(a, &QAction::triggered, this, slot);
    };
    bind("Open Project",      &MainWindow::on_open_project);
    bind("Save Project",      &MainWindow::on_save_project);
    bind("Import CAD",        &MainWindow::on_import_cad);
    bind("Surface Mesh",      &MainWindow::on_generate_surface_mesh);
    bind("Start Solver",      &MainWindow::on_run_solver);
    bind("Preferences",       &MainWindow::on_preferences);
}

void MainWindow::wire_ribbon_actions() { /* reserved for richer bindings */ }

void MainWindow::build_docks() {
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
    docks_->add_panel("Project",   "Project",   workflow_,   gui_core::DockArea::Left);

    properties_ = new PropertyEditor(this);
    docks_->add_panel("Properties","Properties",properties_, gui_core::DockArea::Right);

    // Legacy console kept for backwards compatibility — primary console is
    // now `consolePanel_` (richer filter / level UI).
    console_ = new QPlainTextEdit(this);
    console_->setReadOnly(true);
    console_->setStyleSheet(
        "QPlainTextEdit{background:#1B1B1B;color:#D4D4D4;"
        "font-family:Consolas,monospace;font-size:10pt}");

    consolePanel_ = new ConsolePanel(this);
    diagPanel_    = new DiagnosticsPanel(this);
    meshPanel_    = new MeshStatsPanel(this);
    solverPanel_  = new SolverMonitorPanel(this);
    pythonPanel_  = new PythonConsolePanel(this);
    residuals_    = new ResidualPlot(this);

    // Wire the embedded Python interpreter (Week 16) into the panel.  When
    // SimAll was built without pybind11 the engine returns a friendly
    // "not built in" string, leaving the panel functional but read-only.
    {
        static scripting::PythonEngine s_pyEngine;
        scripting::python::install_output(
            [this](const std::string& s) {
                QString q = QString::fromStdString(s);
                QMetaObject::invokeMethod(pythonPanel_, [this, q]{
                    pythonPanel_->post_output(q);
                }, Qt::QueuedConnection);
            });
        pythonPanel_->set_executor([](const QString& line) -> QString {
            static scripting::Repl s_repl;
            auto r = s_repl.feed(line.toStdString());
            if (r.state == scripting::ReplState::NeedsMore) return QStringLiteral("... ");
            if (r.state == scripting::ReplState::SyntaxError)
                return QString::fromStdString(r.message);
            std::string buf = s_repl.take_buffer();
            s_repl.push_history(buf);
            return QString::fromStdString(
                scripting::python::is_built_in()
                    ? scripting::python::exec(buf)
                    : std::string{"Python scripting not built in."});
        });
    }

    auto* dConsole = docks_->add_panel("Console",      "Console",       consolePanel_, gui_core::DockArea::Bottom);
    auto* dDiag    = docks_->add_panel("Diagnostics",  "Diagnostics",   diagPanel_,    gui_core::DockArea::Bottom);
    auto* dMesh    = docks_->add_panel("MeshStats",    "Mesh stats",    meshPanel_,    gui_core::DockArea::Bottom);
    auto* dSolver  = docks_->add_panel("SolverMonitor","Solver monitor",solverPanel_,  gui_core::DockArea::Bottom);
    auto* dPython  = docks_->add_panel("PythonConsole","Python",        pythonPanel_,  gui_core::DockArea::Bottom);
    auto* dRes     = docks_->add_panel("Residuals",    "Residuals",     residuals_,    gui_core::DockArea::Bottom);

    tabifyDockWidget(dConsole, dDiag);
    tabifyDockWidget(dDiag,    dMesh);
    tabifyDo    if (solverPanel_) {
                    solverPanel_->append_residual("residual_max", e.iteration, e.residualMax);
                    solverPanel_->set_iteration(e.iteration, 0.0, 0.0);
                }
                if (consolePanel_) {
                    consolePanel_->post(2, "Solver",
                        QString("iter %1  r=%2").arg(e.iteration).arg(e.residualMax));
                }
            }, Qt::QueuedConnection);
        });

    core::EventBus::instance().subscribe<core::events::CadImported>(
        [this](const core::events::CadImported& e){
            QMetaObject::invokeMethod(this, [this, e]{
                if (console_) console_->appendPlainText(
                    QString("CAD imported: %1").arg(QString::fromStdString(e.path)));
                if (consolePanel_) consolePanel_->post(
                    2, "CAD",
                    QString("I

void MainWindow::save_default_perspectives() {
    docks_->save_perspective("Default");
    // Mesh-focused: hide solver/python/residuals, show project + properties + mesh stats.
    docks_->show_panel("PythonConsole", false);
    docks_->show_panel("SolverMonitor", false);
    docks_->show_panel("Residuals",     false);
    docks_->save_perspective("Meshing");
    // Solve-focused.
    docks_->show_panel("SolverMonitor", true);
    docks_->show_panel("Residuals",     true);
    docks_->show_panel("MeshStats",     false);
    docks_->save_perspective("Solving");
    // Restore Default afterwards.
    docks_->load_perspective("Default"ribe<core::events::SolverIteration>(
        [this](const core::events::SolverIteration& e){
            QMetaObject::invokeMethod(residuals_, [this, e]{
                residuals_->append(e.iteration, e.residualMax);
            }, Qt::QueuedConnection);
        });

    core::EventBus::instance().subscribe<core::events::CadImported>(
        [this](const core::events::CadImported& e){
            QMetaObject::invokeMethod(this, [this, e]{
                console_->appendPlainText(QString("CAD imported: %1").arg(QString::fromStdString(e.path)));
            }, Qt::QueuedConnection);
        });
}

// ---------------- slots ----------------

void MainWindow::on_open_project() {
    QString f = QFileDialog::getOpenFileName(this, "Open SimAll Project",
                                             {}, "SimAll Project (*.simall)");
    if (f.isEmpty()) return;
    try { app_.open_project(f.toStdString()); } catch (const std::exception& ex) {
        QMessageBox::critical(this, "Open failed", ex.what()); return;
    }
    if (auto* p = dynamic_cast<io::Project*>(app_.project()))
        workflow_->rebuild_from(p->root());
}

void MainWindow::on_save_project() {
    QString f = QFileDialog::getSaveFileName(this, "Save SimAll Project",
                                             {}, "SimAll Project (*.simall)");
    if (f.isEmpty()) return;
    try { app_.save_project(f.toStdString()); }
    catch (const std::exception& ex) {
        QMessageBox::critical(this, "Save failed", ex.what());
    }
}

void MainWindow::on_import_cad() {
    QString f = QFileDialog::getOpenFileName(this, "Import CAD",
        {}, "CAD Files (*.step *.stp *.iges *.igs *.stl)");
    if (f.isEmpty()) return;

    // Phase 2.3 — heavy work off GUI thread.
    auto* watcher = new QFutureWatcher<cad::TriangleMesh>(this);
    connect(watcher, &QFutureWatcher<cad::TriangleMesh>::finished, this,
        [this, watcher]{
            auto tri = watcher->result();
            viewport_->show_cad_mesh(tri);
            watcher->deleteLater();
            statusBar()->showMessage(
                QString("CAD imported: %1 triangles").arg(tri.triangles.size()), 5000);
        });
    const std::string path = f.toStdString();
    QFuture<cad::TriangleMesh> fut = QtConcurrent::run([this, path]{
        auto shape = cad_kernel_->import(path);
        cad_kernel_->heal(shape);
        auto tri   = cad_kernel_->tessellate(shape);
        core::EventBus::instance().publish(core::events::CadImported{path, 1});
        return tri;
    });
    watcher->setFuture(fut);
}

void MainWindow::on_generate_surface_mesh() {
    QMessageBox::information(this, "Mesh",
        "Surface meshing dispatch: wired through meshing::SurfaceMesher.\n"
        "Concrete advancing-front kernel hand-off is Phase 5.2.");
}

void MainWindow::on_run_solver() {
    // Demonstrates non-blocking dispatch (Phase 2.3) — runs the solver in a
    // worker; residual events propagate back to the plot.
    if (solverPanel_) {
        solverPanel_->register_variable("residual_max", QColor(0x00, 0x7A, 0xCC));
        solverPanel_->set_run_state("Running");
        solverPanel_->clear_history();
    }
    QtConcurrent::run([]{
        for (int i = 1; i <= 100; ++i) {
            core::EventBus::instance().publish(
                core::events::SolverIteration{i, std::exp(-i*0.05)});
        }
        core::EventBus::instance().publish(core::events::SolverConverged{100});
    });
}

void MainWindow::on_preferences() {
    PreferencesDialog dlg(this);
    if (dlg.exec() == QDialog::Accepted) {
        if (consolePanel_)
            consolePanel_->post(2, "Preferences", "Preferences updated.");
    }
}

void MainWindow::on_about() {
    AboutDialog dlg(this);
    dlg.exec();
}

}  // namespace simall::gui
