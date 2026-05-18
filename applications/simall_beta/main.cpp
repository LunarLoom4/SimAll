// =============================================================================
// SimAll Beta - Application Entry Point
// File   : applications/simall_beta/main.cpp
//
// Bootstraps subsystems in the correct order:
//   1. Qt application + GPU surface format
//   2. parallel::initialize (MPI if enabled)
//   3. io::Project::install_factory  → register concrete project type
//   4. core::Application::bootstrap  → create default project, logger, etc.
//   5. gui::MainWindow               → show ribbon + viewport
// Shutdown is in strict reverse order via Application::shutdown + RAII.
// =============================================================================
#include "core/Application.hpp"
#include "core/Logger.hpp"
#include "io/Project.hpp"
#include "parallel/Parallel.hpp"
#include "gui/MainWindow.hpp"

#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>

int main(int argc, char* argv[]) {
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

    QApplication qt(argc, argv);
    qt.setApplicationName("SimAll Beta");
    qt.setOrganizationName("SimAll");
    qt.setApplicationVersion("0.1.0");

    simall::core::Logger::instance().add_file_sink("simall_beta.log");
    SIMALL_LOG_INFO("Boot", "SimAll Beta starting…");

    simall::parallel::initialize(argc, argv);
    simall::io::Project::install_factory();

    simall::core::Application app;
    app.bootstrap();

    simall::gui::MainWindow win(app);
    win.show();

    const int rc = qt.exec();

    app.shutdown();
    simall::parallel::finalize();
    return rc;
}
