// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/Application.cpp
// =============================================================================
#include "core/Application.hpp"

#include "core/IProject.hpp"

namespace simall::core
{

Application::Application() = default;
Application::~Application()
{
    shutdown();
}

void Application::bootstrap()
{
    if (booted_)
        return;
    Logger::instance().set_level(LogLevel::Info);
    SIMALL_LOG_INFO("Application", "SimAll Beta bootstrapping…");
    new_project();
    booted_ = true;
}

void Application::shutdown()
{
    if (!booted_)
        return;
    SIMALL_LOG_INFO("Application", "Shutting down…");
    history_.clear();
    project_.reset();
    Logger::instance().flush();
    booted_ = false;
}

static ProjectFactory& factory_slot()
{
    static ProjectFactory f;
    return f;
}

void register_project_factory(ProjectFactory f)
{
    factory_slot() = std::move(f);
}

void Application::new_project()
{
    if (factory_slot())
        project_ = factory_slot()("");
    EventBus::instance().publish(events::ProjectOpened{""});
}

void Application::open_project(const std::string& path)
{
    if (factory_slot())
        project_ = factory_slot()(path);
    EventBus::instance().publish(events::ProjectOpened{path});
}

void Application::save_project(const std::string& path)
{
    if (project_) {
        project_->save(path);
        EventBus::instance().publish(events::ProjectSaved{path});
    }
}

} // namespace simall::core
