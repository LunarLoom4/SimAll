// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/Application.hpp / .cpp companion
// Phase  : 1.3 (APPLICATION CORE) + Phase 21 (PROJECT SYSTEM)
//
// Top-level orchestrator. Owns the long-lived singletons (logger, event bus,
// command history) and the active Project. Invoked by both the GUI shell
// (applications/simall_beta) and the headless scripting front-end.
// =============================================================================
#pragma once

#include "Command.hpp"
#include "EventBus.hpp"
#include "IProject.hpp"
#include "Logger.hpp"
#include "ServiceLocator.hpp"

#include <memory>
#include <string>

namespace simall::core
{

class Application
{
public:
    Application();
    ~Application();

    Application(const Application&) = delete;
    Application& operator=(const Application&) = delete;

    void bootstrap();
    void shutdown();

    CommandHistory& history() noexcept { return history_; }
    IProject* project() noexcept { return project_.get(); }

    void new_project();
    void open_project(const std::string& path);
    void save_project(const std::string& path);

private:
    CommandHistory history_;
    std::unique_ptr<IProject> project_;
    bool booted_ = false;
};

/// Installed by the io subsystem (or any module that owns the project format).
void register_project_factory(ProjectFactory f);

} // namespace simall::core
