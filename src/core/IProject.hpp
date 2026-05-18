// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/IProject.hpp
// Phase  : 21 (PROJECT SYSTEM — abstract interface)
//
// Pure interface decoupling the core Application from the IO subsystem.
// The concrete Project implementation lives in src/io/Project.* and is
// constructed via core::ProjectFactory (registered with ServiceLocator by
// the io module during bootstrap).
// =============================================================================
#pragma once

#include <functional>
#include <memory>
#include <string>

namespace simall::core {

class IProject {
public:
    virtual ~IProject()                              = default;
    virtual void        save(const std::string&)    = 0;
    virtual std::string path() const                = 0;
    virtual bool        is_dirty() const            = 0;
};

using ProjectFactory =
    std::function<std::unique_ptr<IProject>(const std::string& path /* "" = new */)>;

}  // namespace simall::core
