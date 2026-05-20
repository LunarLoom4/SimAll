// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/Project.hpp
// Phase  : 21 (PROJECT SYSTEM), 22 (FILE FORMATS)
//
// Concrete Project. Owns references — NOT copies — to:
//   - imported CAD shapes  (cad::ShapeHandle)
//   - generated meshes     (meshing::MeshHandle)
//   - physics setup        (solver::PhysicsConfig)
//   - boundary conditions  (solver::BoundaryConditions)
//   - postprocessing state (visualization::PostState)
//
// Persistence uses chunked binary format `.simall` (Section 16 of the
// ultra-detailed spec) — implemented in ProjectSerializer.cpp. The serializer
// writes a manifest + content chunks for partial loading.
// =============================================================================
#pragma once

#include "core/IProject.hpp"
#include "utilities/MathTypes.hpp"

#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace simall::io
{

class ProjectNode; // generic dependency-graph node (Phase 21)

class Project final : public core::IProject
{
public:
    Project();
    explicit Project(std::string path);
    ~Project() override;

    // core::IProject
    void save(const std::string& path) override;
    std::string path() const override { return path_; }
    bool is_dirty() const override { return dirty_.load(); }

    static std::unique_ptr<Project> load(const std::string& path);

    // Dependency graph access (Phase 21 — automatic invalidation)
    ProjectNode* root() noexcept { return root_.get(); }
    void mark_dirty() noexcept { dirty_.store(true); }

    /// Register the io::Project factory with core::Application.
    static void install_factory();

private:
    std::string path_;
    std::atomic<bool> dirty_{false};
    std::unique_ptr<ProjectNode> root_;
    mutable std::mutex mtx_;
};

} // namespace simall::io
