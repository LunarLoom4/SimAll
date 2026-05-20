// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/Project.cpp
// =============================================================================
#include "io/Project.hpp"

#include "core/Application.hpp"
#include "core/Logger.hpp"
#include "io/ProjectNode.hpp"

#include <fstream>
#include <stdexcept>

namespace simall::io
{

namespace
{

std::unique_ptr<ProjectNode> build_default_tree()
{
    auto root = std::make_unique<ProjectNode>(NodeKind::Project, "Project");
    auto geom = std::make_unique<ProjectNode>(NodeKind::Geometry, "Geometry");
    auto sel = std::make_unique<ProjectNode>(NodeKind::NamedSelection, "Named Selections");
    auto mat = std::make_unique<ProjectNode>(NodeKind::Material, "Materials");
    auto mesh = std::make_unique<ProjectNode>(NodeKind::Mesh, "Mesh");
    auto phy = std::make_unique<ProjectNode>(NodeKind::Physics, "Physics");
    auto bc = std::make_unique<ProjectNode>(NodeKind::BoundaryCondition, "Boundary Conditions");
    auto init = std::make_unique<ProjectNode>(NodeKind::Initialization, "Initialization");
    auto sol = std::make_unique<ProjectNode>(NodeKind::Solution, "Solution");
    auto res = std::make_unique<ProjectNode>(NodeKind::Result, "Results");

    // wire dataflow dependencies (Phase 21): geometry→mesh→physics→sol→results
    geom->add_downstream(mesh.get());
    mesh->add_downstream(phy.get());
    phy->add_downstream(bc.get());
    bc->add_downstream(init.get());
    init->add_downstream(sol.get());
    sol->add_downstream(res.get());

    root->add_child(std::move(geom));
    root->add_child(std::move(sel));
    root->add_child(std::move(mat));
    root->add_child(std::move(mesh));
    root->add_child(std::move(phy));
    root->add_child(std::move(bc));
    root->add_child(std::move(init));
    root->add_child(std::move(sol));
    root->add_child(std::move(res));
    return root;
}

} // namespace

Project::Project() : root_(build_default_tree()) {}
Project::Project(std::string p) : path_(std::move(p)), root_(build_default_tree()) {}
Project::~Project() = default;

void Project::save(const std::string& p)
{
    std::lock_guard lk(mtx_);
    // Chunked binary serializer goes here (Phase 22 — implementation in
    // ProjectSerializer.cpp once the field schemas stabilize per phase).
    std::ofstream out(p, std::ios::binary);
    if (!out)
        throw std::runtime_error("Cannot open project file for write: " + p);
    const char magic[8] = {'S', 'I', 'M', 'A', 'L', 'L', 0x01, 0x00};
    out.write(magic, sizeof(magic));
    path_ = p;
    dirty_.store(false);
    SIMALL_LOG_INFO("Project", "Saved → ", p);
}

std::unique_ptr<Project> Project::load(const std::string& p)
{
    std::ifstream in(p, std::ios::binary);
    if (!in)
        throw std::runtime_error("Cannot open project file: " + p);
    char magic[8];
    in.read(magic, sizeof(magic));
    if (magic[0] != 'S' || magic[1] != 'I')
        throw std::runtime_error("Bad .simall magic");
    return std::make_unique<Project>(p);
}

void Project::install_factory()
{
    core::register_project_factory([](const std::string& path) -> std::unique_ptr<core::IProject> {
        return path.empty() ? std::make_unique<Project>() : Project::load(path);
    });
}

} // namespace simall::io
