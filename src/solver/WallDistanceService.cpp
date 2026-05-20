// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/WallDistanceService.cpp
// Phase  : 22 Pass 3
// =============================================================================
#include "solver/WallDistanceService.hpp"

#include <utility>

namespace simall::solver {

WallDistanceService::WallDistanceService(const meshing::Mesh& mesh,
                                         const std::vector<BoundarySpec>& bcs,
                                         std::unique_ptr<IWallDistance> backend,
                                         std::string implName)
    : mesh_(mesh), bcs_(bcs),
      backend_(std::move(backend)),
      implName_(std::move(implName)) {}

std::size_t WallDistanceService::prepare() {
    backend_->compute(mesh_, bcs_, field_);
    return field_.size();
}

std::unique_ptr<WallDistanceService>
make_wall_distance_service_exact(const meshing::Mesh& mesh,
                                 const std::vector<BoundarySpec>& bcs) {
    return std::make_unique<WallDistanceService>(
        mesh, bcs, make_wall_distance_exact(), "exact-nearest");
}

std::unique_ptr<WallDistanceService>
make_wall_distance_service_poisson(const meshing::Mesh& mesh,
                                   const std::vector<BoundarySpec>& bcs,
                                   ILinearSolver& linear) {
    return std::make_unique<WallDistanceService>(
        mesh, bcs, make_wall_distance_poisson(linear), "poisson-eikonal");
}

}  // namespace simall::solver
