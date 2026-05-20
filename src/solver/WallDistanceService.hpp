// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/WallDistanceService.hpp
// Phase  : 22 Pass 3 — adapter promoting solver::IWallDistance back-ends to
//                      the public core::IWallDistanceService interface.
//
// Architecture
// ------------
// Two paths into the same field:
//
//   (a) Solver-internal (legacy): call make_wall_distance_exact() or
//       make_wall_distance_poisson() and invoke compute(mesh, bcs, out)
//       directly. Useful when the caller already owns the ScalarField.
//
//   (b) Service (new):  construct WallDistanceService(mesh, bcs, backend)
//       once, register it in core::ServiceLocator, let any consumer call
//       prepare()/data() through the core::IWallDistanceService interface.
//
// The service owns the ScalarField storage and forwards prepare() to the
// wrapped solver back-end. Re-calling prepare() recomputes the field
// (useful after mesh deformation / AMR refinement).
// =============================================================================
#pragma once

#include "core/IWallDistanceService.hpp"
#include "solver/Solver.hpp"
#include "solver/WallDistance.hpp"
#include "meshing/MeshStorage.hpp"

#include <memory>
#include <string>
#include <vector>

namespace simall::solver {

class WallDistanceService final : public core::IWallDistanceService {
public:
    /// Takes ownership of `backend` (an ExactNearest or Poisson instance).
    /// Keeps references to `mesh` and `bcs` — caller must outlive the
    /// service.
    WallDistanceService(const meshing::Mesh& mesh,
                        const std::vector<BoundarySpec>& bcs,
                        std::unique_ptr<IWallDistance> backend,
                        std::string implName);

    std::size_t   prepare() override;
    const double* data() const noexcept override { return field_.data(); }
    std::size_t   size() const noexcept override { return field_.size(); }
    const char*   implementation_name() const noexcept override {
        return implName_.c_str();
    }

    /// Direct access for solver-internal consumers that need the typed
    /// ScalarField (e.g. to alias into a FieldRegistry slot).
    const ScalarField& field() const noexcept { return field_; }

private:
    const meshing::Mesh&                   mesh_;
    const std::vector<BoundarySpec>&       bcs_;
    std::unique_ptr<IWallDistance>         backend_;
    ScalarField                            field_;
    std::string                            implName_;
};

/// Convenience factories.
std::unique_ptr<WallDistanceService>
make_wall_distance_service_exact(const meshing::Mesh& mesh,
                                 const std::vector<BoundarySpec>& bcs);

std::unique_ptr<WallDistanceService>
make_wall_distance_service_poisson(const meshing::Mesh& mesh,
                                   const std::vector<BoundarySpec>& bcs,
                                   ILinearSolver& linear);

}  // namespace simall::solver
