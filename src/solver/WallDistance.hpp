// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/WallDistance.hpp
// Phase  : 8 — wall-distance field required by RANS turbulence models.
//
// Two strategies, both fully 3-D, polyhedral-mesh compatible:
//
//   1. ExactNearestSurface — for every cell centroid, computes the shortest
//      distance to ANY wall face (closest point on each wall triangle).
//      O(N_cells · N_wallFaces). Accurate; recommended for steady RANS where
//      wall distance is computed once at start-up.
//
//   2. PoissonEikonal — solves Spalding/Sondak's modified Poisson approach:
//      ∇·(∇φ) = -1   subject to φ = 0 on walls, ∂φ/∂n = 0 elsewhere.
//      Then  d = -|∇φ| + √(|∇φ|² + 2φ).
//      Fully implicit, sparse, scalable to billions of cells. Recommended for
//      complex geometries where exact search would be too slow.
//
// Both produce a ScalarField named "wallDistance" suitable for k-ω SST.
// =============================================================================
#pragma once

#include "FieldRegistry.hpp"
#include "Solver.hpp"

#include "meshing/MeshStorage.hpp"

#include <memory>
#include <vector>

namespace simall::solver
{

class IWallDistance
{
public:
    virtual ~IWallDistance() = default;
    virtual void compute(const meshing::Mesh& mesh,
                         const std::vector<BoundarySpec>& bcs,
                         ScalarField& out) = 0;
};

std::unique_ptr<IWallDistance> make_wall_distance_exact();
std::unique_ptr<IWallDistance> make_wall_distance_poisson(ILinearSolver& solver);

} // namespace simall::solver
