// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/CartesianMesher.hpp
// Phase  : 5 — production Cartesian hex mesher.
//
// Generates a structured (Nx × Ny × Nz) hex grid over an AABB and emits a
// fully-connected Mesh (owner/neighbour/areas/volumes ready for FV).
//
// Boundary face zones are auto-assigned:
//   zone 1 = xMin (left)    zone 2 = xMax (right)
//   zone 3 = yMin (bottom)  zone 4 = yMax (top)
//   zone 5 = zMin (back)    zone 6 = zMax (front)
// =============================================================================
#pragma once

#include "MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

namespace simall::meshing {

struct CartesianGridSpec {
    util::Vec3d origin{0, 0, 0};
    util::Vec3d extent{1, 1, 1};
    int Nx = 10, Ny = 10, Nz = 10;
};

class CartesianMesher {
public:
    explicit CartesianMesher(CartesianGridSpec s) : spec_(s) {}
    void generate(Mesh& out);
private:
    CartesianGridSpec spec_;
};

}  // namespace simall::meshing
