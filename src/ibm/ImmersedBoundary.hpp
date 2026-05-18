// =============================================================================
// SimAll Beta - Immersed Boundary Subsystem
// File   : src/ibm/ImmersedBoundary.hpp
// Phase  : 5.7 / 12 — Direct-forcing immersed boundary method (Fadlun 2000,
// Uhlmann 2005) on arbitrary 3-D polyhedral background mesh.
//
// Given a triangulated rigid surface (StlSurface), the IBM:
//   1. Classifies each background cell as FLUID / SOLID / IB (band of one
//      layer of cells immediately adjacent to the surface).
//   2. On IB cells, computes per-cell signed distance to the surface and a
//      "target velocity" U^* (zero for stationary walls, rigid-body velocity
//      for moving bodies).
//   3. Adds a body-force source S_ibm = ρ (U^* - U) / dt to the momentum
//      equation through the standard FieldRegistry vector channel
//      "S_ibm" — the solver picks it up via its existing source-injection
//      machinery.
//
// This avoids any boundary-fitted re-meshing for moving / rotating bodies
// (propellers, valves, biological flows, turbomachinery rotor-stator).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "meshing/StlImporter.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <vector>

namespace simall::ibm {

enum class CellTag : std::int8_t { Fluid = 0, IB = 1, Solid = 2 };

struct RigidBodyKinematics {
    util::Vec3d linearVelocity{0, 0, 0};
    util::Vec3d angularVelocity{0, 0, 0};   // rad/s
    util::Vec3d centreOfRotation{0, 0, 0};
};

class ImmersedBoundary {
public:
    /// Configure with the immersed surface and the background fluid mesh.
    /// Builds the per-cell signed-distance field and the cell tag array.
    void initialize(const meshing::Mesh& fluid,
                    const meshing::StlSurface& surface);

    /// Override / update the rigid body motion. Cheap; can be called every
    /// time step for moving / rotating bodies.
    void set_kinematics(const RigidBodyKinematics& k) { kin_ = k; }

    /// Update the IBM momentum source for the current fluid velocity field
    /// and integration time-step dt. Writes into fields["S_ibm"] (vector).
    void apply_forcing(double dt, solver::FieldRegistry& fields);

    /// Zero velocity & enforce solid-body velocity in pure SOLID cells.
    /// Should be called right after each pressure-correction sub-step.
    void enforce_solid(solver::FieldRegistry& fields);

    const std::vector<CellTag>&              tags() const { return tag_; }
    const util::aligned_vector<double>&      signed_distance() const { return sdf_; }

private:
    util::Vec3d body_velocity(const util::Vec3d& x) const;

    const meshing::Mesh* mesh_ = nullptr;
    RigidBodyKinematics kin_{};
    std::vector<CellTag>            tag_;
    util::aligned_vector<double>    sdf_;        // signed distance, < 0 inside
};

}  // namespace simall::ibm
