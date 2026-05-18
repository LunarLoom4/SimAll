// =============================================================================
// SimAll Beta - Rotating-Frame Subsystem
// File   : src/rotating/MultipleReferenceFrame.hpp
// Phase  : 15 — Multiple Reference Frame (MRF) for turbomachinery, pumps,
// mixers, fans, propellers. Steady-state approximation: cells inside an
// "MRF zone" are solved in a rotating frame; the rest of the domain in the
// inertial frame. Momentum is augmented with Coriolis + centrifugal sources:
//
//     S_MRF = -ρ [ 2 ω × u  +  ω × ( ω × r ) ]                  [N/m³]
//
//   where r = x_cell - x_origin, ω = rotationAxis · rotationRate (rad/s).
//
// Outputs a vector source channel "S_MRF" that the SIMPLE momentum
// assembly already integrates as part of its source aggregation pass.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <vector>

namespace simall::rotating {

struct MrfZone {
    meshing::ZoneId zone     = 0;
    util::Vec3d     origin   {0,0,0};         // rotation axis pivot point
    util::Vec3d     axis     {0,0,1};         // unit rotation axis
    double          omega    = 0.0;           // angular speed [rad/s]
};

class MultipleReferenceFrame {
public:
    void add_zone(MrfZone z) { zones_.push_back(z); }
    void clear()             { zones_.clear(); }

    /// Pre-compute per-cell angular-velocity vector and pivot-relative
    /// position; tag cells inside any MRF zone.
    void initialize(const meshing::Mesh& m,
                    const std::vector<meshing::ZoneId>& cellZone);

    /// Update fields["S_MRF"] from current fields["U"] and density ρ.
    void apply(solver::FieldRegistry& fields, double rho) const;

    /// Rotation-frame velocity field at cell c (utility for post-processing).
    util::Vec3d frame_velocity(std::size_t cell) const;

    std::size_t active_count() const noexcept;

private:
    const meshing::Mesh*      mesh_ = nullptr;
    std::vector<MrfZone>      zones_;
    // Per-cell ω vector (= axis * rate) and r = x_c - origin. Stored SoA.
    util::aligned_vector<double> wx_, wy_, wz_;
    util::aligned_vector<double> rx_, ry_, rz_;
    std::vector<std::uint8_t>    active_;
};

}  // namespace simall::rotating
