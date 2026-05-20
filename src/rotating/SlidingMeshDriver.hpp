// =============================================================================
// SimAll Beta - Rotating Subsystem
// File   : src/rotating/SlidingMeshDriver.hpp
// Week   : 18
//
// Sliding-mesh driver for rotating machinery (fans, pumps, compressors).
// The mesh is split into a stationary outer domain and a rotating inner
// domain; the two communicate across a cylindrical or planar *sliding
// interface*.  Each time step the rotor is incremented by Ω · dt; the
// interface re-builds its donor/receiver weighting and the solver gets
// a refreshed map of face-to-face flux exchange.
//
// This module:
//   * stores the cylindrical / planar interface geometry,
//   * tracks the cumulative rotor angle,
//   * computes nearest-neighbour donor weights between rotor-side and
//     stator-side face centroids (cell-centred FV friendly),
//   * exposes `advance_rotor(dt)` and `interface_map()` so the solver
//     can call it once per time step.
// =============================================================================
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace simall::rotating
{

using FaceIdx = std::uint32_t;

struct InterfaceFace
{
    FaceIdx id;
    std::array<double, 3> centroid;
    std::array<double, 3> normal;
    double area = 0.0;
};

struct DonorWeight
{
    FaceIdx receiver;
    FaceIdx donor;
    double w; // ∈ [0,1]
};

struct SlidingInterfaceMap
{
    std::vector<DonorWeight> entries;
};

class SlidingMeshDriver
{
public:
    SlidingMeshDriver(std::array<double, 3> axisOrigin,
                      std::array<double, 3> axisDir,
                      double omega);

    void set_rotor_faces(std::vector<InterfaceFace> faces);
    void set_stator_faces(std::vector<InterfaceFace> faces);

    /// Advance the rotor by Ω·dt; returns the new cumulative angle (rad).
    double advance_rotor(double dt);

    [[nodiscard]] double cumulative_angle() const noexcept { return angle_; }

    /// Build a donor/receiver map.  Each rotor face is matched against the
    /// k=3 closest stator faces (by Euclidean distance of centroids in the
    /// *current* rotor frame) and the weights are inverse-distance.
    [[nodiscard]] SlidingInterfaceMap build_interface_map() const;

private:
    std::array<double, 3> axisOrigin_;
    std::array<double, 3> axisDir_;
    double omega_;
    double angle_ = 0.0;
    std::vector<InterfaceFace> rotor_;
    std::vector<InterfaceFace> stator_;

    [[nodiscard]] InterfaceFace rotated_rotor_face(const InterfaceFace& f) const;
};

} // namespace simall::rotating
