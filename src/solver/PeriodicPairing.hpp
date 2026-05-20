// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/PeriodicPairing.hpp
// Phase  : 6.6 — General periodic-face pairing utility supporting both
//          translational AND rotational periodicity (turbomachinery).
//
// Given two boundary zones and a rigid-body transform T such that
//     x_B = T(x_A) = R · (x_A - c) + c + Δ
//
// builds a per-face twin array compatible with SimpleAlgorithm via
// SimpleAlgorithm::set_periodic_pairs(...).
//
// Translation-only: leave axis=any, angle=0.
// Rotation-only   : leave Δ=(0,0,0).
//
// The matching is exact only up to mesh tolerance; nearest-centroid
// matching with a configurable absolute tolerance is used internally.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

#include <vector>

namespace simall::solver
{

struct PeriodicTransform
{
    meshing::ZoneId zoneA = 0;
    meshing::ZoneId zoneB = 0;
    util::Vec3d rotationAxis{0, 0, 1};
    double rotationAngleRad = 0.0;
    util::Vec3d rotationCentre{0, 0, 0};
    util::Vec3d translation{0, 0, 0};
    double tolerance = 1e-6;
};

/// Build a face-twin array of size mesh.faces().size(). Interior face entries
/// are set to -1; periodic boundary face entries point to their twin.
std::vector<int> build_periodic_pairs(const meshing::Mesh& mesh,
                                      const std::vector<PeriodicTransform>& transforms);

} // namespace simall::solver
