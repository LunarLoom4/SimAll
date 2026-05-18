// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/SlidingInterface.hpp
// Phase  : 6.7 — General Grid Interface (GGI) / sliding-mesh patch coupling.
//
// Couples two non-conformal mesh patches (typically an inner rotating zone
// and an outer stationary zone) by computing area-weighted intersections
// between source faces and target faces in a common projection plane.
//
// For each TARGET face f_T, the interpolation weights to a set of SOURCE
// faces { f_S } are
//
//     w_{T,S} = A( f_T  ∩  Π(f_S) )   /   A(f_T)
//
// where Π(f_S) is the source face projected through the mean interface
// normal onto the target plane (handles both translational and rotational
// sliding by applying the per-side transform before projection).
//
// At assembly time the SIMPLE momentum/pressure solver replaces the
// missing neighbour-cell contributions across the interface with weighted
// linear combinations of cell-centred values on the opposite side. This
// header exposes:
//
//   - A SlidingPair table (target face f_T, source faces f_S[], weights w[])
//   - A build_pairs() routine that constructs the table from the two patch
//     descriptors and an optional rigid-body transform (rotation + offset).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <vector>

namespace simall::solver {

struct SlidingPatch {
    meshing::ZoneId boundaryZone = 0;   // face zone forming this side
    // Rigid transform applied to take this patch into the COMMON frame.
    util::Vec3d rotAxis  {0, 0, 1};
    util::Vec3d rotOrigin{0, 0, 0};
    double      rotAngle = 0.0;          // [rad]
    util::Vec3d translation{0, 0, 0};
};

struct SlidingWeight {
    meshing::FaceId sourceFace;
    double          w;
};

struct SlidingPairTable {
    // For each target face, an offset-range into entries[].
    std::vector<meshing::FaceId>  targetFaces;
    std::vector<int>              offsets;        // size = targetFaces.size()+1
    std::vector<SlidingWeight>    entries;
};

class SlidingInterface {
public:
    /// Build the interpolation table between two non-conformal patches.
    /// `tolerance` is in mesh length units; intersections below it are
    /// dropped.
    static SlidingPairTable build_pairs(const meshing::Mesh& mesh,
                                        const SlidingPatch& target,
                                        const SlidingPatch& source,
                                        double tolerance = 1e-10);

    /// Sample a cell-centred scalar field across the interface: returns the
    /// interpolated face value on the target side from source neighbours.
    static double interpolate_scalar(const SlidingPairTable& table,
                                     std::size_t targetIndex,
                                     const meshing::Mesh& mesh,
                                     const std::vector<double>& cellPhi);
};

}  // namespace simall::solver
