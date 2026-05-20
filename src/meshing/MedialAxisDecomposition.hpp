// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MedialAxisDecomposition.hpp
// Phase  : 6.12 — Medial-axis-based volume decomposition for sweep meshing.
//
// The medial axis of a 3-D solid is the locus of centres of maximal
// inscribed spheres.  Decomposing a complex extrusion body along its
// medial axis yields sweepable sub-regions for HexSweepMesher.
//
// Algorithm sketch (Sheffer-Etzion-Bercovier 1998; Tam-Armstrong 1991):
//   1. From an input STL / triangulated boundary, compute Delaunay tets;
//      circumcentres of Delaunay tets approximate the medial axis.
//   2. Cluster medial points by axial continuity → medial branches.
//   3. For each branch, the swept-disc decomposition yields a sub-region
//      bounded by branch endpoints and projection planes orthogonal to
//      the branch tangent.
//
// The full medial-axis transform is complex.  This module provides:
//   * compute_medial_points(): centroids/radii of inscribed spheres
//     approximated by per-tet circumcenters
//   * decompose(): a greedy axis-walker that returns axially-aligned
//     bounding boxes per sub-region, sufficient downstream for sweeping
//     HexSweepMesher.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <vector>

namespace simall::meshing
{

struct MedialPoint
{
    util::Vec3d position;
    double radius; // inscribed-sphere radius
    int branchId;
};

struct MedialRegion
{
    util::Vec3d minCorner;
    util::Vec3d maxCorner;
    util::Vec3d axisStart;
    util::Vec3d axisEnd;
    int branchId;
};

struct MedialProps
{
    double clusterRadius = 5e-2; // medial-point clustering distance
    double branchMinLength = 0.05;
    std::size_t maxBranches = 1024;
};

class MedialAxisDecomposition
{
public:
    void initialize(MedialProps props);

    /// Build medial points from a tetrahedral mesh (uses cell centroids
    /// and approximate inscribed-sphere radii from face distances).
    std::size_t compute_medial_points(const Mesh& tetMesh);

    /// Decompose into axial regions.  Returns region count.
    std::size_t decompose(std::vector<MedialRegion>& outRegions);

    const std::vector<MedialPoint>& medial_points() const noexcept { return pts_; }
    const MedialProps& props() const noexcept { return p_; }

private:
    MedialProps p_{};
    std::vector<MedialPoint> pts_;
};

} // namespace simall::meshing
