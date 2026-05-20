// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MeshOps.hpp
// Phase  : 23 Pass 11
//
// Mesh-utility operations consolidating the small but recurring tasks that
// canonical CFD pipelines invoke as standalone CLI tools (transformPoints,
// rotateMesh, mirrorMesh, mergeMeshes, stitchMesh, renumberMesh, checkMesh).
// Each operation is exposed as a free function in `simall::meshing::ops`
// so the same routine can be driven by GUI buttons, batch-mode CLIs, and
// in-process pipelines (e.g. SnappyHexMesher's prism-layer stitch).
//
// Operations covered in Pass 11:
//   * Geometry transforms      : transform_points, translate, rotate, mirror
//   * Topology composition     : merge_meshes (raw concat),
//                                 stitch_meshes (merge + node weld +
//                                                duplicate-face dedupe)
//   * Numbering optimisation   : renumber_cells_cuthill_mckee (reverse CM)
//   * Validation               : check_mesh (consolidates MeshQuality
//                                            + topological consistency
//                                            + geometric positivity)
//
// Refinement / splitting (refineMesh, splitMesh) require true remeshing
// and are queued for Pass 11b; this module deliberately ships only the
// operations that can be implemented without changing volume connectivity.
// =============================================================================
#pragma once

#include "MeshStorage.hpp"
#include "MeshQuality.hpp"
#include "utilities/MathTypes.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace simall::meshing::ops {

// =============================================================================
// Geometry transforms
// =============================================================================

/// 3x3 rotation/scale + 3x1 translation packed row-major.  x' = R*x + t.
struct Affine {
    /// Row-major 3x3 matrix (R[0..2] = first row, R[3..5] = second, ...).
    std::array<double, 9> R{1, 0, 0,
                            0, 1, 0,
                            0, 0, 1};
    util::Vec3d           t{0.0, 0.0, 0.0};

    /// Determinant of the rotational part.  Negative values indicate a
    /// reflection and trigger face-orientation reversal in
    /// `transform_points()`.
    [[nodiscard]] double det() const noexcept;
};

/// Apply an affine transform to every node, then recompute geometry.
/// If the determinant of `a.R` is negative the per-face node ordering is
/// reversed so outward normals continue to point outward.
void transform_points(Mesh& m, const Affine& a);

/// Translate by (dx, dy, dz).
void translate(Mesh& m, double dx, double dy, double dz);

/// Rotate by `angleRad` about `axis` (right-hand rule) through `pivot`.
/// `axis` does not need to be unit-length (it is normalised inside).
void rotate(Mesh& m, util::Vec3d axis, double angleRad,
            util::Vec3d pivot = {0, 0, 0});

/// Mirror across the plane defined by point `p0` and unit-length-or-not
/// normal `normal`.  Face orientations are reversed by the underlying
/// `transform_points()` call so the mesh remains right-handed.
void mirror(Mesh& m, util::Vec3d normal, util::Vec3d p0 = {0, 0, 0});

// =============================================================================
// Topology composition
// =============================================================================

/// Append mesh `b` onto mesh `a` with no deduplication and write the
/// result to `out`.  Cell, face, and node indices in `b` are shifted to
/// occupy the contiguous range immediately after `a`'s, boundary faces
/// remain boundary faces, and `out.compute_geometry()` is called at the
/// end.  `a` and `out` may be the same object.
void merge_meshes(const Mesh& a, const Mesh& b, Mesh& out);

struct StitchOptions {
    /// Two nodes are welded when their Euclidean distance is below this
    /// tolerance.  Default scales to ~ machine-eps for unit meshes.
    double weldTolerance = 1.0e-9;
};

struct StitchStats {
    std::size_t weldedNodes      = 0;  // (nA + nB) - (final nNodes)
    std::size_t deduplicatedFaces = 0; // boundary face pairs that became internal
};

/// `merge_meshes` followed by:
///   1. spatial hash weld of coincident nodes (within `weldTolerance`),
///   2. detection of boundary face pairs with identical node sets,
///   3. promotion of each such pair to a single INTERNAL face whose
///      owner/neighbour are the two parent cells (the duplicate face is
///      removed and the surviving cell's face list is rewritten).
///
/// After return, `out` has a fully-connected polyhedral mesh suitable for
/// FV discretisation across the welded interface.
StitchStats stitch_meshes(const Mesh& a, const Mesh& b, Mesh& out,
                          StitchOptions opt = {});

// =============================================================================
// Numbering optimisation
// =============================================================================

struct RenumberStats {
    std::size_t bandwidthBefore = 0;
    std::size_t bandwidthAfter  = 0;
    std::size_t cellsPermuted   = 0;
};

/// Reverse Cuthill-McKee renumbering of cell indices using the face-based
/// adjacency graph (two cells are adjacent iff they share a face).  Faces
/// are kept in place; only `face.owner` and `face.neighbor` are remapped
/// to use the new cell numbering, and the per-cell face/centroid/volume
/// arrays are permuted accordingly.  Bandwidth here is the maximum
/// |owner - neighbor| across interior faces.
RenumberStats renumber_cells_cuthill_mckee(Mesh& m);

// =============================================================================
// Validation
// =============================================================================

struct CheckReport {
    bool                     ok                = true;
    std::size_t              nNodes            = 0;
    std::size_t              nFaces            = 0;
    std::size_t              nCells            = 0;
    std::size_t              nBoundaryFaces    = 0;
    std::size_t              nInteriorFaces    = 0;
    std::size_t              nDegenerateFaces  = 0;  // |area| < eps
    std::size_t              nZeroVolumeCells  = 0;  // |V|    < eps
    std::size_t              nNegativeVolume   = 0;
    std::size_t              nOrphanedFaces    = 0;  // owner out of range
    double                   minFaceArea       = 0.0;
    double                   maxFaceArea       = 0.0;
    double                   minCellVolume     = 0.0;
    double                   maxCellVolume     = 0.0;
    double                   maxNonOrthoDeg    = 0.0;
    double                   maxSkewness       = 0.0;
    double                   maxAspectRatio    = 0.0;
    std::vector<std::string> errors;      // human-readable, populated on !ok

    /// Render a multi-line report (suitable for stdout / log).
    [[nodiscard]] std::string format() const;
};

/// Topological + geometric consistency check.  Calls `MeshQuality::evaluate`
/// internally to obtain skewness / non-orthogonality / aspect-ratio
/// histograms, then layers per-mesh structural checks on top.
CheckReport check_mesh(const Mesh& m);

// =============================================================================
// Sub-mesh extraction & refinement  (Pass 11b)
// =============================================================================

struct SplitStats {
    std::size_t selectedCells = 0;
    std::size_t outputCells   = 0;
    std::size_t outputFaces   = 0;
    std::size_t outputNodes   = 0;
};

/// Extract the sub-mesh consisting of every cell `c` for which
/// `cellMask[c]` is non-zero.  Faces are kept iff at least one of their
/// adjacent cells survives; faces that previously joined a kept cell to a
/// dropped cell become boundary faces in the output (zone 0).  Nodes are
/// renumbered to a contiguous range.  `cellMask.size()` must equal
/// `src.cells().size()`; `&src == &out` is supported.
SplitStats split_mesh(const Mesh& src,
                      const std::vector<std::uint8_t>& cellMask,
                      Mesh& out);

struct RefineStats {
    std::size_t cellsRefined   = 0;
    std::size_t cellsBefore    = 0;
    std::size_t cellsAfter     = 0;
    std::size_t nodesBefore    = 0;
    std::size_t nodesAfter     = 0;
    std::size_t rejectedCells  = 0;   // cells that were not topological hexes
};

/// Uniform 1-to-8 hexahedral refinement.  Every cell of the input mesh
/// MUST be a topological hex (6 quadrilateral faces, 8 unique nodes, each
/// node on exactly 3 faces).  Cells that fail this check are passed
/// through unchanged and counted in `rejectedCells` rather than aborting.
///
/// For each accepted cell the routine generates 12 edge midpoints
/// (deduplicated across cells sharing the edge), 6 face centroids
/// (deduplicated across cells sharing the face), and 1 cell centroid,
/// then emits 8 sub-hexes in standard (i,j,k) octant order.  Topology is
/// detected purely from the per-face node lists - the routine does NOT
/// require any particular face ordering inside the cell.
RefineStats refine_hex(Mesh& m);

// =============================================================================
// Domain decomposition  (Pass 20)
// =============================================================================

struct SubdomainStats {
    std::size_t ownedCells     = 0;   // cells with cellRank[c] == rank
    std::size_t ghostCells     = 0;   // halo cells included via -ghost layer
    std::size_t interfaceFaces = 0;   // faces between owned and non-owned cells
    SplitStats  split;                // wrapped split_mesh result
};

/// Extract the subdomain owned by `rank` from a partitioned global mesh.
/// `cellRank[c]` is the rank to which cell `c` belongs; must have size
/// equal to `src.cells().size()`.  When `includeGhostLayer` is true, the
/// output also includes every non-owned cell that shares a face with at
/// least one owned cell (one-deep halo).  Faces on the partition cut
/// degrade to boundary faces in the output (via the standard `split_mesh`
/// path).  `interfaceFaces` is computed on the input mesh before split.
SubdomainStats extract_subdomain(const Mesh&                       src,
                                 const std::vector<std::int32_t>&  cellRank,
                                 std::int32_t                      rank,
                                 bool                              includeGhostLayer,
                                 Mesh&                             out);

}  // namespace simall::meshing::ops
