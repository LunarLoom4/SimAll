// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/HexSweepMesher.hpp
// Phase  : 6.8 — Hex sweep mesher (Cubit-style sweep) for sweepable volumes.
//
// Given a source surface quadrilateral mesh, a target surface, and a
// linking-side surface, the sweep operation extrudes (or "sweeps") the
// source-quad mesh along the linking-side curve to fill the volume with
// structured hex cells.  This is the dominant technique for sweepable
// extrusion bodies (pipes, ducts, fan blades).
//
// Algorithm:
//   1. Discretise the sweep path P(s) ∈ ℝ³, s ∈ [0,1] in N_layers stations.
//   2. At each station s_k compute a local frame (parallel-transported
//      Frenet-Serret to avoid twist) and map each source node by the
//      frame's affine transform (scale, translate, rotate).
//   3. Build hex cells layer-by-layer:  hex_k = (node_k,*, node_{k+1,*})
//      with 8 nodes (4 source + 4 next-layer images).
//   4. Side faces inherit boundaryZone equal to props.sideZone; the
//      source plane gets sourceZone, target plane gets targetZone.
//
// API consumes the source quad mesh as parallel arrays of node indices
// (4 per quad) and the path as a poly-line of 3-D positions.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace simall::meshing
{

struct SweepProps
{
    std::uint32_t sourceZone = 1;
    std::uint32_t targetZone = 2;
    std::uint32_t sideZone = 3;
    std::size_t nLayers = 10; // path layers (≥2)
    bool scaleAlongPath = false;
    double startScale = 1.0;
    double endScale = 1.0;
};

struct SourceQuadMesh
{
    std::vector<util::Vec3d> nodes;           // node positions
    std::vector<std::array<NodeId, 4>> quads; // CCW node indices
};

class HexSweepMesher
{
public:
    void initialize(SweepProps props);

    /// Sweep `source` along path `pathPts` (≥ 2 control points; linear
    /// interp inside) and emit cells into `outMesh`. Returns hex count.
    std::size_t sweep(const SourceQuadMesh& source,
                      const std::vector<util::Vec3d>& pathPts,
                      Mesh& outMesh);

    const SweepProps& props() const noexcept { return p_; }

private:
    SweepProps p_{};
};

} // namespace simall::meshing
