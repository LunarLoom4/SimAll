// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/Mesher.hpp
// Phase  : 5.2–5.6 (Surface/Volume/Boundary-layer/Octree/Quality)
//
// Abstract mesher API. Concrete implementations:
//   - SurfaceMesher  (Advancing Front + constrained Delaunay)        Phase 5.2
//   - TetMesher      (Constrained Delaunay tetrahedralization)        Phase 5.3
//   - HexMesher      (Sweep + medial-axis decomposition)              Phase 5.3
//   - PolyMesher     (Voronoi agglomeration)                          Phase 5.3
//   - PrismExtruder  (Boundary-layer extrusion, y+ targeting)         Phase 5.4
//   - OctreeMesher   (Cut-cell / immersed boundary)                   Phase 5.5
//   - QualityImprover(Laplacian + optimization smoothing, edge swaps) Phase 5.6
// =============================================================================
#pragma once

#include "MeshStorage.hpp"
#include "cad/CadKernel.hpp"

#include <memory>

namespace simall::meshing {

struct SurfaceMeshParams {
    double targetEdgeLength = 1.0e-2;
    double curvatureFactor  = 0.1;       // refine when h > κ·R
    double minAngleDeg      = 20.0;      // reject triangles below
    double maxAspectRatio   = 5.0;
};

struct VolumeMeshParams {
    double targetEdgeLength = 1.0e-2;
    double growthRate       = 1.2;
    bool   suppressSlivers  = true;
};

struct BoundaryLayerParams {
    double firstLayerHeight = 1.0e-5;    // m (use y+ targeting in helper)
    double growthRate       = 1.2;
    int    numLayers        = 10;
    double totalThickness   = 0.0;       // auto if 0
};

struct QualityReport {
    double minSkewness, maxSkewness, avgSkewness;
    double minOrthogonality, maxAspectRatio, minJacobian;
    std::size_t badCells;
};

class IMesher {
public:
    virtual ~IMesher() = default;
    virtual void execute(Mesh& output) = 0;
};

class SurfaceMesher : public IMesher {
public:
    SurfaceMesher(const cad::ShapeHandle&, SurfaceMeshParams);
    void execute(Mesh&) override;
private:
    const cad::ShapeHandle& shape_; SurfaceMeshParams p_;
};

class TetMesher : public IMesher {
public:
    explicit TetMesher(VolumeMeshParams p) : p_(p) {}
    void execute(Mesh&) override;
private:
    VolumeMeshParams p_;
};

class PrismExtruder {
public:
    explicit PrismExtruder(BoundaryLayerParams p) : p_(p) {}
    void extrude(Mesh& surface_mesh, Mesh& output_with_prisms);

    /// First-cell estimation from y+ target.
    /// Δy = y+ * ν / uτ ; here uτ from skin-friction estimate Cf*0.5*ρU²/ρ.
    static double first_layer_from_yplus(double yPlusTarget, double nu, double uTau);
private:
    BoundaryLayerParams p_;
};

class QualityAnalyzer {
public:
    static QualityReport analyze(const Mesh&);
};

}  // namespace simall::meshing
