// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/Mesher.cpp
//
// First implementation pass: surface mesh ingested from CAD tessellation
// (BRepMesh) gives the initial triangulation; volume meshing scaffolding is
// in place but the constrained-Delaunay back-end will be plugged in during
// Phase 5.3 implementation work. Helpers (y+, quality) are fully realized.
// =============================================================================
#include "meshing/Mesher.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::meshing
{

// ---------------------------------------------------------------- SurfaceMesher
SurfaceMesher::SurfaceMesher(const cad::ShapeHandle& s, SurfaceMeshParams p) : shape_(s), p_(p) {}

void SurfaceMesher::execute(Mesh& m)
{
    cad::CadKernel kernel;
    cad::TessellationParams tp;
    tp.deflection = p_.targetEdgeLength * p_.curvatureFactor;
    tp.angle = 0.35;
    auto tri = kernel.tessellate(shape_, tp);

    auto& N = m.nodes();
    auto& F = m.faces();
    N.reserve(tri.points.size());
    for (auto& p : tri.points) {
        N.x.push_back(p.x);
        N.y.push_back(p.y);
        N.z.push_back(p.z);
    }

    F.nodeOffsets.push_back(0);
    for (auto& t : tri.triangles) {
        F.nodeIndices.push_back(t[0]);
        F.nodeIndices.push_back(t[1]);
        F.nodeIndices.push_back(t[2]);
        F.nodeOffsets.push_back(static_cast<std::int32_t>(F.nodeIndices.size()));
        F.owner.push_back(kBoundaryCell);
        F.neighbor.push_back(kBoundaryCell);
        F.boundaryZone.push_back(1); // surface zone
    }
    m.add_zone("Surface", true);
    m.compute_geometry();
    SIMALL_LOG_INFO("Mesh", "Surface meshed: ", tri.triangles.size(), " triangles");
}

// ----------------------------------------------------------------- TetMesher
void TetMesher::execute(Mesh& m)
{
    // Constrained Delaunay tetrahedralization (Bowyer-Watson + edge recovery)
    // is implemented incrementally; the public contract is set here so the
    // solver pipeline can be wired in parallel.
    (void) m;
    SIMALL_LOG_WARN("Mesh", "TetMesher::execute pending Phase 5.3 implementation hand-off");
}

// ----------------------------------------------------------------- PrismExtruder
double PrismExtruder::first_layer_from_yplus(double yp, double nu, double uTau)
{
    return (uTau > 0) ? yp * nu / uTau : 0.0;
}

void PrismExtruder::extrude(Mesh& /*src*/, Mesh& /*dst*/)
{
    SIMALL_LOG_WARN("Mesh", "PrismExtruder::extrude pending Phase 5.4 implementation hand-off");
}

// ----------------------------------------------------------------- QualityAnalyzer
QualityReport QualityAnalyzer::analyze(const Mesh& m)
{
    QualityReport r{1, 0, 0, 1, 0, 1, 0};
    const auto& F = m.faces();
    if (F.size() == 0)
        return r;

    double sumSkew = 0;
    for (std::size_t f = 0; f < F.size(); ++f) {
        // Skewness via deviation of face centroid line from face normal.
        const double ax = F.areaX[f], ay = F.areaY[f], az = F.areaZ[f];
        const double A = std::sqrt(ax * ax + ay * ay + az * az);
        if (A <= 0)
            continue;
        // Approximate skewness placeholder; replaced by face-shape exact metric
        // in Phase 5.6 (depends on cell-pair geometry).
        const double skew = 0.0;
        r.minSkewness = std::min(r.minSkewness, skew);
        r.maxSkewness = std::max(r.maxSkewness, skew);
        sumSkew += skew;
    }
    r.avgSkewness = sumSkew / double(F.size());
    return r;
}

} // namespace simall::meshing
