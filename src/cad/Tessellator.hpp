// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/Tessellator.hpp
// Phase  : 4.4 (BRep → triangulation)
//
// Promotes the previous in-line implementation from CadKernel to a focused
// module. Two modes:
//   * `tessellate`     — uniform OCC IncrementalMesh with caller-supplied
//                        chordal deflection + angular tolerance.
//   * `adaptive`       — curvature-aware: target edge length
//                          h = min(maxEdge, sqrt(2 * ε * R))
//                        where R is local face curvature radius and ε is
//                        the chordal tolerance. Coarse on planes, fine on
//                        small fillets.
// Outputs preserve the face↔triangle persistent-id mapping required for
// picking and result reporting.
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"

namespace simall::cad {

class Tessellator {
public:
    TriangleMesh tessellate(const ShapeHandle& shape, const TessellationParams& p = {});
    TriangleMesh adaptive  (const ShapeHandle& shape, double maxEdge, double chordalEpsilon);
};

}  // namespace simall::cad
