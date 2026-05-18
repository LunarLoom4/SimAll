// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/Defeaturing.hpp
// Phase  : 4.3 (Geometry simplification — defeaturing)
//
// Removes detail features from a model so a meshable proxy can be produced
// without polluting the simulation domain with sub-mm geometry. Supports:
//   * fillet removal (radius < threshold)
//   * small-hole filling (disc-shaped face whose bounding circle <
//     `minHoleRadius`)
//   * face removal by user-supplied id list (named selections)
// Implemented atop OCC's `BRepAlgoAPI_Defeaturing` (OCC 7.4+).
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"
#include "utilities/MathTypes.hpp"

#include <vector>

namespace simall::cad {

struct DefeatureOptions {
    double maxFilletRadius = 1.0e-3;   // m
    double minHoleRadius   = 0.0;      // m; 0 disables hole-fill
    bool   removeFillets   = true;
    bool   removeHoles     = true;
    /// Explicit face IDs to remove (from TopologyGraph).
    std::vector<util::PersistentId> faceIds;
};

struct DefeatureReport {
    int filletsRemoved = 0;
    int holesFilled    = 0;
    int facesRemoved   = 0;
    bool anyChange     = false;
};

class Defeaturing {
public:
    void             apply(ShapeHandle& shape, const DefeatureOptions& opts);
    DefeatureReport  lastReport() const noexcept { return report_; }

private:
    DefeatureReport report_;
};

}  // namespace simall::cad
