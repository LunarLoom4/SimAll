// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/Healing.hpp
// Phase  : 4.3 (Geometry repair)
//
// Wraps OCC's ShapeFix machinery. The healing pipeline performs:
//   1. tolerance harmonisation
//   2. wireframe closure (stitching gaps)
//   3. small-edge collapse and sliver-face removal
//   4. orientation / shell consistency
//   5. self-intersection cleanup
// Each stage can be disabled per `HealingOptions`. The operation is
// destructive — the ShapeHandle's underlying TopoDS_Shape and topology
// graph are replaced. PersistentIdManager (separate module) can be used to
// re-map IDs across the repair if required by named selections.
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"

#include <vector>
#include <string>

namespace simall::cad {

struct HealingReport {
    int    smallEdgesRemoved   = 0;
    int    sliverFacesRemoved  = 0;
    int    gapsStitched        = 0;
    int    orientationFixed    = 0;
    bool   anyChange           = false;
    std::vector<std::string> messages;
};

class Healing {
public:
    void              repair(ShapeHandle& shape, const HealingOptions& opts = {});
    HealingReport     lastReport() const noexcept { return report_; }

private:
    HealingReport report_;
};

}  // namespace simall::cad
