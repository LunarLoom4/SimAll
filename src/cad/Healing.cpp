// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/Healing.cpp
// =============================================================================
#include "cad/Healing.hpp"
#include "cad/ShapeHandleInternal.hpp"
#include "core/Logger.hpp"

#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <ShapeFix_FixSmallFace.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeUpgrade_RemoveInternalWires.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <ShapeAnalysis_ShapeTolerance.hxx>
#include <TopAbs.hxx>

namespace simall::cad {

void Healing::repair(ShapeHandle& shape, const HealingOptions& opts) {
    report_ = {};
    if (!shape.valid()) return;

    TopoDS_Shape work = ShapeHandleAccess::shape(shape);

    if (opts.harmonizeTolerances) {
        Handle(ShapeFix_Shape) fixer = new ShapeFix_Shape(work);
        fixer->SetPrecision(opts.tolerance);
        fixer->SetMinTolerance(opts.tolerance * 0.01);
        fixer->SetMaxTolerance(opts.tolerance * 100.0);
        fixer->Perform();
        if (!fixer->Shape().IsNull()) {
            work = fixer->Shape();
            report_.anyChange = true;
            report_.messages.emplace_back("ShapeFix_Shape: harmonised tolerances");
        }
    }

    if (opts.stitchGaps || opts.collapseTinyEdges) {
        Handle(ShapeFix_Wireframe) wf = new ShapeFix_Wireframe(work);
        wf->SetPrecision(opts.tolerance);
        wf->SetMaxTolerance(opts.tolerance * 100.0);
        wf->ModeDropSmallEdges() = opts.collapseTinyEdges ? Standard_True : Standard_False;
        if (opts.stitchGaps        && wf->FixWireGaps())   { report_.gapsStitched   = 1; report_.anyChange = true; }
        if (opts.collapseTinyEdges && wf->FixSmallEdges()) { report_.smallEdgesRemoved = 1; report_.anyChange = true; }
        if (!wf->Shape().IsNull()) work = wf->Shape();
    }

    if (opts.removeSlivers) {
        Handle(ShapeFix_FixSmallFace) sff = new ShapeFix_FixSmallFace();
        sff->Init(work);
        sff->SetPrecision(opts.tolerance);
        sff->Perform();
        TopoDS_Shape after = sff->FixShape();
        if (!after.IsNull() && !after.IsSame(work)) {
            work = after;
            report_.sliverFacesRemoved = 1;
            report_.anyChange = true;
            report_.messages.emplace_back("ShapeFix_FixSmallFace: removed slivers");
        }
    }

    // Final validity sanity-check; messages but never aborts (caller may
    // accept partial repairs).
    BRepCheck_Analyzer ana(work);
    if (!ana.IsValid())
        report_.messages.emplace_back("BRepCheck_Analyzer: residual issues remain");

    // Reseat the shape and rebuild the topology graph in place.
    ShapeHandleAccess::shape(shape) = work;
    rebuildTopologyGraph(work, ShapeHandleAccess::impl(shape).graph);

    SIMALL_LOG_INFO("CAD/Heal",
                    "tol=", opts.tolerance,
                    " changed=", (report_.anyChange ? "yes" : "no"),
                    " stitched=", report_.gapsStitched,
                    " smallEdges=", report_.smallEdgesRemoved,
                    " slivers=", report_.sliverFacesRemoved);
}

}  // namespace simall::cad
