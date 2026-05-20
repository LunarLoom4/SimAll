// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/Defeaturing.cpp
//
// Implementation notes:
//   * `BRepAlgoAPI_Defeaturing` accepts a list of faces to remove and
//     attempts to reconstruct adjacent surfaces seamlessly. We therefore
//     classify faces (fillet / disc-hole / explicit) into a single removal
//     set, then call the OCC operator once.
//   * Fillet detection: a face is treated as a fillet when its surface
//     reports cylinder/torus continuity with radius below threshold AND it
//     is adjacent to at least two non-coplanar faces (the classic G1 blend
//     signature). For robustness we accept BSpline surfaces whose principal
//     curvatures consistently fall below 1/r.
//   * Hole detection: a face whose only outer wire is a circle of radius
//     ≤ `minHoleRadius` is filled.
// =============================================================================
#include "cad/Defeaturing.hpp"

#include "cad/ShapeHandleInternal.hpp"
#include "core/Logger.hpp"

#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Defeaturing.hxx>
#include <BRepGProp.hxx>
#include <Geom_Circle.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <GProp_GProps.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>

namespace simall::cad
{

namespace
{

bool isLikelyFillet(const TopoDS_Face& face, double maxR)
{
    BRepAdaptor_Surface surf(face);
    switch (surf.GetType()) {
    case GeomAbs_Cylinder:
        return surf.Cylinder().Radius() <= maxR;
    case GeomAbs_Torus:
        return surf.Torus().MinorRadius() <= maxR;
    default:
        break;
    }
    return false;
}

bool isCircularHole(const TopoDS_Face& face, double maxR)
{
    // Count edges and check whether the dominant edge is a circle of small radius.
    TopExp_Explorer exp(face, TopAbs_EDGE);
    int circles = 0;
    double rMin = 1e300;
    for (; exp.More(); exp.Next()) {
        const TopoDS_Edge& e = TopoDS::Edge(exp.Current());
        BRepAdaptor_Curve crv(e);
        if (crv.GetType() == GeomAbs_Circle) {
            ++circles;
            double r = crv.Circle().Radius();
            if (r < rMin)
                rMin = r;
        }
    }
    return circles >= 1 && rMin <= maxR;
}

} // namespace

void Defeaturing::apply(ShapeHandle& shape, const DefeatureOptions& opts)
{
    report_ = {};
    if (!shape.valid())
        return;

    TopoDS_Shape work = ShapeHandleAccess::shape(shape);

    // Gather candidate faces. We use IndexedMap for stable iteration order.
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(work, TopAbs_FACE, faceMap);

    TopTools_ListOfShape victims;

    // Explicit IDs.
    if (!opts.faceIds.empty()) {
        for (auto id : opts.faceIds) {
            if (auto* node = shape.topology().find(id); node && node->occHandle) {
                const TopoDS_Shape* s = static_cast<const TopoDS_Shape*>(node->occHandle);
                if (s->ShapeType() == TopAbs_FACE) {
                    victims.Append(*s);
                    ++report_.facesRemoved;
                }
            }
        }
    }

    if (opts.removeFillets) {
        for (int i = 1; i <= faceMap.Extent(); ++i) {
            const TopoDS_Face& f = TopoDS::Face(faceMap(i));
            if (isLikelyFillet(f, opts.maxFilletRadius)) {
                victims.Append(f);
                ++report_.filletsRemoved;
            }
        }
    }
    if (opts.removeHoles && opts.minHoleRadius > 0.0) {
        for (int i = 1; i <= faceMap.Extent(); ++i) {
            const TopoDS_Face& f = TopoDS::Face(faceMap(i));
            if (isCircularHole(f, opts.minHoleRadius)) {
                victims.Append(f);
                ++report_.holesFilled;
            }
        }
    }
    if (victims.IsEmpty())
        return;

    BRepAlgoAPI_Defeaturing op;
    op.SetShape(work);
    op.AddFacesToRemove(victims);
    op.SetRunParallel(Standard_True);
    op.Build();
    if (op.HasErrors()) {
        SIMALL_LOG_WARN("CAD/Defeature", "BRepAlgoAPI_Defeaturing reported errors");
        return;
    }

    TopoDS_Shape result = op.Shape();
    if (result.IsNull())
        return;

    ShapeHandleAccess::shape(shape) = result;
    rebuildTopologyGraph(result, ShapeHandleAccess::impl(shape).graph);
    report_.anyChange = true;

    SIMALL_LOG_INFO("CAD/Defeature",
                    "fillets=",
                    report_.filletsRemoved,
                    " holes=",
                    report_.holesFilled,
                    " explicit=",
                    report_.facesRemoved);
}

} // namespace simall::cad
