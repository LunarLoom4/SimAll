// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/FeatureEdges.cpp
//
// Algorithm:
//   1. Build the EDGE→list-of-FACES adjacency map via TopExp::MapShapesAndAncestors.
//   2. Classify each edge:
//      * 0 incident faces      → ignore
//      * 1 incident face       → boundary
//      * ≥ 3 incident faces    → non-manifold
//      * 2 faces, dihedral ≥ τ → sharp
//   3. For each retained edge, sample the underlying curve uniformly into
//      a polyline. We use BRepAdaptor_Curve so analytic curves are sampled
//      coarsely (lines/circles) and BSplines finely.
// =============================================================================
#include "cad/FeatureEdges.hpp"
#include "cad/ShapeHandleInternal.hpp"
#include "core/Logger.hpp"

#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepGProp_Face.hxx>
#include <GeomLProp_SLProps.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>

#include <algorithm>
#include <cmath>

namespace simall::cad {

namespace {

util::PersistentId resolveEdgeId(const TopologyGraph& g, const TopoDS_Shape& edge) {
    for (auto const& [id, n] : g.nodes()) {
        if (n.type == TopologyType::Edge && n.occHandle &&
            static_cast<const TopoDS_Shape*>(n.occHandle)->IsSame(edge))
            return id;
    }
    return util::kInvalidId;
}

// Sample a face's outward normal at the midpoint of an edge.
gp_Dir faceNormalAt(const TopoDS_Face& face, const TopoDS_Edge& edge, double t) {
    Standard_Real f, l;
    Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(edge, face, f, l);
    if (pcurve.IsNull()) return gp_Dir(0, 0, 1);
    gp_Pnt2d uv = pcurve->Value(f + (l - f) * t);
    BRepAdaptor_Surface surf(face);
    GeomLProp_SLProps p(surf.Surface().Surface(), uv.X(), uv.Y(), 1, 1e-9);
    if (!p.IsNormalDefined()) return gp_Dir(0, 0, 1);
    gp_Dir n = p.Normal();
    if (face.Orientation() == TopAbs_REVERSED) n.Reverse();
    return n;
}

double dihedralDeg(const TopoDS_Face& f1, const TopoDS_Face& f2, const TopoDS_Edge& e) {
    gp_Dir n1 = faceNormalAt(f1, e, 0.5);
    gp_Dir n2 = faceNormalAt(f2, e, 0.5);
    double c  = std::clamp(n1.Dot(n2), -1.0, 1.0);
    return std::acos(c) * 180.0 / 3.14159265358979323846;
}

void sampleEdgePolyline(const TopoDS_Edge& edge,
                        std::vector<util::Vec3d>& out,
                        int samples) {
    BRepAdaptor_Curve crv(edge);
    Standard_Real u0 = crv.FirstParameter();
    Standard_Real u1 = crv.LastParameter();
    if (!std::isfinite(u0) || !std::isfinite(u1) || u1 <= u0) return;

    // Analytic curves don't benefit from heavy sampling; line: 2 points.
    if (crv.GetType() == GeomAbs_Line) samples = 2;
    else if (crv.GetType() == GeomAbs_Circle || crv.GetType() == GeomAbs_Ellipse)
        samples = std::max(samples, 24);

    out.reserve(out.size() + std::size_t(samples));
    for (int i = 0; i < samples; ++i) {
        double t = double(i) / double(samples - 1);
        gp_Pnt p = crv.Value(u0 + (u1 - u0) * t);
        out.push_back({p.X(), p.Y(), p.Z()});
    }
}

}  // namespace

std::vector<FeatureEdgePolyline> FeatureEdges::extract(const ShapeHandle& shape,
                                                       const FeatureEdgeOptions& opts) {
    std::vector<FeatureEdgePolyline> out;
    if (!shape.valid()) return out;

    const TopoDS_Shape& s = ShapeHandleAccess::shape(shape);

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(s, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    const int nEdges = edgeFaceMap.Extent();
    out.reserve(std::size_t(nEdges));

    for (int i = 1; i <= nEdges; ++i) {
        const TopoDS_Edge& edge = TopoDS::Edge(edgeFaceMap.FindKey(i));
        const TopTools_ListOfShape& faces = edgeFaceMap.FindFromIndex(i);
        int nf = faces.Extent();

        FeatureEdgePolyline poly;
        poly.edgeId = resolveEdgeId(shape.topology(), edge);

        if (nf == 0) continue;
        if (nf == 1 && opts.includeBoundary) {
            poly.isBoundary = true;
            sampleEdgePolyline(edge, poly.points, 16);
            if (!poly.points.empty()) out.push_back(std::move(poly));
            continue;
        }
        if (nf >= 3 && opts.includeNonManifold) {
            poly.isNonManifold = true;
            sampleEdgePolyline(edge, poly.points, 16);
            if (!poly.points.empty()) out.push_back(std::move(poly));
            continue;
        }
        if (nf == 2 && opts.includeSharp) {
            TopTools_ListIteratorOfListOfShape lit(faces);
            const TopoDS_Face& f1 = TopoDS::Face(lit.Value()); lit.Next();
            const TopoDS_Face& f2 = TopoDS::Face(lit.Value());
            double ang = dihedralDeg(f1, f2, edge);
            // dihedral measured between normals; "sharp" means deviation from
            // parallel (i.e. angle far from 0). Use 180 - ang for convex/concave.
            double sharpness = std::min(ang, 180.0 - ang);
            if (sharpness >= opts.dihedralAngleDeg) {
                poly.isSharp = true;
                sampleEdgePolyline(edge, poly.points, 16);
                if (!poly.points.empty()) out.push_back(std::move(poly));
            }
        }
    }
    SIMALL_LOG_INFO("CAD/FeatureEdges",
                    "Extracted ", out.size(), " feature edges (",
                    nEdges, " total)");
    return out;
}

}  // namespace simall::cad
