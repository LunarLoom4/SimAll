// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/Tessellator.cpp
// =============================================================================
#include "cad/Tessellator.hpp"

#include "cad/ShapeHandleInternal.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <GeomLProp_SLProps.hxx>
#include <gp_Pnt.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

namespace simall::cad
{

namespace
{

util::PersistentId resolveFaceId(const TopologyGraph& g, const TopoDS_Shape& face)
{
    for (auto const& [id, n] : g.nodes()) {
        if (n.type == TopologyType::Face && n.occHandle
            && static_cast<const TopoDS_Shape*>(n.occHandle)->IsSame(face))
            return id;
    }
    return util::kInvalidId;
}

double estimateMinCurvatureRadius(const TopoDS_Face& face)
{
    BRepAdaptor_Surface surf(face);
    Standard_Real u0 = surf.FirstUParameter(), u1 = surf.LastUParameter();
    Standard_Real v0 = surf.FirstVParameter(), v1 = surf.LastVParameter();
    if (!std::isfinite(u0) || !std::isfinite(u1) || !std::isfinite(v0) || !std::isfinite(v1))
        return 1e30;

    double rMin = 1e30;
    constexpr int N = 5;
    auto h = surf.Surface().Surface();
    if (h.IsNull())
        return rMin;
    for (int iu = 0; iu < N; ++iu) {
        for (int iv = 0; iv < N; ++iv) {
            double u = u0 + (u1 - u0) * (iu + 0.5) / N;
            double v = v0 + (v1 - v0) * (iv + 0.5) / N;
            GeomLProp_SLProps p(h, u, v, 2, 1e-9);
            if (!p.IsCurvatureDefined())
                continue;
            double k1 = std::abs(p.MaxCurvature());
            double k2 = std::abs(p.MinCurvature());
            double kMax = std::max(k1, k2);
            if (kMax > 1e-12)
                rMin = std::min(rMin, 1.0 / kMax);
        }
    }
    return rMin;
}

void emitFace(const TopoDS_Face& face, const TopologyGraph& g, TriangleMesh& out)
{
    TopLoc_Location loc;
    Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
    if (tri.IsNull())
        return;

    util::PersistentId faceId = resolveFaceId(g, face);

    const std::uint32_t base = static_cast<std::uint32_t>(out.points.size());
    const auto& trsf = loc.Transformation();
    const Standard_Integer nNodes = tri->NbNodes();
    out.points.reserve(out.points.size() + std::size_t(nNodes));
    for (Standard_Integer i = 1; i <= nNodes; ++i) {
        gp_Pnt pt = tri->Node(i).Transformed(trsf);
        out.points.push_back({pt.X(), pt.Y(), pt.Z()});
    }
    const Standard_Integer nTris = tri->NbTriangles();
    out.triangles.reserve(out.triangles.size() + std::size_t(nTris));
    out.triangleFaceId.reserve(out.triangleFaceId.size() + std::size_t(nTris));

    // Respect the parent face orientation when emitting triangle winding.
    bool reversed = (face.Orientation() == TopAbs_REVERSED);
    for (Standard_Integer i = 1; i <= nTris; ++i) {
        Standard_Integer a, b, c;
        tri->Triangle(i).Get(a, b, c);
        if (reversed)
            std::swap(b, c);
        out.triangles.push_back({base + std::uint32_t(a - 1),
                                 base + std::uint32_t(b - 1),
                                 base + std::uint32_t(c - 1)});
        out.triangleFaceId.push_back(faceId);
    }
}

} // namespace

TriangleMesh Tessellator::tessellate(const ShapeHandle& shape, const TessellationParams& p)
{
    TriangleMesh out;
    if (!shape.valid())
        return out;

    const TopoDS_Shape& s = ShapeHandleAccess::shape(shape);

    BRepMesh_IncrementalMesh mesher(
        s, p.deflection, p.relative ? Standard_True : Standard_False, p.angle, Standard_True);
    mesher.Perform();

    for (TopExp_Explorer it(s, TopAbs_FACE); it.More(); it.Next()) {
        emitFace(TopoDS::Face(it.Current()), shape.topology(), out);
    }
    SIMALL_LOG_INFO("CAD/Tess",
                    "uniform deflection=",
                    p.deflection,
                    " angle=",
                    p.angle,
                    " → ",
                    out.triangles.size(),
                    " tris");
    return out;
}

TriangleMesh Tessellator::adaptive(const ShapeHandle& shape, double maxEdge, double chordalEpsilon)
{
    TriangleMesh out;
    if (!shape.valid())
        return out;

    const TopoDS_Shape& s = ShapeHandleAccess::shape(shape);

    // Per-face IncrementalMesh: pick per-face deflection from local curvature.
    for (TopExp_Explorer it(s, TopAbs_FACE); it.More(); it.Next()) {
        const TopoDS_Face& f = TopoDS::Face(it.Current());
        double R = estimateMinCurvatureRadius(f);
        double hCurv = std::sqrt(std::max(0.0, 2.0 * chordalEpsilon * R));
        double h = std::min(maxEdge, hCurv > 0 ? hCurv : maxEdge);
        double defl = std::max(chordalEpsilon, 0.5 * h);

        BRepMesh_IncrementalMesh m(f, defl, Standard_False, 0.35, Standard_True);
        m.Perform();
        emitFace(f, shape.topology(), out);
    }
    SIMALL_LOG_INFO("CAD/Tess",
                    "adaptive maxEdge=",
                    maxEdge,
                    " ε=",
                    chordalEpsilon,
                    " → ",
                    out.triangles.size(),
                    " tris");
    return out;
}

} // namespace simall::cad
