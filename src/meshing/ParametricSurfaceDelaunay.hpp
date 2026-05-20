// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/ParametricSurfaceDelaunay.hpp
// Phase  : 5.2 — surface meshing helper (Bowyer-Watson in CAD (u,v) space).
//
// SCOPE: this class operates in the 2-D PARAMETER space of a CAD face. The
// 3-D surface mesh is produced by the SurfaceMesher (Mesher.cpp) which:
//   1. samples each CAD face in (u,v),
//   2. runs ParametricSurfaceDelaunay to triangulate the parameter domain
//      while honouring trimming-loop edges as constraints,
//   3. back-projects every (u,v) vertex through Geom_Surface::D0 to obtain
//      the corresponding (x,y,z) point on the NURBS surface.
//
// This is the standard pipeline used by ANSYS Meshing, cfMesh, Pointwise
// and Siemens Simcenter. It is NOT a 2-D solver path.
// =============================================================================
#pragma once

#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace simall::meshing
{

struct ParamPoint
{
    double u, v;
};
struct ParamTriangle
{
    std::array<std::uint32_t, 3> v;
};
struct ParamEdge
{
    std::uint32_t a, b;
    bool operator==(const ParamEdge&) const = default;
};

class ParametricSurfaceDelaunay
{
public:
    /// Triangulate the input parameter-space point cloud.
    /// `constraints` are oriented edges from CAD trimming loops that MUST
    /// appear in the final mesh.
    void triangulate(const std::vector<ParamPoint>& points,
                     const std::vector<ParamEdge>& constraints = {});

    const std::vector<ParamTriangle>& triangles() const noexcept { return tris_; }
    const std::vector<ParamPoint>& points() const noexcept { return pts_; }

private:
    std::vector<ParamPoint> pts_;
    std::vector<ParamTriangle> tris_;

    void insert_point(std::uint32_t p);
    bool in_circle(std::uint32_t a, std::uint32_t b, std::uint32_t c, std::uint32_t d) const;
    double orient2d(std::uint32_t a, std::uint32_t b, std::uint32_t c) const;
    void recover_edges(const std::vector<ParamEdge>& constraints);
};

} // namespace simall::meshing
