// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/TetrahedralMesher.hpp
// Phase  : 5.3 — full 3-D Delaunay tetrahedralization (Bowyer-Watson).
//
// Produces a conforming tetrahedral volume mesh from:
//   - a triangulated surface boundary (closed, watertight, outward-oriented)
//   - optional internal Steiner points (for size-field control)
//
// Algorithm:
//   1. Construct an enclosing super-tetrahedron containing all boundary nodes.
//   2. Insert every boundary node via Bowyer-Watson cavity rebuild
//      (4D in-sphere predicate; orient3d positivity check on cavity faces).
//   3. Insert size-driven Steiner points until target edge length / quality
//      criteria are satisfied (longest-edge bisection of large tets, plus
//      circumcentre insertion for poorly-shaped tets — "Delaunay refinement"
//      à la Shewchuk/Ruppert in 3-D).
//   4. Recover boundary faces via tetrahedral flipping (small-polyhedron
//      recovery — Si & Gärtner).
//   5. Tag exterior tets (outside the boundary surface) by flood-fill from
//      the super-tet vertices, then strip them so the result is a watertight
//      tetrahedral volume mesh of the input shell's interior.
//   6. Materialise the result as a polyhedral simall::meshing::Mesh via
//      ConnectivityBuilder.
// =============================================================================
#pragma once

#include "MeshStorage.hpp"

#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace simall::meshing
{

struct TetMeshOptions
{
    double targetEdgeLength = 0.1;     // target h
    double minDihedralAngleDeg = 10.0; // quality filter
    int maxRefinementPasses = 6;
    bool recoverBoundary = true;
};

class TetrahedralMesher
{
public:
    /// Boundary surface: triangulated, watertight, outward-oriented.
    /// `points` are the (x,y,z) coordinates; `boundaryTris` indexes into points.
    void mesh(const std::vector<util::Vec3d>& points,
              const std::vector<std::array<std::uint32_t, 3>>& boundaryTris,
              TetMeshOptions opt,
              Mesh& out);

private:
    struct Tet
    {
        std::array<std::uint32_t, 4> v;
    };
    struct Face3
    {
        std::array<std::uint32_t, 3> v;
    };

    std::vector<util::Vec3d> pts_;
    std::vector<Tet> tets_;
    std::uint32_t superStart_ = 0;

    // ---- predicates (robust enough for typical CAD meshes) ----------------
    double orient3d(std::uint32_t a, std::uint32_t b, std::uint32_t c, std::uint32_t d) const;
    bool in_sphere(
        std::uint32_t a, std::uint32_t b, std::uint32_t c, std::uint32_t d, std::uint32_t e) const;

    // ---- core operations --------------------------------------------------
    void insert_point(std::uint32_t p);
    void initialize_super_tet(const std::vector<util::Vec3d>& nodes);
    void strip_exterior(const std::vector<std::array<std::uint32_t, 3>>& boundaryTris);
    void refine(TetMeshOptions opt);

    double tet_quality_ar(const Tet& t) const;
    util::Vec3d circumcentre(const Tet& t) const;
};

} // namespace simall::meshing
