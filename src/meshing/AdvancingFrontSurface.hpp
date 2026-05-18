// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/AdvancingFrontSurface.hpp
// Phase  : 6.7 — Advancing-front surface mesh generator (Lo 1985).
//
// Generates a triangular surface mesh on a planar polygon (or projected
// parametric patch) from its boundary loop using the Lo advancing-front
// scheme:
//
//   1. Seed the active front with the input boundary edges.
//   2. For each front edge e = (A, B), choose the "ideal" apex P such that
//      |AP| = |BP| = h(midpoint) and P lies on the active-domain side.
//   3. Search candidate apex nodes within radius R for any existing front
//      vertex; if one is closer than the ideal apex, snap to it.
//   4. Validate the new triangle (no self-intersection with existing
//      front edges, positive orientation, non-degenerate aspect).
//   5. Insert the triangle, replace edge e on the front with the two new
//      open edges, remove pairs that close on themselves.
//   6. Continue until the front is empty.
//
// Mesh sizing h(x) is provided via a callable (constant or analytic).
// Output is a triangle list (NodeId triples) attached to a meshing::Mesh
// face list with boundary zone equal to props.outputZone.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <functional>
#include <vector>

namespace simall::meshing {

struct AfsProps {
    double  baseSize        = 0.05;   // default characteristic edge h
    double  searchRadiusMul = 1.6;    // R = mul · h
    double  minQuality      = 0.2;    // min normalised triangle quality
    double  planeNormalX    = 0.0;    // surface normal for 2-D plane case
    double  planeNormalY    = 0.0;
    double  planeNormalZ    = 1.0;
    std::uint32_t outputZone= 1;
    std::size_t maxIters    = 200000;
};

class AdvancingFrontSurface {
public:
    using SizingFn = std::function<double(double x, double y, double z)>;

    void initialize(AfsProps props, SizingFn sizing = nullptr);

    /// Seed the front with an ordered polygon loop (closed). Multiple loops
    /// allowed (call seed_boundary repeatedly).  Inner loops define holes.
    void seed_boundary(const std::vector<util::Vec3d>& loop);

    /// Run the advancing front until completion. Returns triangle count.
    std::size_t generate(Mesh& outMesh);

    std::size_t triangle_count() const noexcept { return tris_.size(); }

    const AfsProps& props() const noexcept { return p_; }

private:
    struct Edge { NodeId a, b; };
    struct Tri  { NodeId a, b, c; };

    bool triangle_valid(const util::Vec3d& A, const util::Vec3d& B,
                        const util::Vec3d& P, double h_ref) const;
    NodeId find_or_create_node(const util::Vec3d& p, double tol);

    AfsProps             p_{};
    SizingFn             sizing_;
    std::vector<util::Vec3d> nodes_;
    std::vector<Edge>    front_;
    std::vector<Tri>     tris_;
};

}  // namespace simall::meshing
