// =============================================================================
// SimAll Beta - Morphing Subsystem
// File   : src/morphing/RbfMorpher.hpp
// Phase  : 14 — Radial Basis Function mesh morphing for FSI / moving-mesh.
//
// Given prescribed displacements at "control points" (typically wall nodes
// or surface point sets), interpolate a smooth volume displacement onto
// every interior mesh node:
//
//     s(x) = Σ_i α_i φ(‖x - x_i‖) + p(x)
//
// where φ is a compactly-supported / global radial basis (here we use the
// Wendland C2 compactly-supported function as the default — sparse system,
// scalable to large meshes) and p(x) is a linear polynomial (1, x, y, z)
// added to recover exact rigid motions.
//
// Wendland C2:  φ(r) = (1 - r/R)_+^4 (4 r/R + 1)
//
// Solve (Φ + λ I)·α + P·β = d ;   P^T α = 0
// via dense LDL^T (control-point counts typically O(1e3-1e5); larger sets
// would need a parallel sparse path which is reserved for Pass 7).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

#include <vector>

namespace simall::morphing {

struct RbfOptions {
    double supportRadius   = 0.0;   // 0 ⇒ auto-set to bbox diag
    double regularization  = 1e-12; // ridge to keep system SPD
    bool   addPolynomial   = true;  // append linear poly for affine recovery
};

class RbfMorpher {
public:
    /// Set control points + their prescribed displacements (3-vector each).
    /// Sizes must match.
    void set_controls(const std::vector<util::Vec3d>& positions,
                      const std::vector<util::Vec3d>& displacements);

    /// Train the interpolant (factorise the dense system once per FSI step).
    void train(RbfOptions opt = {});

    /// Apply displacement to every node of the supplied mesh in place. Faces
    /// & cells' geometric quantities are recomputed afterwards via
    /// Mesh::compute_geometry().
    void apply(meshing::Mesh& mesh) const;

    /// Standalone evaluation at an arbitrary point (used by particle injection
    /// or post-processing).
    util::Vec3d evaluate(const util::Vec3d& x) const;

private:
    std::vector<util::Vec3d> ctrl_;
    std::vector<util::Vec3d> disp_;
    // Coefficient vectors per Cartesian component (αx, αy, αz).
    std::vector<double> ax_, ay_, az_;
    // Polynomial coefficients [c0, cx, cy, cz] per component.
    double bx_[4]{0,0,0,0}, by_[4]{0,0,0,0}, bz_[4]{0,0,0,0};
    double R_ = 1.0;
    bool   havePoly_ = true;
    bool   trained_  = false;
};

}  // namespace simall::morphing
