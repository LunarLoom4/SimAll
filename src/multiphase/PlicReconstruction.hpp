// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/PlicReconstruction.hpp
// Phase  : 12.2 — Piecewise-Linear Interface Calculation (PLIC) using
// Youngs' (1982) gradient-based normal estimation.
//
// For every "interface cell" (0 < α < 1) we approximate the immiscible
// fluid interface by an oriented plane through the cell:
//
//     n · (x - x_c) = δ                with n = -∇α / |∇α|
//
// The signed offset δ is chosen so the half-space {n·(x-x_c) < δ} occupies
// exactly volume α · V_cell. Solved iteratively by bisection on δ against
// the cut-volume of the cell's polyhedron sampled via a tetrahedral
// decomposition of each face fan.
//
// Output (per interface cell):
//   - normal (n_x, n_y, n_z)
//   - offset δ
//   - centroid of the planar interface patch
//   - patch area
//
// These quantities feed advection schemes (geometric MULES, VoFLib-style
// flux operators) and CSF surface-tension models requiring a sharper κ
// estimate than the smoothed-α curvature.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"
#include "utilities/MathTypes.hpp"

#include <vector>

namespace simall::multiphase {

struct PlicPatch {
    util::Vec3d normal   {0,0,0};
    util::Vec3d centroid {0,0,0};
    double      offset   = 0.0;
    double      area     = 0.0;
    bool        valid    = false;
};

class PlicReconstruction {
public:
    /// Compute a PLIC patch for every cell with 0 < α < 1. Cells outside
    /// that range get `valid = false`.
    void reconstruct(const meshing::Mesh& mesh,
                     const solver::ScalarField& alpha);

    const std::vector<PlicPatch>& patches() const noexcept { return patches_; }

private:
    /// Cut-volume of the half-space {n·x < δ} inside cell c.
    double cell_cut_volume(const meshing::Mesh& m, meshing::CellId c,
                           const util::Vec3d& n, double delta) const;

    std::vector<PlicPatch> patches_;
};

}  // namespace simall::multiphase
