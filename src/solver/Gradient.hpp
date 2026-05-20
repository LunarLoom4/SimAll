// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/Gradient.hpp
// Phase  : 6.3 — least-squares cell-centred gradient reconstruction.
//
// For each cell c with neighbours n (across face f), solve the over-
// determined system
//        d_n · ∇φ_c ≈ (φ_n - φ_c)
// in weighted least-squares sense (weight = 1/|d|). The 3×3 normal matrix
// A^T W A is precomputed once and stored per-cell; gradient evaluation is
// O(neighbours).
// =============================================================================
#pragma once

#include "FieldRegistry.hpp"

#include "meshing/MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

namespace simall::solver
{

class LeastSquaresGradient
{
public:
    explicit LeastSquaresGradient(const meshing::Mesh& mesh);

    /// Compute ∇φ for every interior cell. Boundary contribution is via the
    /// neighbour-cell projection across the face (Dirichlet-style closure).
    void evaluate(const ScalarField& phi, VectorField& grad) const;

private:
    void precompute();
    const meshing::Mesh& mesh_;
    // Symmetric inverse of (A^T W A) stored as 6 doubles per cell:
    // [xx, xy, xz, yy, yz, zz].
    util::aligned_vector<double> mInv_;
};

} // namespace simall::solver
