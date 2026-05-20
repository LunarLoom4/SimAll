// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/GradientLimiter.hpp
// Phase  : 6.6 — Per-cell scalar gradient limiters for high-resolution
// schemes near discontinuities. Implements Barth-Jespersen (1989) and
// Venkatakrishnan (1993). Both produce a per-cell ψ ∈ [0,1] so that the
// limited reconstruction
//
//     φ_face = φ_c + ψ_c · (∇φ_c · (x_face - x_c))
//
// remains bounded by the cell's neighbourhood φ_min, φ_max — i.e. no new
// extrema are introduced. Venkatakrishnan additionally relaxes the strict
// monotonicity in smooth regions via a tunable threshold K · h^{3/2}, which
// avoids the convergence stalls characteristic of Barth-Jespersen.
//
// Output: scalar field "limiter_<name>" with the per-cell ψ value.
// =============================================================================
#pragma once

#include "FieldRegistry.hpp"

#include "meshing/MeshStorage.hpp"

namespace simall::solver
{

enum class LimiterKind
{
    None,
    BarthJespersen,
    Venkatakrishnan
};

class GradientLimiter
{
public:
    explicit GradientLimiter(const meshing::Mesh& mesh);

    /// Compute ψ_c for every cell from the cell-centred field φ and its
    /// (unlimited) gradient ∇φ. The result is written into psi.
    void evaluate(const ScalarField& phi,
                  const VectorField& grad,
                  ScalarField& psi,
                  LimiterKind kind,
                  double venkatK = 5.0) const;

private:
    const meshing::Mesh& mesh_;
};

} // namespace simall::solver
