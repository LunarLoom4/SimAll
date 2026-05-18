// =============================================================================
// SimAll Beta - AMR Subsystem
// File   : src/amr/AdaptiveRefinement.hpp
// Phase  : 16 — Solution-adaptive mesh refinement controller.
//
// Computes a per-cell scalar error indicator from a chosen field, then
// produces refinement/coarsening flags using a fraction-based strategy
// (refine top X% of cells, coarsen bottom Y%), with hard min/max level
// constraints to avoid runaway refinement.
//
// Available indicators:
//   - GradientMagnitude:  |∇φ| · h    (h = V^{1/3})
//   - SecondDerivative :  |∇²φ| · h²  (estimated via gradient-of-gradient)
//   - Jump            :   max neighbour jump |φ_n - φ_c|
//   - Vorticity       :   |∇×U|       (vector field "U")
//
// The flags vector is intended to drive a downstream conformal/octree
// re-meshing pass; this controller does NOT mutate the mesh in place.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace simall::amr {

enum class Indicator { GradientMagnitude, SecondDerivative, Jump, Vorticity };

enum class RefineFlag : std::int8_t { Coarsen = -1, Keep = 0, Refine = +1 };

struct AmrOptions {
    Indicator indicator       = Indicator::GradientMagnitude;
    std::string fieldName     = "p";       // scalar field (Vorticity uses "U")
    double refineFraction     = 0.10;      // top 10 % flagged refine
    double coarsenFraction    = 0.30;      // bottom 30 % flagged coarsen
    int    minLevel           = 0;
    int    maxLevel           = 5;
};

class AdaptiveRefinement {
public:
    /// Compute per-cell indicators ε_c for the chosen field on the mesh.
    void compute_indicator(const meshing::Mesh& mesh,
                           const solver::FieldRegistry& fields,
                           AmrOptions opt);

    /// Mark cells. `level[c]` is the current refinement level per cell
    /// (provided by the mesh manager). Flags are written into `flags`.
    void mark(const std::vector<std::int32_t>& level,
              std::vector<RefineFlag>& flags,
              AmrOptions opt) const;

    const std::vector<double>& indicator() const noexcept { return ind_; }

private:
    std::vector<double> ind_;
};

}  // namespace simall::amr
