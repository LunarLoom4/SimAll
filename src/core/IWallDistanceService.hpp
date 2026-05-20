// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/IWallDistanceService.hpp
// Phase  : 22 Pass 3 — public wall-distance service interface.
//
// Why this lives in core/ rather than solver/
// -------------------------------------------
// The wall-distance field d(x) (shortest distance from each cell centroid
// to the nearest wall, satisfying the Eikonal equation |∇d| = 1 with d=0
// on walls) is consumed by every subsystem that models near-wall physics:
//
//   - turbulence : k-ω SST (F1/F2 blending), Spalart-Allmaras (f_v1, f_w),
//                  all DES/IDDES hybrid lengths, low-Re k-ε damping.
//   - heat       : wall-function thermal y+, near-wall conjugate coupling.
//   - emag       : skin-depth boundary-layer resolution in conductors.
//   - chem       : Damköhler-number wall-quench models, surface reaction
//                  near-wall species accumulation.
//   - acoustics  : wall-shear-driven sound source localisation.
//
// To avoid every consumer subsystem pulling in solver/ (and the resulting
// dependency cycle), this interface lives in core/ and is implemented by an
// adapter in solver/WallDistanceService.cpp that wraps the existing
// solver::IWallDistance back-ends (ExactNearest / PoissonWallDistance).
//
// Lifecycle:
//   1. Application constructs the solver-side adapter (which captures
//      references to mesh + boundaries) and registers it under
//      `simall::core::ServiceLocator::set<IWallDistanceService>(...)`.
//   2. Consumers retrieve the service from the locator and call
//      prepare() once after mesh changes (or at every solver iteration if
//      using deforming meshes), then read data()[c] / operator[](c).
//   3. The pointer returned by data() remains valid until the next
//      prepare() call or until the service is destroyed.
// =============================================================================
#pragma once

#include <cstddef>

namespace simall::core {

class IWallDistanceService {
public:
    virtual ~IWallDistanceService() = default;

    /// Recompute (or compute for the first time) the wall-distance field.
    /// Returns the number of cells in the underlying mesh. Idempotent —
    /// calling twice with no mesh change must produce the same result.
    virtual std::size_t prepare() = 0;

    /// Read-only pointer to the per-cell distance array. Size equals the
    /// last prepare() return value (zero if prepare() has not been called).
    /// The pointer is invalidated by the next prepare() call.
    virtual const double* data() const noexcept = 0;

    /// Convenience accessor — undefined behaviour for c >= size().
    double operator[](std::size_t c) const noexcept { return data()[c]; }

    /// Cached size (equals the most recent prepare() return value, or 0).
    virtual std::size_t size() const noexcept = 0;

    /// Implementation tag for diagnostics: "exact-nearest" / "poisson-eikonal"
    /// / "fast-marching" / etc.
    virtual const char* implementation_name() const noexcept = 0;
};

}  // namespace simall::core
