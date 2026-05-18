// =============================================================================
// SimAll Beta - Acoustics Subsystem
// File   : src/acoustics/FwhSurface.hpp
// Phase  : 18 — Ffowcs Williams-Hawkings (FW-H) far-field acoustic
// integral over a permeable / impermeable data surface. Used to predict
// fan noise, rotor noise, jet noise, etc. from LES/DES surface data.
//
// Formulation (Farassat 1A retarded-time form, permeable surface):
//
//   p'(x,t) = p'_T (thickness) + p'_L (loading)
//
//   4π p'_T(x,t) = ∫_S [ ρ_0 (U_n + u_n)  / (r (1-M_r)²) ]_τ dS
//                + ∫_S [ ρ_0 (U_n + u_n) (rṀ_r + a_0(M_r - M²)) / (r²(1-M_r)³) ]_τ dS
//
//   4π p'_L(x,t) = (1/a_0) ∫_S [ L̇_r / (r(1-M_r)²) ]_τ dS
//                + ∫_S [ (L_r - L_M) / (r²(1-M_r)²) ]_τ dS
//                + (1/a_0) ∫_S [ L_r (r Ṁ_r + a_0(M_r - M²)) / (r²(1-M_r)³) ]_τ dS
//
// Here, () indicates evaluation at emission (retarded) time τ such that
//   t - τ = |x - y(τ)| / a_0     (Newton iteration per source panel)
//
// Symbols
//   S     : data surface (set of boundary faces)
//   ρ_0   : far-field mean density          [kg/m³]
//   a_0   : far-field speed of sound        [m/s]
//   U     : surface velocity (impermeable: rigid body motion; permeable: 0)
//   u     : fluid velocity at surface       [m/s]
//   L     : loading vector  L_i = (p δ_ij + ρ u_i (u_j - U_j)) n_j  [Pa]
//
// The class is **storage-aware**: callers stream per-iteration surface
// samples to `record_sample(t, ...)`, and `emit(observers, tStart, tEnd, dt)`
// returns p'(t) sampled at the observer points.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <vector>

namespace simall::acoustics {

struct FwhObserver {
    util::Vec3d position{0, 0, 0};
};

struct FwhSurfaceProps {
    double rho0 = 1.225;             // far-field density [kg/m³]
    double a0   = 340.0;             // far-field sound speed [m/s]
    bool   permeable = true;         // false ⇒ skip thickness term entirely
};

class FwhSurface {
public:
    /// Define the integration surface from a set of boundary face zones.
    void initialize(const meshing::Mesh& mesh,
                    const std::vector<meshing::ZoneId>& zones,
                    FwhSurfaceProps props);

    /// Record one time-history sample from the running solver. Samples must
    /// be supplied at uniform Δt; the integrator interpolates linearly in
    /// retarded time between adjacent samples.
    void record_sample(double t, const solver::FieldRegistry& fields,
                       const util::Vec3d& surfaceVelocity = {0, 0, 0});

    /// Compute the acoustic pressure p'(t_obs) for each observer at the
    /// requested observer-time stamps (uniform sampling [tStart, tEnd] in
    /// step dt). Returns a flat array indexed [obsIdx * Nt + it].
    std::vector<double> emit(const std::vector<FwhObserver>& observers,
                             double tStart, double tEnd, double dt) const;

    std::size_t panel_count() const noexcept { return panels_.size(); }
    std::size_t sample_count() const noexcept { return times_.size(); }

private:
    struct Panel {
        util::Vec3d centroid;
        util::Vec3d normal;       // unit outward normal
        double      area = 0;
    };
    struct Snapshot {
        // Per-panel cached scalars at one time sample.
        std::vector<util::Vec3d> u;          // fluid velocity at panel
        std::vector<double>      p;          // gauge pressure
        std::vector<double>      rho;        // local density (or ρ_0)
        util::Vec3d              surfaceVel; // rigid-body velocity of S
    };

    const meshing::Mesh*  mesh_ = nullptr;
    FwhSurfaceProps       props_{};
    std::vector<Panel>    panels_;
    std::vector<double>   times_;
    std::vector<Snapshot> snaps_;
};

}  // namespace simall::acoustics
