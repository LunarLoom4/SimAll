// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/LagrangianTracker.hpp
// Phase  : 13 — Discrete Phase Model (Lagrangian particle tracking).
//
//   m_p dv_p/dt = F_drag + F_grav + F_pressure + F_user
//   dx_p/dt     = v_p
//
//   F_drag = (3/4) (μ ρ / d_p) C_d Re_p (u - v_p)        with
//   Re_p   = ρ |u - v_p| d_p / μ
//   C_d    = Schiller-Naumann (1933):
//             24/Re_p (1 + 0.15 Re_p^0.687)   if Re_p < 1000
//             0.44                            else
//
// 3-D polyhedral mesh, fully unstructured. Particle location uses
// face-walking from previous cell. Time integration via classical RK4.
//
// One-way coupling out of the box; two-way momentum coupling exposed via
// `momentum_source(cellId)` for the SIMPLE outer loop to subtract from
// the fluid momentum equation.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <vector>

namespace simall::particles {

struct ParticleSpec {
    double diameter = 1e-5;     // m
    double density  = 2500.0;   // kg/m³ (e.g. silica)
    util::Vec3d gravity{0, 0, -9.81};
    bool   trackRotation = false;
    bool   saffmanLift   = false;   ///< Saffman shear-induced lift (Mei 1992)
    bool   magnusLift    = false;   ///< Magnus rotation-induced lift
};

struct ParticleState {
    util::Vec3d x{0,0,0};
    util::Vec3d v{0,0,0};
    util::Vec3d omega{0,0,0};   ///< particle angular velocity (if trackRotation)
    double      mass     = 0.0;
    double      inertia  = 0.0; ///< moment of inertia = (1/10) m d_p²
    meshing::CellId cell = static_cast<meshing::CellId>(-1);
    bool        active   = true;
};

class LagrangianTracker {
public:
    void initialize(const meshing::Mesh& mesh,
                    double fluidDensity, double fluidViscosity,
                    const ParticleSpec& spec);

    /// Inject N particles uniformly from a point with a velocity.
    void inject_point(const util::Vec3d& location,
                      const util::Vec3d& velocity,
                      std::size_t count);

    /// Advance all particles by dt using fluid cell-centred velocity field
    /// "U". Updates positions, velocities, and current host cell.
    void advance(double dt, const solver::FieldRegistry& fields);

    /// Two-way coupling: per-cell aggregated momentum exchange source (N/m³)
    /// accumulated over the most recent `advance()` call.
    void momentum_source(solver::VectorField& Smom) const;

    std::size_t live() const noexcept;
    const std::vector<ParticleState>& particles() const noexcept { return parts_; }
    std::vector<ParticleState>&       mutable_particles()       noexcept { return parts_; }

private:
    /// Locate cell containing point p starting from cell hint. Returns -1 if
    /// the particle leaves the domain.
    meshing::CellId locate(const util::Vec3d& p, meshing::CellId hint) const;

    const meshing::Mesh* mesh_ = nullptr;
    double rhoF_ = 1.0, muF_ = 1.0e-3;
    ParticleSpec spec_{};

    std::vector<ParticleState> parts_;
    // Per-cell aggregated drag force F_d (last advance call), divided by V_cell.
    util::aligned_vector<double> srcX_, srcY_, srcZ_;
};

}  // namespace simall::particles
