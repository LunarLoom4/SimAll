// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/DemCollision.hpp
// Phase  : 13.5 — Soft-sphere Discrete Element Method (DEM) collisions on
// the Lagrangian particle cloud. Implements the standard Hertz-Mindlin
// non-linear spring-dashpot contact law with Coulomb tangential friction.
//
//   Normal   : F_n = -k_n δ_n^{3/2} n  - γ_n v_{rel,n}
//   Tangent  : F_t = -k_t δ_t          - γ_t v_{rel,t}     (capped by μ|F_n|)
//
// k_n, k_t derive from Hertz contact mechanics with material
// Young's modulus E, Poisson ν, and effective radius / mass.
//
// Particle-particle neighbour search uses a uniform spatial hash grid
// sized to the largest particle diameter (the canonical O(N) DEM lookup).
//
// Particle-wall collisions use the host cell's face normals and
// face-particle distance for the same contact law.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "particles/LagrangianTracker.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace simall::particles
{

struct DemMaterial
{
    double youngsModulus = 1.0e9; // Pa
    double poissonRatio = 0.30;
    double restitution = 0.85; // coefficient of restitution
    double friction = 0.30;    // Coulomb friction coefficient
};

class DemCollision
{
public:
    void initialize(const meshing::Mesh& mesh, const DemMaterial& mat);

    /// Apply contact forces by mutating particle velocities/angular
    /// velocities in-place over time step `dt`. Pure substep — call after
    /// LagrangianTracker::advance(dt, ...).
    void apply(double dt, std::vector<ParticleState>& particles, const ParticleSpec& spec);

    /// Number of contacts processed during the most recent `apply()`.
    std::size_t last_contact_count() const noexcept { return lastContacts_; }

private:
    struct HashKey
    {
        std::int32_t i, j, k;
        bool operator==(const HashKey& o) const { return i == o.i && j == o.j && k == o.k; }
    };
    struct HashHasher
    {
        std::size_t operator()(const HashKey& k) const noexcept
        {
            return (std::size_t) (k.i * 73856093u ^ k.j * 19349663u ^ k.k * 83492791u);
        }
    };

    void contact_pair(
        double dt, ParticleState& a, ParticleState& b, double dA, double dB, double mA, double mB);
    void contact_wall(double dt,
                      ParticleState& p,
                      double dP,
                      double mP,
                      const util::Vec3d& faceNormal,
                      double penetration);

    const meshing::Mesh* mesh_ = nullptr;
    DemMaterial mat_{};
    std::size_t lastContacts_ = 0;
};

} // namespace simall::particles
