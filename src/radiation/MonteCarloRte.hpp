// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/MonteCarloRte.hpp
// Phase  : 11.5 — Reverse Monte-Carlo Ray Tracing of the RTE.
//
// Algorithm (Howell 1968, Modest 2013 §20):
//   1. Per cell, emit N_p photon bundles with energy
//        E_p = 4 κ σ T⁴ V_cell / N_p
//   2. Sample isotropic direction (sphere); travel distance to next
//      collision is L = -ln(U) / (κ + σ_s) with U ~ U(0,1).
//   3. At collision, the bundle is absorbed (prob κ/β) or scattered
//      (prob σ_s/β with phase Φ — isotropic here).
//   4. Track ray-segment endpoints by walking face graph; on hitting a
//      boundary wall, deposit (1-ρ_w)·E_p into the wall, reflect rest
//      diffusely.
//   5. Aggregate per-cell absorbed energy → ΔE_abs[c].
//      S_rad[c] = (ΔE_abs[c] - 4 κ σ T⁴ V_cell) / V_cell .
//
// Mesh walk: rays advance one cell at a time via cell-to-face adjacency.
// Cell-segment exit face is found by minimum positive t = (face·n -
// origin·n)/(s·n) — orthogonal projection (centroid-based approximation;
// adequate for moderate-skewness meshes — flagged TODO for stricter test).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <cstdint>
#include <random>
#include <vector>

namespace simall::radiation
{

struct McWallSpec
{
    meshing::ZoneId zone;
    double emissivity = 1.0;
    double temperature = 300.0;
};

struct MonteCarloProps
{
    double absorption = 0.5; // κ
    double scattering = 0.0; // σ_s
    double refractiveN = 1.0;
    std::size_t nPhotonsPerCell = 200;
    std::size_t maxSegments = 64;
    std::uint64_t rngSeed = 0xCAFEFEED'BABE5EE5ULL;
};

class MonteCarloRte
{
public:
    bool initialize(const meshing::Mesh& mesh,
                    solver::FieldRegistry& fields,
                    const MonteCarloProps& props);

    void add_wall(McWallSpec w) { walls_.push_back(w); }

    /// Performs one MC pass.  Writes additive "S_rad" (W/m³) and
    /// stores per-wall-face absorbed flux internally (queryable).
    /// Returns the total radiative energy emitted this pass.
    double run();

    double wall_face_flux(std::size_t faceIdx) const;
    std::size_t face_count() const noexcept { return q_wall_.size(); }

    const MonteCarloProps& props() const noexcept { return p_; }

private:
    const McWallSpec* find_wall(meshing::ZoneId z) const;

    const meshing::Mesh* mesh_ = nullptr;
    solver::FieldRegistry* F_ = nullptr;
    MonteCarloProps p_{};
    std::vector<McWallSpec> walls_;
    std::vector<double> q_wall_;
    std::mt19937_64 rng_;
};

} // namespace simall::radiation
