// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ParcelCoalescence.hpp
// Phase  : 13.15 — Stochastic parcel collision / coalescence (O'Rourke 1981).
//
// For each pair of parcels (i, j) sharing a host cell, evaluate
//
//   ν_coll = π (r_i + r_j)² · |v_i - v_j| · N_j / V_cell
//   P_coll = 1 - exp(-ν_coll · dt)
//
// where N_j is the *number of physical droplets* in parcel j.  Sampling a
// uniform U ∈ [0,1]; if U < P_coll, decide outcome via Weber criterion:
//
//   We = ρ_l |v_i - v_j|² r_smaller / σ
//   coalesce  if We < We_crit (≈ 5)
//   bounce / separate otherwise (momentum conserved, sizes unchanged)
//
// On coalescence: r_new = (r_i³ + r_j³)^{1/3}, mass merged into the larger
// parcel; smaller parcel deactivated.  Mass-weighted velocity averaging.
//
// Pairing strategy: for each cell with k parcels, evaluate the ⌈k/2⌉
// random pair sweeps (O'Rourke's algorithm, O(k·log k) average).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "particles/LagrangianTracker.hpp"

#include <cstdint>
#include <random>

namespace simall::particles {

struct CoalescenceProps {
    double sigma   = 0.072;     // [N/m] liquid surface tension
    double rho_l   = 998.2;
    double We_crit = 5.0;
    std::uint64_t rngSeed = 0xC0AE'CE99;
};

class ParcelCoalescence {
public:
    void initialize(const meshing::Mesh& mesh, CoalescenceProps props);

    /// Walks all cells, runs the stochastic O'Rourke sweep per cell.
    /// Returns the number of coalescence events that occurred.
    std::size_t apply(double dt, LagrangianTracker& tracker);

    const CoalescenceProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    CoalescenceProps     p_{};
    std::mt19937_64      rng_;
};

}  // namespace simall::particles
