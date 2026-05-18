// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/DiscreteOrdinates.hpp
// Phase  : 11 — Discrete Ordinates Method (DOM) with S_N quadrature.
//
// Solves the Radiative Transfer Equation (RTE) for grey participating media:
//
//   s · ∇I(r,s) + (κ + σ_s) I = κ I_b + (σ_s / 4π) ∫ I(r,s') Φ(s,s') dΩ'
//
// I_b = σ T^4 / π  (blackbody intensity).  With isotropic scattering and a
// grey assumption the in-scatter integral reduces to (σ_s / 4π) Σ_m w_m I_m.
//
// Quadrature: Level-Symmetric LSn (S4 → 24 directions, S6 → 48 directions).
// Spatial discretisation: cell-centred finite volume with first-order upwind
// along each ordinate (Fiveland 1984). Iterative source-iteration / Jacobi
// sweeping; for a single grey gas this converges in O(10) sweeps.
//
// Boundary handling: diffusely emitting/reflecting walls
//   I(s) = ε σ T_w^4 / π + (1-ε)/π · ∫_{n·s'<0} |n·s'| I(s') dΩ'
//
// Outputs:
//   - cell-centred mean intensity G = Σ_m w_m I_m
//   - radiative heat-source S_rad = κ (G - 4 σ T^4)  (added to energy eq)
//
// Fully 3-D, arbitrary polyhedral mesh.
// =============================================================================
#pragma once

#include "solver/FieldRegistry.hpp"
#include "solver/Solver.hpp"
#include "meshing/MeshStorage.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <array>
#include <vector>

namespace simall::radiation {

enum class QuadratureOrder { S4 = 4, S6 = 6 };

struct WallEmission {
    meshing::ZoneId zone;
    double          emissivity = 1.0;
    double          temperature = 300.0;
};

class DiscreteOrdinates {
public:
    DiscreteOrdinates(meshing::Mesh& mesh,
                      solver::FieldRegistry& fields,
                      QuadratureOrder order = QuadratureOrder::S4);

    void set_absorption(double kappa)    { kappa_ = kappa; }
    void set_scattering(double sigmaS)   { sigmaS_ = sigmaS; }
    void set_refractive_index(double n)  { nRef_  = n; }
    void add_wall(WallEmission w)        { walls_.push_back(w); }

    /// One DOM sweep over all ordinates. Returns max change in I.
    double sweep();

    /// Compute radiative source term S_rad per cell (J/m³/s).
    void   compute_source(util::aligned_vector<double>& Srad) const;

    int num_ordinates() const { return static_cast<int>(omega_.size()); }

private:
    meshing::Mesh&             mesh_;
    solver::FieldRegistry&     F_;
    QuadratureOrder            order_;
    double kappa_  = 0.0;
    double sigmaS_ = 0.0;
    double nRef_   = 1.0;
    std::vector<WallEmission>  walls_;

    // Ordinate directions and weights.
    struct Ordinate { double sx, sy, sz, w; };
    std::vector<Ordinate>      omega_;
    // I[m * nCells + c]  — intensity for ordinate m at cell c
    std::vector<util::aligned_vector<double>> I_;

    void build_quadrature();
    const WallEmission* find_wall(meshing::ZoneId z) const;
    void sweep_ordinate(int m, double& maxDelta);
};

}  // namespace simall::radiation
