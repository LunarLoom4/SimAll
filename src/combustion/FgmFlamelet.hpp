// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/FgmFlamelet.hpp
// Phase  : 11.12 — Flamelet-Generated Manifold (FGM) model.
//
// Tabulates 1-D laminar premixed flamelet solutions, parameterised by:
//   Z (mixture fraction) and c (progress variable, e.g. Y_CO2 + Y_H2O).
//
// Look-up provides ψ(Z, c) for ψ ∈ {ρ, T, ω̇_c, μ, Y_k(species), ...} so
// the 3-D CFD solver only transports Z, Z'', c, c'' (Favre-mean and
// variance) and queries the manifold for everything else.
//
// Reference:
//   van Oijen & de Goey, "Modelling of premixed laminar flames using
//   flamelet-generated manifolds", Combust. Sci. Tech. 161, 113-137 (2000).
//
// Beta-PDF for Z subgrid distribution; delta-PDF for c.  Variance closure
// uses standard k-ε model with C_Z = 2, C_d = 2.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace simall::combustion {

struct FgmTableAxis {
    std::vector<double> Z;     // mixture fraction axis (sorted ascending)
    std::vector<double> c;     // progress variable axis (sorted ascending)
};

struct FgmTable {
    FgmTableAxis axis;
    // table[var][iZ * Nc + ic] flat storage
    std::unordered_map<std::string, std::vector<double>> var;
};

class FgmFlamelet {
public:
    /// Build a uniform analytic FGM table for a simple methane-air flame
    /// (Burke-Schumann limit with Damköhler-style progress variable).
    /// Useful for verification and unit testing without external data.
    void build_analytic(std::size_t NZ = 21, std::size_t Nc = 21,
                        double Z_st = 0.055, double T_un = 300.0,
                        double T_ad = 2225.0, double rho_un = 1.18,
                        double rho_ad = 0.18);

    /// Replace the table with externally-supplied data (e.g. CANTERA output).
    void set_table(FgmTable table) { table_ = std::move(table); }
    const FgmTable& table() const noexcept { return table_; }

    /// Bilinear (Z, c) interpolation with end-point clamping.
    double lookup(const std::string& var, double Z, double c) const;

    /// Initialise mesh-bound progress / mixture fields and reaction-rate
    /// scratch field. Call once.
    void initialize(const meshing::Mesh& mesh, solver::FieldRegistry& fields);

    /// Per-cell evaluation: reads "Z" and "c"; writes "T", "rho",
    /// "S_progress" (= ω̇_c) and "S_combustion" (Q from latent ΔH).
    void apply(solver::FieldRegistry& fields);

private:
    FgmTable table_;
    const meshing::Mesh* mesh_ = nullptr;
};

}  // namespace simall::combustion
