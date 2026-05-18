// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/SlfmFlamelet.hpp
// Phase  : 11.13 — Steady Laminar Flamelet Model (SLFM) of Peters (1984).
//
// Diffusion flames in the flamelet regime obey, in mixture-fraction space:
//
//   ρ ∂Y_k/∂t  =  (ρ χ / 2) · ∂²Y_k/∂Z²   +   ω̇_k
//
// where χ = 2 D |∇Z|² is the scalar dissipation rate.  At steady state and
// for given χ, the flamelet structure ψ(Z; χ) is tabulated and queried by
// the 3-D CFD solver:  ψ̄ = ∫∫ ψ(Z, χ) P̃(Z) P̃(χ) dZ dχ.
//
// Subgrid Z distribution: β-PDF parameterised by (Z̄, Z̃").
// Subgrid χ distribution: log-normal PDF (σ = 1, default).
//
// The SLFM table maps (Z, log10 χ) → {T, ρ, Y_k, ω̇_k} and is queried via
// bilinear interpolation.  An analytic Burke-Schumann + exponential-decay
// flamelet is provided to make the module standalone for verification.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace simall::combustion {

struct SlfmAxis {
    std::vector<double> Z;          // mixture fraction
    std::vector<double> log10Chi;   // log10 of scalar dissipation rate (1/s)
};

struct SlfmTable {
    SlfmAxis axis;
    std::unordered_map<std::string, std::vector<double>> var;
};

class SlfmFlamelet {
public:
    /// Build an analytic SLFM table: Burke-Schumann triangular T(Z)
    /// modulated by exp(-χ/χ_q) extinction factor.
    void build_analytic(std::size_t NZ      = 41,
                        std::size_t NlogChi = 11,
                        double Z_st = 0.055, double T_un = 300.0,
                        double T_ad = 2225.0,
                        double chi_q = 100.0,
                        double log10chi_min = -3.0,
                        double log10chi_max =  3.0);

    void set_table(SlfmTable t) { table_ = std::move(t); }
    const SlfmTable& table() const noexcept { return table_; }

    /// Bilinear lookup in (Z, log10 χ).
    double lookup(const std::string& var, double Z, double chi) const;

    void initialize(const meshing::Mesh& mesh, solver::FieldRegistry& fields);

    /// Reads "Z", "chi"; writes "T", "rho", "S_combustion".
    void apply(solver::FieldRegistry& fields);

private:
    SlfmTable           table_;
    const meshing::Mesh* mesh_ = nullptr;
};

}  // namespace simall::combustion
