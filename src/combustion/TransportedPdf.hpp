// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/TransportedPdf.hpp
// Phase  : 11.14 — Transported composition PDF (T-PDF) model.
//
// Solves a Monte Carlo realisation of the transport equation for the
// Favre-mass-density-weighted joint PDF f̃(ψ; x, t) of N_s composition
// scalars Y_k (and optionally h_s).  Notional particles are advanced by:
//
//   1. Convection by the Favre-mean velocity Ũ(x).
//   2. Random walk by isotropic turbulent diffusion D_t = μ_t/(ρ Sc_t).
//   3. Mixing model (Curl 1963 or IEM/LMSE - Villermaux 1972).
//   4. Chemistry step (Arrhenius source ω̇_k integrated by EE-implicit).
//
// Each cell c keeps Np particles ψ^{(c,p)}; cell-averaged quantities are
// reconstructed by simple ensemble averages.
//
// References:
//   Pope, "PDF methods for turbulent reactive flows", Prog. Energy Combust.
//   Sci. 11, 119–192 (1985).
//   Subramaniam & Pope, Combust. Flame 115, 487-514 (1998) — IEM/Curl.
// =============================================================================
#pragma once

#include "combustion/ChemkinParser.hpp"
#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <random>
#include <vector>

namespace simall::combustion
{

enum class PdfMixingModel
{
    IEM,
    ModifiedCurl
};

struct TransportedPdfParams
{
    std::size_t particlesPerCell = 50;
    PdfMixingModel mixing = PdfMixingModel::IEM;
    double C_phi = 2.0;          // mixing-time constant
    double Sc_t = 0.7;           // turbulent Schmidt
    double T_reference = 1500.0; // K  (chemistry linearisation)
    std::uint64_t rngSeed = 0xC0DE'BEEFULL;
};

struct PdfParticle
{
    std::vector<double> Y; // mass fractions (size = nSpecies)
    double h = 0.0;        // sensible enthalpy (optional channel)
    double age = 0.0;      // residence time, diagnostic
};

class TransportedPdf
{
public:
    void initialize(const meshing::Mesh& mesh,
                    const std::vector<std::string>& speciesNames,
                    TransportedPdfParams params = {});

    /// Advance Np particles in every cell by dt:
    ///   convection (deterministic), random walk (Wiener), mixing (IEM/MC),
    ///   chemistry (one-step LFR with provided reactions).
    /// Writes ensemble-averaged Y_<name> + T to the field registry.
    void step(double dt,
              solver::FieldRegistry& fields,
              const std::vector<ChemkinReaction>& reactions = {});

    std::size_t num_species() const noexcept { return names_.size(); }
    std::size_t num_particles_per_cell() const noexcept { return p_.particlesPerCell; }

private:
    void apply_iem(double dt, std::vector<PdfParticle>& particles, const std::vector<double>& Ybar);
    void apply_modified_curl(double dt, std::vector<PdfParticle>& particles);
    void apply_chemistry(double dt,
                         std::vector<PdfParticle>& particles,
                         const std::vector<ChemkinReaction>& reactions);

    const meshing::Mesh* mesh_ = nullptr;
    TransportedPdfParams p_{};
    std::vector<std::string> names_;
    std::vector<std::vector<PdfParticle>> cellParticles_; // per cell
    std::mt19937_64 rng_;
};

} // namespace simall::combustion
