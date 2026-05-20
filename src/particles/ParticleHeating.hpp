// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ParticleHeating.hpp
// Phase  : 13.14 — Convective + radiative particle heating.
//
//   dT_p/dt = (Q_conv + Q_rad) / (m_p c_pl)
//
// where:
//   Q_conv = h_c · A_p · (T_g - T_p)
//   h_c    = Nu k_g / d_p,  Nu = 2 + 0.6 Re^{1/2} Pr^{1/3}   (Ranz-Marshall)
//   Q_rad  = ε σ A_p (T_R^4 - T_p^4)        (optional, requires radiation T_R)
//
// Energy is removed from the surrounding fluid via the cell field
// "S_particle_energy" (W/m³).  Per-parcel temperature is stored as a
// field "T_parcel" sized to track-particle indexing (matches the
// LagrangianTracker particle ordering).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "particles/LagrangianTracker.hpp"
#include "solver/FieldRegistry.hpp"

#include <vector>

namespace simall::particles
{

struct HeatingProps
{
    double cp_l = 4186.0; // [J/(kg·K)] water default
    double k_g = 2.6e-2;  // gas thermal conductivity W/(m·K)
    double cp_g = 1006.0;
    double mu_g = 1.81e-5;
    double rho_g_ref = 1.225;
    double T_p_init = 300.0;
    double T_p_min = 200.0;
    double T_p_max = 5000.0;
    double emissivity = 0.85;
    bool includeRadiation = false;
};

class ParticleHeating
{
public:
    void initialize(const meshing::Mesh& mesh, HeatingProps props);

    /// Reads gas "T" and (optionally) "T_radiation" cell fields; writes
    /// additive cell field "S_particle_energy" (W/m³).  Updates parcel
    /// temperatures stored internally (vector aligned with tracker).
    double apply(double dt, LagrangianTracker& tracker, solver::FieldRegistry& fields);

    double parcel_temperature(std::size_t i) const;
    void set_parcel_temperature(std::size_t i, double T);
    void resize_to(std::size_t nParticles);

    const HeatingProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    HeatingProps p_{};
    std::vector<double> T_p_;
};

} // namespace simall::particles
