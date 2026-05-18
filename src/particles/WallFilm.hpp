// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/WallFilm.hpp
// Phase  : 13.16 — Wall-film model (Stanton & Rutland 1996,
//          Bai-Gosman 1995 impingement regimes).
//
// Thin liquid film is tracked per *wall face* via:
//
//   ∂(ρ_l h_f)/∂t + ∇·(ρ_l h_f U_f) = ṁ_imp - ṁ_evap
//
// Impingement regime determined by impact Weber number:
//
//   We = ρ_l |v_n|² d_p / σ
//
// Bai-Gosman thresholds (cold/wet):
//   We < 5   : stick
//   We < 10  : spread (full deposition into film, no rebound)
//   We < 18  : rebound (no film mass added)
//   else     : splash (Wang/Mundo: ε_s ≈ min(1, k(We-We_c)/We) re-injected
//                       as secondary parcels; the rest sticks)
//
// Film advection uses a simple explicit upwind on face neighbours (faces
// sharing nodes).  Evaporation is delegated to a host model — here we only
// track film thickness, mass and momentum.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "particles/LagrangianTracker.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <vector>

namespace simall::particles {

struct WallFilmProps {
    double rho_l   = 998.2;
    double mu_l    = 1.0e-3;
    double sigma   = 0.072;
    double h_min   = 1.0e-9;
    double h_max   = 5.0e-3;
    double Weks    = 5.0;
    double Wespd   = 10.0;
    double Wereb   = 18.0;
    double splash_ratio = 0.6;       // mass fraction splashed when We > Wereb
    std::vector<std::int32_t> wallZones;  // boundary zone IDs treated as walls
};

class WallFilm {
public:
    void initialize(const meshing::Mesh& mesh, WallFilmProps props);

    /// Accept parcels which crossed into a wall face this step.
    /// Updates film mass/thickness and (optionally) deactivates absorbed
    /// parcels.  Returns total deposited mass.
    double deposit(double dt, LagrangianTracker& tracker);

    /// Apply explicit face-based advection of film mass under shear.
    void   advect(double dt);

    double film_thickness(std::size_t faceIdx) const;
    double film_mass(std::size_t faceIdx) const;
    std::size_t face_count() const noexcept { return h_.size(); }

    const WallFilmProps& props() const noexcept { return p_; }

private:
    bool is_wall_face(std::int32_t zone) const;

    const meshing::Mesh* mesh_ = nullptr;
    WallFilmProps        p_{};
    std::vector<double>  h_;        // film thickness per face
    std::vector<double>  m_;        // film mass per face
    std::vector<util::Vec3d> u_;    // film velocity per face
};

}  // namespace simall::particles
