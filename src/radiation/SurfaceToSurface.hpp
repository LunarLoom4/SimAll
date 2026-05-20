// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/SurfaceToSurface.hpp
// Phase  : 11.6 — Surface-to-Surface (S2S) view-factor radiation exchange.
//
// Models radiative exchange between diffuse-grey opaque surfaces of an
// enclosure with NON-participating medium.  Per the radiosity equation:
//
//   J_i = ε_i σ T_i⁴ + (1-ε_i) Σ_j F_ij J_j
//   q_i = ε_i / (1-ε_i) · (σ T_i⁴ - J_i)        (for ε_i < 1)
//   q_i = σ T_i⁴ - Σ_j F_ij J_j                 (for ε_i = 1, black)
//
// View factors F_ij are computed by hemisphere ray-sampling between
// face centroids using the formula
//
//   F_ij = (1 / (π A_i)) ∫∫ (cos θ_i · cos θ_j / r²) dA_i dA_j
//
// approximated with stratified Monte-Carlo sampling (Nray rays per face).
// Shadowing accounted for via a simple all-pair LOS test against the
// mesh-wall-face plane (cheap-but-conservative; flagged TODO for BVH).
//
// Output: per-face heat flux q_face (W/m²) written into "q_rad_face"
// scratch vector accessible via face_flux(i); volumetric S_rad sourcing
// is delegated to a host coupler when needed.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <cstdint>
#include <vector>

namespace simall::radiation
{

struct S2SFaceSpec
{
    meshing::ZoneId zone;
    double emissivity = 0.9;
    double temperature = 300.0;
};

struct S2SProps
{
    std::size_t nRaysPerFace = 256;
    bool checkOcclusion = false;
    std::uint64_t rngSeed = 0xD15E'A5E5ULL;
};

class SurfaceToSurface
{
public:
    bool initialize(const meshing::Mesh& mesh, const S2SProps& props);

    void add_zone(S2SFaceSpec z) { zones_.push_back(z); }

    /// Build the F_ij matrix for every face belonging to a registered zone.
    /// Returns the number of participating faces.
    std::size_t build_view_factors();

    /// Solve the radiosity system and write per-face fluxes (W/m²).
    /// Returns max |q_i|.
    double solve();

    double face_flux(std::size_t faceIdx) const;
    double view_factor(std::size_t i, std::size_t j) const;
    std::size_t participating_face_count() const noexcept { return faceList_.size(); }

    const S2SProps& props() const noexcept { return p_; }

private:
    const S2SFaceSpec* spec_for_zone(meshing::ZoneId z) const;

    const meshing::Mesh* mesh_ = nullptr;
    S2SProps p_{};
    std::vector<S2SFaceSpec> zones_;
    std::vector<std::size_t> faceList_; // index → mesh-face id
    std::vector<double> F_;             // row-major nF·nF view factor
    std::vector<double> q_;             // per face (W/m²)
};

} // namespace simall::radiation
