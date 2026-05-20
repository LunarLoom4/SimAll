// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/SolarLoad.hpp
// Phase  : 11.7 — Solar load model (direct beam + diffuse) for exterior
//          surfaces (building envelopes, vehicle thermal management).
//
// Direct beam:
//   q_dir,i = I_dir · max( -ŝ · n̂_i, 0 ) · (1 - α_i)
//
// Sky-diffuse (isotropic Liu-Jordan):
//   q_dif,i = I_dif · F_sky,i · (1 - α_i)
// where F_sky,i = (1 + n̂_i·ẑ) / 2  (face-to-sky view factor).
//
// Ground-reflected:
//   q_grd,i = ρ_grd · (I_dir·cosθ_z + I_dif) · F_grd,i · (1 - α_i)
//   F_grd,i = (1 - n̂_i·ẑ) / 2
//
// Sun direction ŝ is provided either explicitly (props.sunDir, must point
// FROM the sun) or computed via ASHRAE solar geometry from latitude,
// day-of-year, hour (ŝ recomputed each call from `time_seconds`).
//
// Output: per-face absorbed flux written to "q_solar_face" (W/m²) and
// optionally added as a boundary heat flux on the named zones.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

#include <vector>

namespace simall::radiation
{

struct SolarFaceSpec
{
    meshing::ZoneId zone;
    double absorptivity = 0.6;
};

struct SolarProps
{
    util::Vec3d sunDir{-0.5, -0.5, -0.7071}; // FROM sun toward surface
    double I_direct = 900.0;                 // W/m² beam normal
    double I_diffuse = 100.0;                // W/m² isotropic diffuse
    double I_total_max = 1367.0;             // solar constant cap
    double groundReflectance = 0.2;
    bool useAshrae = false; // if true, override sunDir from time
    double latitude_deg = 0.0;
    double longitude_deg = 0.0;
    double timezone_hours = 0.0;
};

class SolarLoad
{
public:
    bool initialize(const meshing::Mesh& mesh,
                    solver::FieldRegistry& fields,
                    const SolarProps& props);

    void add_zone(SolarFaceSpec z) { zones_.push_back(z); }

    /// Apply solar load: for each registered exterior face compute the
    /// absorbed direct + diffuse + ground-reflected flux.  `time_seconds`
    /// is used only when props.useAshrae=true.  Returns sum-flux (W).
    double apply(double time_seconds = 0.0);

    double face_flux(std::size_t faceIdx) const;
    std::size_t face_count() const noexcept { return q_face_.size(); }

    const SolarProps& props() const noexcept { return p_; }

private:
    const SolarFaceSpec* spec_for_zone(meshing::ZoneId z) const;
    util::Vec3d ashrae_sun_dir(double t) const;

    const meshing::Mesh* mesh_ = nullptr;
    solver::FieldRegistry* F_ = nullptr;
    SolarProps p_{};
    std::vector<SolarFaceSpec> zones_;
    std::vector<double> q_face_;
};

} // namespace simall::radiation
