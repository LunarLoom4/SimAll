// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/BlobInjector.hpp
// Phase  : 13.11 — Blob (large-droplet) injection model (Reitz 1987).
//
// Diesel-type sprays where atomisation occurs *downstream* of the nozzle
// are seeded with parcels whose initial diameter equals the nozzle
// diameter ("blobs").  Each blob then undergoes secondary breakup (KH-RT
// in src/particles/SprayBreakup.hpp).
//
// Parcel velocity magnitude derives from the discharge coefficient:
//
//   U_inj = C_d · √(2 (p_inj - p_amb) / ρ_l)
//
// Mass flow rate ṁ = ρ_l · U_inj · A_nozzle · C_d / sqrt-correction is
// distributed uniformly over the injected parcels per time step.
//
// Reference:
//   Reitz, R.D., "Modeling Atomization Processes in High-Pressure Vaporizing
//   Sprays", Atomisation and Spray Technology, 3, 309-337 (1987).
// =============================================================================
#pragma once

#include "particles/LagrangianTracker.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <random>

namespace simall::particles {

struct BlobInjectorProps {
    double      nozzle_diameter = 2.0e-4;   // [m]
    double      Cd              = 0.8;      // discharge coefficient
    double      p_inj           = 1.0e8;    // injection pressure [Pa]
    double      p_amb           = 4.0e6;    // ambient (chamber) pressure [Pa]
    double      rho_l           = 830.0;    // diesel-like density
    double      duration        = 1.5e-3;   // injection duration window [s]
    util::Vec3d origin{0,0,0};
    util::Vec3d axis  {1,0,0};
    std::uint64_t rngSeed       = 0xBEEFCAFE;
};

class BlobInjector {
public:
    void initialize(BlobInjectorProps props);

    /// Inject `Nparcels` blob parcels — but only while elapsed time
    /// `elapsed` ∈ [0, duration].
    /// Returns the number of parcels actually injected (0 if window closed).
    std::size_t inject(double dt, std::size_t Nparcels,
                       double elapsed, LagrangianTracker& tracker);

    double injection_velocity() const noexcept;

    const BlobInjectorProps& props() const noexcept { return p_; }

private:
    BlobInjectorProps p_{};
    std::mt19937_64   rng_;
};

}  // namespace simall::particles
