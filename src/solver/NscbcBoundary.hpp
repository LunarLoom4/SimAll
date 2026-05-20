// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/NscbcBoundary.hpp
// Phase  : 7.6 — Navier-Stokes Characteristic Boundary Conditions for the
// compressible Euler/Navier-Stokes system (Poinsot & Lele JCP 1992).
//
// The five characteristic wave amplitudes at a boundary face are
//
//   L1 = (u_n - c) (∂p/∂x_n - ρ c ∂u_n/∂x_n)
//   L2 = u_n      (c² ∂ρ/∂x_n - ∂p/∂x_n)
//   L3 = u_n      (∂u_t1/∂x_n)
//   L4 = u_n      (∂u_t2/∂x_n)
//   L5 = (u_n + c) (∂p/∂x_n + ρ c ∂u_n/∂x_n)
//
// At a subsonic non-reflecting outflow we replace L1 (incoming) by the
// LODI relaxation L1 = K (p − p_∞),  K = σ (1 − M²) c / L.
// At a subsonic inflow we prescribe T, u, u_t and back out L5.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

#include <vector>

namespace simall::solver
{

enum class NscbcType
{
    SubsonicOutflow,   ///< relax pressure to p∞
    SubsonicInflow,    ///< prescribe T, u, v, w
    NonReflectingWall, ///< slip with characteristic update
};

struct NscbcZoneSpec
{
    meshing::ZoneId zone = 0;
    NscbcType type = NscbcType::SubsonicOutflow;
    double pInf = 101325.0;
    double TInf = 300.0;
    util::Vec3d uInf{0, 0, 0};
    double sigma = 0.25; ///< relaxation factor (0.15-0.3 typical)
    double length = 1.0; ///< reference domain length L
};

class NscbcBoundary
{
public:
    void initialize(const meshing::Mesh& mesh,
                    std::vector<NscbcZoneSpec> zones,
                    double gamma = 1.4,
                    double Rgas = 287.0);

    /// Apply LODI updates by writing time-derivative source contributions
    /// into the compressible solver's "rho", "rhoU", "rhoE" fields.
    /// Caller is responsible for time-integrating those sources.
    /// Returns the number of boundary faces processed.
    std::size_t apply(FieldRegistry& F);

private:
    const meshing::Mesh* mesh_ = nullptr;
    std::vector<NscbcZoneSpec> zones_;
    double gamma_ = 1.4;
    double Rgas_ = 287.0;
};

} // namespace simall::solver
