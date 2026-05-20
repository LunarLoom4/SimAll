// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/SootMomMethod.hpp
// Phase  : 11.19 — Method of Moments with Interpolative Closure (MOMIC).
//
// Tracks the first K moments M_r of the soot particle size distribution:
//
//   M_r = Σ_i N_i v_i^r          r = 0, 1, 2, ..., K-1
//
// Population balance after Smoluchowski coagulation + Lindstedt-style
// inception + HACA surface growth + OH/O2 oxidation:
//
//   dM_r/dt  =  R_inc,r  +  R_coag,r  +  R_sg,r  -  R_ox,r
//
// MOMIC closure: fractional moments M_{r+1/2}, M_{r-1/2} approximated by
// log-linear interpolation between adjacent integer moments (Frenklach &
// Harris 1987, Frenklach 2002).
//
// Coagulation kernel uses the free-molecular harmonic mean with continuum
// regime (Kazakov & Frenklach 1998).
//
// Reference:
//   Frenklach, "Method of moments with interpolative closure",
//   Chem. Eng. Sci. 57, 2229-2239 (2002).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <array>
#include <vector>

namespace simall::combustion
{

inline constexpr int MOMIC_K = 4; // number of moments tracked (0..3)

struct MomicProps
{
    double rho_soot = 1800.0;
    double Mw_C2H2 = 26.04e-3; // HACA growth precursor
    double C_inc = 1.0e5;      // nucleation prefactor
    double T_inc = 21000.0;
    double C_sg = 6000.0; // HACA surface growth
    double T_sg = 12100.0;
    double C_ox_O2 = 1.0e3;
    double T_ox_O2 = 19680.0;
    double C_ox_OH = 0.36; // Neoh collision efficiency
};

class SootMomMethod
{
public:
    void initialize(const meshing::Mesh& mesh, MomicProps props = {});

    /// Reads "T", "rho_mix", "Y_C2H2", "Y_O2", "Y_OH", and per-cell
    /// moment fields "M0_soot" .. "M3_soot"; writes "S_M0_soot" ..
    /// "S_M3_soot" additive sources.  Returns volume-integrated dM_0/dt.
    double apply(solver::FieldRegistry& fields);

    /// Interpolative closure: log-linear estimate of fractional moment M_r.
    static double fractional_moment(const std::array<double, MOMIC_K>& M, double r);

    const MomicProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    MomicProps p_{};
};

} // namespace simall::combustion
