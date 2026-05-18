// =============================================================================
// SimAll Beta - Electromagnetics Subsystem
// File   : src/emag/Mhd.hpp
// Phase  : 19 — Single-fluid resistive magnetohydrodynamics in the
// low-magnetic-Reynolds approximation suitable for liquid-metal stirring,
// MHD pumps, induction furnaces.
//
// Governing equations augment the incompressible Navier-Stokes:
//
//   Momentum  : ρ Du/Dt = -∇p + ∇·τ + J × B          (Lorentz force)
//   Charge    : ∇·J = 0,           J = σ ( -∇φ + u × B )
//
// where φ is the electric scalar potential (Poisson problem) and B is the
// applied (externally specified) magnetic field. The induced magnetic
// field is neglected (small Rm).
//
// Outputs written into FieldRegistry:
//   - "S_Lorentz" (vector) — body force per unit volume = J × B
//   - "phi_e"     (scalar) — electric potential
//   - "J_e"       (vector) — current density
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/CSRMatrix.hpp"
#include "solver/FieldRegistry.hpp"
#include "solver/LinearSolvers.hpp"
#include "utilities/MathTypes.hpp"

#include <functional>
#include <string>

namespace simall::emag {

struct MhdProps {
    double sigma     = 1.0e6;             // electrical conductivity [S/m]
    // Externally specified B(x) as a callable, evaluated per cell centroid.
    std::function<util::Vec3d(const util::Vec3d&)> Bfield;
    // Voltage boundary conditions on insulator/electrode zones.
    struct PhiBC { meshing::ZoneId zone; double value; bool isDirichlet; };
    std::vector<PhiBC> bcs;
};

class Mhd {
public:
    Mhd(meshing::Mesh& mesh, solver::FieldRegistry& fields,
        solver::ILinearSolver& linear);

    void initialize(MhdProps props);

    /// One solve of (i) potential Poisson then (ii) J = σ(-∇φ + u×B) then
    /// (iii) Lorentz force into "S_Lorentz". Returns the L2 residual of
    /// the potential solve.
    double solve_iteration();

private:
    void build_sparsity();

    meshing::Mesh&         mesh_;
    solver::FieldRegistry& F_;
    solver::ILinearSolver& lin_;
    MhdProps               props_{};
    solver::CSRMatrix      A_;
    util::aligned_vector<double> rhs_;
    bool                   sparsity_built_ = false;
};

}  // namespace simall::emag
