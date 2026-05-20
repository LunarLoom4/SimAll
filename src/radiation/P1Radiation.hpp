// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/P1Radiation.hpp
// Phase  : 11 — P1 (first-order spherical harmonics) approximation for the
// Radiative Transfer Equation in grey participating media.
//
// PDE solved (Modest, "Radiative Heat Transfer", 3rd ed., §16.4):
//
//      −∇·(Γ ∇G) + a G = 4 π a I_b ,    Γ = 1 / (3(a + σ_s − A σ_s))
//
// with the Marshak diffuse-grey wall boundary condition
//
//      −Γ ∂G/∂n = ε_w / (2 (2 − ε_w)) · (4 σ T_w⁴ − G)
//
// I_b = σ T⁴ / π (n=1). Radiative heat source added to energy equation:
//
//      S_rad = a (G − 4 σ T⁴)
//
// O(N) cost per nonlinear iteration (single Helmholtz solve). Much cheaper
// than DO for optically-thick media.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "solver/Solver.hpp"

#include <memory>
#include <vector>

namespace simall::solver
{
class ILinearSolver;
}

namespace simall::radiation
{

struct P1WallBC
{
    meshing::ZoneId zone;
    double emissivity = 1.0;
    double temperature = 300.0;
};

struct P1Props
{
    double absorption = 0.5;  // a   [1/m]
    double scattering = 0.0;  // σ_s [1/m]
    double asymmetry = 0.0;   // A ∈ [-1,1], linear-aniso scattering
    double refractiveN = 1.0; // (kept for future spectral extensions)
};

class P1Radiation
{
public:
    bool initialize(const meshing::Mesh& mesh, solver::FieldRegistry& fields, const P1Props& props);

    void add_wall(P1WallBC w) { walls_.push_back(w); }

    /// One Helmholtz solve. Reads "T"; writes "G" and "S_rad".
    /// Returns L2 residual reported by the linear solver.
    double solve();

private:
    const meshing::Mesh* mesh_ = nullptr;
    solver::FieldRegistry* F_ = nullptr;
    P1Props p_{};
    std::vector<P1WallBC> walls_;
    std::unique_ptr<solver::ILinearSolver> lin_;
};

} // namespace simall::radiation
