// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/ScalarTransport.hpp
// Phase  : 6.4 / 9 — generic implicit scalar transport equation.
//
//   ∂(ρφ)/∂t + ∇·(ρ U φ) = ∇·(Γ ∇φ) + S_φ
//
// Used by:
//   - Energy equation       φ = T  Γ = k/cp        S = viscous + radiation
//   - Species mass fraction φ = Yk Γ = ρ D_k
//   - Passive scalars       any φ
//
// Operates on the same polyhedral Mesh, FieldRegistry, ILinearSolver and
// boundary-condition table as the SIMPLE momentum/pressure assembly. Fully
// 3-D, supports arbitrary cells. Convection scheme is upwind by default;
// second-order schemes plug in through the SpatialScheme enum.
// =============================================================================
#pragma once

#include "CSRMatrix.hpp"
#include "FieldRegistry.hpp"
#include "Solver.hpp"

#include "meshing/MeshStorage.hpp"

#include <string>
#include <vector>

namespace simall::solver
{

struct ScalarBC
{
    meshing::ZoneId zone;
    enum class Kind
    {
        Dirichlet,
        Neumann,
        Robin
    } kind = Kind::Dirichlet;
    double value = 0.0;  // φ for Dirichlet, flux for Neumann, h_ref for Robin
    double valueB = 0.0; // φ_∞ for Robin
};

class ScalarTransport
{
public:
    ScalarTransport(meshing::Mesh& mesh, FieldRegistry& fields, ILinearSolver& linear);

    void set_field(std::string name) { name_ = std::move(name); }
    void set_diffusivity(double gamma) { gamma_ = gamma; }
    void set_density(double rho) { rho_ = rho; }
    void set_source(const std::string& srcField) { sourceField_ = srcField; }
    void add_bc(ScalarBC bc) { bcs_.push_back(bc); }
    /// Replace an existing zone's BC in place (or append if not present).
    /// Used by coupling drivers (CHT, FSI) to swap interface conditions
    /// between Dirichlet and Neumann between outer iterations.
    void set_zone_bc(ScalarBC bc);
    void set_scheme(SpatialScheme s) { scheme_ = s; }
    void set_urf(double urf) { urf_ = urf; }

    /// Returns the L2 residual norm of (A x - b).
    double solve_iteration();

private:
    void build_sparsity();

    meshing::Mesh& mesh_;
    FieldRegistry& F_;
    ILinearSolver& lin_;
    std::vector<ScalarBC> bcs_;
    std::string name_ = "T";
    std::string sourceField_;
    double gamma_ = 1.0;
    double rho_ = 1.0;
    double urf_ = 0.7;
    SpatialScheme scheme_ = SpatialScheme::FirstOrderUpwind;
    CSRMatrix A_;
    bool sparsity_built_ = false;
    util::aligned_vector<double> rhs_;
};

} // namespace simall::solver
