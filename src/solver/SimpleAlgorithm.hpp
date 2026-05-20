// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/SimpleAlgorithm.hpp
// Phase  : 6.3 — SIMPLE (Semi-Implicit Method for Pressure-Linked Equations).
//
// Cell-centred collocated finite-volume implementation for incompressible
// flow. Rhie-Chow interpolation prevents the checker-board mode on the
// collocated grid. The implementation supports arbitrary polyhedral meshes
// produced by the meshing subsystem (face owner/neighbour with CSR cell→face).
//
// Public entry points:
//   - SimpleAlgorithm::assemble_momentum(component)
//   - SimpleAlgorithm::solve_momentum()
//   - SimpleAlgorithm::assemble_pressure_correction()
//   - SimpleAlgorithm::solve_and_correct()
//   - SimpleAlgorithm::iterate() → wraps the above into one outer iteration
//
// Boundary conditions are consumed from a BoundaryConditionTable that is
// populated by Solver::add_boundary().
// =============================================================================
#pragma once

#include "CSRMatrix.hpp"
#include "FieldRegistry.hpp"
#include "Solver.hpp"

#include "meshing/MeshStorage.hpp"
#include "parallel/IFieldSynchronizer.hpp"

#include <vector>

namespace simall::solver
{

/// Pressure-velocity coupling variant used by SimpleAlgorithm and its
/// subclasses.  Set in SimpleOptions::algorithm.
enum class PvCouplingVariant
{
    SIMPLE,  // Patankar (1980)              — under-relax p, urfP ≈ 0.3
    SIMPLEC, // Van Doormaal-Raithby (1984)  — denom = (a_P - Σ a_NB), urfP ≈ 1.0
};

struct SimpleOptions
{
    double urfU = 0.7;  // momentum under-relaxation
    double urfP = 0.3;  // pressure under-relaxation
    double rho = 1.0;   // density (incompressible)
    double mu = 1.0e-3; // dynamic viscosity
    SpatialScheme scheme = SpatialScheme::FirstOrderUpwind;
    // Transient (Phase 6.2). dt > 0 enables time-accurate integration. The
    // first physical step degrades to implicit Euler automatically.
    double dt = 0.0;
    TemporalScheme timeScheme = TemporalScheme::ImplicitEuler;
    /// Coupling-equation variant. PISO is a separate class (PisoAlgorithm)
    /// because it changes the iteration *structure*, not just a coefficient.
    PvCouplingVariant algorithm = PvCouplingVariant::SIMPLE;
};

struct SimpleResiduals
{
    double mom[3] = {0, 0, 0};
    double cont = 0.0;
};

class SimpleAlgorithm
{
public:
    SimpleAlgorithm(meshing::Mesh& mesh,
                    FieldRegistry& fields,
                    const std::vector<BoundarySpec>& boundaries,
                    ILinearSolver& momentumSolver,
                    ILinearSolver& pressureSolver,
                    SimpleOptions opts);
    virtual ~SimpleAlgorithm() = default;

    /// Perform one full SIMPLE outer iteration.
    /// Returns residuals (normalized L2 of the unscaled system).
    virtual SimpleResiduals iterate();

    /// Advance one physical time step (calls iterate() N times then rolls
    /// the U^{n}, U^{n-1} history forward). Use when opt.dt > 0.
    SimpleResiduals advance_time_step(int nInnerIters = 5);

    /// Inject an externally-computed periodic-face twin array (built e.g.
    /// by build_periodic_pairs() from PeriodicPairing.hpp). When non-empty
    /// this overrides the built-in translational pairing heuristic and is
    /// the only path that supports rotational periodicity.
    void set_periodic_pairs(std::vector<int> twin)
    {
        externalTwin_ = std::move(twin);
        externalTwinSet_ = true;
    }

    /// Install a halo-exchange driver (e.g. parallel::GhostExchange).  When
    /// non-null, ghost cells of U and p are updated after every momentum
    /// solve and pressure correction so that face-flux computations see
    /// consistent neighbour values across MPI rank boundaries.  Passing
    /// nullptr disables the hook and reverts to single-rank behaviour.
    void set_field_synchronizer(parallel::IFieldSynchronizer* sync) noexcept { sync_ = sync; }
    parallel::IFieldSynchronizer* field_synchronizer() const noexcept { return sync_; }

    /// Direct read-access for diagnostics and for subclasses / drivers.
    const SimpleOptions& options() const noexcept { return opt_; }
    SimpleOptions& options() noexcept { return opt_; }

protected:
    // --- internal hooks (protected so PisoAlgorithm / SimplecAlgorithm can reuse) ---
    void assemble_momentum(int component);
    void assemble_pressure_correction();
    void apply_velocity_bcs(int component, util::aligned_vector<double>& rhs);
    void apply_pressure_bcs(util::aligned_vector<double>& rhs);
    void compute_face_fluxes();
    void correct_fields();
    void build_periodic_pairs();
    void build_sparsity();

    /// Solve the most recently assembled momentum component and write the
    /// solution back into U.<component>. Returns L2 residual norm.
    double solve_momentum_component(int component);

    /// Solve the most recently assembled pressure-correction equation and
    /// apply the field correction.  Returns the continuity residual.
    double solve_pressure_correction_and_correct();

    meshing::Mesh& mesh_;
    FieldRegistry& F_;
    const std::vector<BoundarySpec>& bcs_;
    ILinearSolver& linMom_;
    ILinearSolver& linP_;
    SimpleOptions opt_;

    // CSR storage shared between momentum components (sparsity is identical).
    CSRMatrix A_mom_;
    CSRMatrix A_p_;
    // Per-cell scratch
    util::aligned_vector<double> aP_;      // momentum diagonal (for Rhie-Chow)
    util::aligned_vector<double> phiFlux_; // mass flux per face
    util::aligned_vector<double> rhs_;
    util::aligned_vector<double> sol_;
    util::aligned_vector<double> pPrime_;
    // BDF2 history (Phase 6.2)
    util::aligned_vector<double> Ux_n_, Uy_n_, Uz_n_;       // current step start
    util::aligned_vector<double> Ux_nm1_, Uy_nm1_, Uz_nm1_; // previous step
    int timeStep_ = 0;
    // Periodic-pair table: for each boundary face f tagged Periodic, the
    // index of its paired face (or -1 if not periodic).
    std::vector<int> periodicTwin_;
    std::vector<int> externalTwin_;
    bool externalTwinSet_ = false;
    bool sparsity_built_ = false;
    // Optional ghost-cell exchange driver for MPI runs (W12).  Non-owning.
    parallel::IFieldSynchronizer* sync_ = nullptr;
};

} // namespace simall::solver
