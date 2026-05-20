// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/Solver.hpp
// Phase  : 6 / 7 / 16  (finite-volume + linear solvers + boundary conditions)
//
// The Solver framework is PHYSICS-AGNOSTIC. Concrete schemes
// (SIMPLE / SIMPLEC / PISO / Coupled) are pluggable via ISolverScheme.
// Linear systems are solved by ILinearSolver. Both interfaces live here so
// turbulence/combustion/multiphase plugins can target them.
// =============================================================================
#pragma once

#include "CSRMatrix.hpp"
#include "FieldRegistry.hpp"

#include "materials/Material.hpp"
#include "meshing/MeshStorage.hpp"

#include <memory>
#include <string>
#include <vector>

namespace simall::turbulence
{
class ITurbulenceModel;
}

namespace simall::solver
{

// ---------------- Discretization schemes (Phase 6.2) -----------------------
enum class SpatialScheme
{
    FirstOrderUpwind,
    SecondOrderUpwind,
    QUICK,
    MUSCL,
    Central,
    BoundedCentral,
    SuperBeeTVD,
    VanLeerTVD
};
enum class TemporalScheme
{
    ExplicitEuler,
    ImplicitEuler,
    CrankNicolson,
    BDF2,
    DualTime
};

// ---------------- Pressure-velocity coupling (Phase 6.3 + 22.1) ------------
// PIMPLE = transient SIMPLE-outer / PISO-inner hybrid (Issa+Patankar).
//   nOuterCorrectors == 1  → degrades to PISO (no momentum/pressure URF).
//   nOuterCorrectors  > 1  → SIMPLE-style outer loop with URF on each step,
//                            PISO-style inner pressure correctors per outer.
enum class CouplingAlgorithm
{
    SIMPLE,
    SIMPLEC,
    PISO,
    PIMPLE,
    Coupled
};

// ---------------- Linear solver interface (Phase 7) ------------------------
enum class LinearSolverKind
{
    GMRES,
    BiCGSTAB,
    CG,
    TFQMR
};
enum class PreconditionerKind
{
    Jacobi,
    ILU,
    SSOR,
    AMG,
    None
};

struct LinearSolverConfig
{
    LinearSolverKind kind = LinearSolverKind::GMRES;
    PreconditionerKind preconditioner = PreconditionerKind::ILU;
    double tolerance = 1.0e-8;
    int maxIterations = 1000;
    int restart = 30; // GMRES
};

class ILinearSolver
{
public:
    virtual ~ILinearSolver() = default;
    virtual int solve(const CSRMatrix& A,
                      const util::aligned_vector<double>& b,
                      util::aligned_vector<double>& x) = 0;
    virtual double last_residual() const = 0;
};

std::unique_ptr<ILinearSolver> make_linear_solver(LinearSolverConfig cfg);

// ---------------- Boundary conditions (Phase 16) ---------------------------
enum class BCType
{
    Wall,
    NoSlipWall,
    MovingWall,
    SlipWall,
    VelocityInlet,
    PressureInlet,
    MassFlowInlet,
    PressureOutlet,
    OutFlow,
    Symmetry,
    Axisymmetric,
    Periodic,
    Interface,
    PorousJump,
    Fan,
    Overset
};

struct BoundarySpec
{
    meshing::ZoneId zone;
    BCType type;
    double scalarValue = 0.0;
    double vectorValue[3]{0, 0, 0};
    std::string turbulenceProfile;
};

// ---------------- Physics configuration ------------------------------------
struct PhysicsConfig
{
    CouplingAlgorithm coupling = CouplingAlgorithm::SIMPLE;
    SpatialScheme momentumScheme = SpatialScheme::SecondOrderUpwind;
    SpatialScheme scalarScheme = SpatialScheme::SecondOrderUpwind;
    TemporalScheme timeScheme = TemporalScheme::ImplicitEuler;
    bool compressible = false;
    bool energyEquation = true;
    std::string turbulenceModel = "kOmegaSST";
};

// ---------------- Convergence monitor --------------------------------------
struct ResidualSnapshot
{
    int iteration;
    double continuity;
    double momentumX, momentumY, momentumZ;
    double energy;
    double turbulence_k, turbulence_eps;
};

class IResidualMonitor
{
public:
    virtual ~IResidualMonitor() = default;
    virtual void on_residual(const ResidualSnapshot&) = 0;
};

// ---------------- Top-level driver -----------------------------------------
class Solver
{
public:
    Solver(meshing::Mesh& mesh, materials::MaterialDatabase& materials);
    ~Solver();

    void configure(PhysicsConfig cfg);
    void add_boundary(BoundarySpec bc);
    void set_linear_solver(LinearSolverConfig cfg);
    void add_monitor(std::shared_ptr<IResidualMonitor> mon);
    /// Inject a turbulence model. The application owns/creates the model
    /// (e.g. via TurbulenceRegistry::create) and configures its density,
    /// viscosity, and boundary list before passing it here. The Solver then
    /// calls model->initialize / solve on each outer iteration.
    void set_turbulence_model(std::shared_ptr<::simall::turbulence::ITurbulenceModel> m);

    /// Initialize fields (zero / hybrid / user).
    void initialize();

    /// Single outer iteration (SIMPLE/PISO sweep).
    bool step();

    /// Run until convergence or max iterations reached.
    void run(int maxOuterIter, double residualTarget);

    FieldRegistry& fields() noexcept { return fields_; }
    const PhysicsConfig& physics() const noexcept { return physics_; }

private:
    meshing::Mesh& mesh_;
    materials::MaterialDatabase& materials_;
    PhysicsConfig physics_;
    LinearSolverConfig linearCfg_;
    std::vector<BoundarySpec> boundaries_;
    std::vector<std::shared_ptr<IResidualMonitor>> monitors_;
    FieldRegistry fields_;
    std::unique_ptr<ILinearSolver> linear_;
    std::unique_ptr<ILinearSolver> pressureLinear_;
    std::unique_ptr<ILinearSolver> scalarLinear_;
    std::unique_ptr<class SimpleAlgorithm> simple_;
    std::unique_ptr<class EnergyEquation> energy_;
    std::shared_ptr<::simall::turbulence::ITurbulenceModel> turbulence_;
    int iteration_ = 0;
};

} // namespace simall::solver
