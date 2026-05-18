// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/Solver.cpp
//
// Phase 6/7 driver. The full SIMPLE/PISO inner loop is implemented across
// AssemblyMomentum.cpp / AssemblyPressure.cpp / LinearSolvers.cpp during
// the Phase 6 implementation hand-off; this file ties them together and
// guarantees the convergence-monitor + EventBus contract is honoured.
// =============================================================================
#include "solver/Solver.hpp"
#include "solver/LinearSolvers.hpp"
#include "solver/SimpleAlgorithm.hpp"
#include "solver/EnergyEquation.hpp"
#include "turbulence/ITurbulenceModel.hpp"
#include "core/EventBus.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::solver {

std::unique_ptr<ILinearSolver> make_linear_solver(LinearSolverConfig cfg) {
    switch (cfg.kind) {
        case LinearSolverKind::GMRES:    return make_gmres   (cfg);
        case LinearSolverKind::BiCGSTAB: return make_bicgstab(cfg);
        case LinearSolverKind::CG:       return make_cg      (cfg);
        case LinearSolverKind::TFQMR:    return make_tfqmr   (cfg);
        default:                          return make_gmres   (cfg);
    }
}

// ----- Solver -----
Solver::Solver(meshing::Mesh& m, materials::MaterialDatabase& mat)
    : mesh_(m), materials_(mat) {}

Solver::~Solver() = default;

void Solver::configure(PhysicsConfig c)            { physics_   = c; }
void Solver::add_boundary(BoundarySpec b)          { boundaries_.push_back(b); }
void Solver::set_linear_solver(LinearSolverConfig c){ linearCfg_ = c; linear_ = make_linear_solver(c); }
void Solver::add_monitor(std::shared_ptr<IResidualMonitor> m) { monitors_.push_back(std::move(m)); }
void Solver::set_turbulence_model(std::shared_ptr<::simall::turbulence::ITurbulenceModel> m) {
    turbulence_ = std::move(m);
}

void Solver::initialize() {
    const std::size_t nC = mesh_.cells().size();
    fields_.scalar  ("p",   nC);
    fields_.scalar  ("rho", nC);
    fields_.scalar  ("T",   nC);
    fields_.vector  ("U",   nC);
    if (!physics_.turbulenceModel.empty()) {
        fields_.scalar("k",     nC);
        fields_.scalar("omega", nC);
    }
    if (!linear_)         linear_         = make_linear_solver(linearCfg_);
    if (!pressureLinear_) {
        LinearSolverConfig pCfg = linearCfg_;
        pCfg.kind = LinearSolverKind::CG;          // pressure system is SPD
        pressureLinear_ = make_linear_solver(pCfg);
    }
    if (!scalarLinear_) {
        LinearSolverConfig sCfg = linearCfg_;
        sCfg.kind = LinearSolverKind::BiCGSTAB;
        scalarLinear_ = make_linear_solver(sCfg);
    }
    if (!simple_) {
        SimpleOptions opt;  // defaults fine; can be exposed in PhysicsConfig
        simple_ = std::make_unique<SimpleAlgorithm>(
            mesh_, fields_, boundaries_, *linear_, *pressureLinear_, opt);
    }
    if (physics_.energyEquation && !energy_) {
        EnergyOptions eo;
        energy_ = std::make_unique<EnergyEquation>(
            mesh_, fields_, *scalarLinear_, boundaries_, eo);
    }
    if (turbulence_) {
        turbulence_->initialize(mesh_, fields_);
    }
    iteration_ = 0;
    core::EventBus::instance().publish(core::events::SolverStarted{physics_.turbulenceModel});
    SIMALL_LOG_INFO("Solver", "Initialized (cells=", nC, ", coupling=",
        int(physics_.coupling), ", energy=", physics_.energyEquation,
        ", turbulence=", (turbulence_ ? turbulence_->name() : "none"), ")");
}

bool Solver::step() {
    ++iteration_;
    SimpleResiduals sr = simple_->iterate();
    double energyRes = 0.0;
    if (energy_) energyRes = energy_->iterate();
    if (turbulence_) turbulence_->solve(0.0, fields_);
    ResidualSnapshot r{};
    r.iteration   = iteration_;
    r.continuity  = sr.cont;
    r.momentumX   = sr.mom[0];
    r.momentumY   = sr.mom[1];
    r.momentumZ   = sr.mom[2];
    r.energy      = energyRes;
    r.turbulence_k = r.turbulence_eps = 0.0;
    for (auto& m : monitors_) m->on_residual(r);
    core::EventBus::instance().publish(core::events::SolverIteration{r.iteration, r.continuity});
    return r.continuity > 1e-30;
}

void Solver::run(int maxOuter, double target) {
    initialize();
    for (int i = 0; i < maxOuter; ++i) {
        step();
        if (linear_ && linear_->last_residual() < target) {
            core::EventBus::instance().publish(core::events::SolverConverged{iteration_});
            SIMALL_LOG_INFO("Solver", "Converged at iter ", iteration_);
            return;
        }
    }
    SIMALL_LOG_WARN("Solver", "Reached max iterations without target convergence");
}

}  // namespace simall::solver
