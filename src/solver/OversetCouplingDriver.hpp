// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/OversetCouplingDriver.hpp
// Phase  : 6.5 — Multi-mesh Chimera coupling orchestrator.
//
// Wraps the per-domain solvers (SIMPLE / PISO / Coupled) and the
// meshing::OversetInterpolation engine so that on every outer iteration:
//
//   for each (donor, receptor) pair:
//        interp.interpolate_scalar/vector(donor.fields[X], rcvBuffer)
//        bc.OversetBc::setDonorValues(rcvBuffer at receptor cells)
//   for each domain:
//        domainSolver.iterate()
//
// Convergence of the overset coupling is monitored by the max change in
// donor values between successive interpolation rounds.  We stop when
// either the iteration count or the relative change tolerance is met.
//
// Reference: Steger, Dougherty & Benek (1983); Chan (2009) "Overset Grid
// Methods for the Aerospace Industry"; Henshaw & Schwendeman, J. Comp.
// Phys. 216 (2006).
// =============================================================================
#pragma once

#include "meshing/OversetInterpolation.hpp"
#include "solver/bc/OversetBc.hpp"
#include "solver/FieldRegistry.hpp"
#include "solver/SimpleAlgorithm.hpp"

#include <memory>
#include <string>
#include <vector>

namespace simall::solver
{

struct OversetDomain
{
    meshing::Mesh* mesh = nullptr;
    FieldRegistry* fields = nullptr;
    SimpleAlgorithm* solver = nullptr; // SIMPLE/PISO/SIMPLEC subclass
    /// Receptor cells in this mesh that need donor values each iter.
    std::vector<meshing::CellId> receptorCells;
    /// Overset BC objects keyed by zone whose donor values we will refresh
    /// every outer iteration (one per scalar/vector channel).
    std::vector<bc::OversetBc*> oversetBcs;
};

struct OversetLink
{
    /// Index of the donor and receptor domains in the domains vector.
    int donor = -1;
    int receptor = -1;
    std::unique_ptr<meshing::OversetInterpolation> interp;
    /// Scalar field names to transfer (e.g., "p"); vector names go in vecVars.
    std::vector<std::string> scalarVars;
    std::vector<std::string> vectorVars;
    std::string holeMaskField; // optional
};

struct OversetDriverOptions
{
    int maxOuterIters = 25;
    double relTol = 1.0e-4; // relative change of donor-side variable across iter
    bool verbose = false;
};

struct OversetReport
{
    int iters = 0;
    double finalChange = 0.0;
    bool converged = false;
};

class OversetCouplingDriver
{
public:
    OversetCouplingDriver(std::vector<OversetDomain> domains,
                          std::vector<OversetLink> links,
                          OversetDriverOptions opts = {});

    /// Build interpolation stencils for all links (called once after
    /// meshes are immutable). Returns false if any link fails to build.
    bool build();

    /// Run the outer coupling loop. Each outer iter:
    ///   1) interp every link (donor→receptor field samples)
    ///   2) call each domain's solver.iterate()
    ///   3) measure max relative change of monitored variables across the
    ///      bridge; stop if below tol.
    OversetReport iterate_to_convergence();

    /// Single outer iteration (for use inside an outer time-step loop).
    /// Returns max relative change observed across all links.
    double iterate_once();

private:
    void interpolate_all_links();
    void apply_donor_values_to_bcs();
    double measure_max_change(); // compares donor field snapshots

    std::vector<OversetDomain> domains_;
    std::vector<OversetLink> links_;
    OversetDriverOptions opt_;
    bool built_ = false;

    // Snapshot of last-iteration donor field values (per link, per receptor)
    // used for change measurement. Each entry is a flat array of scalars
    // accumulated across all scalarVars+vectorVars (vectors flattened to 3).
    std::vector<std::vector<double>> lastSnapshot_;
};

} // namespace simall::solver
