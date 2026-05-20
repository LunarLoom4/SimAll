// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/SimplecAlgorithm.hpp
// Phase  : 6.5 — SIMPLE-Consistent (Van Doormaal & Raithby 1984)
//
// SIMPLEC differs from SIMPLE only in the pressure-correction equation
// coefficient. SIMPLE uses (V/aP)_f; SIMPLEC uses (V/(aP - Σ a_NB))_f.
// Because Σa_NB is typically of comparable magnitude to aP, the SIMPLEC
// denominator is much smaller, yielding stronger pressure corrections per
// iteration. The practical consequence is that pressure under-relaxation
// urfP ≈ 1.0 is stable (vs ~0.3 for SIMPLE), giving 2-3× fewer outer iters.
//
// Implemented as a thin façade over SimpleAlgorithm; the branch lives in
// assemble_pressure_correction() guarded by algorithm == SIMPLEC.
//
// References:
//   Van Doormaal & Raithby, Num. Heat Transfer 7, 147-163 (1984).
//   Versteeg & Malalasekera, ch 6.6.
// =============================================================================
#pragma once

#include "solver/SimpleAlgorithm.hpp"

namespace simall::solver
{

struct SimplecOptions
{
    double urfU = 0.9; // can be much closer to 1 than SIMPLE
    double urfP = 1.0; // SIMPLEC's whole point
    double rho = 1.0;
    double mu = 1.0e-3;
    double dt = 0.0;
    TemporalScheme timeScheme = TemporalScheme::ImplicitEuler;
};

class SimplecAlgorithm : public SimpleAlgorithm
{
public:
    SimplecAlgorithm(meshing::Mesh& mesh,
                     FieldRegistry& fields,
                     const std::vector<BoundarySpec>& boundaries,
                     ILinearSolver& momentumSolver,
                     ILinearSolver& pressureSolver,
                     SimplecOptions opts)
        : SimpleAlgorithm(
              mesh, fields, boundaries, momentumSolver, pressureSolver, to_simple_opts(opts))
    {
    }

private:
    static SimpleOptions to_simple_opts(const SimplecOptions& s)
    {
        SimpleOptions o;
        o.urfU = s.urfU;
        o.urfP = s.urfP;
        o.rho = s.rho;
        o.mu = s.mu;
        o.dt = s.dt;
        o.timeScheme = s.timeScheme;
        o.algorithm = PvCouplingVariant::SIMPLEC;
        return o;
    }
};

} // namespace simall::solver
