// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/SimplecAlgorithm.cpp
//
// SIMPLEC is a SimpleAlgorithm with PvCouplingVariant::SIMPLEC set in
// SimpleOptions::algorithm. The SIMPLEC denominator branch lives in
// SimpleAlgorithm::assemble_pressure_correction(). This translation unit
// exists to give the subclass a non-header anchor for symbol export and
// to host any future SIMPLEC-only diagnostics.
// =============================================================================
#include "solver/SimplecAlgorithm.hpp"

namespace simall::solver
{
// (currently no SIMPLEC-only out-of-line members)
} // namespace simall::solver
