// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/InterfaceBc.hpp
// Phase  : 16.10 — Conjugate / fluid–solid interface.
//
// The interface enforces:
//   T_fluid  = T_solid                       (Dirichlet match)
//   q_fluid  = -q_solid                      (heat-flux conservation)
//
// Implemented as a paired-zone BC: each side of the interface holds a
// pointer to a shared InterfaceState that the solver advances each outer
// iteration with the most recent (T, q) trace from the other side.
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"

#include <memory>
#include <vector>

namespace simall::solver::bc
{

/// Trace stored per face of the interface zone — populated by the solver
/// after each outer iteration on the OTHER side of the interface.
struct InterfaceTrace
{
    std::vector<double> T_other; // size = nFaces on this zone
    std::vector<double> q_other;
};

class InterfaceBc : public IBoundaryCondition
{
public:
    explicit InterfaceBc(std::shared_ptr<InterfaceTrace> shared) : trace_(std::move(shared)) {}

    BcKind kind() const noexcept override { return BcKind::Interface; }
    const char* name() const noexcept override { return "Interface"; }

    std::size_t apply(BcContext& ctx) override;

    std::unique_ptr<IBoundaryCondition> clone() const override
    {
        auto c = std::make_unique<InterfaceBc>(trace_);
        c->setZone(zone());
        return c;
    }

    InterfaceTrace* trace() noexcept { return trace_.get(); }
    const InterfaceTrace* trace() const noexcept { return trace_.get(); }

private:
    std::shared_ptr<InterfaceTrace> trace_;
};

} // namespace simall::solver::bc
