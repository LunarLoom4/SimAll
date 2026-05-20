// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/PeriodicBc.hpp
// Phase  : 16.6 — Periodic boundary using PeriodicPairing twin table.
//
// PeriodicBc is bound to **one** of the two zones in a periodic pair.
// At apply time it copies the field value from the twin cell into the
// boundary cell via an off-diagonal coupling (penalty-symmetric).
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"
#include "solver/PeriodicPairing.hpp"

#include <vector>

namespace simall::solver::bc
{

class PeriodicBc : public IBoundaryCondition
{
public:
    explicit PeriodicBc(PeriodicTransform t = {}) : xform_(t) {}

    BcKind kind() const noexcept override { return BcKind::Periodic; }
    const char* name() const noexcept override { return "Periodic"; }

    void initialize(const meshing::Mesh& mesh) override;
    std::size_t apply(BcContext& ctx) override;

    std::unique_ptr<IBoundaryCondition> clone() const override
    {
        auto c = std::make_unique<PeriodicBc>(xform_);
        c->setZone(zone());
        c->pairs_ = pairs_;
        return c;
    }

    const PeriodicTransform& transform() const noexcept { return xform_; }

private:
    PeriodicTransform xform_;
    std::vector<int> pairs_; // face -> twin face (or -1)
};

} // namespace simall::solver::bc
