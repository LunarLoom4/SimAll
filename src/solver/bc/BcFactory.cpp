// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/BcFactory.cpp
// =============================================================================
#include "solver/bc/BcFactory.hpp"

#include "solver/bc/AxisymmetricBc.hpp"
#include "solver/bc/FanBc.hpp"
#include "solver/bc/InletBc.hpp"
#include "solver/bc/InterfaceBc.hpp"
#include "solver/bc/OutletBc.hpp"
#include "solver/bc/OversetBc.hpp"
#include "solver/bc/PeriodicBc.hpp"
#include "solver/bc/PorousJumpBc.hpp"
#include "solver/bc/SymmetryBc.hpp"
#include "solver/bc/WallBc.hpp"

namespace simall::solver::bc
{

BcFactory& BcFactory::instance()
{
    static BcFactory inst;
    return inst;
}

bool BcFactory::registerBuilder(BcKind kind, Builder b)
{
    return builders_.emplace(static_cast<int>(kind), std::move(b)).second;
}

std::unique_ptr<IBoundaryCondition> BcFactory::create(BcKind kind) const
{
    auto it = builders_.find(static_cast<int>(kind));
    return it == builders_.end() ? nullptr : it->second();
}

std::unique_ptr<IBoundaryCondition> BcFactory::create(const std::string& name) const
{
    for (auto const& [_, b] : builders_) {
        auto inst = b();
        if (inst && name == inst->name())
            return inst;
    }
    return nullptr;
}

void BcFactory::registerBuiltins()
{
    auto& f = instance();
    f.registerBuilder(BcKind::Wall, [] { return std::make_unique<WallBc>(); });
    f.registerBuilder(BcKind::Inlet, [] { return std::make_unique<InletBc>(); });
    f.registerBuilder(BcKind::Outlet, [] { return std::make_unique<OutletBc>(); });
    f.registerBuilder(BcKind::Symmetry, [] { return std::make_unique<SymmetryBc>(); });
    f.registerBuilder(BcKind::Axisymmetric, [] { return std::make_unique<AxisymmetricBc>(); });
    f.registerBuilder(BcKind::Periodic, [] { return std::make_unique<PeriodicBc>(); });
    f.registerBuilder(BcKind::PorousJump, [] { return std::make_unique<PorousJumpBc>(); });
    f.registerBuilder(BcKind::Fan, [] {
        return std::make_unique<FanBc>(materials::PolynomialFit{{0.0}, 0.0, 1.0e9, true});
    });
    f.registerBuilder(BcKind::Overset, [] { return std::make_unique<OversetBc>(); });
    f.registerBuilder(BcKind::Interface, [] {
        return std::make_unique<InterfaceBc>(std::make_shared<InterfaceTrace>());
    });
}

} // namespace simall::solver::bc
