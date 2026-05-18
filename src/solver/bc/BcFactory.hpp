// =============================================================================
// SimAll Beta - Solver/BC
// File   : src/solver/bc/BcFactory.hpp
// Phase  : 16.11 — Registry / factory for boundary conditions.
//
// Concrete BC classes register a builder under their BcKind enum.  Code that
// constructs a BC by kind (e.g. project loaders, the UI, the scripting API)
// uses BcFactory::create(kind) to instantiate a default instance, which the
// caller then configures via the typed params() accessor.
// =============================================================================
#pragma once

#include "solver/bc/Bc.hpp"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

namespace simall::solver::bc {

class BcFactory {
public:
    using Builder = std::function<std::unique_ptr<IBoundaryCondition>()>;

    /// Singleton access.
    static BcFactory& instance();

    /// Register a builder under a given kind. Returns true if newly added.
    bool registerBuilder(BcKind kind, Builder b);

    /// Create a default-configured instance of the requested kind, or
    /// nullptr if no builder is registered.
    std::unique_ptr<IBoundaryCondition> create(BcKind kind) const;

    /// Lookup by string name (matches IBoundaryCondition::name()).
    std::unique_ptr<IBoundaryCondition> create(const std::string& name) const;

    /// Static convenience: instantiates each of the built-in BCs.  Idempotent.
    static void registerBuiltins();

private:
    BcFactory() = default;
    std::unordered_map<int, Builder> builders_;
};

}  // namespace simall::solver::bc
