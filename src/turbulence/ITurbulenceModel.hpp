// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/ITurbulenceModel.hpp
// Phase  : 8 (TURBULENCE MODELS)
//
// All RANS/LES/Hybrid models implement this interface and register with the
// TurbulenceRegistry. The solver framework treats them as opaque
// "turbulence plugins" per Master Blueprint Principle 4.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <memory>
#include <string>
#include <unordered_map>

namespace simall::turbulence
{

class ITurbulenceModel
{
public:
    virtual ~ITurbulenceModel() = default;
    virtual std::string name() const = 0;
    virtual void initialize(meshing::Mesh&, solver::FieldRegistry&) = 0;
    virtual void solve(double dt, solver::FieldRegistry&) = 0;
    virtual double turbulent_viscosity(std::size_t cellId) const = 0;
};

class TurbulenceRegistry
{
public:
    static TurbulenceRegistry& instance()
    {
        static TurbulenceRegistry r;
        return r;
    }

    using Factory = std::unique_ptr<ITurbulenceModel> (*)();
    void register_model(const std::string& key, Factory f) { factories_[key] = f; }
    std::unique_ptr<ITurbulenceModel> create(const std::string& key) const
    {
        auto it = factories_.find(key);
        return it == factories_.end() ? nullptr : it->second();
    }

private:
    std::unordered_map<std::string, Factory> factories_;
};

} // namespace simall::turbulence
