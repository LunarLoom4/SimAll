// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/CellAdapterRegistry.hpp
// Phase  : 22 Pass 22.5
//
// CellAdapterRegistry -- a name-keyed registry of ICellAdapter factory
// functions.  See CellAdapter.hpp for the broader design.
//
// Threading
// ---------
// The registry is intended to be populated once at application boot and
// read-many afterwards.  We still take an internal mutex on every API
// call so that test rigs (and a future plugin-reload feature) can swap
// factories at runtime safely.  The default Qt main-thread usage
// pattern incurs uncontested-mutex cost only, which is negligible
// compared to a CAD import or a solver run.
//
// Identity
// --------
// Adapter IDs are arbitrary dotted strings.  Convention (Pass 22.5):
//   * "cad.import.<format>"     -- CAD readers
//   * "cad.heal"                -- single-shape healing
//   * "cad.tessellate"          -- BREP -> triangle mesh
//   * "mesh.surface.<algo>"     -- surface meshers
//   * "mesh.volume.<algo>"      -- volume meshers
//   * "setup.boundary.<scheme>" -- BC application
//   * "solver.run.<algo>"       -- pressure-velocity coupled drivers
//   * "results.load"            -- read a results file into the viewport
//   * "results.export.<format>" -- write results out
// =============================================================================
#pragma once

#include "workbench/CellAdapter.hpp"

#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::workbench
{

class CellAdapterRegistry
{
public:
    // A factory returns a freshly-constructed adapter every call so the
    // engine can run several cells concurrently without aliasing
    // `last_error_` between them.
    using Factory = std::function<std::unique_ptr<ICellAdapter>()>;

    CellAdapterRegistry() = default;
    CellAdapterRegistry(const CellAdapterRegistry&) = delete;
    CellAdapterRegistry& operator=(const CellAdapterRegistry&) = delete;

    // Register (or overwrite) a factory under the given id.  Returns
    // true if the id was new, false if it overwrote an existing entry.
    bool register_factory(std::string id, Factory f);

    // Remove a registration.  Returns true if something was actually
    // removed.  Idempotent.
    bool unregister(const std::string& id);

    [[nodiscard]] bool has(const std::string& id) const;

    // Construct a fresh adapter.  Returns nullptr when the id is not
    // registered.
    [[nodiscard]] std::unique_ptr<ICellAdapter> make(const std::string& id) const;

    // Snapshot of all registered ids -- useful for the GUI "Available
    // adapters" picker.  Sorted lexicographically for deterministic UI.
    [[nodiscard]] std::vector<std::string> keys() const;

    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] bool empty() const;

    void clear();

private:
    mutable std::mutex mu_;
    std::unordered_map<std::string, Factory> factories_;
};

} // namespace simall::workbench
