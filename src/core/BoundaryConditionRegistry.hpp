// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/BoundaryConditionRegistry.hpp
// Phase  : 22 Pass 2 — central, physics-agnostic BC catalogue.
//
// Why this exists in core/ rather than solver/bc/
// -----------------------------------------------
// The pre-existing src/solver/bc/BcFactory is a *builder* (BcKind → default
// instance constructor) tightly bound to fluid-flow boundary conditions
// (Wall / Inlet / Outlet / Symmetry / Periodic / PorousJump / Fan / Overset
// / Interface). It depends on solver/CSRMatrix, FieldRegistry, MeshStorage.
//
// Other subsystems also own boundary-condition concepts:
//   - electromagnetics : MagneticInsulation, PerfectConductor, Floquet
//   - structural       : FixedConstraint, Roller, Spring, Bolt
//   - thermal          : ConvectiveCooling, Radiosity, Thin Layer
//   - acoustics        : ImpedanceWall, SoundHard, PML
//   - electrochemistry : ElectrodeBoundary, Insulating, Concentration
//
// These cannot live under solver::bc without inverting the dependency layer
// (em / structural / chemistry must not transitively pull in solver). This
// registry holds them all via a type-erased opaque handle plus metadata
// (subsystem, kind id, zone, variable). Drivers iterate the registry; each
// subsystem casts back to its own concrete interface inside its own apply
// routines.
//
// Layering: core depends only on utilities. We therefore mirror meshing's
// ZoneId underlying type (uint32) here as BcZoneId rather than including
// meshing/MeshStorage.hpp from core. solver::bc adapters convert at the
// boundary (one-line static_cast).
//
// Thread safety: register/remove/clear are synchronised by an internal
// mutex. find() returning a pointer is also synchronised but the pointer
// is only safe to dereference under the assumption that no concurrent
// remove() is in flight for that handle (typical SimAll lifecycle: register
// during setup, iterate during solve, clear at shutdown).
// =============================================================================
#pragma once

#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace simall::core
{

/// Generic zone identifier. Matches meshing::ZoneId underlying type so
/// callers can static_cast at the boundary without lossy conversion.
using BcZoneId = std::uint32_t;

/// Per-registration record. The opaque pointer carries a type-erased BC
/// instance (e.g. std::shared_ptr<solver::bc::IBoundaryCondition>) and is
/// retrieved by the owning subsystem via BoundaryConditionRegistry::as<T>().
struct BoundaryConditionEntry
{
    int kindId = 0;       // subsystem-defined enum value
    std::string kindName; // e.g. "Wall", "MagneticInsulation"
    BcZoneId zone = 0;
    std::string variable;         // "U", "p", "T", "E", "B", "phi", ...
    std::string subsystem;        // "solver", "em", "structural", ...
    std::shared_ptr<void> opaque; // type-erased concrete BC instance
};

class BoundaryConditionRegistry
{
public:
    using Handle = std::uint64_t;
    static constexpr Handle kInvalid = 0;

    BoundaryConditionRegistry() = default;
    BoundaryConditionRegistry(const BoundaryConditionRegistry&) = delete;
    BoundaryConditionRegistry& operator=(const BoundaryConditionRegistry&) = delete;

    /// Insert a new entry. Returns a non-zero handle on success.
    Handle add(BoundaryConditionEntry entry);

    /// Remove an entry by handle. Returns true if anything was removed.
    bool remove(Handle h);

    /// Drop all entries. Iteration order resets.
    void clear() noexcept;

    /// Read-only lookup by handle. Returns nullptr if not present.
    const BoundaryConditionEntry* find(Handle h) const noexcept;

    /// Lookup helpers — return handles in insertion order.
    std::vector<Handle> by_zone(BcZoneId z) const;
    std::vector<Handle> by_variable(std::string_view v) const;
    std::vector<Handle> by_subsystem(std::string_view s) const;

    /// First match for a (zone, variable) pair. Returns kInvalid on miss.
    Handle first(BcZoneId z, std::string_view variable) const noexcept;

    /// Iterate all entries in insertion order, invoking fn(handle, entry).
    /// fn must not call add/remove/clear on this registry (the internal
    /// mutex is held for the whole walk).
    template <class Fn> void for_each(Fn&& fn) const
    {
        std::lock_guard<std::mutex> lk(mu_);
        for (auto h : order_) {
            auto it = entries_.find(h);
            if (it != entries_.end())
                fn(h, it->second);
        }
    }

    std::size_t size() const noexcept;
    bool empty() const noexcept;

    /// Typed accessor: visits opaque as a T* (caller-owned cast). Returns
    /// nullptr if the handle is unknown or opaque is null.
    template <class T> T* as(Handle h) const noexcept
    {
        const auto* e = find(h);
        return e ? static_cast<T*>(e->opaque.get()) : nullptr;
    }

    /// Process-wide singleton. Subsystems that own their own Project may
    /// also instantiate a private registry — they are independent.
    static BoundaryConditionRegistry& instance();

private:
    mutable std::mutex mu_;
    Handle next_ = 1;
    std::unordered_map<Handle, BoundaryConditionEntry> entries_;
    std::vector<Handle> order_; // insertion order
};

} // namespace simall::core
