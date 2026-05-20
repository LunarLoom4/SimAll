// =============================================================================
// SimAll Beta - Zones / Selection Subsystem
// File   : src/zones/SelectionManager.hpp
// Phase  : 25
//
// Project-wide registry of NamedSelections.  Provides lookup, creation,
// merging, set algebra (union/intersect/subtract), and stable JSON
// serialization so that selections persist alongside the project file.
// =============================================================================
#pragma once

#include "zones/NamedSelection.hpp"

#include <optional>
#include <string>
#include <unordered_map>

namespace simall::zones
{

class SelectionManager
{
public:
    /// Insert (or replace) a selection.  Returns true if a previous
    /// definition with the same name was overwritten.
    bool define(const NamedSelection& s);

    /// Erase by name.  Returns true if a selection was removed.
    bool remove(const std::string& name);

    /// Lookup by name.  Returns `nullptr` if missing.
    const NamedSelection* find(const std::string& name) const noexcept;
    NamedSelection* find(const std::string& name) noexcept;

    /// Rename.  Fails if `newName` already exists.
    bool rename(const std::string& oldName, const std::string& newName);

    /// Set algebra — results are returned as fresh NamedSelections.
    std::optional<NamedSelection> set_union(const std::string& a,
                                            const std::string& b,
                                            const std::string& result) const;
    std::optional<NamedSelection> set_intersect(const std::string& a,
                                                const std::string& b,
                                                const std::string& result) const;
    std::optional<NamedSelection> set_subtract(const std::string& a,
                                               const std::string& b,
                                               const std::string& result) const;

    /// Bulk access.
    std::vector<std::string> names() const;
    std::size_t size() const noexcept { return store_.size(); }
    void clear() noexcept { store_.clear(); }

    /// Stable round-trip JSON serialization (used by the project I/O layer).
    std::string to_json() const;
    bool from_json(const std::string& json);

private:
    std::unordered_map<std::string, NamedSelection> store_;
};

} // namespace simall::zones
