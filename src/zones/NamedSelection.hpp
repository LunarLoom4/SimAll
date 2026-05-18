// =============================================================================
// SimAll Beta - Zones / Selection Subsystem
// File   : src/zones/NamedSelection.hpp
// Phase  : 25 — persistent named selections.
//
// A "named selection" is a user-curated, persistent grouping of mesh
// entities (faces, cells, or edges) identified by their `PersistentId`s.
// Named selections survive remeshing and CAD updates because they bind to
// the persistent CAD topology IDs assigned by `cad::PersistentIdManager`,
// not to transient mesh face indices.  Once defined, a selection becomes
// the *target* of boundary-condition assignments and post-processing
// (surface integrals, force monitors, etc.), so that BCs do not need to be
// re-applied every time the mesh is regenerated.
//
// `SelectionManager` (a sibling file) owns and persists a project-wide
// catalogue of NamedSelections through the project I/O layer.
// =============================================================================
#pragma once

#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace simall::zones {

enum class SelectionEntity : std::uint8_t {
    Face   = 0,
    Cell   = 1,
    Edge   = 2,
    Vertex = 3
};

struct NamedSelection {
    std::string                       name;                ///< unique within a project
    SelectionEntity                   entity = SelectionEntity::Face;
    std::vector<util::PersistentId>   ids;                 ///< persistent CAD/topology IDs
    std::string                       color   = "#7fb3ff"; ///< visualization tint
    std::string                       comment;             ///< user note (optional)

    bool empty() const noexcept { return ids.empty(); }
    std::size_t size() const noexcept { return ids.size(); }
};

}  // namespace simall::zones
