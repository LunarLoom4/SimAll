// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/MaterialAssignment.hpp
// Phase  : 15.7 — Zone → Material mapping.
//
// Bridges the meshing layer (ZoneId per cell / face zone) with the materials
// database.  At solver assemble time the solver queries this assignment per
// cell to fetch the active Material (and hence its property functions).
//
//   * Cell zone assignment   : (ZoneId  → MaterialKey)   for fluid/solid bulk
//   * Face zone assignment   : (ZoneId  → MaterialKey)   for thin walls
//   * Default fallback       : when a zone is not explicitly assigned
//
// MaterialKey is the std::string name registered in MaterialDatabase.
// =============================================================================
#pragma once

#include "materials/Material.hpp"

#include <cstdint>
#include <string>
#include <unordered_map>

namespace simall::materials
{

using ZoneKey = std::uint32_t;

class MaterialAssignment
{
public:
    explicit MaterialAssignment(MaterialDatabase& db) : db_(db) {}

    void setCellZone(ZoneKey zone, std::string materialName)
    {
        cellAssign_[zone] = std::move(materialName);
    }
    void setFaceZone(ZoneKey zone, std::string materialName)
    {
        faceAssign_[zone] = std::move(materialName);
    }
    void setDefault(std::string materialName) { default_ = std::move(materialName); }

    /// Look up the material assigned to a cell-zone; returns nullptr if
    /// neither an explicit nor default assignment resolves.
    Material* materialForCellZone(ZoneKey zone) const
    {
        auto it = cellAssign_.find(zone);
        const std::string& key = (it != cellAssign_.end()) ? it->second : default_;
        if (key.empty())
            return nullptr;
        return db_.find(key);
    }

    Material* materialForFaceZone(ZoneKey zone) const
    {
        auto it = faceAssign_.find(zone);
        if (it == faceAssign_.end())
            return nullptr;
        return db_.find(it->second);
    }

    const std::unordered_map<ZoneKey, std::string>& cellAssignments() const noexcept
    {
        return cellAssign_;
    }
    const std::unordered_map<ZoneKey, std::string>& faceAssignments() const noexcept
    {
        return faceAssign_;
    }
    const std::string& defaultMaterial() const noexcept { return default_; }

private:
    MaterialDatabase& db_;
    std::unordered_map<ZoneKey, std::string> cellAssign_;
    std::unordered_map<ZoneKey, std::string> faceAssign_;
    std::string default_;
};

} // namespace simall::materials
