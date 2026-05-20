// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/StlImporter.hpp
// Phase  : 5.5 — STL (ASCII or binary) triangulated surface import.
//
// STL has no topology — vertices are duplicated per triangle. The importer
// welds coincident vertices using a tolerance-based hash grid (Akenine-Möller
// 2008 §17.6) so the resulting triangle soup is suitable as a tight closed
// surface for inside/outside tests by the octree volume mesher.
//
// Both ASCII ("solid ... endsolid") and the IEEE binary format
// (80-byte header, uint32 count, 50-byte triangles) are auto-detected.
//
// Fully 3-D. The y- coordinate is NOT special; this is a general geometry
// importer with no symmetry assumption.
// =============================================================================
#pragma once

#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace simall::meshing
{

struct StlSurface
{
    std::vector<util::Vec3d> vertices;
    std::vector<std::array<std::uint32_t, 3>> triangles;
    std::vector<util::Vec3d> normals; // per triangle
};

class StlImporter
{
public:
    /// Returns true on success. Both ASCII and binary STL are auto-detected.
    static bool load(const std::string& path, StlSurface& out, double weldTolerance = 1.0e-9);
};

} // namespace simall::meshing
