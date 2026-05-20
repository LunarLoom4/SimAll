// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MeshStorage.hpp
// Phase  : 5.1 / 6 / Section 6 of ultra-detailed spec.
//
// Structure-of-Arrays (SoA) computational mesh storage. AVX-512-friendly
// 64-byte alignment via util::aligned_vector. Layout EXACTLY matches the
// spec (Sections 6.2–6.4 of the ultra-detailed blueprint):
//
//   NodeStorage : x[], y[], z[]
//   FaceStorage : owner[], neighbor[], area{X,Y,Z}[], centroid{X,Y,Z}[]
//   CellStorage : volume[], centroid{X,Y,Z}[]
//
// Field variables (ρ, u, p, k, ε, …) live in solver::FieldRegistry — NOT
// here, because the mesh is solver-agnostic.
// =============================================================================
#pragma once

#include "utilities/AlignedAllocator.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::meshing
{

using CellId = std::uint64_t;
using FaceId = std::uint64_t;
using NodeId = std::uint64_t;
using ZoneId = std::uint32_t;

inline constexpr CellId kBoundaryCell = static_cast<CellId>(-1);

struct NodeStorage
{
    util::aligned_vector<double> x, y, z;
    std::size_t size() const noexcept { return x.size(); }
    void reserve(std::size_t n)
    {
        x.reserve(n);
        y.reserve(n);
        z.reserve(n);
    }
};

struct FaceStorage
{
    util::aligned_vector<CellId> owner;
    util::aligned_vector<CellId> neighbor; // kBoundaryCell if boundary
    util::aligned_vector<double> areaX, areaY, areaZ;
    util::aligned_vector<double> centroidX, centroidY, centroidZ;
    util::aligned_vector<ZoneId> boundaryZone; // 0 = interior

    // Variable-stride node connectivity in CSR.
    util::aligned_vector<std::int32_t> nodeOffsets; // size = nFaces+1
    util::aligned_vector<NodeId> nodeIndices;

    std::size_t size() const noexcept { return owner.size(); }
};

struct CellStorage
{
    util::aligned_vector<double> volume;
    util::aligned_vector<double> centroidX, centroidY, centroidZ;

    // CSR face connectivity per cell.
    util::aligned_vector<std::int32_t> faceOffsets; // size = nCells+1
    util::aligned_vector<FaceId> faceIndices;

    std::size_t size() const noexcept { return volume.size(); }
};

struct ZoneInfo
{
    ZoneId id;
    std::string name;
    bool isBoundary;
};

class Mesh
{
public:
    NodeStorage& nodes() noexcept { return nodes_; }
    const NodeStorage& nodes() const noexcept { return nodes_; }
    FaceStorage& faces() noexcept { return faces_; }
    const FaceStorage& faces() const noexcept { return faces_; }
    CellStorage& cells() noexcept { return cells_; }
    const CellStorage& cells() const noexcept { return cells_; }

    ZoneId add_zone(std::string name, bool isBoundary);
    ZoneInfo* find_zone(ZoneId id);

    /// Recompute face areas, centroids, and cell volumes from node coords.
    void compute_geometry();

private:
    NodeStorage nodes_;
    FaceStorage faces_;
    CellStorage cells_;
    std::unordered_map<ZoneId, ZoneInfo> zones_;
    ZoneId next_zone_ = 1;
};

} // namespace simall::meshing
