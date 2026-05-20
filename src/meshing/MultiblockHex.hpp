// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MultiblockHex.hpp
// Phase  : 6.9 — Multi-block structured hex assembler.
//
// Each block is a topologically-cubical region (i, j, k) with explicit
// corner coordinates and per-direction divisions; multiple blocks are
// joined along shared faces.  Connectivity is reconciled by matching node
// coordinates within tolerance — internal faces have neighbour cells set,
// boundary faces carry the user-specified zone tag (per block face).
//
// Block face tagging convention (per Hexa/Plot3D community):
//   imin (-x): face index 0           iWall0
//   imax (+x): face index 1           iWall1
//   jmin (-y): face index 2           jWall0
//   jmax (+y): face index 3           jWall1
//   kmin (-z): face index 4           kWall0
//   kmax (+z): face index 5           kWall1
//
// Use cases: turbomachinery cascades, square-cross-section ducts with
// branch fittings, automotive HVAC, conjugate-heat-transfer enclosures.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace simall::meshing {

struct HexBlock {
    std::array<util::Vec3d, 8> corners;       // standard hex CCW ordering
    std::array<std::uint32_t, 3> divisions{8,8,8};
    std::array<double, 3>        grading{1.0, 1.0, 1.0};   // simpleGrading
                                                            // expansion ratio
                                                            // (last/first cell)
                                                            // 1.0 == uniform
    std::array<std::uint32_t, 6> faceZones{0,0,0,0,0,0};
    std::string name;
};

struct MultiblockProps {
    double matchTol = 1e-6;
    bool   reorder  = true;          // RCM ordering after assembly
};

class MultiblockHex {
public:
    void initialize(MultiblockProps props);

    void add_block(HexBlock block) { blocks_.push_back(std::move(block)); }

    /// Assemble all blocks into a single connected mesh.  Returns total
    /// hex count.
    std::size_t assemble(Mesh& outMesh);

    std::size_t block_count() const noexcept { return blocks_.size(); }

    const MultiblockProps& props() const noexcept { return p_; }

private:
    MultiblockProps         p_{};
    std::vector<HexBlock>   blocks_;
};

}  // namespace simall::meshing
