// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/PrismLayerExtruder.cpp
// =============================================================================
#include "meshing/PrismLayerExtruder.hpp"

#include "core/Logger.hpp"
#include "meshing/Connectivity.hpp"

namespace simall::meshing
{

void PrismLayerExtruder::extrude(const std::vector<util::Vec3d>& wallNodes,
                                 const std::vector<std::array<std::uint32_t, 3>>& wallTris,
                                 const std::vector<util::Vec3d>& normals,
                                 PrismLayerOptions opt,
                                 Mesh& out)
{
    const std::size_t nW = wallNodes.size();
    const int N = std::max(1, opt.nLayers);

    // Compute per-node cumulative offsets (negative along outward normal so
    // the layers grow INTO the fluid domain).
    std::vector<double> heights(N + 1, 0.0);
    double t = opt.firstLayerHeight;
    for (int i = 1; i <= N; ++i) {
        heights[i] = heights[i - 1] + t;
        t *= opt.growthRatio;
    }

    NodeStorage ns;
    ns.reserve(nW * (N + 1));
    for (int i = 0; i <= N; ++i) {
        for (std::size_t k = 0; k < nW; ++k) {
            ns.x.push_back(wallNodes[k].x - normals[k].x * heights[i]);
            ns.y.push_back(wallNodes[k].y - normals[k].y * heights[i]);
            ns.z.push_back(wallNodes[k].z - normals[k].z * heights[i]);
        }
    }

    auto pid = [&](int layer, std::uint32_t wallNode) -> NodeId {
        return NodeId(layer * nW + wallNode);
    };

    std::vector<CellDescriptor> cells;
    cells.reserve(wallTris.size() * N);

    for (int i = 0; i < N; ++i) {
        for (const auto& tri : wallTris) {
            // Wedge / triangular prism with bottom on layer i, top on i+1.
            NodeId b0 = pid(i, tri[0]);
            NodeId b1 = pid(i, tri[1]);
            NodeId b2 = pid(i, tri[2]);
            NodeId t0 = pid(i + 1, tri[0]);
            NodeId t1 = pid(i + 1, tri[1]);
            NodeId t2 = pid(i + 1, tri[2]);
            CellDescriptor cd;
            cd.faces = {
                {b0, b2, b1},     // bottom triangle (outward-CCW)
                {t0, t1, t2},     // top triangle
                {b0, b1, t1, t0}, // side quad 0-1
                {b1, b2, t2, t1}, // side quad 1-2
                {b2, b0, t0, t2}  // side quad 2-0
            };
            cells.push_back(std::move(cd));
        }
    }

    ConnectivityBuilder::build(out, ns, cells);
    SIMALL_LOG_INFO("Mesh",
                    "Prism layers extruded: ",
                    N,
                    " layers x ",
                    wallTris.size(),
                    " wall tris = ",
                    cells.size(),
                    " prisms");
}

} // namespace simall::meshing
