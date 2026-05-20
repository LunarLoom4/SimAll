// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/StlSurfaceReader.cpp
// =============================================================================
#include "cad/StlSurfaceReader.hpp"

#include "cad/ShapeHandleInternal.hpp"
#include "core/Logger.hpp"

#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <unordered_map>

#include <BRep_Builder.hxx>
#include <Poly_Triangulation.hxx>
#include <RWStl.hxx>
#include <StlAPI_Reader.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>

namespace simall::cad
{

namespace
{

struct VertexKey
{
    long long ix, iy, iz;
    bool operator==(const VertexKey& o) const noexcept
    {
        return ix == o.ix && iy == o.iy && iz == o.iz;
    }
};
struct VertexKeyHash
{
    std::size_t operator()(const VertexKey& k) const noexcept
    {
        std::size_t h = std::hash<long long>{}(k.ix);
        h ^= std::hash<long long>{}(k.iy) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<long long>{}(k.iz) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

} // namespace

TriangleMesh StlSurfaceReader::readTriangles(const std::string& path, const StlReadOptions& opts)
{
    if (!std::filesystem::exists(path))
        throw std::runtime_error("StlSurfaceReader: file not found: " + path);

    Handle(Poly_Triangulation) poly = RWStl::ReadFile(path.c_str());
    if (poly.IsNull())
        throw std::runtime_error("StlSurfaceReader: failed to read " + path);

    TriangleMesh out;
    const Standard_Integer nNodes = poly->NbNodes();
    const Standard_Integer nTris = poly->NbTriangles();
    out.points.reserve(nNodes);
    out.triangles.reserve(nTris);
    out.triangleFaceId.assign(nTris, util::PersistentId{1});

    std::vector<std::uint32_t> remap(static_cast<std::size_t>(nNodes), 0);

    if (opts.mergeCoincidentVertices && opts.mergeTolerance > 0.0) {
        std::unordered_map<VertexKey, std::uint32_t, VertexKeyHash> table;
        const double inv = 1.0 / opts.mergeTolerance;
        for (Standard_Integer i = 1; i <= nNodes; ++i) {
            gp_Pnt p = poly->Node(i);
            VertexKey k{static_cast<long long>(std::llround(p.X() * inv)),
                        static_cast<long long>(std::llround(p.Y() * inv)),
                        static_cast<long long>(std::llround(p.Z() * inv))};
            auto it = table.find(k);
            if (it == table.end()) {
                std::uint32_t idx = static_cast<std::uint32_t>(out.points.size());
                table.emplace(k, idx);
                out.points.push_back({p.X(), p.Y(), p.Z()});
                remap[static_cast<std::size_t>(i - 1)] = idx;
            } else {
                remap[static_cast<std::size_t>(i - 1)] = it->second;
            }
        }
    } else {
        for (Standard_Integer i = 1; i <= nNodes; ++i) {
            gp_Pnt p = poly->Node(i);
            out.points.push_back({p.X(), p.Y(), p.Z()});
            remap[static_cast<std::size_t>(i - 1)] = static_cast<std::uint32_t>(i - 1);
        }
    }

    for (Standard_Integer i = 1; i <= nTris; ++i) {
        Standard_Integer a, b, c;
        poly->Triangle(i).Get(a, b, c);
        out.triangles.push_back({remap[static_cast<std::size_t>(a - 1)],
                                 remap[static_cast<std::size_t>(b - 1)],
                                 remap[static_cast<std::size_t>(c - 1)]});
    }
    SIMALL_LOG_INFO("CAD/STL",
                    "Read ",
                    path,
                    " — ",
                    out.points.size(),
                    " verts, ",
                    out.triangles.size(),
                    " tris");
    return out;
}

ShapeHandle StlSurfaceReader::read(const std::string& path, const StlReadOptions& /*opts*/)
{
    if (!std::filesystem::exists(path))
        throw std::runtime_error("StlSurfaceReader: file not found: " + path);

    TopoDS_Shape shape;
    StlAPI_Reader reader;
    if (!reader.Read(shape, path.c_str()))
        throw std::runtime_error("StlSurfaceReader: failed to read " + path);
    if (shape.IsNull())
        throw std::runtime_error("StlSurfaceReader: empty shape after read");

    auto h = makeHandle(std::move(shape));
    SIMALL_LOG_INFO("CAD/STL", "Imported ", path, " (", h.topology().size(), " topo nodes)");
    return h;
}

} // namespace simall::cad
