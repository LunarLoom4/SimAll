// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/ObjReader.cpp
// Phase  : 4.1
// =============================================================================
#include "cad/ObjReader.hpp"
#include "cad/ShapeHandleInternal.hpp"
#include "core/Logger.hpp"

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace simall::cad {

namespace {

struct VertexKey {
    long long ix, iy, iz;
    bool operator==(const VertexKey& o) const noexcept {
        return ix == o.ix && iy == o.iy && iz == o.iz;
    }
};
struct VertexKeyHash {
    std::size_t operator()(const VertexKey& k) const noexcept {
        std::size_t h = std::hash<long long>{}(k.ix);
        h ^= std::hash<long long>{}(k.iy) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<long long>{}(k.iz) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

inline void ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
        [](unsigned char c){ return !std::isspace(c); }));
}

// Parse an OBJ "f" token of forms "v", "v/vt", "v//vn", or "v/vt/vn" and
// return the 1-based vertex index (negative indices wrap from the end).
int parse_face_index(const std::string& tok, int nVerts) {
    if (tok.empty()) return 0;
    int i = 0;
    try { i = std::stoi(tok); } catch (...) { return 0; }
    if (i < 0) i = nVerts + 1 + i;
    return i;
}

}  // namespace

TriangleMesh ObjReader::readTriangles(const std::string& path,
                                      const ObjReadOptions& opts) {
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("ObjReader: cannot open " + path);

    std::vector<util::Vec3d>             raw;
    std::vector<std::array<int, 3>>      tris;
    raw.reserve(4096);
    tris.reserve(8192);

    std::string line;
    while (std::getline(in, line)) {
        ltrim(line);
        if (line.empty() || line[0] == '#') continue;

        if (line[0] == 'v' && line.size() > 1 &&
            std::isspace(static_cast<unsigned char>(line[1]))) {
            std::istringstream is(line.substr(1));
            double x = 0, y = 0, z = 0;
            if (is >> x >> y >> z) raw.emplace_back(x, y, z);
            continue;
        }
        if (line[0] == 'f' && line.size() > 1 &&
            std::isspace(static_cast<unsigned char>(line[1]))) {
            std::istringstream is(line.substr(1));
            std::vector<int> idx;
            std::string tok;
            while (is >> tok) {
                const auto slash = tok.find('/');
                const std::string vTok =
                    (slash == std::string::npos) ? tok : tok.substr(0, slash);
                int i = parse_face_index(vTok, static_cast<int>(raw.size()));
                if (i > 0 && static_cast<std::size_t>(i) <= raw.size())
                    idx.push_back(i - 1);
            }
            if (idx.size() < 3) continue;
            if (idx.size() == 3) {
                tris.push_back({idx[0], idx[1], idx[2]});
            } else if (opts.triangulateNGons) {
                for (std::size_t k = 1; k + 1 < idx.size(); ++k)
                    tris.push_back({idx[0], idx[k], idx[k + 1]});
            }
        }
    }

    TriangleMesh mesh;
    if (!opts.mergeCoincidentVertices || opts.mergeTolerance <= 0.0) {
        mesh.points = std::move(raw);
        mesh.triangles.reserve(tris.size());
        for (const auto& t : tris) {
            mesh.triangles.push_back({std::uint32_t(t[0]),
                                      std::uint32_t(t[1]),
                                      std::uint32_t(t[2])});
        }
    } else {
        std::unordered_map<VertexKey, std::uint32_t, VertexKeyHash> table;
        table.reserve(raw.size() * 2);
        const double inv = 1.0 / opts.mergeTolerance;
        std::vector<std::uint32_t> remap(raw.size(), 0);
        for (std::size_t i = 0; i < raw.size(); ++i) {
            VertexKey k{static_cast<long long>(std::llround(raw[i].x * inv)),
                        static_cast<long long>(std::llround(raw[i].y * inv)),
                        static_cast<long long>(std::llround(raw[i].z * inv))};
            auto it = table.find(k);
            if (it == table.end()) {
                const auto idx = static_cast<std::uint32_t>(mesh.points.size());
                table.emplace(k, idx);
                mesh.points.push_back(raw[i]);
                remap[i] = idx;
            } else {
                remap[i] = it->second;
            }
        }
        for (const auto& t : tris) {
            const std::uint32_t a = remap[std::size_t(t[0])];
            const std::uint32_t b = remap[std::size_t(t[1])];
            const std::uint32_t c = remap[std::size_t(t[2])];
            if (a == b || b == c || a == c) continue;        // degenerate
            mesh.triangles.push_back({a, b, c});
        }
    }
    mesh.triangleFaceId.assign(mesh.triangles.size(), util::PersistentId{1});

    (void)opts.recomputeNormals;  // OBJ normals not stored in TriangleMesh
    SIMALL_LOG_INFO("CAD/OBJ", "Read ", path, " — ", mesh.points.size(),
                    " verts, ", mesh.triangles.size(), " tris");
    return mesh;
}

ShapeHandle ObjReader::read(const std::string& path, const ObjReadOptions& opts) {
    TriangleMesh tri = readTriangles(path, opts);

    BRep_Builder builder;
    TopoDS_Compound comp;
    builder.MakeCompound(comp);

    for (const auto& t : tri.triangles) {
        const auto& p0 = tri.points[t[0]];
        const auto& p1 = tri.points[t[1]];
        const auto& p2 = tri.points[t[2]];
        try {
            BRepBuilderAPI_MakePolygon poly(
                gp_Pnt(p0.x, p0.y, p0.z),
                gp_Pnt(p1.x, p1.y, p1.z),
                gp_Pnt(p2.x, p2.y, p2.z),
                /*close=*/true);
            if (!poly.IsDone()) continue;
            BRepBuilderAPI_MakeFace face(poly.Wire(), /*onlyPlane=*/true);
            if (face.IsDone()) builder.Add(comp, face.Face());
        } catch (...) {
            // ignore one bad triangle, keep going
        }
    }

    auto h = makeHandle(std::move(comp));
    SIMALL_LOG_INFO("CAD/OBJ", "Imported ", path, " (",
                    h.topology().size(), " topo nodes)");
    return h;
}

}  // namespace simall::cad
