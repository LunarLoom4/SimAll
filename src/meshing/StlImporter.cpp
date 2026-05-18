// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/StlImporter.cpp
// =============================================================================
#include "meshing/StlImporter.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace simall::meshing {

namespace {

struct VertexKey {
    std::int64_t ix, iy, iz;
    bool operator==(const VertexKey& o) const noexcept
        { return ix == o.ix && iy == o.iy && iz == o.iz; }
};
struct VertexKeyHash {
    std::size_t operator()(const VertexKey& k) const noexcept {
        std::size_t h = static_cast<std::size_t>(k.ix) * 73856093ULL;
        h ^= static_cast<std::size_t>(k.iy) * 19349663ULL;
        h ^= static_cast<std::size_t>(k.iz) * 83492791ULL;
        return h;
    }
};

std::uint32_t weld(std::unordered_map<VertexKey, std::uint32_t, VertexKeyHash>& map,
                   std::vector<util::Vec3d>& verts,
                   double tol, double x, double y, double z) {
    const double inv = 1.0 / std::max(tol, 1e-30);
    VertexKey k{ static_cast<std::int64_t>(std::llround(x * inv)),
                 static_cast<std::int64_t>(std::llround(y * inv)),
                 static_cast<std::int64_t>(std::llround(z * inv)) };
    auto it = map.find(k);
    if (it != map.end()) return it->second;
    const std::uint32_t id = static_cast<std::uint32_t>(verts.size());
    verts.push_back({x, y, z});
    map.emplace(k, id);
    return id;
}

bool load_binary(std::ifstream& f, StlSurface& out, double tol) {
    char header[80];
    f.read(header, 80);
    std::uint32_t nTri = 0;
    f.read(reinterpret_cast<char*>(&nTri), 4);
    out.vertices.clear(); out.triangles.clear(); out.normals.clear();
    out.triangles.reserve(nTri); out.normals.reserve(nTri);
    std::unordered_map<VertexKey, std::uint32_t, VertexKeyHash> map;
    map.reserve(nTri * 3);
    for (std::uint32_t t = 0; t < nTri; ++t) {
        float buf[12];
        f.read(reinterpret_cast<char*>(buf), 48);
        f.ignore(2);                              // attribute byte count
        if (!f) return false;
        out.normals.push_back({buf[0], buf[1], buf[2]});
        const std::uint32_t a = weld(map, out.vertices, tol, buf[3],  buf[4],  buf[5]);
        const std::uint32_t b = weld(map, out.vertices, tol, buf[6],  buf[7],  buf[8]);
        const std::uint32_t c = weld(map, out.vertices, tol, buf[9],  buf[10], buf[11]);
        out.triangles.push_back({a, b, c});
    }
    return true;
}

bool load_ascii(std::ifstream& f, StlSurface& out, double tol) {
    out.vertices.clear(); out.triangles.clear(); out.normals.clear();
    std::unordered_map<VertexKey, std::uint32_t, VertexKeyHash> map;
    std::string tok;
    util::Vec3d normal{0,0,0};
    std::vector<std::uint32_t> currentTri;
    while (f >> tok) {
        if (tok == "facet") {
            f >> tok;       // "normal"
            f >> normal.x >> normal.y >> normal.z;
            currentTri.clear();
        } else if (tok == "vertex") {
            double x,y,z; f >> x >> y >> z;
            currentTri.push_back(weld(map, out.vertices, tol, x, y, z));
        } else if (tok == "endfacet") {
            if (currentTri.size() == 3) {
                out.triangles.push_back({currentTri[0], currentTri[1], currentTri[2]});
                out.normals.push_back(normal);
            }
        }
    }
    return !out.triangles.empty();
}

}  // namespace

bool StlImporter::load(const std::string& path, StlSurface& out, double tol) {
    std::ifstream test(path, std::ios::binary);
    if (!test) { SIMALL_LOG_ERROR("STL", "open failed: ", path); return false; }
    // Auto-detect: ASCII starts with "solid " followed by printable text on
    // the first 80 bytes. Binary headers may also start with "solid" (bad
    // exporters!), so check whether file size matches 80 + 4 + 50·n.
    test.seekg(0, std::ios::end);
    const std::streamsize sz = test.tellg();
    test.seekg(0, std::ios::beg);
    bool isBinary = false;
    if (sz >= 84) {
        char hdr[80]; test.read(hdr, 80);
        std::uint32_t n = 0; test.read(reinterpret_cast<char*>(&n), 4);
        if (sz == static_cast<std::streamsize>(80 + 4 + 50 * std::int64_t(n))) isBinary = true;
        test.seekg(0, std::ios::beg);
    }
    bool ok = isBinary ? load_binary(test, out, tol) : load_ascii(test, out, tol);
    if (ok)
        SIMALL_LOG_INFO("STL", path, ": ", out.triangles.size(), " tris, ",
            out.vertices.size(), " welded verts (", isBinary ? "binary" : "ascii", ")");
    return ok;
}

}  // namespace simall::meshing
