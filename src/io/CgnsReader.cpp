// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/CgnsReader.cpp
//
// Companion reader for CgnsWriter.cpp's "CGNS-native" chunked format.
// =============================================================================
#include "io/CgnsReader.hpp"

#include <cstring>
#include <fstream>
#include <sstream>

namespace simall::io {

bool cgns_libcgns_available() noexcept {
#ifdef SIMALL_HAVE_CGNS
    return true;
#else
    return false;
#endif
}

namespace {

struct Cursor {
    const std::string& buf;
    std::size_t        i = 0;

    bool need(std::size_t n) const { return i + n <= buf.size(); }
    template<typename T> bool read(T& out) {
        if (!need(sizeof(T))) return false;
        std::memcpy(&out, buf.data() + i, sizeof(T));
        i += sizeof(T);
        return true;
    }
    bool read_bytes(void* dst, std::size_t n) {
        if (!need(n)) return false;
        std::memcpy(dst, buf.data() + i, n);
        i += n;
        return true;
    }
    bool read_str(std::string& out) {
        std::uint32_t n = 0;
        if (!read(n)) return false;
        if (!need(n)) return false;
        out.assign(buf.data() + i, n);
        i += n;
        return true;
    }
};

bool read_zone(Cursor& c, UnstructuredZone& z) {
    return c.read_str(z.name);
}
bool read_coords(Cursor& c, UnstructuredZone& z) {
    std::uint64_t n = 0;
    if (!c.read(n)) return false;
    z.x.resize(n); z.y.resize(n); z.z.resize(n);
    return c.read_bytes(z.x.data(), n * sizeof(double))
        && c.read_bytes(z.y.data(), n * sizeof(double))
        && c.read_bytes(z.z.data(), n * sizeof(double));
}
bool read_section(Cursor& c, ElementSection& s) {
    if (!c.read_str(s.name)) return false;
    std::uint8_t et = 0;
    if (!c.read(et)) return false;
    s.type = static_cast<ElementType>(et);
    std::uint64_t ne = 0;
    if (!c.read(ne)) return false;
    if (s.type == ElementType::Poly) {
        std::uint32_t no = 0;
        if (!c.read(no)) return false;
        s.polyOffsets.resize(no);
        if (!c.read_bytes(s.polyOffsets.data(), no * sizeof(std::uint32_t))) return false;
        std::uint32_t nn = 0;
        if (!c.read(nn)) return false;
        s.polyNodes.resize(nn);
        return c.read_bytes(s.polyNodes.data(), nn * sizeof(NodeIdx));
    }
    std::uint32_t nn = 0;
    if (!c.read(nn)) return false;
    s.nodes.resize(nn);
    return c.read_bytes(s.nodes.data(), nn * sizeof(NodeIdx));
}
bool read_bc(Cursor& c, BoundaryPatch& bp) {
    if (!c.read_str(bp.name)) return false;
    if (!c.read_str(bp.bcType)) return false;
    std::uint32_t n = 0;
    if (!c.read(n)) return false;
    bp.faceElementIndices.resize(n);
    return c.read_bytes(bp.faceElementIndices.data(), n * sizeof(std::uint32_t));
}

}  // namespace

CgnsReadResult read_cgns(const std::string& path) {
    CgnsReadResult r;
    r.backend = "cgns_native";
    std::ifstream f(path, std::ios::binary);
    if (!f) { r.error = "Cannot open: " + path; return r; }
    std::ostringstream ss; ss << f.rdbuf();
    const std::string bytes = ss.str();

    if (bytes.size() < 16
        || std::memcmp(bytes.data(), "CGNS-NATIVE\n", 12) != 0) {
        r.error = "Not a CGNS-native file (and libcgns backend not built in)";
        return r;
    }
    std::uint32_t version = 0;
    std::memcpy(&version, bytes.data() + 12, 4);
    if (version != 1) {
        r.error = "Unsupported CGNS-native version: " + std::to_string(version);
        return r;
    }

    Cursor c{bytes};
    c.i = 16;
    UnstructuredZone* current = nullptr;
    while (c.i < bytes.size()) {
        std::uint8_t  kind = 0;
        std::uint32_t plen = 0;
        if (!c.read(kind)) { r.error = "Truncated chunk header"; return r; }
        if (!c.read(plen)) { r.error = "Truncated chunk length"; return r; }
        if (!c.need(plen))  { r.error = "Truncated chunk payload"; return r; }
        std::string payload(bytes.data() + c.i, plen);
        c.i += plen;
        Cursor pc{payload};
        switch (kind) {
            case 1: {
                r.mesh.zones.emplace_back();
                current = &r.mesh.zones.back();
                if (!read_zone(pc, *current)) { r.error = "Bad zone chunk"; return r; }
                break;
            }
            case 2: {
                if (!current) { r.error = "Coords before zone"; return r; }
                if (!read_coords(pc, *current)) { r.error = "Bad coords chunk"; return r; }
                break;
            }
            case 3: {
                if (!current) { r.error = "Section before zone"; return r; }
                ElementSection s;
                if (!read_section(pc, s)) { r.error = "Bad section chunk"; return r; }
                current->sections.push_back(std::move(s));
                break;
            }
            case 4: {
                if (!current) { r.error = "BC before zone"; return r; }
                BoundaryPatch bp;
                if (!read_bc(pc, bp)) { r.error = "Bad BC chunk"; return r; }
                current->boundaries.push_back(std::move(bp));
                break;
            }
            case 255: c.i = bytes.size(); break;
            default:  break;        // forward-compatible skip
        }
    }
    r.mesh.sourceFormat = "cgns_native";
    r.mesh.sourcePath   = path;
    r.ok = true;
    return r;
}

}  // namespace simall::io
