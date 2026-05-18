// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/CgnsWriter.cpp
//
// "CGNS-native" chunked binary format used by the fallback reader.  Layout:
//
//   bytes  0.. 11   magic   "CGNS-NATIVE\n"
//   bytes 12.. 15   version (uint32 LE, currently 1)
//   bytes 16..      sequence of chunks:
//
//   chunk:
//     uint8     kind                  (1=zone,2=coords,3=section,4=bc,255=end)
//     uint32    payloadBytes
//     payload (kind-dependent)
//
//   zone payload    : u32 nameLen | char[nameLen]
//   coords payload  : u64 nNodes  | f64 x[nNodes] | f64 y[]  | f64 z[]
//   section payload : u32 nameLen | char[nameLen] | u8 elemType
//                     u64 nElems  | u32 nodes[nElems * verts] OR
//                     u64 nElems  | u32 offsets[nElems+1]    | u32 nodes[]
//   bc payload      : u32 nameLen | char[nameLen] | u32 typeLen | char[typeLen]
//                     u32 count   | u32 faceIdx[count]
// =============================================================================
#include "io/CgnsWriter.hpp"
#include "io/CgnsReader.hpp"

#include <cstring>
#include <fstream>

namespace simall::io {

namespace {

void put_u8 (std::ostream& o, std::uint8_t v)  { o.write(reinterpret_cast<const char*>(&v), 1); }
void put_u32(std::ostream& o, std::uint32_t v) { o.write(reinterpret_cast<const char*>(&v), 4); }
void put_u64(std::ostream& o, std::uint64_t v) { o.write(reinterpret_cast<const char*>(&v), 8); }
void put_str(std::ostream& o, const std::string& s) {
    put_u32(o, std::uint32_t(s.size()));
    o.write(s.data(), std::streamsize(s.size()));
}

void put_chunk(std::ostream& o, std::uint8_t kind, const std::string& payload) {
    put_u8 (o, kind);
    put_u32(o, std::uint32_t(payload.size()));
    o.write(payload.data(), std::streamsize(payload.size()));
}

std::string serialize_zone(const UnstructuredZone& z) {
    std::ostringstream os(std::ios::binary);
    put_str(os, z.name);
    return os.str();
}

std::string serialize_coords(const UnstructuredZone& z) {
    std::ostringstream os(std::ios::binary);
    const std::uint64_t n = z.x.size();
    put_u64(os, n);
    os.write(reinterpret_cast<const char*>(z.x.data()), std::streamsize(n * sizeof(double)));
    os.write(reinterpret_cast<const char*>(z.y.data()), std::streamsize(n * sizeof(double)));
    os.write(reinterpret_cast<const char*>(z.z.data()), std::streamsize(n * sizeof(double)));
    return os.str();
}

std::string serialize_section(const ElementSection& s) {
    std::ostringstream os(std::ios::binary);
    put_str(os, s.name);
    put_u8 (os, static_cast<std::uint8_t>(s.type));
    if (s.type == ElementType::Poly) {
        const std::uint64_t ne = s.polyOffsets.empty() ? 0 : s.polyOffsets.size() - 1;
        put_u64(os, ne);
        put_u32(os, std::uint32_t(s.polyOffsets.size()));
        os.write(reinterpret_cast<const char*>(s.polyOffsets.data()),
                 std::streamsize(s.polyOffsets.size() * sizeof(std::uint32_t)));
        put_u32(os, std::uint32_t(s.polyNodes.size()));
        os.write(reinterpret_cast<const char*>(s.polyNodes.data()),
                 std::streamsize(s.polyNodes.size() * sizeof(NodeIdx)));
    } else {
        put_u64(os, s.element_count());
        put_u32(os, std::uint32_t(s.nodes.size()));
        os.write(reinterpret_cast<const char*>(s.nodes.data()),
                 std::streamsize(s.nodes.size() * sizeof(NodeIdx)));
    }
    return os.str();
}

std::string serialize_bc(const BoundaryPatch& bp) {
    std::ostringstream os(std::ios::binary);
    put_str(os, bp.name);
    put_str(os, bp.bcType);
    put_u32(os, std::uint32_t(bp.faceElementIndices.size()));
    os.write(reinterpret_cast<const char*>(bp.faceElementIndices.data()),
             std::streamsize(bp.faceElementIndices.size() * sizeof(std::uint32_t)));
    return os.str();
}

}  // namespace

CgnsWriteResult write_cgns_native(const std::string& path,
                                   const ImportedMesh& mesh) {
    CgnsWriteResult r;
    r.backend = "cgns_native";
    std::ofstream f(path, std::ios::binary);
    if (!f) { r.error = "Cannot open: " + path; return r; }
    const char magic[12] = {'C','G','N','S','-','N','A','T','I','V','E','\n'};
    f.write(magic, 12);
    put_u32(f, 1);                                 // version
    for (const auto& z : mesh.zones) {
        put_chunk(f, 1, serialize_zone(z));
        put_chunk(f, 2, serialize_coords(z));
        for (const auto& s : z.sections)   put_chunk(f, 3, serialize_section(s));
        for (const auto& b : z.boundaries) put_chunk(f, 4, serialize_bc(b));
    }
    put_chunk(f, 255, {});
    r.ok = f.good();
    if (!r.ok) r.error = "I/O error writing " + path;
    return r;
}

CgnsWriteResult write_cgns(const std::string& path, const ImportedMesh& mesh) {
    // libcgns backend would slot in here when SIMALL_HAVE_CGNS is defined.
    return write_cgns_native(path, mesh);
}

}  // namespace simall::io
