// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/ProjectSerializer.hpp
// Week   : 17
//
// `.simall` project format (Section 16 of the ultra-detailed spec).
// Chunked binary container — magic + manifest + arbitrary chunks.  The
// project owns references to imported mesh / physics / BCs / etc.; the
// serializer just round-trips a name-keyed chunk store.
//
// Chunk layout:
//   bytes  0..12   magic        "SIMALL-PROJ\n"
//   bytes 13..16   version      uint32 LE  (currently 1)
//   bytes 17..20   chunkCount   uint32 LE
//   <chunk>{chunkCount}
//
//   chunk: u32 nameLen | char[nameLen] | u32 mimeLen | char[mimeLen]
//          u64 payloadBytes | u8[payloadBytes]
// =============================================================================
#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace simall::io
{

struct ProjectChunk
{
    std::string mime; // "application/vnd.simall.mesh.cgns" etc.
    std::vector<std::uint8_t> bytes;
};

class ProjectSerializer
{
public:
    void put(std::string name, ProjectChunk chunk);
    [[nodiscard]] const ProjectChunk* get(std::string_view name) const noexcept;
    [[nodiscard]] std::vector<std::string> names() const;
    [[nodiscard]] std::size_t size() const noexcept { return chunks_.size(); }
    void clear();

    bool save(const std::string& path) const;
    bool load(const std::string& path);

private:
    std::map<std::string, ProjectChunk> chunks_;
};

} // namespace simall::io
