// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/Checkpoint.hpp
// Phase  : 17.3 — Binary checkpoint / restart of a running simulation.
//
// Captures the FULL solver state needed to resume bit-identical execution:
//   - Mesh (nodes, faces, cells, connectivity)
//   - All scalar and vector fields in a FieldRegistry
//   - Iteration counter and simulation time
//   - Solver settings hash (consistency check on load)
//
// Format: little-endian binary with a 16-byte magic header, version,
// section table, then payload blocks. Each block is prefixed with its
// byte length so unknown blocks can be skipped (forward compatibility).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <cstdint>
#include <string>

namespace simall::io {

struct CheckpointHeader {
    char        magic[8] = {'S','I','M','A','L','L','C','P'};
    std::uint32_t version = 1;
    std::uint32_t flags   = 0;
    std::uint64_t timestamp = 0;
};

struct CheckpointMeta {
    int     iteration = 0;
    double  time      = 0.0;
    std::uint64_t settingsHash = 0;
};

class Checkpoint {
public:
    /// Atomic write: first to `.tmp`, then rename over `path`.
    static bool write(const std::string& path,
                      const meshing::Mesh& mesh,
                      const solver::FieldRegistry& fields,
                      const CheckpointMeta& meta);

    /// Restore mesh & fields from a checkpoint file. Returns false on any
    /// structural error (magic mismatch, version too new, truncated file).
    static bool read(const std::string& path,
                     meshing::Mesh& mesh,
                     solver::FieldRegistry& fields,
                     CheckpointMeta& meta);
};

}  // namespace simall::io
