// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/Hdf5ResultStore.hpp
// Week   : 17
//
// Hierarchical solver-result store.  Models the HDF5 group / dataset /
// attribute trio so application code reads the same way regardless of
// which backend is compiled in:
//
//   * SIMALL_HAVE_HDF5=1 — `read`/`write` round-trip a real `.h5` file
//     using libhdf5 (deferred plumbing; not part of the W17 scope but the
//     interface is HDF5-shaped so the slot-in is mechanical).
//
//   * Native fallback — a self-contained "SimAll Result Store" (SRS)
//     chunked binary file (magic `SIMALL-SRS\n`).  Hierarchical paths
//     ("/run/iter_000123/cellPressure"), little-endian, no compression,
//     fully round-trippable.  All tests target this backend.
//
// The store is intended for end-of-step writes and ad-hoc reads (post-
// processing, restart).  Streaming output uses Checkpoint.cpp instead.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace simall::io {

enum class StoreDtype : std::uint8_t { F64 = 1, F32 = 2, I32 = 3, I64 = 4, U8 = 5 };

struct Dataset {
    StoreDtype                 dtype = StoreDtype::F64;
    std::vector<std::uint64_t> shape;       // length = rank
    std::vector<std::uint8_t>  bytes;       // packed raw values (LE)
    [[nodiscard]] std::uint64_t element_count() const noexcept;
    [[nodiscard]] std::size_t  element_size()  const noexcept;
};

struct Group {
    std::map<std::string, std::string>             attributes;
    std::map<std::string, Dataset>                 datasets;
    std::map<std::string, std::unique_ptr<Group>>  children;

    Group*       create_group(const std::string& name);
    Group*       find_group  (const std::string& path);          // "/run/iter_000123"
    Dataset*     find_dataset(const std::string& path);
    void         set_dataset (const std::string& path, Dataset ds);
    void         set_attribute(const std::string& path,
                                const std::string& key,
                                const std::string& value);
};

class Hdf5ResultStore {
public:
    Hdf5ResultStore();
    ~Hdf5ResultStore();

    [[nodiscard]] Group& root() { return *root_; }

    // -- typed convenience helpers ------------------------------------------
    void write_field(const std::string& groupPath,
                     const std::string& name,
                     const std::vector<double>& values,
                     FieldLocation loc = FieldLocation::Cell);
    [[nodiscard]] std::vector<double> read_field_f64(
                     const std::string& groupPath,
                     const std::string& name) const;

    // -- file IO -------------------------------------------------------------
    bool save_native(const std::string& path) const;
    bool load_native(const std::string& path);

    // Returns true when the libhdf5 backend is compiled in.
    [[nodiscard]] static bool has_hdf5_backend() noexcept;

private:
    std::unique_ptr<Group> root_;
};

}  // namespace simall::io
