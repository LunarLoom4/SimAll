// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/VtkXmlWriter.hpp
// Phase  : 22.1 — Native VTK XML output (.vtu, UnstructuredGrid).
//
// Writes the SimAll polyhedral mesh + arbitrary scalar/vector cell-data into
// VTK's modern XML format (Paraview-readable). Supports arbitrary cells via
// VTK_POLYHEDRON (cell type 42) so curved + skewed polyhedra round-trip
// correctly. Ascii encoding for portability; switch to base64-raw for big
// runs (future pass).
//
// API:
//   VtkXmlWriter w(path);
//   w.add_cell_scalar("p", &p[0]);
//   w.add_cell_vector("U", &Ux[0], &Uy[0], &Uz[0]);
//   w.write(mesh);
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"

#include <string>
#include <vector>

namespace simall::io {

class VtkXmlWriter {
public:
    explicit VtkXmlWriter(std::string path) : path_(std::move(path)) {}

    void add_cell_scalar(const std::string& name, const double* data) {
        scalars_.push_back({name, data}); }
    void add_cell_vector(const std::string& name,
                         const double* x, const double* y, const double* z) {
        vectors_.push_back({name, x, y, z}); }

    /// Serialize mesh and registered cell-data to disk. Returns true on success.
    bool write(const meshing::Mesh& mesh) const;

private:
    struct ScalarRef { std::string name; const double* data; };
    struct VectorRef { std::string name; const double* x; const double* y; const double* z; };
    std::string path_;
    std::vector<ScalarRef> scalars_;
    std::vector<VectorRef> vectors_;
};

}  // namespace simall::io
