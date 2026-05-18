// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/VtkXmlWriter.cpp
// =============================================================================
#include "io/VtkXmlWriter.hpp"
#include "core/Logger.hpp"

#include <fstream>
#include <set>
#include <unordered_map>

namespace simall::io {

namespace {

// VTK_POLYHEDRON descriptor format per cell:
//   nFaces, [nFaceNodes, n0,n1,...], [nFaceNodes, n0,n1,...], ...
// We emit polyhedra to keep arbitrary topologies (hex/tet/prism/pyr/poly)
// round-trippable.
struct PolyDesc {
    std::vector<int> faceCounts;     // nodes per face for this cell
    std::vector<int> faceNodes;      // flat node IDs
};

void build_poly_desc(const meshing::Mesh& mesh,
                     std::size_t c, PolyDesc& d) {
    const auto& C = mesh.cells();
    const auto& F = mesh.faces();
    d.faceCounts.clear();
    d.faceNodes.clear();
    const int beg = C.faceOffsets[c], end = C.faceOffsets[c + 1];
    for (int k = beg; k < end; ++k) {
        const meshing::FaceId fid = C.faceIndices[k];
        const int nb = F.nodeOffsets[fid];
        const int ne = F.nodeOffsets[fid + 1];
        d.faceCounts.push_back(ne - nb);
        for (int q = nb; q < ne; ++q)
            d.faceNodes.push_back(static_cast<int>(F.nodeIndices[q]));
    }
}

}  // namespace

bool VtkXmlWriter::write(const meshing::Mesh& mesh) const {
    std::ofstream f(path_);
    if (!f) { SIMALL_LOG_ERROR("VTK", "open failed: ", path_); return false; }
    const auto& N = mesh.nodes();
    const auto& C = mesh.cells();
    const std::size_t nN = N.size();
    const std::size_t nC = C.size();

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    f << "<UnstructuredGrid>\n";
    f << "<Piece NumberOfPoints=\"" << nN << "\" NumberOfCells=\"" << nC << "\">\n";

    // ---- Points
    f << "<Points>\n";
    f << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < nN; ++i)
        f << N.x[i] << ' ' << N.y[i] << ' ' << N.z[i] << '\n';
    f << "</DataArray>\n</Points>\n";

    // ---- Cells (VTK_POLYHEDRON)
    f << "<Cells>\n";
    // Connectivity and offsets: VTK demands per-cell unique-point lists.
    std::vector<int> connectivity, cellOffsets;
    std::vector<int> faces, faceOffsets;
    cellOffsets.reserve(nC);
    faceOffsets.reserve(nC);
    int totalConn = 0;
    int totalFaceEntries = 0;
    PolyDesc d;
    for (std::size_t c = 0; c < nC; ++c) {
        build_poly_desc(mesh, c, d);
        // Unique node set for this cell.
        std::set<int> unique(d.faceNodes.begin(), d.faceNodes.end());
        for (int n : unique) connectivity.push_back(n);
        totalConn += static_cast<int>(unique.size());
        cellOffsets.push_back(totalConn);

        // Face stream:  nFaces, [nNodesFace, n0,n1,...], ...
        faces.push_back(static_cast<int>(d.faceCounts.size()));
        ++totalFaceEntries;
        std::size_t cur = 0;
        for (int nFn : d.faceCounts) {
            faces.push_back(nFn);
            ++totalFaceEntries;
            for (int q = 0; q < nFn; ++q) {
                faces.push_back(d.faceNodes[cur++]);
                ++totalFaceEntries;
            }
        }
        faceOffsets.push_back(totalFaceEntries);
    }

    f << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int v : connectivity) f << v << ' ';
    f << "\n</DataArray>\n";

    f << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int v : cellOffsets) f << v << ' ';
    f << "\n</DataArray>\n";

    f << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t c = 0; c < nC; ++c) f << 42 << ' ';   // VTK_POLYHEDRON
    f << "\n</DataArray>\n";

    f << "<DataArray type=\"Int32\" Name=\"faces\" format=\"ascii\">\n";
    for (int v : faces) f << v << ' ';
    f << "\n</DataArray>\n";

    f << "<DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"ascii\">\n";
    for (int v : faceOffsets) f << v << ' ';
    f << "\n</DataArray>\n";
    f << "</Cells>\n";

    // ---- Cell data
    if (!scalars_.empty() || !vectors_.empty()) {
        f << "<CellData>\n";
        for (const auto& s : scalars_) {
            f << "<DataArray type=\"Float64\" Name=\"" << s.name
              << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (std::size_t c = 0; c < nC; ++c) f << s.data[c] << ' ';
            f << "\n</DataArray>\n";
        }
        for (const auto& v : vectors_) {
            f << "<DataArray type=\"Float64\" Name=\"" << v.name
              << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (std::size_t c = 0; c < nC; ++c)
                f << v.x[c] << ' ' << v.y[c] << ' ' << v.z[c] << '\n';
            f << "</DataArray>\n";
        }
        f << "</CellData>\n";
    }

    f << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    SIMALL_LOG_INFO("VTK", "wrote ", path_, " (", nC, " cells, ", nN, " points)");
    return true;
}

}  // namespace simall::io
