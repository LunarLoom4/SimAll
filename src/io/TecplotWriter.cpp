// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/TecplotWriter.cpp
// =============================================================================
#include "io/TecplotWriter.hpp"

#include <fstream>
#include <sstream>

namespace simall::io {

namespace {

const char* zone_etype(ElementType t) {
    switch (t) {
        case ElementType::Tri3:   return "FETRIANGLE";
        case ElementType::Quad4:  return "FEQUADRILATERAL";
        case ElementType::Tetra4: return "FETETRAHEDRON";
        case ElementType::Hexa8:  return "FEBRICK";
        case ElementType::Penta6: return "FEBRICK";        // Tecplot collapses prisms onto a hex
        case ElementType::Pyra5:  return "FEBRICK";
        case ElementType::Poly:   return "FEPOLYHEDRON";
        default:                  return "FETRIANGLE";
    }
}

void write_block_doubles(std::ofstream& f, const std::vector<double>& v) {
    int col = 0;
    for (double x : v) {
        f << x;
        if (++col == 5) { f << '\n'; col = 0; } else { f << ' '; }
    }
    if (col) f << '\n';
}

}  // namespace

TecplotWriteResult write_tecplot_ascii(const std::string& path,
                                        const ImportedMesh& mesh,
                                        const FieldFrame* fields,
                                        const std::string& title) {
    TecplotWriteResult r;
    std::ofstream f(path);
    if (!f) { r.error = "Cannot open: " + path; return r; }

    f << "TITLE = \"" << title << "\"\n";
    f << "VARIABLES = \"X\" \"Y\" \"Z\"";
    if (fields) for (const auto& c : fields->components) f << " \"" << c.name << "\"";
    f << '\n';

    for (std::size_t zi = 0; zi < mesh.zones.size(); ++zi) {
        const auto& z = mesh.zones[zi];

        // Determine dominant element type from the largest volume section.
        const ElementSection* dominant = nullptr;
        std::size_t           bestCount = 0;
        for (const auto& s : z.sections) {
            const auto ne = s.element_count();
            if (ne > bestCount) { bestCount = ne; dominant = &s; }
        }
        if (!dominant) continue;

        f << "ZONE T=\"" << z.name << "\"\n";
        f << " NODES=" << z.x.size() << " ELEMENTS=" << dominant->element_count() << "\n";
        f << " DATAPACKING=BLOCK ZONETYPE=" << zone_etype(dominant->type) << "\n";

        // Mark every field as CELLCENTERED (X/Y/Z always at vertices).
        if (fields && !fields->components.empty()) {
            f << " VARLOCATION=(";
            for (std::size_t k = 0; k < fields->components.size(); ++k) {
                if (k) f << ",";
                f << "[" << (4 + k) << "]=CELLCENTERED";
            }
            f << ")\n";
        }
        f << " SOLUTIONTIME=" << (fields ? fields->time : 0.0)
          << " STRANDID="     << (zi + 1) << "\n";

        write_block_doubles(f, z.x);
        write_block_doubles(f, z.y);
        write_block_doubles(f, z.z);
        if (fields) for (const auto& c : fields->components) write_block_doubles(f, c.values);

        // Connectivity (1-based for Tecplot).
        const std::uint8_t v = vertices_per_element(dominant->type);
        if (v > 0) {
            for (std::size_t e = 0; e < dominant->element_count(); ++e) {
                for (std::uint8_t k = 0; k < v; ++k)
                    f << (dominant->nodes[e * v + k] + 1) << ' ';
                f << '\n';
            }
        }
    }
    r.ok = f.good();
    if (!r.ok) r.error = "I/O error writing " + path;
    return r;
}

}  // namespace simall::io
