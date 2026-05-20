// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/EnSightGoldWriter.cpp
//
// EnSight Gold C-Binary uses 80-byte string blocks, little-endian 32-bit
// ints and floats. Per the CEI specification, geometry begins with two
// 80-byte description lines, the "node id assign" / "element id assign"
// flags, "extents" line, then per-part: "part" header, "coordinates" block
// with x[],y[],z[] as float32, then element blocks (we use "nfaced" for
// fully polyhedral generality).
// =============================================================================
#include "io/EnSightGoldWriter.hpp"

#include "core/Logger.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>

namespace simall::io
{

namespace
{
void write_string80(std::ostream& s, const std::string& str)
{
    std::array<char, 80> buf{};
    std::memcpy(buf.data(), str.c_str(), std::min<std::size_t>(79, str.size()));
    s.write(buf.data(), 80);
}
void write_i32(std::ostream& s, std::int32_t v)
{
    s.write(reinterpret_cast<const char*>(&v), 4);
}
void write_f32(std::ostream& s, float v)
{
    s.write(reinterpret_cast<const char*>(&v), 4);
}
std::string step_suffix(int step)
{
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%05d", step);
    return std::string(buf);
}
} // namespace

bool EnSightGoldWriter::open(const std::string& caseDir, const std::string& caseName)
{
    caseDir_ = caseDir;
    caseName_ = caseName;
    std::error_code ec;
    std::filesystem::create_directories(caseDir_, ec);
    if (ec) {
        SIMALL_LOG_ERROR("EnSight", "create_directories failed: ", ec.message());
        return false;
    }
    times_.clear();
    scalarVars_.clear();
    vectorVars_.clear();
    geometryWritten_ = false;
    isOpen_ = true;
    return true;
}

bool EnSightGoldWriter::write_geometry(const meshing::Mesh& mesh)
{
    const std::string path = caseDir_ + "/" + caseName_ + ".geo";
    std::ofstream s(path, std::ios::binary);
    if (!s) {
        SIMALL_LOG_ERROR("EnSight", "open geometry: ", path);
        return false;
    }
    write_string80(s, "C Binary");
    write_string80(s, "SimAll Beta export");
    write_string80(s, caseName_);
    write_string80(s, "node id off");
    write_string80(s, "element id off");
    write_string80(s, "part");
    write_i32(s, 1);
    write_string80(s, "fluid");
    write_string80(s, "coordinates");
    const auto& N = mesh.nodes();
    const auto& F = mesh.faces();
    const auto& C = mesh.cells();
    write_i32(s, static_cast<std::int32_t>(N.size()));
    for (std::size_t i = 0; i < N.size(); ++i)
        write_f32(s, static_cast<float>(N.x[i]));
    for (std::size_t i = 0; i < N.size(); ++i)
        write_f32(s, static_cast<float>(N.y[i]));
    for (std::size_t i = 0; i < N.size(); ++i)
        write_f32(s, static_cast<float>(N.z[i]));

    // Polyhedral element block: "nfaced".
    write_string80(s, "nfaced");
    write_i32(s, static_cast<std::int32_t>(C.size()));
    // For each cell: number of faces.
    for (std::size_t c = 0; c < C.size(); ++c) {
        const int nF = C.faceOffsets[c + 1] - C.faceOffsets[c];
        write_i32(s, nF);
    }
    // For each face of each cell: number of nodes.
    for (std::size_t c = 0; c < C.size(); ++c) {
        for (int k = C.faceOffsets[c]; k < C.faceOffsets[c + 1]; ++k) {
            const auto fid = C.faceIndices[k];
            const int nN = F.nodeOffsets[fid + 1] - F.nodeOffsets[fid];
            write_i32(s, nN);
        }
    }
    // For each face of each cell: node ids (1-based per EnSight convention).
    for (std::size_t c = 0; c < C.size(); ++c) {
        for (int k = C.faceOffsets[c]; k < C.faceOffsets[c + 1]; ++k) {
            const auto fid = C.faceIndices[k];
            for (int n = F.nodeOffsets[fid]; n < F.nodeOffsets[fid + 1]; ++n) {
                write_i32(s, static_cast<std::int32_t>(F.nodeIndices[n] + 1));
            }
        }
    }
    geometryWritten_ = true;
    return static_cast<bool>(s);
}

bool EnSightGoldWriter::write_variable_scalar(const std::string& name,
                                              const solver::ScalarField& f,
                                              int step)
{
    const std::string path = caseDir_ + "/" + caseName_ + "." + name + "." + step_suffix(step);
    std::ofstream s(path, std::ios::binary);
    if (!s)
        return false;
    write_string80(s, "Per element scalar: " + name);
    write_string80(s, "part");
    write_i32(s, 1);
    write_string80(s, "nfaced");
    for (double v : f)
        write_f32(s, static_cast<float>(v));
    return static_cast<bool>(s);
}

bool EnSightGoldWriter::write_variable_vector(const std::string& name,
                                              const solver::VectorField& f,
                                              int step)
{
    const std::string path = caseDir_ + "/" + caseName_ + "." + name + "." + step_suffix(step);
    std::ofstream s(path, std::ios::binary);
    if (!s)
        return false;
    write_string80(s, "Per element vector: " + name);
    write_string80(s, "part");
    write_i32(s, 1);
    write_string80(s, "nfaced");
    for (double v : f.x)
        write_f32(s, static_cast<float>(v));
    for (double v : f.y)
        write_f32(s, static_cast<float>(v));
    for (double v : f.z)
        write_f32(s, static_cast<float>(v));
    return static_cast<bool>(s);
}

bool EnSightGoldWriter::flush_case_file()
{
    const std::string path = caseDir_ + "/" + caseName_ + ".case";
    std::ofstream s(path);
    if (!s)
        return false;
    s << "FORMAT\n type: ensight gold\n\n";
    s << "GEOMETRY\n model: 1 " << caseName_ << ".geo\n\n";
    s << "VARIABLE\n";
    for (const auto& v : scalarVars_)
        s << " scalar per element: 1 " << v << " " << caseName_ << "." << v << ".*****\n";
    for (const auto& v : vectorVars_)
        s << " vector per element: 1 " << v << " " << caseName_ << "." << v << ".*****\n";
    s << "\nTIME\n";
    s << " time set: 1\n";
    s << " number of steps: " << times_.size() << "\n";
    s << " filename start number: 0\n";
    s << " filename increment: 1\n";
    s << " time values:\n";
    for (double t : times_)
        s << "  " << t << "\n";
    return static_cast<bool>(s);
}

bool EnSightGoldWriter::write_step(const meshing::Mesh& mesh,
                                   const solver::FieldRegistry& fields,
                                   double t)
{
    if (!isOpen_)
        return false;
    if (!geometryWritten_) {
        if (!write_geometry(mesh))
            return false;
    }
    const int step = static_cast<int>(times_.size());
    for (const auto& [name, fld] : fields.scalars()) {
        if (!write_variable_scalar(name, fld, step))
            return false;
        if (std::find(scalarVars_.begin(), scalarVars_.end(), name) == scalarVars_.end())
            scalarVars_.push_back(name);
    }
    for (const auto& [name, fld] : fields.vectors()) {
        if (!write_variable_vector(name, fld, step))
            return false;
        if (std::find(vectorVars_.begin(), vectorVars_.end(), name) == vectorVars_.end())
            vectorVars_.push_back(name);
    }
    times_.push_back(t);
    return flush_case_file();
}

bool EnSightGoldWriter::close()
{
    if (!isOpen_)
        return false;
    const bool ok = flush_case_file();
    isOpen_ = false;
    return ok;
}

} // namespace simall::io
