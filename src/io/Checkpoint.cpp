// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/Checkpoint.cpp
// =============================================================================
#include "io/Checkpoint.hpp"

#include "core/Logger.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>

namespace simall::io
{

namespace
{

template <class T> void wb(std::ostream& s, const T& v)
{
    s.write(reinterpret_cast<const char*>(&v), sizeof(T));
}
template <class T> bool rb(std::istream& s, T& v)
{
    s.read(reinterpret_cast<char*>(&v), sizeof(T));
    return s.good();
}
void w_str(std::ostream& s, const std::string& v)
{
    std::uint32_t n = static_cast<std::uint32_t>(v.size());
    wb(s, n);
    s.write(v.data(), n);
}
bool r_str(std::istream& s, std::string& v)
{
    std::uint32_t n = 0;
    if (!rb(s, n))
        return false;
    v.resize(n);
    s.read(v.data(), n);
    return s.good();
}
template <class V> void w_vec(std::ostream& s, const V& v)
{
    std::uint64_t n = v.size();
    wb(s, n);
    if (n)
        s.write(reinterpret_cast<const char*>(v.data()), n * sizeof(typename V::value_type));
}
template <class V> bool r_vec(std::istream& s, V& v)
{
    std::uint64_t n = 0;
    if (!rb(s, n))
        return false;
    v.assign(n, typename V::value_type{});
    if (n)
        s.read(reinterpret_cast<char*>(v.data()), n * sizeof(typename V::value_type));
    return s.good();
}

} // namespace

bool Checkpoint::write(const std::string& path,
                       const meshing::Mesh& m,
                       const solver::FieldRegistry& F,
                       const CheckpointMeta& meta)
{
    const std::string tmp = path + ".tmp";
    std::ofstream s(tmp, std::ios::binary | std::ios::trunc);
    if (!s) {
        SIMALL_LOG_WARN("Checkpoint", "cannot open ", tmp);
        return false;
    }

    CheckpointHeader h;
    h.timestamp =
        static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::seconds>(
                                       std::chrono::system_clock::now().time_since_epoch())
                                       .count());
    wb(s, h);
    wb(s, meta);

    // ---- Section 1: Mesh ----
    const auto& N = m.nodes();
    const auto& Ff = m.faces();
    const auto& C = m.cells();
    w_vec(s, N.x);
    w_vec(s, N.y);
    w_vec(s, N.z);
    w_vec(s, Ff.owner);
    w_vec(s, Ff.neighbor);
    w_vec(s, Ff.areaX);
    w_vec(s, Ff.areaY);
    w_vec(s, Ff.areaZ);
    w_vec(s, Ff.centroidX);
    w_vec(s, Ff.centroidY);
    w_vec(s, Ff.centroidZ);
    w_vec(s, Ff.boundaryZone);
    w_vec(s, Ff.nodeOffsets);
    w_vec(s, Ff.nodeIndices);
    w_vec(s, C.volume);
    w_vec(s, C.centroidX);
    w_vec(s, C.centroidY);
    w_vec(s, C.centroidZ);
    w_vec(s, C.faceOffsets);
    w_vec(s, C.faceIndices);

    // ---- Section 2: Field Registry ----
    const auto& sMap = F.scalars();
    const auto& vMap = F.vectors();
    std::uint32_t nScalars = static_cast<std::uint32_t>(sMap.size());
    std::uint32_t nVectors = static_cast<std::uint32_t>(vMap.size());
    wb(s, nScalars);
    wb(s, nVectors);
    for (const auto& [name, fld] : sMap) {
        w_str(s, name);
        w_vec(s, fld);
    }
    for (const auto& [name, fld] : vMap) {
        w_str(s, name);
        w_vec(s, fld.x);
        w_vec(s, fld.y);
        w_vec(s, fld.z);
    }

    s.flush();
    s.close();
    std::error_code ec;
    std::filesystem::rename(tmp, path, ec);
    if (ec) {
        SIMALL_LOG_WARN("Checkpoint", "rename failed: ", ec.message());
        return false;
    }
    SIMALL_LOG_INFO(
        "Checkpoint", "wrote ", path, "  cells=", m.cells().size(), "  iter=", meta.iteration);
    return true;
}

bool Checkpoint::read(const std::string& path,
                      meshing::Mesh& m,
                      solver::FieldRegistry& F,
                      CheckpointMeta& meta)
{
    std::ifstream s(path, std::ios::binary);
    if (!s) {
        SIMALL_LOG_WARN("Checkpoint", "cannot open ", path);
        return false;
    }
    CheckpointHeader h;
    if (!rb(s, h))
        return false;
    if (std::memcmp(h.magic, "SIMALLCP", 8) != 0) {
        SIMALL_LOG_WARN("Checkpoint", "bad magic");
        return false;
    }
    if (h.version > 1) {
        SIMALL_LOG_WARN("Checkpoint", "version ", h.version, " not supported");
        return false;
    }
    if (!rb(s, meta))
        return false;

    auto& N = m.nodes();
    auto& Ff = m.faces();
    auto& C = m.cells();
    if (!r_vec(s, N.x) || !r_vec(s, N.y) || !r_vec(s, N.z))
        return false;
    if (!r_vec(s, Ff.owner) || !r_vec(s, Ff.neighbor))
        return false;
    if (!r_vec(s, Ff.areaX) || !r_vec(s, Ff.areaY) || !r_vec(s, Ff.areaZ))
        return false;
    if (!r_vec(s, Ff.centroidX) || !r_vec(s, Ff.centroidY) || !r_vec(s, Ff.centroidZ))
        return false;
    if (!r_vec(s, Ff.boundaryZone))
        return false;
    if (!r_vec(s, Ff.nodeOffsets) || !r_vec(s, Ff.nodeIndices))
        return false;
    if (!r_vec(s, C.volume))
        return false;
    if (!r_vec(s, C.centroidX) || !r_vec(s, C.centroidY) || !r_vec(s, C.centroidZ))
        return false;
    if (!r_vec(s, C.faceOffsets) || !r_vec(s, C.faceIndices))
        return false;

    std::uint32_t nScalars = 0, nVectors = 0;
    if (!rb(s, nScalars) || !rb(s, nVectors))
        return false;
    for (std::uint32_t i = 0; i < nScalars; ++i) {
        std::string name;
        if (!r_str(s, name))
            return false;
        std::uint64_t n = 0;
        if (!rb(s, n))
            return false;
        auto& fld = F.scalar(name, static_cast<std::size_t>(n));
        if (n)
            s.read(reinterpret_cast<char*>(fld.data()), n * sizeof(double));
        if (!s.good())
            return false;
    }
    for (std::uint32_t i = 0; i < nVectors; ++i) {
        std::string name;
        if (!r_str(s, name))
            return false;
        std::uint64_t nx = 0;
        if (!rb(s, nx))
            return false;
        auto& fld = F.vector(name, static_cast<std::size_t>(nx));
        if (nx)
            s.read(reinterpret_cast<char*>(fld.x.data()), nx * sizeof(double));
        std::uint64_t ny = 0;
        if (!rb(s, ny))
            return false;
        if (ny)
            s.read(reinterpret_cast<char*>(fld.y.data()), ny * sizeof(double));
        std::uint64_t nz = 0;
        if (!rb(s, nz))
            return false;
        if (nz)
            s.read(reinterpret_cast<char*>(fld.z.data()), nz * sizeof(double));
        if (!s.good())
            return false;
    }
    SIMALL_LOG_INFO(
        "Checkpoint", "loaded ", path, "  cells=", m.cells().size(), "  iter=", meta.iteration);
    return true;
}

} // namespace simall::io
