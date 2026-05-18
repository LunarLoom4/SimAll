// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/Plot3dReader.cpp
// =============================================================================
#include "io/Plot3dReader.hpp"

#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>

namespace simall::io {

namespace {

bool sniff_is_formatted(const std::string& bytes) {
    const std::size_t n = std::min<std::size_t>(bytes.size(), 128);
    for (std::size_t i = 0; i < n; ++i) {
        const unsigned char c = static_cast<unsigned char>(bytes[i]);
        if (c == '\n' || c == '\r' || c == '\t' || c == ' ') continue;
        if (c < 32 || c > 126) return false;
    }
    return n > 0;
}

}  // namespace

Plot3dReadResult parse_plot3d_formatted(const std::string& text,
                                         std::string sourceHint) {
    Plot3dReadResult r;
    r.grid.sourceFormat = "plot3d_formatted";
    r.grid.sourcePath   = std::move(sourceHint);

    std::istringstream is(text);
    int nblocks = 0;
    if (!(is >> nblocks) || nblocks <= 0) {
        r.error = "Plot3D: expected block count as first token";
        return r;
    }
    r.grid.blocks.resize(std::size_t(nblocks));
    for (int b = 0; b < nblocks; ++b) {
        if (!(is >> r.grid.blocks[std::size_t(b)].ni
                  >> r.grid.blocks[std::size_t(b)].nj
                  >> r.grid.blocks[std::size_t(b)].nk)) {
            r.error = "Plot3D: missing ni/nj/nk for block " + std::to_string(b);
            return r;
        }
        r.grid.blocks[std::size_t(b)].name = "block_" + std::to_string(b);
    }
    for (auto& bk : r.grid.blocks) {
        const std::size_t N = bk.point_count();
        bk.x.resize(N); bk.y.resize(N); bk.z.resize(N);
        for (std::size_t i = 0; i < N; ++i) if (!(is >> bk.x[i])) {
            r.error = "Plot3D: short read for X in " + bk.name;
            return r;
        }
        for (std::size_t i = 0; i < N; ++i) if (!(is >> bk.y[i])) {
            r.error = "Plot3D: short read for Y in " + bk.name;
            return r;
        }
        for (std::size_t i = 0; i < N; ++i) if (!(is >> bk.z[i])) {
            r.error = "Plot3D: short read for Z in " + bk.name;
            return r;
        }
    }
    r.ok = true;
    return r;
}

namespace {

// Skip the 4-byte Fortran record header/footer; tolerant to both endians
// (we just trust the recorded length if it matches what we expect to read).
bool skip_record_marker(std::istream& is, std::uint32_t& outLen) {
    char buf[4];
    if (!is.read(buf, 4)) return false;
    std::memcpy(&outLen, buf, 4);
    return true;
}

template <typename T>
bool read_record(std::istream& is, T* dst, std::size_t count) {
    std::uint32_t hdr = 0, ftr = 0;
    if (!skip_record_marker(is, hdr)) return false;
    if (!is.read(reinterpret_cast<char*>(dst),
                 std::streamsize(count * sizeof(T)))) return false;
    if (!skip_record_marker(is, ftr)) return false;
    return hdr == ftr;
}

Plot3dReadResult read_plot3d_binary(const std::string& path,
                                     const std::string& bytes,
                                     bool dp) {
    Plot3dReadResult r;
    r.grid.sourceFormat = "plot3d_binary";
    r.grid.sourcePath   = path;
    std::istringstream is(bytes, std::ios::binary);

    std::int32_t nblocks = 0;
    if (!read_record(is, &nblocks, 1) || nblocks <= 0) {
        r.error = "Plot3D binary: bad block-count record";
        return r;
    }
    r.grid.blocks.resize(std::size_t(nblocks));
    std::vector<std::int32_t> dims(std::size_t(nblocks) * 3);
    if (!read_record(is, dims.data(), dims.size())) {
        r.error = "Plot3D binary: bad dimensions record";
        return r;
    }
    for (std::int32_t b = 0; b < nblocks; ++b) {
        auto& bk = r.grid.blocks[std::size_t(b)];
        bk.ni = std::uint32_t(dims[std::size_t(b)*3 + 0]);
        bk.nj = std::uint32_t(dims[std::size_t(b)*3 + 1]);
        bk.nk = std::uint32_t(dims[std::size_t(b)*3 + 2]);
        bk.name = "block_" + std::to_string(b);
    }
    for (auto& bk : r.grid.blocks) {
        const std::size_t N = bk.point_count();
        bk.x.resize(N); bk.y.resize(N); bk.z.resize(N);
        if (dp) {
            if (!read_record(is, bk.x.data(), N) ||
                !read_record(is, bk.y.data(), N) ||
                !read_record(is, bk.z.data(), N)) {
                r.error = "Plot3D binary: short coord record for " + bk.name;
                return r;
            }
        } else {
            std::vector<float> tmp(N);
            if (!read_record(is, tmp.data(), N)) {
                r.error = "Plot3D binary: short X record for " + bk.name; return r;
            }
            for (std::size_t i = 0; i < N; ++i) bk.x[i] = double(tmp[i]);
            if (!read_record(is, tmp.data(), N)) {
                r.error = "Plot3D binary: short Y record for " + bk.name; return r;
            }
            for (std::size_t i = 0; i < N; ++i) bk.y[i] = double(tmp[i]);
            if (!read_record(is, tmp.data(), N)) {
                r.error = "Plot3D binary: short Z record for " + bk.name; return r;
            }
            for (std::size_t i = 0; i < N; ++i) bk.z[i] = double(tmp[i]);
        }
    }
    r.ok = true;
    return r;
}

}  // namespace

Plot3dReadResult read_plot3d(const std::string& path, Plot3dOptions opts) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        Plot3dReadResult r;
        r.error = "Cannot open: " + path;
        return r;
    }
    std::ostringstream ss;
    ss << f.rdbuf();
    const std::string bytes = ss.str();

    bool formatted;
    if (opts.forceFormatted) formatted = true;
    else if (opts.forceBinary) formatted = false;
    else formatted = sniff_is_formatted(bytes);

    if (formatted) {
        auto r = parse_plot3d_formatted(bytes, path);
        r.grid.sourceFormat = "plot3d_formatted";
        return r;
    }
    return read_plot3d_binary(path, bytes, opts.doublePrecision);
}

}  // namespace simall::io
