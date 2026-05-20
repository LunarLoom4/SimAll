// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/ProjectSerializer.cpp
// =============================================================================
#include "io/ProjectSerializer.hpp"

#include <cstring>
#include <fstream>
#include <sstream>

namespace simall::io
{

namespace
{
void put_u32(std::ostream& o, std::uint32_t v)
{
    o.write(reinterpret_cast<const char*>(&v), 4);
}
void put_u64(std::ostream& o, std::uint64_t v)
{
    o.write(reinterpret_cast<const char*>(&v), 8);
}
void put_str(std::ostream& o, const std::string& s)
{
    put_u32(o, std::uint32_t(s.size()));
    o.write(s.data(), std::streamsize(s.size()));
}
} // namespace

void ProjectSerializer::put(std::string name, ProjectChunk c)
{
    chunks_[std::move(name)] = std::move(c);
}

const ProjectChunk* ProjectSerializer::get(std::string_view name) const noexcept
{
    auto it = chunks_.find(std::string(name));
    return it == chunks_.end() ? nullptr : &it->second;
}

std::vector<std::string> ProjectSerializer::names() const
{
    std::vector<std::string> out;
    out.reserve(chunks_.size());
    for (auto& [n, _] : chunks_)
        out.push_back(n);
    return out;
}

void ProjectSerializer::clear()
{
    chunks_.clear();
}

bool ProjectSerializer::save(const std::string& path) const
{
    std::ofstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.write("SIMALL-PROJ\n", 12);
    put_u32(f, 1);
    put_u32(f, std::uint32_t(chunks_.size()));
    for (auto& [name, c] : chunks_) {
        put_str(f, name);
        put_str(f, c.mime);
        put_u64(f, c.bytes.size());
        if (!c.bytes.empty())
            f.write(reinterpret_cast<const char*>(c.bytes.data()), std::streamsize(c.bytes.size()));
    }
    return f.good();
}

bool ProjectSerializer::load(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    std::ostringstream ss;
    ss << f.rdbuf();
    const std::string b = ss.str();
    if (b.size() < 20 || std::memcmp(b.data(), "SIMALL-PROJ\n", 12) != 0)
        return false;
    std::uint32_t version = 0, count = 0;
    std::memcpy(&version, b.data() + 12, 4);
    std::memcpy(&count, b.data() + 16, 4);
    if (version != 1)
        return false;
    chunks_.clear();
    std::size_t i = 20;

    auto read_str = [&](std::string& s) -> bool {
        if (i + 4 > b.size())
            return false;
        std::uint32_t n = 0;
        std::memcpy(&n, b.data() + i, 4);
        i += 4;
        if (i + n > b.size())
            return false;
        s.assign(b.data() + i, n);
        i += n;
        return true;
    };
    auto read_u64 = [&](std::uint64_t& v) -> bool {
        if (i + 8 > b.size())
            return false;
        std::memcpy(&v, b.data() + i, 8);
        i += 8;
        return true;
    };

    for (std::uint32_t k = 0; k < count; ++k) {
        std::string name, mime;
        std::uint64_t nb = 0;
        if (!read_str(name) || !read_str(mime) || !read_u64(nb))
            return false;
        if (i + nb > b.size())
            return false;
        ProjectChunk c;
        c.mime = std::move(mime);
        c.bytes.assign(b.data() + i, b.data() + i + nb);
        i += nb;
        chunks_[std::move(name)] = std::move(c);
    }
    return true;
}

} // namespace simall::io
