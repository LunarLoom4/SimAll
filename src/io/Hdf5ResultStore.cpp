// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/Hdf5ResultStore.cpp
// =============================================================================
#include "io/Hdf5ResultStore.hpp"

#include <cstring>
#include <fstream>
#include <sstream>

namespace simall::io {

// ---------------------------------------------------------------------------
// Dataset helpers
// ---------------------------------------------------------------------------
std::size_t Dataset::element_size() const noexcept {
    switch (dtype) {
        case StoreDtype::F64: return 8;
        case StoreDtype::F32: return 4;
        case StoreDtype::I32: return 4;
        case StoreDtype::I64: return 8;
        case StoreDtype::U8:  return 1;
    }
    return 0;
}
std::uint64_t Dataset::element_count() const noexcept {
    std::uint64_t n = 1;
    for (auto s : shape) n *= s;
    return n;
}

// ---------------------------------------------------------------------------
// Path resolution (HDF5-style "/a/b/c")
// ---------------------------------------------------------------------------
namespace {
std::vector<std::string> split_path(const std::string& path) {
    std::vector<std::string> out;
    std::string              cur;
    for (char c : path) {
        if (c == '/') { if (!cur.empty()) { out.push_back(std::move(cur)); cur.clear(); } }
        else cur.push_back(c);
    }
    if (!cur.empty()) out.push_back(std::move(cur));
    return out;
}
}  // namespace

Group* Group::create_group(const std::string& name) {
    auto& slot = children[name];
    if (!slot) slot = std::make_unique<Group>();
    return slot.get();
}

Group* Group::find_group(const std::string& path) {
    auto parts = split_path(path);
    Group* g = this;
    for (auto& p : parts) g = g->create_group(p);
    return g;
}

Dataset* Group::find_dataset(const std::string& path) {
    auto parts = split_path(path);
    if (parts.empty()) return nullptr;
    auto leaf = parts.back();
    parts.pop_back();
    Group* g = this;
    for (auto& p : parts) g = g->create_group(p);
    auto it = g->datasets.find(leaf);
    return it == g->datasets.end() ? nullptr : &it->second;
}

void Group::set_dataset(const std::string& path, Dataset ds) {
    auto parts = split_path(path);
    if (parts.empty()) return;
    auto leaf = parts.back();
    parts.pop_back();
    Group* g = this;
    for (auto& p : parts) g = g->create_group(p);
    g->datasets[leaf] = std::move(ds);
}

void Group::set_attribute(const std::string& path, const std::string& key,
                           const std::string& value) {
    Group* g = find_group(path);
    g->attributes[key] = value;
}

// ---------------------------------------------------------------------------
// Store
// ---------------------------------------------------------------------------
Hdf5ResultStore::Hdf5ResultStore() : root_(std::make_unique<Group>()) {}
Hdf5ResultStore::~Hdf5ResultStore() = default;

bool Hdf5ResultStore::has_hdf5_backend() noexcept {
#ifdef SIMALL_HAVE_HDF5
    return true;
#else
    return false;
#endif
}

void Hdf5ResultStore::write_field(const std::string& groupPath,
                                   const std::string& name,
                                   const std::vector<double>& values,
                                   FieldLocation loc) {
    Group* g = root_->find_group(groupPath);
    Dataset ds;
    ds.dtype = StoreDtype::F64;
    ds.shape = { std::uint64_t(values.size()) };
    ds.bytes.resize(values.size() * sizeof(double));
    if (!values.empty())
        std::memcpy(ds.bytes.data(), values.data(), values.size() * sizeof(double));
    g->datasets[name] = std::move(ds);
    const char* locStr = (loc == FieldLocation::Cell) ? "cell"
                        : (loc == FieldLocation::Node) ? "node" : "face";
    g->attributes[name + ".location"] = locStr;
}

std::vector<double> Hdf5ResultStore::read_field_f64(
        const std::string& groupPath, const std::string& name) const {
    auto* g = root_.get();
    for (auto& p : split_path(groupPath)) {
        auto it = g->children.find(p);
        if (it == g->children.end()) return {};
        g = it->second.get();
    }
    auto it = g->datasets.find(name);
    if (it == g->datasets.end()) return {};
    const auto& ds = it->second;
    if (ds.dtype != StoreDtype::F64) return {};
    std::vector<double> out(ds.element_count());
    if (!out.empty())
        std::memcpy(out.data(), ds.bytes.data(), out.size() * sizeof(double));
    return out;
}

// ---------------------------------------------------------------------------
// Native chunked binary IO
// ---------------------------------------------------------------------------
namespace {

void put_u8 (std::ostream& o, std::uint8_t v)  { o.write(reinterpret_cast<const char*>(&v), 1); }
void put_u32(std::ostream& o, std::uint32_t v) { o.write(reinterpret_cast<const char*>(&v), 4); }
void put_u64(std::ostream& o, std::uint64_t v) { o.write(reinterpret_cast<const char*>(&v), 8); }
void put_str(std::ostream& o, const std::string& s) {
    put_u32(o, std::uint32_t(s.size()));
    o.write(s.data(), std::streamsize(s.size()));
}

void serialize_group(std::ostream& o, const Group& g) {
    // attributes
    put_u32(o, std::uint32_t(g.attributes.size()));
    for (auto& [k, v] : g.attributes) { put_str(o, k); put_str(o, v); }
    // datasets
    put_u32(o, std::uint32_t(g.datasets.size()));
    for (auto& [name, ds] : g.datasets) {
        put_str(o, name);
        put_u8 (o, static_cast<std::uint8_t>(ds.dtype));
        put_u32(o, std::uint32_t(ds.shape.size()));
        for (auto s : ds.shape) put_u64(o, s);
        put_u64(o, ds.bytes.size());
        o.write(reinterpret_cast<const char*>(ds.bytes.data()),
                std::streamsize(ds.bytes.size()));
    }
    // children
    put_u32(o, std::uint32_t(g.children.size()));
    for (auto& [name, child] : g.children) {
        put_str(o, name);
        serialize_group(o, *child);
    }
}

struct Cur {
    const std::string& b;
    std::size_t        i = 0;
    bool need(std::size_t n) const { return i + n <= b.size(); }
    template<typename T> bool read(T& v) {
        if (!need(sizeof(T))) return false;
        std::memcpy(&v, b.data() + i, sizeof(T));
        i += sizeof(T);
        return true;
    }
    bool read_bytes(void* dst, std::size_t n) {
        if (!need(n)) return false;
        std::memcpy(dst, b.data() + i, n);
        i += n;
        return true;
    }
    bool read_str(std::string& s) {
        std::uint32_t n = 0;
        if (!read(n)) return false;
        if (!need(n)) return false;
        s.assign(b.data() + i, n);
        i += n;
        return true;
    }
};

bool deserialize_group(Cur& c, Group& g) {
    std::uint32_t na = 0;
    if (!c.read(na)) return false;
    for (std::uint32_t i = 0; i < na; ++i) {
        std::string k, v;
        if (!c.read_str(k) || !c.read_str(v)) return false;
        g.attributes[k] = v;
    }
    std::uint32_t nd = 0;
    if (!c.read(nd)) return false;
    for (std::uint32_t i = 0; i < nd; ++i) {
        std::string name;
        std::uint8_t dt = 0;
        std::uint32_t rank = 0;
        if (!c.read_str(name) || !c.read(dt) || !c.read(rank)) return false;
        Dataset ds;
        ds.dtype = static_cast<StoreDtype>(dt);
        ds.shape.resize(rank);
        for (std::uint32_t k = 0; k < rank; ++k) if (!c.read(ds.shape[k])) return false;
        std::uint64_t nb = 0;
        if (!c.read(nb)) return false;
        ds.bytes.resize(nb);
        if (nb && !c.read_bytes(ds.bytes.data(), std::size_t(nb))) return false;
        g.datasets[name] = std::move(ds);
    }
    std::uint32_t nc = 0;
    if (!c.read(nc)) return false;
    for (std::uint32_t i = 0; i < nc; ++i) {
        std::string name;
        if (!c.read_str(name)) return false;
        auto& slot = g.children[name];
        slot = std::make_unique<Group>();
        if (!deserialize_group(c, *slot)) return false;
    }
    return true;
}

}  // namespace

bool Hdf5ResultStore::save_native(const std::string& path) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    f.write("SIMALL-SRS\n", 11);
    put_u32(f, 1);
    serialize_group(f, *root_);
    return f.good();
}

bool Hdf5ResultStore::load_native(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::ostringstream ss;
    ss << f.rdbuf();
    const std::string bytes = ss.str();
    if (bytes.size() < 15 || std::memcmp(bytes.data(), "SIMALL-SRS\n", 11) != 0) return false;
    std::uint32_t version = 0;
    std::memcpy(&version, bytes.data() + 11, 4);
    if (version != 1) return false;
    Cur c{bytes};
    c.i = 15;
    root_ = std::make_unique<Group>();
    return deserialize_group(c, *root_);
}

}  // namespace simall::io
