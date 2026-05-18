// =============================================================================
// SimAll Beta - Zones / Selection Subsystem
// File   : src/zones/SelectionManager.cpp
// Phase  : 25
// =============================================================================
#include "zones/SelectionManager.hpp"

#include <algorithm>
#include <set>
#include <sstream>

namespace simall::zones {

namespace {

const char* entity_to_str(SelectionEntity e) {
    switch (e) {
        case SelectionEntity::Face:   return "face";
        case SelectionEntity::Cell:   return "cell";
        case SelectionEntity::Edge:   return "edge";
        case SelectionEntity::Vertex: return "vertex";
    }
    return "face";
}
SelectionEntity entity_from_str(const std::string& s) {
    if (s == "cell")   return SelectionEntity::Cell;
    if (s == "edge")   return SelectionEntity::Edge;
    if (s == "vertex") return SelectionEntity::Vertex;
    return SelectionEntity::Face;
}

std::string escape_json(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 4);
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:   out += c;
        }
    }
    return out;
}

// Minimalist JSON object scanner sufficient for round-tripping our own
// `to_json()` output. Not a general JSON parser; that lives in
// utilities/JsonIo if a richer parser is ever needed.
struct Cursor { const std::string& s; std::size_t i = 0; };
void skip_ws(Cursor& c) { while (c.i < c.s.size() && std::isspace(static_cast<unsigned char>(c.s[c.i]))) ++c.i; }
bool eat(Cursor& c, char ch) { skip_ws(c); if (c.i < c.s.size() && c.s[c.i] == ch) { ++c.i; return true; } return false; }
std::string read_string(Cursor& c) {
    skip_ws(c);
    if (c.i >= c.s.size() || c.s[c.i] != '"') return {};
    ++c.i;
    std::string out;
    while (c.i < c.s.size() && c.s[c.i] != '"') {
        if (c.s[c.i] == '\\' && c.i + 1 < c.s.size()) {
            char n = c.s[c.i + 1];
            if      (n == 'n')  { out += '\n'; c.i += 2; }
            else if (n == 't')  { out += '\t'; c.i += 2; }
            else if (n == 'r')  { out += '\r'; c.i += 2; }
            else                { out += n;    c.i += 2; }
        } else {
            out += c.s[c.i++];
        }
    }
    if (c.i < c.s.size()) ++c.i;  // closing quote
    return out;
}
long long read_int(Cursor& c) {
    skip_ws(c);
    long long v = 0;
    bool neg = false;
    if (c.i < c.s.size() && (c.s[c.i] == '-' || c.s[c.i] == '+')) { neg = (c.s[c.i] == '-'); ++c.i; }
    while (c.i < c.s.size() && std::isdigit(static_cast<unsigned char>(c.s[c.i]))) {
        v = v * 10 + (c.s[c.i] - '0');
        ++c.i;
    }
    return neg ? -v : v;
}

}  // namespace

bool SelectionManager::define(const NamedSelection& s) {
    if (s.name.empty()) return false;
    const bool overwrite = store_.find(s.name) != store_.end();
    store_[s.name] = s;
    return overwrite;
}

bool SelectionManager::remove(const std::string& name) {
    return store_.erase(name) > 0;
}

const NamedSelection* SelectionManager::find(const std::string& name) const noexcept {
    auto it = store_.find(name);
    return it == store_.end() ? nullptr : &it->second;
}
NamedSelection* SelectionManager::find(const std::string& name) noexcept {
    auto it = store_.find(name);
    return it == store_.end() ? nullptr : &it->second;
}

bool SelectionManager::rename(const std::string& oldName, const std::string& newName) {
    if (newName.empty() || store_.find(newName) != store_.end()) return false;
    auto it = store_.find(oldName);
    if (it == store_.end()) return false;
    NamedSelection sel = std::move(it->second);
    store_.erase(it);
    sel.name = newName;
    store_.emplace(newName, std::move(sel));
    return true;
}

namespace {
std::optional<NamedSelection> combine(const SelectionManager& mgr,
                                       const std::string& a, const std::string& b,
                                       const std::string& result,
                                       int op /*0=union,1=intersect,2=subtract*/) {
    const NamedSelection* sa = mgr.find(a);
    const NamedSelection* sb = mgr.find(b);
    if (!sa || !sb || sa->entity != sb->entity) return std::nullopt;

    std::set<util::PersistentId> A(sa->ids.begin(), sa->ids.end());
    std::set<util::PersistentId> B(sb->ids.begin(), sb->ids.end());
    NamedSelection out;
    out.name   = result;
    out.entity = sa->entity;
    if (op == 0) {
        out.ids.assign(A.begin(), A.end());
        for (auto id : B) if (!A.count(id)) out.ids.push_back(id);
    } else if (op == 1) {
        for (auto id : A) if (B.count(id)) out.ids.push_back(id);
    } else {
        for (auto id : A) if (!B.count(id)) out.ids.push_back(id);
    }
    return out;
}
}  // namespace

std::optional<NamedSelection>
SelectionManager::set_union(const std::string& a, const std::string& b,
                            const std::string& result) const { return combine(*this, a, b, result, 0); }
std::optional<NamedSelection>
SelectionManager::set_intersect(const std::string& a, const std::string& b,
                                const std::string& result) const { return combine(*this, a, b, result, 1); }
std::optional<NamedSelection>
SelectionManager::set_subtract(const std::string& a, const std::string& b,
                               const std::string& result) const { return combine(*this, a, b, result, 2); }

std::vector<std::string> SelectionManager::names() const {
    std::vector<std::string> out;
    out.reserve(store_.size());
    for (const auto& kv : store_) out.push_back(kv.first);
    std::sort(out.begin(), out.end());
    return out;
}

std::string SelectionManager::to_json() const {
    std::ostringstream os;
    os << "{\"selections\":[";
    bool first = true;
    auto sortedNames = names();
    for (const auto& n : sortedNames) {
        const auto& s = store_.at(n);
        if (!first) os << ',';
        first = false;
        os << "{\"name\":\""    << escape_json(s.name)
           << "\",\"entity\":\"" << entity_to_str(s.entity)
           << "\",\"color\":\""  << escape_json(s.color)
           << "\",\"comment\":\"" << escape_json(s.comment)
           << "\",\"ids\":[";
        for (std::size_t i = 0; i < s.ids.size(); ++i) {
            if (i) os << ',';
            os << s.ids[i];
        }
        os << "]}";
    }
    os << "]}";
    return os.str();
}

bool SelectionManager::from_json(const std::string& json) {
    store_.clear();
    Cursor c{json};
    if (!eat(c, '{')) return false;
    // Expect "selections": [ ... ]
    skip_ws(c);
    std::string key = read_string(c);
    if (key != "selections") return false;
    if (!eat(c, ':')) return false;
    if (!eat(c, '[')) return false;

    while (true) {
        skip_ws(c);
        if (c.i < c.s.size() && c.s[c.i] == ']') { ++c.i; break; }
        if (!eat(c, '{')) return false;
        NamedSelection sel;
        while (true) {
            skip_ws(c);
            if (c.i < c.s.size() && c.s[c.i] == '}') { ++c.i; break; }
            const std::string k = read_string(c);
            if (!eat(c, ':')) return false;
            if (k == "name")        sel.name    = read_string(c);
            else if (k == "color")  sel.color   = read_string(c);
            else if (k == "comment") sel.comment = read_string(c);
            else if (k == "entity") sel.entity  = entity_from_str(read_string(c));
            else if (k == "ids") {
                if (!eat(c, '[')) return false;
                while (true) {
                    skip_ws(c);
                    if (c.i < c.s.size() && c.s[c.i] == ']') { ++c.i; break; }
                    sel.ids.push_back(static_cast<util::PersistentId>(read_int(c)));
                    skip_ws(c);
                    if (c.i < c.s.size() && c.s[c.i] == ',') ++c.i;
                }
            } else {
                // skip unknown value (string or number) until comma/brace.
                skip_ws(c);
                if (c.i < c.s.size() && c.s[c.i] == '"') (void)read_string(c);
                else                                      (void)read_int(c);
            }
            skip_ws(c);
            if (c.i < c.s.size() && c.s[c.i] == ',') ++c.i;
        }
        if (!sel.name.empty()) store_[sel.name] = std::move(sel);
        skip_ws(c);
        if (c.i < c.s.size() && c.s[c.i] == ',') ++c.i;
    }
    return true;
}

}  // namespace simall::zones
