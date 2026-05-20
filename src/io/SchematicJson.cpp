// =============================================================================
// SimAll Beta — IO Subsystem
// File   : src/io/SchematicJson.cpp
// Phase  : 22 Pass 22.3
// =============================================================================
#include "io/SchematicJson.hpp"

#include "workbench/Cell.hpp"
#include "workbench/CellLink.hpp"
#include "workbench/CellPort.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/Workbench.hpp"

#include <cctype>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace simall::io {

using simall::workbench::Cell;
using simall::workbench::CellId;
using simall::workbench::CellKind;
using simall::workbench::CellLink;
using simall::workbench::CellPort;
using simall::workbench::CellState;
using simall::workbench::PortDirection;
using simall::workbench::PortId;
using simall::workbench::Schematic;
using simall::workbench::StateMachine;

// ---------------------------------------------------------------------------
// Enum parsing helpers (string-view -> enum).  to_string() is already in
// workbench::; we only need the reverse direction here.
// ---------------------------------------------------------------------------
namespace {

CellKind kind_from_str(std::string_view s) {
    if (s == "Geometry") return CellKind::Geometry;
    if (s == "Mesh")     return CellKind::Mesh;
    if (s == "Setup")    return CellKind::Setup;
    if (s == "Solution") return CellKind::Solution;
    if (s == "Results")  return CellKind::Results;
    return CellKind::Custom;
}

CellState state_from_str(std::string_view s) {
    if (s == "UpToDate")        return CellState::UpToDate;
    if (s == "RefreshRequired") return CellState::RefreshRequired;
    if (s == "Failed")          return CellState::Failed;
    return CellState::Unfulfilled;
}

PortDirection dir_from_str(std::string_view s) {
    return (s == "Output") ? PortDirection::Output : PortDirection::Input;
}

// ---------------------------------------------------------------------------
// JSON emit helpers (same style as zones::SelectionManager::to_json).
// ---------------------------------------------------------------------------
std::string escape_json(std::string_view s) {
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

// ---------------------------------------------------------------------------
// Minimalist JSON scanner, sufficient for round-tripping our own output.
// Mirrors the pattern in zones::SelectionManager.
// ---------------------------------------------------------------------------
struct Cursor { const std::string& s; std::size_t i = 0; };

void skip_ws(Cursor& c) {
    while (c.i < c.s.size() && std::isspace(static_cast<unsigned char>(c.s[c.i]))) ++c.i;
}
bool eat(Cursor& c, char ch) {
    skip_ws(c);
    if (c.i < c.s.size() && c.s[c.i] == ch) { ++c.i; return true; }
    return false;
}
bool peek(Cursor& c, char ch) {
    skip_ws(c);
    return c.i < c.s.size() && c.s[c.i] == ch;
}
std::string read_string(Cursor& c) {
    skip_ws(c);
    if (c.i >= c.s.size() || c.s[c.i] != '"') return {};
    ++c.i;
    std::string out;
    while (c.i < c.s.size() && c.s[c.i] != '"') {
        if (c.s[c.i] == '\\' && c.i + 1 < c.s.size()) {
            const char n = c.s[c.i + 1];
            if      (n == 'n') { out += '\n'; c.i += 2; }
            else if (n == 't') { out += '\t'; c.i += 2; }
            else if (n == 'r') { out += '\r'; c.i += 2; }
            else               { out += n;    c.i += 2; }
        } else {
            out += c.s[c.i++];
        }
    }
    if (c.i < c.s.size()) ++c.i;   // closing quote
    return out;
}
long long read_int(Cursor& c) {
    skip_ws(c);
    long long v = 0;
    bool neg = false;
    if (c.i < c.s.size() && (c.s[c.i] == '-' || c.s[c.i] == '+')) {
        neg = (c.s[c.i] == '-'); ++c.i;
    }
    while (c.i < c.s.size() && std::isdigit(static_cast<unsigned char>(c.s[c.i]))) {
        v = v * 10 + (c.s[c.i] - '0');
        ++c.i;
    }
    return neg ? -v : v;
}
bool read_bool(Cursor& c) {
    skip_ws(c);
    if (c.i + 4 <= c.s.size() && c.s.compare(c.i, 4, "true") == 0)  { c.i += 4; return true;  }
    if (c.i + 5 <= c.s.size() && c.s.compare(c.i, 5, "false") == 0) { c.i += 5; return false; }
    return false;
}

}  // namespace


// ---------------------------------------------------------------------------
// schematic_to_json
// ---------------------------------------------------------------------------
std::string schematic_to_json(const Schematic& s, const StateMachine* /*sm*/) {
    std::ostringstream os;
    os << R"({"schema":"simall.workbench.schematic/v1","cells":[)";
    bool first_cell = true;
    for (const Cell& cell : s.cells()) {
        if (!first_cell) os << ',';
        first_cell = false;
        os << R"({"id":)" << cell.id()
           << R"(,"kind":")"  << to_string(cell.kind())
           << R"(","label":")"<< escape_json(cell.label())
           << R"(","adapter":")"<< escape_json(cell.adapter_id())
           << R"(","state":")"<< to_string(cell.state())
           << R"(","ports":[)";
        bool first_port = true;
        for (const CellPort& p : cell.ports()) {
            if (!first_port) os << ',';
            first_port = false;
            os << R"({"id":)" << p.id
               << R"(,"dir":")"  << to_string(p.direction)
               << R"(","name":")"<< escape_json(p.name)
               << R"(","type":")"<< escape_json(p.data_type)
               << R"(","required":)" << (p.required ? "true" : "false") << '}';
        }
        os << "]}";
    }
    os << R"(],"links":[)";
    bool first_link = true;
    for (const CellLink& l : s.links()) {
        if (!first_link) os << ',';
        first_link = false;
        os << R"({"from_cell":)" << l.from_cell
           << R"(,"from_port":)" << l.from_port
           << R"(,"to_cell":)"   << l.to_cell
           << R"(,"to_port":)"   << l.to_port << '}';
    }
    os << "]}";
    return os.str();
}


// ---------------------------------------------------------------------------
// schematic_from_json
// ---------------------------------------------------------------------------
namespace {

// One row of the cells[] array, with per-cell port id remap embedded so we
// can rewire links to the freshly-minted ids that Schematic / Cell will
// hand out.
struct ParsedPort {
    PortId        old_id{};
    PortDirection direction{PortDirection::Input};
    std::string   name;
    std::string   data_type;
    bool          required{true};
};
struct ParsedCell {
    CellId      old_id{};
    CellKind    kind{CellKind::Custom};
    std::string label;
    std::string adapter_id;
    CellState   state{CellState::Unfulfilled};
    std::vector<ParsedPort> ports;
};
struct ParsedLink {
    CellId from_cell{}, to_cell{};
    PortId from_port{}, to_port{};
};

// Skip a JSON value of arbitrary shape (used for unknown keys).
void skip_value(Cursor& c) {
    skip_ws(c);
    if (c.i >= c.s.size()) return;
    const char ch = c.s[c.i];
    if (ch == '"')                   { (void)read_string(c); return; }
    if (ch == 't' || ch == 'f')      { (void)read_bool(c);   return; }
    if (ch == '-' || std::isdigit(static_cast<unsigned char>(ch))) {
        (void)read_int(c); return;
    }
    if (ch == '{' || ch == '[') {
        const char open  = ch;
        const char close = (open == '{') ? '}' : ']';
        int depth = 0;
        while (c.i < c.s.size()) {
            const char d = c.s[c.i++];
            if      (d == open)  ++depth;
            else if (d == close) { --depth; if (depth == 0) return; }
            else if (d == '"') {
                // skip until matching unescaped quote
                while (c.i < c.s.size() && c.s[c.i] != '"') {
                    if (c.s[c.i] == '\\' && c.i + 1 < c.s.size()) c.i += 2;
                    else                                            ++c.i;
                }
                if (c.i < c.s.size()) ++c.i;
            }
        }
    }
}

bool parse_port(Cursor& c, ParsedPort& p) {
    if (!eat(c, '{')) return false;
    while (true) {
        skip_ws(c);
        if (peek(c, '}')) { ++c.i; return true; }
        const std::string key = read_string(c);
        if (key.empty() || !eat(c, ':')) return false;
        if      (key == "id")       p.old_id    = static_cast<PortId>(read_int(c));
        else if (key == "dir")      p.direction = dir_from_str(read_string(c));
        else if (key == "name")     p.name      = read_string(c);
        else if (key == "type")     p.data_type = read_string(c);
        else if (key == "required") p.required  = read_bool(c);
        else                        skip_value(c);
        skip_ws(c);
        if (peek(c, ',')) ++c.i;
    }
}

bool parse_cell(Cursor& c, ParsedCell& pc) {
    if (!eat(c, '{')) return false;
    while (true) {
        skip_ws(c);
        if (peek(c, '}')) { ++c.i; return true; }
        const std::string key = read_string(c);
        if (key.empty() || !eat(c, ':')) return false;
        if      (key == "id")    pc.old_id = static_cast<CellId>(read_int(c));
        else if (key == "kind")  pc.kind   = kind_from_str(read_string(c));
        else if (key == "label") pc.label  = read_string(c);
        else if (key == "adapter") pc.adapter_id = read_string(c);
        else if (key == "state") pc.state  = state_from_str(read_string(c));
        else if (key == "ports") {
            if (!eat(c, '[')) return false;
            while (!peek(c, ']')) {
                ParsedPort p;
                if (!parse_port(c, p)) return false;
                pc.ports.push_back(std::move(p));
                skip_ws(c);
                if (peek(c, ',')) ++c.i;
            }
            ++c.i;  // consume ']'
        } else {
            skip_value(c);
        }
        skip_ws(c);
        if (peek(c, ',')) ++c.i;
    }
}

bool parse_link(Cursor& c, ParsedLink& pl) {
    if (!eat(c, '{')) return false;
    while (true) {
        skip_ws(c);
        if (peek(c, '}')) { ++c.i; return true; }
        const std::string key = read_string(c);
        if (key.empty() || !eat(c, ':')) return false;
        if      (key == "from_cell") pl.from_cell = static_cast<CellId>(read_int(c));
        else if (key == "from_port") pl.from_port = static_cast<PortId>(read_int(c));
        else if (key == "to_cell")   pl.to_cell   = static_cast<CellId>(read_int(c));
        else if (key == "to_port")   pl.to_port   = static_cast<PortId>(read_int(c));
        else                         skip_value(c);
        skip_ws(c);
        if (peek(c, ',')) ++c.i;
    }
}

}  // namespace

bool schematic_from_json(const std::string& json, Schematic& out, StateMachine* sm) {
    // Reset target. Schematic has no `clear()`; assign a default-constructed
    // instance instead.  StateMachine holds a Schematic& and re-derives on
    // demand, so it does not need to be re-bound.
    out = Schematic{};

    Cursor c{json};
    if (!eat(c, '{')) return false;

    std::vector<ParsedCell> cells;
    std::vector<ParsedLink> links;

    while (true) {
        skip_ws(c);
        if (peek(c, '}')) { ++c.i; break; }
        const std::string key = read_string(c);
        if (key.empty() || !eat(c, ':')) return false;
        if (key == "cells") {
            if (!eat(c, '[')) return false;
            while (!peek(c, ']')) {
                ParsedCell pc;
                if (!parse_cell(c, pc)) return false;
                cells.push_back(std::move(pc));
                skip_ws(c);
                if (peek(c, ',')) ++c.i;
            }
            ++c.i;
        } else if (key == "links") {
            if (!eat(c, '[')) return false;
            while (!peek(c, ']')) {
                ParsedLink pl;
                if (!parse_link(c, pl)) return false;
                links.push_back(pl);
                skip_ws(c);
                if (peek(c, ',')) ++c.i;
            }
            ++c.i;
        } else {
            skip_value(c);
        }
        skip_ws(c);
        if (peek(c, ',')) ++c.i;
    }

    // Build remap tables as we recreate cells / ports through the public
    // mutators so Schematic / Cell mint fresh ids.
    std::unordered_map<CellId, CellId>                             cell_remap;
    std::unordered_map<CellId, std::unordered_map<PortId, PortId>> port_remap;

    for (const ParsedCell& pc : cells) {
        const CellId new_id = out.add_cell(pc.kind, pc.label);
        cell_remap[pc.old_id] = new_id;
        Cell* c2 = out.cell(new_id);
        if (!c2) return false;
        auto& pm = port_remap[pc.old_id];
        for (const ParsedPort& pp : pc.ports) {
            PortId minted = simall::workbench::kInvalidPortId;
            if (pp.direction == PortDirection::Output) {
                minted = c2->add_output(pp.name, pp.data_type);
            } else {
                minted = c2->add_input(pp.name, pp.data_type, pp.required);
            }
            pm[pp.old_id] = minted;
        }
        // Apply persisted state directly; StateMachine, if provided, will
        // honor it on the next recompute pass.
        c2->set_state(pc.state);
        c2->set_adapter_id(pc.adapter_id);
    }

    for (const ParsedLink& pl : links) {
        auto fc = cell_remap.find(pl.from_cell);
        auto tc = cell_remap.find(pl.to_cell);
        if (fc == cell_remap.end() || tc == cell_remap.end()) return false;
        auto fpit = port_remap[pl.from_cell].find(pl.from_port);
        auto tpit = port_remap[pl.to_cell].find(pl.to_port);
        if (fpit == port_remap[pl.from_cell].end() ||
            tpit == port_remap[pl.to_cell].end()) return false;
        const CellLink link{fc->second, fpit->second, tc->second, tpit->second};
        if (!out.add_link(link)) return false;
    }

    // If the caller provided a StateMachine, leave the persisted per-cell
    // states alone -- the user has already set them via Cell::set_state().
    // Touch the pointer to avoid unused-parameter warnings.
    (void)sm;
    return true;
}


// ---------------------------------------------------------------------------
// File round-trip
// ---------------------------------------------------------------------------
bool save_schematic_json(const std::string& path, const Schematic& s, const StateMachine* sm) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    const std::string text = schematic_to_json(s, sm);
    f.write(text.data(), static_cast<std::streamsize>(text.size()));
    return f.good();
}

bool load_schematic_json(const std::string& path, Schematic& out, StateMachine* sm) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::ostringstream ss; ss << f.rdbuf();
    return schematic_from_json(ss.str(), out, sm);
}

}  // namespace simall::io
