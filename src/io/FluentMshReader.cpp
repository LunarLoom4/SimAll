// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/FluentMshReader.cpp
//
// Recursive-descent parser for Fluent ASCII `.msh`.  The grammar is a
// nested s-expression: `(index data)` where `data` may itself be a list.
// We tokenise hex numbers ([0-9a-fA-F]+) and floats while preserving
// parentheses so the section handlers see one balanced group at a time.
// =============================================================================
#include "io/FluentMshReader.hpp"

#include <cctype>
#include <charconv>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace simall::io {

namespace {

struct Lexer {
    const std::string& src;
    std::size_t        i = 0;

    bool eof() const { return i >= src.size(); }
    void skip_ws() {
        while (i < src.size()) {
            const char c = src[i];
            if (std::isspace(static_cast<unsigned char>(c))) { ++i; continue; }
            // Fluent files can carry `\n`-terminated `! comment` lines.
            if (c == '!') {
                while (i < src.size() && src[i] != '\n') ++i;
                continue;
            }
            break;
        }
    }
    char peek() { skip_ws(); return eof() ? 0 : src[i]; }
    char consume() { skip_ws(); return eof() ? 0 : src[i++]; }

    // Read one token (atom).  Stops at whitespace or paren.
    std::string atom() {
        skip_ws();
        std::string out;
        while (i < src.size()) {
            const char c = src[i];
            if (c == '(' || c == ')' || std::isspace(static_cast<unsigned char>(c))) break;
            out.push_back(c);
            ++i;
        }
        return out;
    }

    bool expect(char c) {
        if (consume() != c) return false;
        return true;
    }
};

std::uint64_t hex_u64(std::string_view s) {
    std::uint64_t v = 0;
    for (char c : s) {
        v <<= 4;
        if      (c >= '0' && c <= '9') v |= std::uint64_t(c - '0');
        else if (c >= 'a' && c <= 'f') v |= std::uint64_t(c - 'a' + 10);
        else if (c >= 'A' && c <= 'F') v |= std::uint64_t(c - 'A' + 10);
        else return std::uint64_t(-1);
    }
    return v;
}

// Read a balanced parenthesised body (excluding the outer parens) into `out`.
bool read_balanced(Lexer& lx, std::string& out) {
    if (lx.consume() != '(') return false;
    int depth = 1;
    while (!lx.eof()) {
        const char c = lx.src[lx.i];
        if (c == '(') ++depth;
        else if (c == ')') {
            if (--depth == 0) { ++lx.i; return true; }
        }
        out.push_back(c);
        ++lx.i;
    }
    return false;
}

// Split a flat numeric body into doubles.
std::vector<double> split_doubles(std::string_view body) {
    std::vector<double> out;
    std::string token;
    auto flush = [&]() {
        if (token.empty()) return;
        try { out.push_back(std::stod(token)); } catch (...) {}
        token.clear();
    };
    for (char c : body) {
        if (std::isspace(static_cast<unsigned char>(c))) flush();
        else token.push_back(c);
    }
    flush();
    return out;
}

ElementType fluent_cell_type(int t) {
    // Fluent cell-type indices (Table 6.2 of the Fluent file format spec).
    switch (t) {
        case 1: return ElementType::Tri3;       // triangular (2-D)
        case 2: return ElementType::Tetra4;
        case 3: return ElementType::Quad4;      // quadrilateral (2-D)
        case 4: return ElementType::Hexa8;
        case 5: return ElementType::Pyra5;
        case 6: return ElementType::Penta6;     // wedge
        case 7: return ElementType::Poly;       // polyhedral
        default: return ElementType::Unknown;
    }
}

const char* fluent_bc_type(int code) {
    switch (code) {
        case  2: return "interior";
        case  3: return "wall";
        case  4: return "pressure-inlet";
        case  5: return "pressure-outlet";
        case  7: return "symmetry";
        case  8: return "periodic-shadow";
        case  9: return "pressure-far-field";
        case 10: return "velocity-inlet";
        case 12: return "periodic";
        case 14: return "fan";
        case 20: return "mass-flow-inlet";
        case 24: return "interface";
        case 31: return "porous-jump";
        case 36: return "outflow";
        default: return "unspecified";
    }
}

void parse_section(int index, std::string_view header, std::string_view body,
                   FluentReadResult& res) {
    auto& mesh = res.mesh;
    if (mesh.zones.empty()) mesh.zones.emplace_back();
    auto& zone = mesh.zones.front();

    if (index == 2) {
        // (2 dim)
        const auto trimmed = std::string(header);
        try { res.dimension = static_cast<std::uint8_t>(std::stoi(trimmed)); } catch (...) {}
        return;
    }
    if (index == 10) {
        // Header: zone-id first last type ND   (all hex)
        std::istringstream h{std::string(header)};
        std::string zid, first, last, type, nd;
        h >> zid >> first >> last >> type >> nd;
        const auto firstId = hex_u64(first);
        const auto lastId  = hex_u64(last);
        const auto zoneId  = hex_u64(zid);
        if (zoneId == 0) {
            // Declaration (no body), record final node count.
            return;
        }
        const auto count = (lastId >= firstId) ? (lastId - firstId + 1) : 0;
        zone.x.reserve(zone.x.size() + count);
        zone.y.reserve(zone.y.size() + count);
        zone.z.reserve(zone.z.size() + count);
        const auto doubles = split_doubles(body);
        const int ND = (nd.empty() ? 3 : int(hex_u64(nd)));
        for (std::size_t k = 0; k + (ND - 1) < doubles.size(); k += std::size_t(ND)) {
            zone.x.push_back(doubles[k]);
            zone.y.push_back(doubles[k + 1]);
            zone.z.push_back(ND == 3 ? doubles[k + 2] : 0.0);
        }
        return;
    }
    if (index == 12) {
        // (12 (zone-id first last type elem-type)) — when zone-id==0 this is
        // a declaration with no body.  Body (if present) lists per-element
        // type codes for mixed-element zones; we currently only honour the
        // declared elem-type and emit a contiguous ElementSection.
        std::istringstream h{std::string(header)};
        std::string zid, first, last, type, et;
        h >> zid >> first >> last >> type >> et;
        const auto zoneId = hex_u64(zid);
        if (zoneId == 0) return;
        const auto firstId = hex_u64(first);
        const auto lastId  = hex_u64(last);
        const auto count   = (lastId >= firstId) ? (lastId - firstId + 1) : 0;
        ElementType et_enum = fluent_cell_type(int(hex_u64(et)));
        ElementSection sec;
        sec.name = "cells_" + std::to_string(zoneId);
        sec.type = et_enum;
        // Connectivity is in the (13 ...) face sections; for now we record the
        // expected element count via empty `nodes` and let the caller build
        // cell-node connectivity from face-owner pairs.
        sec.nodes.resize(count * vertices_per_element(et_enum), 0);
        zone.sections.push_back(std::move(sec));
        return;
    }
    if (index == 13) {
        // (13 (zone-id first last bc-type face-type)(...))
        std::istringstream h{std::string(header)};
        std::string zid, first, last, bc, ft;
        h >> zid >> first >> last >> bc >> ft;
        const auto zoneId = hex_u64(zid);
        if (zoneId == 0) return;
        const auto firstId = hex_u64(first);
        const auto lastId  = hex_u64(last);
        const auto count   = (lastId >= firstId) ? (lastId - firstId + 1) : 0;
        const int faceType = int(hex_u64(ft));   // 0=mixed,2=lin,3=tri,4=quad
        const int bcCode   = int(hex_u64(bc));

        ElementType etype = ElementType::Unknown;
        switch (faceType) {
            case 2: etype = ElementType::Bar2;  break;
            case 3: etype = ElementType::Tri3;  break;
            case 4: etype = ElementType::Quad4; break;
            default: etype = ElementType::Poly; break;
        }
        ElementSection sec;
        sec.name = "faces_" + std::to_string(zoneId);
        sec.type = etype;

        // Parse face connectivity: each face has N node ids + owner-cell + neighbour-cell.
        std::istringstream b{std::string(body)};
        std::string tok;
        std::vector<std::uint64_t> all;
        while (b >> tok) all.push_back(hex_u64(tok));
        const std::uint8_t v = (etype == ElementType::Poly) ? 0 : vertices_per_element(etype);
        if (v != 0 && all.size() >= count * (v + 2)) {
            sec.nodes.reserve(count * v);
            for (std::uint64_t f = 0; f < count; ++f) {
                const auto base = f * (v + 2);
                for (std::uint8_t k = 0; k < v; ++k) {
                    sec.nodes.push_back(NodeIdx(all[base + k] - 1));
                }
            }
        } else if (etype == ElementType::Poly) {
            // n-gon format: first entry is the vertex count
            std::size_t i = 0;
            std::uint32_t off = 0;
            sec.polyOffsets.push_back(0);
            while (i < all.size()) {
                const std::uint64_t n = all[i++];
                if (i + n + 2 > all.size()) break;
                for (std::uint64_t k = 0; k < n; ++k) {
                    sec.polyNodes.push_back(NodeIdx(all[i + k] - 1));
                }
                i  += n + 2;          // skip owner / neighbour
                off += std::uint32_t(n);
                sec.polyOffsets.push_back(off);
            }
        }
        zone.sections.push_back(std::move(sec));

        BoundaryPatch bp;
        bp.name   = "boundary_" + std::to_string(zoneId);
        bp.bcType = fluent_bc_type(bcCode);
        bp.faceElementIndices.reserve(count);
        for (std::uint64_t f = 0; f < count; ++f) bp.faceElementIndices.push_back(std::uint32_t(f));
        zone.boundaries.push_back(std::move(bp));
        return;
    }
    if (index == 39 || index == 45) {
        // Zone label : (39 (zone-id zone-type zone-name) ())
        std::istringstream h{std::string(header)};
        std::string zid, ztype, zname;
        h >> zid >> ztype >> zname;
        for (auto& b : zone.boundaries) {
            if (b.name == "boundary_" + std::to_string(hex_u64(zid))) {
                b.name = zname.empty() ? b.name : zname;
                if (!ztype.empty()) b.bcType = ztype;
            }
        }
        return;
    }
    // 2010 / 3010 / 2012 / etc → binary sections, not supported in ASCII reader.
    if (index >= 2000) {
        res.error = "Binary Fluent .msh sections are not supported";
    }
}

}  // namespace

FluentReadResult parse_fluent_msh_string(const std::string& text,
                                          std::string sourceHint) {
    FluentReadResult res;
    res.mesh.sourceFormat = "fluent_msh";
    res.mesh.sourcePath   = std::move(sourceHint);
    res.mesh.zones.emplace_back();
    res.mesh.zones.front().name = "default";

    Lexer lx{text};
    while (!lx.eof()) {
        lx.skip_ws();
        if (lx.eof()) break;
        if (lx.peek() != '(') { ++lx.i; continue; }
        ++lx.i;                                     // consume '('
        // Read section index as the first token.
        std::string idxTok = lx.atom();
        int index = 0;
        try { index = std::stoi(idxTok); } catch (...) {}

        // The remainder up to the matching ')' is the section payload.
        // It may itself contain nested parens, so we collect a balanced body.
        std::string payload;
        int depth = 1;
        while (lx.i < text.size()) {
            const char c = text[lx.i];
            if (c == '(') ++depth;
            else if (c == ')') { if (--depth == 0) { ++lx.i; break; } }
            payload.push_back(c);
            ++lx.i;
        }
        // Split payload into header (first paren group) + body (second).
        std::string header, body;
        std::size_t p = 0;
        while (p < payload.size() && std::isspace((unsigned char)payload[p])) ++p;
        if (p < payload.size() && payload[p] == '(') {
            int d = 1; ++p;
            while (p < payload.size() && d > 0) {
                const char c = payload[p];
                if (c == '(') ++d;
                else if (c == ')') { --d; if (d == 0) { ++p; break; } }
                header.push_back(c);
                ++p;
            }
            while (p < payload.size() && std::isspace((unsigned char)payload[p])) ++p;
            if (p < payload.size() && payload[p] == '(') {
                int d2 = 1; ++p;
                while (p < payload.size() && d2 > 0) {
                    const char c = payload[p];
                    if (c == '(') ++d2;
                    else if (c == ')') { --d2; if (d2 == 0) { ++p; break; } }
                    body.push_back(c);
                    ++p;
                }
            }
        } else {
            header = std::move(payload);
        }
        parse_section(index, header, body, res);
        if (!res.error.empty()) return res;
    }
    res.ok = true;
    return res;
}

FluentReadResult read_fluent_msh(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        FluentReadResult r;
        r.error = "Cannot open: " + path;
        return r;
    }
    std::ostringstream ss;
    ss << f.rdbuf();
    return parse_fluent_msh_string(ss.str(), path);
}

}  // namespace simall::io
