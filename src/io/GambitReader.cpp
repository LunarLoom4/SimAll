// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/GambitReader.cpp
// Phase  : 23 Pass 7
//
// Implementation of the GAMBIT .neu reader.  See the header for the
// grammar reference and the supported element / BC types.
// =============================================================================
#include "io/GambitReader.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string_view>
#include <unordered_map>

namespace simall::io {

// -- element-type translator --------------------------------------------------

ElementType gambit_element_type(int code, int ndp) noexcept {
    switch (code) {
        case 1: return (ndp == 2) ? ElementType::Bar2   : ElementType::Unknown;
        case 2: return (ndp == 4) ? ElementType::Quad4  : ElementType::Unknown;
        case 3: return (ndp == 3) ? ElementType::Tri3   : ElementType::Unknown;
        case 4: return (ndp == 8) ? ElementType::Hexa8  : ElementType::Unknown;
        case 5: return (ndp == 6) ? ElementType::Penta6 : ElementType::Unknown;
        case 6: return (ndp == 4) ? ElementType::Tetra4 : ElementType::Unknown;
        case 7: return (ndp == 5) ? ElementType::Pyra5  : ElementType::Unknown;
        default: return ElementType::Unknown;
    }
}

namespace {

// -- lexical helpers ---------------------------------------------------------

[[nodiscard]] std::string read_whole_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return {};
    std::ostringstream oss;
    oss << f.rdbuf();
    return oss.str();
}

void rtrim(std::string& s) {
    while (!s.empty()
           && (s.back() == '\r' || s.back() == '\n' || s.back() == ' '
               || s.back() == '\t')) {
        s.pop_back();
    }
}
void ltrim(std::string& s) {
    std::size_t i = 0;
    while (i < s.size()
           && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n')) {
        ++i;
    }
    if (i > 0) s.erase(0, i);
}
void trim(std::string& s) { rtrim(s); ltrim(s); }

[[nodiscard]] std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> out;
    std::string              cur;
    for (char c : line) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            if (!cur.empty()) { out.push_back(std::move(cur)); cur.clear(); }
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) out.push_back(std::move(cur));
    return out;
}

[[nodiscard]] bool to_long(const std::string& s, long long& out) {
    if (s.empty()) return false;
    char* end = nullptr;
    out       = std::strtoll(s.c_str(), &end, 10);
    return end != s.c_str() && (end == nullptr || *end == '\0');
}
[[nodiscard]] bool to_double(const std::string& s, double& out) {
    if (s.empty()) return false;
    char* end = nullptr;
    out       = std::strtod(s.c_str(), &end);
    return end != s.c_str() && (end == nullptr || *end == '\0');
}

[[nodiscard]] bool contains(const std::string& s, std::string_view needle) {
    return s.find(needle) != std::string::npos;
}

// -- line source -------------------------------------------------------------

struct LineSource {
    std::vector<std::string> lines;
    std::size_t              pos = 0;

    explicit LineSource(const std::string& src) {
        std::string cur;
        for (char c : src) {
            if (c == '\n') {
                rtrim(cur);
                lines.push_back(std::move(cur));
                cur.clear();
            } else {
                cur.push_back(c);
            }
        }
        if (!cur.empty()) { rtrim(cur); lines.push_back(std::move(cur)); }
    }

    bool eof() const { return pos >= lines.size(); }

    bool next_raw(std::string& out, std::size_t& lineNo) {
        if (pos >= lines.size()) return false;
        out    = lines[pos++];
        lineNo = pos;
        return true;
    }

    // Non-blank next-line.
    bool next(std::string& out, std::size_t& lineNo) {
        while (pos < lines.size()) {
            std::string s = lines[pos++];
            std::string t = s;
            ltrim(t);
            if (t.empty()) continue;
            out    = t;
            lineNo = pos;
            return true;
        }
        return false;
    }
};

// -- token stream for free-format section bodies -----------------------------
//
// GAMBIT element rows and group element-id lists wrap freely across lines;
// the only structural marker inside a section body is the literal
// `ENDOFSECTION` line.  We therefore flatten the body into a token stream
// and consume tokens until ENDOFSECTION is hit.

struct TokenStream {
    std::vector<std::string> toks;
    std::size_t              pos = 0;

    bool empty() const noexcept { return pos >= toks.size(); }
    const std::string& peek() const { return toks[pos]; }
    const std::string& take() { return toks[pos++]; }

    bool take_long(long long& v) {
        if (empty()) return false;
        return to_long(take(), v);
    }
    bool take_double(double& v) {
        if (empty()) return false;
        return to_double(take(), v);
    }
};

// Collect tokens from successive non-blank lines until ENDOFSECTION is
// encountered.  ENDOFSECTION itself is consumed.  Returns false on EOF.
[[nodiscard]] bool collect_section_body(LineSource&  src,
                                        TokenStream& ts,
                                        std::size_t& lastLine) {
    std::string line;
    std::size_t lineNo = 0;
    while (src.next(line, lineNo)) {
        lastLine = lineNo;
        std::string trimmed = line;
        trim(trimmed);
        if (trimmed == "ENDOFSECTION") return true;
        auto tk = split_ws(line);
        for (auto& t : tk) ts.toks.push_back(std::move(t));
    }
    return false;
}

// -- GAMBIT face conventions -------------------------------------------------
//
// Local face -> local node list (1-based GAMBIT node indices within the
// cell, 1-based face id).  Returns the face element type and the local
// node indices (size 2/3/4); returns empty vector for invalid (cellType,
// faceId) pairs.

struct FaceLocalNodes {
    ElementType            faceType = ElementType::Unknown;
    std::array<int, 4>     local    = {0, 0, 0, 0};
    int                    n        = 0;            // 2, 3 or 4
};

[[nodiscard]] FaceLocalNodes face_local_nodes(int cellGambitType, int faceId) noexcept {
    FaceLocalNodes f;
    switch (cellGambitType) {
        // Edge cell has no faces, skip.
        case 2: {  // Quad4 -- 4 edges (Bar2)
            static const int E[4][2] = {{1,2},{2,3},{3,4},{4,1}};
            if (faceId < 1 || faceId > 4) return f;
            f.faceType = ElementType::Bar2;
            f.n = 2;
            f.local[0] = E[faceId-1][0]; f.local[1] = E[faceId-1][1];
            return f;
        }
        case 3: {  // Tri3 -- 3 edges
            static const int E[3][2] = {{1,2},{2,3},{3,1}};
            if (faceId < 1 || faceId > 3) return f;
            f.faceType = ElementType::Bar2;
            f.n = 2;
            f.local[0] = E[faceId-1][0]; f.local[1] = E[faceId-1][1];
            return f;
        }
        case 4: {  // Hexa8 -- 6 quad faces
            static const int Q[6][4] = {
                {1,2,6,5}, {2,3,7,6}, {3,4,8,7},
                {1,5,8,4}, {1,4,3,2}, {5,6,7,8}};
            if (faceId < 1 || faceId > 6) return f;
            f.faceType = ElementType::Quad4;
            f.n = 4;
            for (int k = 0; k < 4; ++k) f.local[k] = Q[faceId-1][k];
            return f;
        }
        case 5: {  // Penta6 (wedge) -- 3 quads + 2 tris (face ids 1..3 quad, 4..5 tri)
            if (faceId == 1) {
                f.faceType = ElementType::Quad4; f.n = 4;
                f.local = {1, 2, 5, 4};
            } else if (faceId == 2) {
                f.faceType = ElementType::Quad4; f.n = 4;
                f.local = {2, 3, 6, 5};
            } else if (faceId == 3) {
                f.faceType = ElementType::Quad4; f.n = 4;
                f.local = {1, 4, 6, 3};
            } else if (faceId == 4) {
                f.faceType = ElementType::Tri3; f.n = 3;
                f.local = {1, 3, 2, 0};
            } else if (faceId == 5) {
                f.faceType = ElementType::Tri3; f.n = 3;
                f.local = {4, 5, 6, 0};
            }
            return f;
        }
        case 6: {  // Tetra4 -- 4 tri faces
            static const int T[4][3] = {{1,2,3},{1,2,4},{2,3,4},{1,3,4}};
            if (faceId < 1 || faceId > 4) return f;
            f.faceType = ElementType::Tri3;
            f.n = 3;
            for (int k = 0; k < 3; ++k) f.local[k] = T[faceId-1][k];
            return f;
        }
        case 7: {  // Pyra5 -- 1 quad base + 4 tris
            if (faceId == 1) {
                f.faceType = ElementType::Quad4; f.n = 4;
                f.local = {1, 2, 3, 4};
            } else if (faceId >= 2 && faceId <= 5) {
                static const int T[4][3] = {{1,2,5},{2,3,5},{3,4,5},{4,1,5}};
                f.faceType = ElementType::Tri3; f.n = 3;
                for (int k = 0; k < 3; ++k) f.local[k] = T[faceId-2][k];
            }
            return f;
        }
        default: return f;
    }
}

// -- parsed-element table ----------------------------------------------------

struct ParsedElement {
    int                       gambitType = 0;     // 1..7 (raw GAMBIT code)
    int                       ndp        = 0;
    std::vector<long long>    nodeTags;           // 1-based GAMBIT node tags
    long long                 groupId    = 0;     // 0 = ungrouped
};

struct ParsedGroup {
    long long                 groupId = 0;
    std::string               name;
    std::vector<long long>    elementIds;         // 1-based
};

struct ParsedBc {
    std::string                              name;
    int                                      itype  = 0;
    // For itype==1: list of (cellId, cellType, faceId)
    std::vector<std::array<long long, 3>>    faceEntries;
    // For itype==0: list of node ids (values floats discarded)
    std::vector<long long>                   nodeIds;
};

}  // namespace

// -- main parser --------------------------------------------------------------

GambitReadResult parse_gambit_neu_string(const std::string& text,
                                          std::string       sourceHint) {
    GambitReadResult r;
    r.mesh.sourceFormat = "gambit_neu";
    r.mesh.sourcePath   = sourceHint;

    LineSource  src(text);
    std::string line;
    std::size_t lineNo = 0;

    auto fail = [&](const std::string& msg) {
        r.ok    = false;
        r.error = "[" + sourceHint + ":" + std::to_string(lineNo) + "] " + msg;
        return r;
    };

    // -- CONTROL INFO ---------------------------------------------------------
    if (!src.next(line, lineNo) || !contains(line, "CONTROL INFO")) {
        return fail("expected `CONTROL INFO` as first section header");
    }
    // Skim down to the integer-counts line.  The header rows between
    // CONTROL INFO and the counts vary across GAMBIT versions but always
    // include a column-header line containing `NUMNP`.
    long long numNp = 0, numEl = 0, numGrps = 0, numBSets = 0;
    long long ndfcd = 0, ndfvl = 0;
    bool      gotCounts = false;
    while (src.next(line, lineNo)) {
        std::string t = line; trim(t);
        if (t == "ENDOFSECTION") break;
        if (contains(line, "NUMNP")) {
            // Next non-blank line must be the 6 integer values.
            if (!src.next(line, lineNo))
                return fail("EOF after NUMNP column-header line");
            auto toks = split_ws(line);
            if (toks.size() < 6)
                return fail("CONTROL INFO counts line needs 6 integers "
                            "(NUMNP NELEM NGRPS NBSETS NDFCD NDFVL)");
            if (!to_long(toks[0], numNp)    || !to_long(toks[1], numEl)
                || !to_long(toks[2], numGrps) || !to_long(toks[3], numBSets)
                || !to_long(toks[4], ndfcd)   || !to_long(toks[5], ndfvl))
                return fail("CONTROL INFO counts line: malformed integers");
            gotCounts = true;
        }
    }
    if (!gotCounts)
        return fail("CONTROL INFO did not contain a NUMNP counts line");
    if (numNp < 0 || numEl < 0 || numGrps < 0 || numBSets < 0)
        return fail("CONTROL INFO has negative counts");
    if (ndfcd != 2 && ndfcd != 3)
        return fail("CONTROL INFO: NDFCD must be 2 or 3 (got "
                    + std::to_string(ndfcd) + ")");

    UnstructuredZone zone;
    zone.name = "gambit_zone";

    // -- pass over remaining sections ----------------------------------------
    std::unordered_map<long long, NodeIdx>     tagToIdx;
    std::vector<ParsedElement>                 elems;        // index = 0-based gambit id - 1
    std::vector<ParsedGroup>                   groups;
    std::vector<ParsedBc>                      bcs;
    bool                                       sawNodes    = false;
    bool                                       sawElements = false;

    while (src.next(line, lineNo)) {
        std::string hdr = line; trim(hdr);
        if (hdr == "ENDOFSECTION") continue;   // tolerate stray markers

        // -- NODAL COORDINATES ----------------------------------------------
        if (contains(line, "NODAL COORDINATES")) {
            sawNodes = true;
            zone.x.reserve(static_cast<std::size_t>(numNp));
            zone.y.reserve(static_cast<std::size_t>(numNp));
            zone.z.reserve(static_cast<std::size_t>(numNp));
            tagToIdx.reserve(static_cast<std::size_t>(numNp));

            TokenStream ts;
            if (!collect_section_body(src, ts, lineNo))
                return fail("EOF inside NODAL COORDINATES (no ENDOFSECTION)");

            const int coordsPerNode = static_cast<int>(ndfcd);
            for (long long i = 0; i < numNp; ++i) {
                long long tag = 0;
                if (!ts.take_long(tag))
                    return fail("NODAL COORDINATES: missing node tag at row "
                                + std::to_string(i + 1));
                double xv = 0.0, yv = 0.0, zv = 0.0;
                if (!ts.take_double(xv) || !ts.take_double(yv))
                    return fail("NODAL COORDINATES: missing x/y at row "
                                + std::to_string(i + 1));
                if (coordsPerNode == 3) {
                    if (!ts.take_double(zv))
                        return fail("NODAL COORDINATES: missing z at row "
                                    + std::to_string(i + 1));
                }
                const NodeIdx idx = static_cast<NodeIdx>(zone.x.size());
                zone.x.push_back(xv);
                zone.y.push_back(yv);
                zone.z.push_back(zv);
                if (!tagToIdx.emplace(tag, idx).second)
                    return fail("NODAL COORDINATES: duplicate node tag "
                                + std::to_string(tag));
            }
            continue;
        }

        // -- ELEMENTS/CELLS --------------------------------------------------
        if (contains(line, "ELEMENTS/CELLS")) {
            sawElements = true;
            elems.assign(static_cast<std::size_t>(numEl), {});
            TokenStream ts;
            if (!collect_section_body(src, ts, lineNo))
                return fail("EOF inside ELEMENTS/CELLS");
            for (long long i = 0; i < numEl; ++i) {
                long long id = 0, type = 0, ndp = 0;
                if (!ts.take_long(id) || !ts.take_long(type) || !ts.take_long(ndp))
                    return fail("ELEMENTS/CELLS: missing id/type/ndp at row "
                                + std::to_string(i + 1));
                if (id < 1 || id > numEl)
                    return fail("ELEMENTS/CELLS: element id "
                                + std::to_string(id) + " out of range");
                if (ndp <= 0 || ndp > 27)
                    return fail("ELEMENTS/CELLS: implausible ndp "
                                + std::to_string(ndp));
                const ElementType simType = gambit_element_type(
                    static_cast<int>(type), static_cast<int>(ndp));
                if (simType == ElementType::Unknown)
                    return fail("ELEMENTS/CELLS: unsupported element "
                                "(type=" + std::to_string(type)
                                + ", ndp=" + std::to_string(ndp)
                                + "); high-order / 27-node bricks are not "
                                "yet supported");
                ParsedElement pe;
                pe.gambitType = static_cast<int>(type);
                pe.ndp        = static_cast<int>(ndp);
                pe.nodeTags.reserve(static_cast<std::size_t>(ndp));
                for (long long k = 0; k < ndp; ++k) {
                    long long nt = 0;
                    if (!ts.take_long(nt))
                        return fail("ELEMENTS/CELLS: missing node tag at "
                                    "element " + std::to_string(id));
                    pe.nodeTags.push_back(nt);
                }
                elems[static_cast<std::size_t>(id - 1)] = std::move(pe);
            }
            continue;
        }

        // -- ELEMENT GROUP ---------------------------------------------------
        if (contains(line, "ELEMENT GROUP")) {
            // Header: "GROUP:  <id>  ELEMENTS:  <m>  MATERIAL: <mat>  NFLAGS: <nf>"
            if (!src.next(line, lineNo))
                return fail("EOF inside ELEMENT GROUP (header)");
            auto htoks = split_ws(line);
            long long gid = 0, mCount = 0, mat = 0, nFlags = 0;
            // Parse by scanning for keywords.
            for (std::size_t k = 0; k + 1 < htoks.size(); ++k) {
                if (htoks[k].rfind("GROUP", 0) == 0)    to_long(htoks[k+1], gid);
                if (htoks[k].rfind("ELEMENTS", 0) == 0) to_long(htoks[k+1], mCount);
                if (htoks[k].rfind("MATERIAL", 0) == 0) to_long(htoks[k+1], mat);
                if (htoks[k].rfind("NFLAGS", 0) == 0)   to_long(htoks[k+1], nFlags);
            }
            if (mCount < 0)
                return fail("ELEMENT GROUP: negative element count");
            // Next line: group name.
            std::string gname;
            if (!src.next_raw(gname, lineNo))
                return fail("EOF inside ELEMENT GROUP (name)");
            trim(gname);
            if (gname.empty()) gname = "group_" + std::to_string(gid);
            // Body tokens: nFlags ints, then mCount element ids, then ENDOFSECTION.
            TokenStream ts;
            if (!collect_section_body(src, ts, lineNo))
                return fail("EOF inside ELEMENT GROUP body");
            for (long long k = 0; k < nFlags; ++k) {
                long long dummy = 0;
                if (!ts.take_long(dummy))
                    return fail("ELEMENT GROUP: truncated NFLAGS row");
            }
            ParsedGroup pg;
            pg.groupId = gid;
            pg.name    = gname;
            pg.elementIds.reserve(static_cast<std::size_t>(mCount));
            for (long long k = 0; k < mCount; ++k) {
                long long eid = 0;
                if (!ts.take_long(eid))
                    return fail("ELEMENT GROUP: truncated element-id list");
                if (eid < 1 || eid > numEl)
                    return fail("ELEMENT GROUP: element id "
                                + std::to_string(eid) + " out of range");
                pg.elementIds.push_back(eid);
                elems[static_cast<std::size_t>(eid - 1)].groupId = gid;
            }
            groups.push_back(std::move(pg));
            continue;
        }

        // -- BOUNDARY CONDITIONS --------------------------------------------
        if (contains(line, "BOUNDARY CONDITIONS")) {
            if (!src.next(line, lineNo))
                return fail("EOF inside BOUNDARY CONDITIONS (header)");
            auto htoks = split_ws(line);
            if (htoks.size() < 4)
                return fail("BOUNDARY CONDITIONS header needs "
                            "<name> <itype> <nentry> <nvalues> [<ibcodes>]");
            ParsedBc bc;
            bc.name = htoks[0];
            long long itype = 0, nentry = 0, nvalues = 0;
            if (!to_long(htoks[1], itype) || !to_long(htoks[2], nentry)
                || !to_long(htoks[3], nvalues))
                return fail("BOUNDARY CONDITIONS header: malformed integers");
            if (itype != 0 && itype != 1)
                return fail("BOUNDARY CONDITIONS: unsupported itype "
                            + std::to_string(itype) + " (only 0 and 1)");
            if (nentry < 0 || nvalues < 0)
                return fail("BOUNDARY CONDITIONS: negative nentry/nvalues");
            bc.itype = static_cast<int>(itype);

            TokenStream ts;
            if (!collect_section_body(src, ts, lineNo))
                return fail("EOF inside BOUNDARY CONDITIONS body");
            if (itype == 0) {
                for (long long k = 0; k < nentry; ++k) {
                    long long nid = 0;
                    if (!ts.take_long(nid))
                        return fail("BC (node): truncated entry list");
                    bc.nodeIds.push_back(nid);
                    // Discard nvalues floats per entry.
                    for (long long v = 0; v < nvalues; ++v) {
                        double dv = 0.0;
                        if (!ts.take_double(dv))
                            return fail("BC (node): truncated values block");
                    }
                }
            } else {
                for (long long k = 0; k < nentry; ++k) {
                    long long cid = 0, ctype = 0, fid = 0;
                    if (!ts.take_long(cid) || !ts.take_long(ctype)
                        || !ts.take_long(fid))
                        return fail("BC (face): truncated (cell,type,face) entry");
                    bc.faceEntries.push_back({cid, ctype, fid});
                    for (long long v = 0; v < nvalues; ++v) {
                        double dv = 0.0;
                        if (!ts.take_double(dv))
                            return fail("BC (face): truncated values block");
                    }
                }
            }
            bcs.push_back(std::move(bc));
            continue;
        }

        // Unknown header: skip to its ENDOFSECTION.
        TokenStream junk;
        (void)collect_section_body(src, junk, lineNo);
    }

    if (!sawNodes)    return fail("file did not contain a NODAL COORDINATES section");
    if (!sawElements) return fail("file did not contain an ELEMENTS/CELLS section");

    // -- build volume element sections, grouped by (groupId, elemType) -------
    struct VolSecKey {
        long long   groupId;
        ElementType type;
        bool operator==(const VolSecKey& o) const noexcept {
            return groupId == o.groupId && type == o.type;
        }
    };
    struct VolSecHash {
        std::size_t operator()(const VolSecKey& k) const noexcept {
            return std::hash<long long>{}(
                (k.groupId << 4)
                ^ static_cast<std::size_t>(k.type));
        }
    };
    std::unordered_map<VolSecKey, std::size_t, VolSecHash> volSecIdx;
    auto find_group_name = [&](long long gid) -> std::string {
        for (const auto& g : groups)
            if (g.groupId == gid) return g.name;
        return "ungrouped";
    };

    // Use 0-based element index -> {sectionIdx, elementIndexInSection}
    std::vector<std::pair<std::size_t, std::size_t>> elemSectionMap;
    elemSectionMap.assign(static_cast<std::size_t>(numEl), {static_cast<std::size_t>(-1), 0});

    std::uint8_t maxDim = 0;
    for (std::size_t i = 0; i < elems.size(); ++i) {
        const auto& pe = elems[i];
        if (pe.gambitType == 0) continue;   // never populated (gap in id range)
        const ElementType simType =
            gambit_element_type(pe.gambitType, pe.ndp);
        if (simType == ElementType::Unknown) continue;
        VolSecKey key{pe.groupId, simType};
        std::size_t sIdx;
        if (auto it = volSecIdx.find(key); it != volSecIdx.end()) {
            sIdx = it->second;
        } else {
            ElementSection sec;
            sec.type = simType;
            sec.name = find_group_name(pe.groupId) + "_"
                     + element_type_name(simType);
            sIdx = zone.sections.size();
            zone.sections.push_back(std::move(sec));
            volSecIdx.emplace(key, sIdx);
        }
        ElementSection& sec = zone.sections[sIdx];
        const std::size_t localElIdx = sec.element_count();
        for (long long nt : pe.nodeTags) {
            auto it = tagToIdx.find(nt);
            if (it == tagToIdx.end())
                return fail("ELEMENTS/CELLS: element references unknown "
                            "node tag " + std::to_string(nt));
            sec.nodes.push_back(it->second);
        }
        elemSectionMap[i] = {sIdx, localElIdx};
        const std::uint8_t vpe = vertices_per_element(simType);
        std::uint8_t d = 0;
        switch (simType) {
            case ElementType::Bar2:                                  d = 1; break;
            case ElementType::Tri3:  case ElementType::Quad4:        d = 2; break;
            case ElementType::Tetra4:case ElementType::Hexa8:
            case ElementType::Penta6:case ElementType::Pyra5:        d = 3; break;
            default: break;
        }
        (void)vpe;
        if (d > maxDim) maxDim = d;
    }

    // -- build face sections + boundary patches from face-BCs -----------------
    //
    // For each BC patch, walk its (cellId, cellType, faceId) entries,
    // extract face nodes from the parent volume element, and append to a
    // BC-local face section (one per face-element-type, since a single
    // patch may mix tris and quads, e.g. on a wedge mesh).

    for (const auto& bc : bcs) {
        BoundaryPatch bp;
        bp.name   = bc.name;
        bp.bcType = (bc.itype == 1) ? "wall" : "node_set";

        if (bc.itype == 0) {
            // Node-based BC: we keep the patch but it carries no face elements.
            zone.boundaries.push_back(std::move(bp));
            continue;
        }

        // (faceType -> section index inside zone.sections, created on demand)
        std::unordered_map<int, std::size_t> faceSecIdx;
        auto get_face_sec = [&](ElementType ft) -> ElementSection& {
            const int key = static_cast<int>(ft);
            if (auto it = faceSecIdx.find(key); it != faceSecIdx.end())
                return zone.sections[it->second];
            ElementSection sec;
            sec.type = ft;
            sec.name = bc.name + "_" + element_type_name(ft);
            const std::size_t sIdx = zone.sections.size();
            zone.sections.push_back(std::move(sec));
            faceSecIdx.emplace(key, sIdx);
            bp.faceElementIndices.push_back(static_cast<std::uint32_t>(sIdx));
            return zone.sections[sIdx];
        };

        for (const auto& fe : bc.faceEntries) {
            const long long cellId   = fe[0];
            const long long cellType = fe[1];
            const long long faceId   = fe[2];
            if (cellId < 1 || cellId > numEl)
                return fail("BC '" + bc.name + "': cell id "
                            + std::to_string(cellId) + " out of range");
            const ParsedElement& pe =
                elems[static_cast<std::size_t>(cellId - 1)];
            if (pe.gambitType == 0)
                return fail("BC '" + bc.name + "': cell id "
                            + std::to_string(cellId) + " was never defined");
            if (static_cast<long long>(pe.gambitType) != cellType)
                return fail("BC '" + bc.name + "': cell-type mismatch "
                            "for cell " + std::to_string(cellId)
                            + " (expected " + std::to_string(pe.gambitType)
                            + ", got " + std::to_string(cellType) + ")");
            FaceLocalNodes fln = face_local_nodes(pe.gambitType,
                                                  static_cast<int>(faceId));
            if (fln.n == 0)
                return fail("BC '" + bc.name + "': invalid face id "
                            + std::to_string(faceId)
                            + " for cell type " + std::to_string(cellType));
            ElementSection& fsec = get_face_sec(fln.faceType);
            for (int k = 0; k < fln.n; ++k) {
                const int localOneBased = fln.local[k];
                if (localOneBased < 1
                    || localOneBased > static_cast<int>(pe.nodeTags.size()))
                    return fail("internal: face_local_nodes returned "
                                "out-of-range local index");
                const long long nt = pe.nodeTags[
                    static_cast<std::size_t>(localOneBased - 1)];
                auto it = tagToIdx.find(nt);
                if (it == tagToIdx.end())
                    return fail("BC '" + bc.name + "': cell "
                                + std::to_string(cellId)
                                + " references unknown node tag "
                                + std::to_string(nt));
                fsec.nodes.push_back(it->second);
            }
        }
        zone.boundaries.push_back(std::move(bp));
    }

    // Deterministic boundary order by name.
    std::sort(zone.boundaries.begin(), zone.boundaries.end(),
              [](const BoundaryPatch& a, const BoundaryPatch& b) {
                  return a.name < b.name;
              });

    if (zone.x.empty())
        return fail("file contained no nodes");
    if (zone.sections.empty())
        return fail("file contained no supported elements");

    r.mesh.zones.push_back(std::move(zone));
    r.dimension = maxDim ? maxDim : static_cast<std::uint8_t>(ndfcd);
    r.ok        = true;
    return r;
}

GambitReadResult read_gambit_neu(const std::string& path) {
    GambitReadResult r;
    std::string      text = read_whole_file(path);
    if (text.empty()) {
        r.ok    = false;
        r.error = "could not open or read '" + path + "'";
        return r;
    }
    return parse_gambit_neu_string(text, path);
}

}  // namespace simall::io
