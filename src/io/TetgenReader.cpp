// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/TetgenReader.cpp
// Phase  : 23 Pass 8
//
// Implementation of the TetGen `.node`/`.ele`/`.face` triple-file reader.
// See the header for the grammar reference.
// =============================================================================
#include "io/TetgenReader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string_view>
#include <unordered_map>

namespace simall::io
{

namespace
{

// -- lexical helpers ---------------------------------------------------------

[[nodiscard]] std::string read_whole_file(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return {};
    std::ostringstream oss;
    oss << f.rdbuf();
    return oss.str();
}

void rtrim(std::string& s)
{
    while (!s.empty()
           && (s.back() == '\r' || s.back() == '\n' || s.back() == ' ' || s.back() == '\t')) {
        s.pop_back();
    }
}
void ltrim(std::string& s)
{
    std::size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n')) {
        ++i;
    }
    if (i > 0)
        s.erase(0, i);
}
void trim(std::string& s)
{
    rtrim(s);
    ltrim(s);
}

// Strip an in-line `#` comment (everything from the first `#` to EOL).
void strip_comment(std::string& s)
{
    auto p = s.find('#');
    if (p != std::string::npos)
        s.erase(p);
}

[[nodiscard]] std::vector<std::string> split_ws(const std::string& line)
{
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            if (!cur.empty()) {
                out.push_back(std::move(cur));
                cur.clear();
            }
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty())
        out.push_back(std::move(cur));
    return out;
}

[[nodiscard]] bool to_long(const std::string& s, long long& out)
{
    if (s.empty())
        return false;
    char* end = nullptr;
    out = std::strtoll(s.c_str(), &end, 10);
    return end != s.c_str() && (end == nullptr || *end == '\0');
}
[[nodiscard]] bool to_double(const std::string& s, double& out)
{
    if (s.empty())
        return false;
    char* end = nullptr;
    out = std::strtod(s.c_str(), &end);
    return end != s.c_str() && (end == nullptr || *end == '\0');
}

// -- line source that auto-skips comments / blanks ---------------------------

struct DataLineSource
{
    std::vector<std::string> lines;
    std::size_t pos = 0;

    explicit DataLineSource(const std::string& src)
    {
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
        if (!cur.empty()) {
            rtrim(cur);
            lines.push_back(std::move(cur));
        }
    }

    /// Yield the next line that contains at least one non-comment token.
    /// `toks` is populated with the whitespace-split tokens (comment
    /// stripped); `lineNo` is the 1-based line number.
    bool next(std::vector<std::string>& toks, std::size_t& lineNo)
    {
        while (pos < lines.size()) {
            std::string s = lines[pos++];
            strip_comment(s);
            trim(s);
            if (s.empty())
                continue;
            toks = split_ws(s);
            lineNo = pos;
            return true;
        }
        return false;
    }
};

} // namespace

// -- main parser --------------------------------------------------------------

TetgenReadResult parse_tetgen_strings(const std::string& nodeText,
                                      const std::string& eleText,
                                      const std::string& faceText,
                                      std::string sourceHint)
{
    TetgenReadResult r;
    r.mesh.sourceFormat = "tetgen";
    r.mesh.sourcePath = sourceHint;

    std::size_t lineNo = 0;
    auto fail = [&](const std::string& which, const std::string& msg) {
        r.ok = false;
        r.error = "[" + sourceHint + ":" + which + ":" + std::to_string(lineNo) + "] " + msg;
        return r;
    };

    UnstructuredZone zone;
    zone.name = "tetgen_zone";

    // -- .node ----------------------------------------------------------------
    std::unordered_map<long long, NodeIdx> tagToIdx;
    long long firstNodeId = 0;
    bool sawFirstNode = false;
    {
        DataLineSource src(nodeText);
        std::vector<std::string> toks;
        if (!src.next(toks, lineNo))
            return fail(".node", "file is empty");
        if (toks.size() < 4)
            return fail(".node",
                        "header needs 4 ints "
                        "(nNodes dim nAttrs hasMarker)");
        long long nNodes = 0, dim = 0, nAttrs = 0, hasMarker = 0;
        if (!to_long(toks[0], nNodes) || !to_long(toks[1], dim) || !to_long(toks[2], nAttrs)
            || !to_long(toks[3], hasMarker))
            return fail(".node", "header has malformed integers");
        if (nNodes < 0)
            return fail(".node", "negative node count");
        if (dim != 2 && dim != 3)
            return fail(".node", "dim must be 2 or 3 (got " + std::to_string(dim) + ")");
        if (nAttrs < 0 || hasMarker < 0 || hasMarker > 1)
            return fail(".node", "negative nAttrs or out-of-range marker flag");

        zone.x.reserve(static_cast<std::size_t>(nNodes));
        zone.y.reserve(static_cast<std::size_t>(nNodes));
        zone.z.reserve(static_cast<std::size_t>(nNodes));
        tagToIdx.reserve(static_cast<std::size_t>(nNodes));

        const std::size_t expected = 1 + static_cast<std::size_t>(dim)
                                     + static_cast<std::size_t>(nAttrs)
                                     + static_cast<std::size_t>(hasMarker);
        for (long long i = 0; i < nNodes; ++i) {
            if (!src.next(toks, lineNo))
                return fail(".node", "unexpected EOF inside node table");
            if (toks.size() < expected)
                return fail(".node",
                            "row too short (need " + std::to_string(expected) + " tokens, got "
                                + std::to_string(toks.size()) + ")");
            long long tag = 0;
            if (!to_long(toks[0], tag))
                return fail(".node", "malformed node tag");
            if (!sawFirstNode) {
                firstNodeId = tag;
                sawFirstNode = true;
            }
            double xv = 0.0, yv = 0.0, zv = 0.0;
            if (!to_double(toks[1], xv) || !to_double(toks[2], yv))
                return fail(".node", "malformed x/y");
            if (dim == 3) {
                if (!to_double(toks[3], zv))
                    return fail(".node", "malformed z");
            }
            const NodeIdx idx = static_cast<NodeIdx>(zone.x.size());
            zone.x.push_back(xv);
            zone.y.push_back(yv);
            zone.z.push_back(zv);
            if (!tagToIdx.emplace(tag, idx).second)
                return fail(".node", "duplicate node tag " + std::to_string(tag));
        }
    }

    // The detected first-index is informational only; tagToIdx already
    // handles arbitrary tag ranges.

    // -- .ele -----------------------------------------------------------------
    {
        DataLineSource src(eleText);
        std::vector<std::string> toks;
        if (!src.next(toks, lineNo))
            return fail(".ele", "file is empty");
        if (toks.size() < 3)
            return fail(".ele",
                        "header needs 3 ints "
                        "(nTets nodesPerTet nAttrs)");
        long long nTets = 0, npe = 0, nAttrs = 0;
        if (!to_long(toks[0], nTets) || !to_long(toks[1], npe) || !to_long(toks[2], nAttrs))
            return fail(".ele", "header has malformed integers");
        if (nTets < 0)
            return fail(".ele", "negative tet count");
        if (npe != 4)
            return fail(".ele",
                        "only linear tetrahedra are supported "
                        "(nodesPerTet="
                            + std::to_string(npe)
                            + "; 10-node parabolic tets / TetGen `-o2` output "
                              "are not yet implemented)");
        if (nAttrs < 0)
            return fail(".ele", "negative attribute count");

        // Sections keyed by attribute value (or by 0 if no attributes).
        std::unordered_map<long long, std::size_t> attrToSec;
        auto get_section = [&](long long attr) -> ElementSection& {
            if (auto it = attrToSec.find(attr); it != attrToSec.end())
                return zone.sections[it->second];
            ElementSection sec;
            sec.type = ElementType::Tetra4;
            sec.name =
                (nAttrs > 0) ? ("region_" + std::to_string(attr)) : std::string("tetgen_volume");
            const std::size_t sIdx = zone.sections.size();
            zone.sections.push_back(std::move(sec));
            attrToSec.emplace(attr, sIdx);
            return zone.sections[sIdx];
        };

        const std::size_t expected = 1 + 4 + static_cast<std::size_t>(nAttrs);
        for (long long i = 0; i < nTets; ++i) {
            if (!src.next(toks, lineNo))
                return fail(".ele", "unexpected EOF inside tet table");
            if (toks.size() < expected)
                return fail(".ele",
                            "row too short (need " + std::to_string(expected) + " tokens, got "
                                + std::to_string(toks.size()) + ")");
            long long n[4] = {0, 0, 0, 0};
            for (int k = 0; k < 4; ++k) {
                if (!to_long(toks[1 + k], n[k]))
                    return fail(".ele", "malformed node tag in tet row");
            }
            long long attr = 0;
            if (nAttrs > 0) {
                if (!to_long(toks[5], attr))
                    return fail(".ele", "malformed attribute in tet row");
            }
            ElementSection& sec = get_section(attr);
            for (int k = 0; k < 4; ++k) {
                auto it = tagToIdx.find(n[k]);
                if (it == tagToIdx.end())
                    return fail(".ele", "tet references unknown node tag " + std::to_string(n[k]));
                sec.nodes.push_back(it->second);
            }
        }
    }

    // -- .face (optional) -----------------------------------------------------
    if (!faceText.empty()) {
        DataLineSource src(faceText);
        std::vector<std::string> toks;
        if (!src.next(toks, lineNo))
            return fail(".face", "file is empty (use empty string to skip)");
        if (toks.size() < 2)
            return fail(".face", "header needs 2 ints (nFaces hasMarker)");
        long long nFaces = 0, hasMarker = 0;
        if (!to_long(toks[0], nFaces) || !to_long(toks[1], hasMarker))
            return fail(".face", "header has malformed integers");
        if (nFaces < 0)
            return fail(".face", "negative face count");
        if (hasMarker < 0 || hasMarker > 1)
            return fail(".face", "marker flag must be 0 or 1");

        // markerValue -> (faceSecIdx, patchIdx)
        std::unordered_map<long long, std::pair<std::size_t, std::size_t>> markerToBuckets;
        auto get_buckets = [&](long long marker) -> std::pair<std::size_t, std::size_t> {
            if (auto it = markerToBuckets.find(marker); it != markerToBuckets.end())
                return it->second;
            ElementSection sec;
            sec.type = ElementType::Tri3;
            sec.name =
                (hasMarker == 1) ? ("marker_" + std::to_string(marker)) : std::string("all_faces");
            const std::size_t sIdx = zone.sections.size();
            zone.sections.push_back(std::move(sec));

            BoundaryPatch bp;
            bp.name = sec.name;
            bp.bcType = "wall";
            bp.faceElementIndices.push_back(static_cast<std::uint32_t>(sIdx));
            const std::size_t pIdx = zone.boundaries.size();
            zone.boundaries.push_back(std::move(bp));

            markerToBuckets.emplace(marker, std::make_pair(sIdx, pIdx));
            return {sIdx, pIdx};
        };

        const std::size_t expected = 1 + 3 + static_cast<std::size_t>(hasMarker);
        for (long long i = 0; i < nFaces; ++i) {
            if (!src.next(toks, lineNo))
                return fail(".face", "unexpected EOF inside face table");
            if (toks.size() < expected)
                return fail(".face",
                            "row too short (need " + std::to_string(expected) + " tokens, got "
                                + std::to_string(toks.size()) + ")");
            long long n[3] = {0, 0, 0};
            for (int k = 0; k < 3; ++k) {
                if (!to_long(toks[1 + k], n[k]))
                    return fail(".face", "malformed node tag in face row");
            }
            long long marker = 0;
            if (hasMarker == 1) {
                if (!to_long(toks[4], marker))
                    return fail(".face", "malformed marker in face row");
            }
            auto buckets = get_buckets(marker);
            ElementSection& fsec = zone.sections[buckets.first];
            for (int k = 0; k < 3; ++k) {
                auto it = tagToIdx.find(n[k]);
                if (it == tagToIdx.end())
                    return fail(".face",
                                "face references unknown node tag " + std::to_string(n[k]));
                fsec.nodes.push_back(it->second);
            }
        }

        // Deterministic patch order: by name.
        std::sort(zone.boundaries.begin(),
                  zone.boundaries.end(),
                  [](const BoundaryPatch& a, const BoundaryPatch& b) { return a.name < b.name; });
    }

    if (zone.x.empty())
        return fail("", "file contained no nodes");
    if (zone.sections.empty())
        return fail("", "file contained no supported elements");

    r.mesh.zones.push_back(std::move(zone));
    r.dimension = 3;
    r.ok = true;
    return r;
}

TetgenReadResult read_tetgen(const std::string& stemPath)
{
    TetgenReadResult r;
    const std::string nodePath = stemPath + ".node";
    const std::string elePath = stemPath + ".ele";
    const std::string facePath = stemPath + ".face";

    const std::string nodeText = read_whole_file(nodePath);
    if (nodeText.empty()) {
        r.ok = false;
        r.error = "could not open or read '" + nodePath + "'";
        return r;
    }
    const std::string eleText = read_whole_file(elePath);
    if (eleText.empty()) {
        r.ok = false;
        r.error = "could not open or read '" + elePath + "'";
        return r;
    }
    // .face is optional -- silently treat as missing if absent.
    const std::string faceText = read_whole_file(facePath);
    return parse_tetgen_strings(nodeText, eleText, faceText, stemPath);
}

} // namespace simall::io
