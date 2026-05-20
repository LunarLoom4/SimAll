// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/GmshReader.cpp
// Phase  : 23 Pass 1 (v4.1 path); Pass 6 added legacy v2.2 path
//
// Line-oriented parser for Gmsh ASCII `.msh`.  The grammar is strictly
// section-based ($Foo ... $EndFoo) so the implementation is a small state
// machine: read a section header, dispatch to a per-section handler, expect
// the matching $End marker, repeat to EOF.  The $Nodes and $Elements
// handlers branch on the format version detected in $MeshFormat (v2 uses
// a flat list, v4 uses per-entity blocks).
// =============================================================================
#include "io/GmshReader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

namespace simall::io
{

namespace
{

// -- small lexical helpers ---------------------------------------------------

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

[[nodiscard]] std::vector<std::string> split_ws(const std::string& line)
{
    std::vector<std::string> out;
    std::string cur;
    bool inQuote = false;
    for (char c : line) {
        if (c == '"') {
            inQuote = !inQuote;
            cur.push_back(c);
            continue;
        }
        if (!inQuote && std::isspace(static_cast<unsigned char>(c))) {
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

[[nodiscard]] bool to_int(const std::string& s, long long& out)
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

// Strip surrounding ASCII double quotes from a token if present.
[[nodiscard]] std::string strip_quotes(std::string s)
{
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
        return s.substr(1, s.size() - 2);
    }
    return s;
}

// -- line-based source --------------------------------------------------------
//
// We slurp the entire file once and walk it line-by-line. Errors are emitted
// with 1-based line numbers for diagnostics.

struct LineSource
{
    std::vector<std::string> lines;
    std::size_t pos = 0;

    explicit LineSource(const std::string& src)
    {
        lines.reserve(256);
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

    bool eof() const { return pos >= lines.size(); }

    // Return next non-blank line (advancing pos).  Returns false at EOF.
    bool next(std::string& out, std::size_t& lineNo)
    {
        while (pos < lines.size()) {
            std::string& s = lines[pos++];
            std::string trimmed = s;
            ltrim(trimmed);
            if (trimmed.empty())
                continue;
            out = trimmed;
            lineNo = pos; // 1-based: pos has already advanced
            return true;
        }
        return false;
    }
};

// -- physical-name entry ------------------------------------------------------

struct PhysicalName
{
    int dim = 0;
    long long tag = 0;
    std::string name;
};

// One Gmsh "entity" (carries physical-tag membership).
struct Entity
{
    int dim = 0;
    long long tag = 0;
    std::vector<long long> physicalTags;
};

} // namespace

// -- public element-type translator ------------------------------------------

ElementType gmsh_element_type(int gmshCode) noexcept
{
    switch (gmshCode) {
    case 1:
        return ElementType::Bar2;
    case 2:
        return ElementType::Tri3;
    case 3:
        return ElementType::Quad4;
    case 4:
        return ElementType::Tetra4;
    case 5:
        return ElementType::Hexa8;
    case 6:
        return ElementType::Penta6;
    case 7:
        return ElementType::Pyra5;
    default:
        return ElementType::Unknown;
    }
}

namespace
{

// Topological dimension of a Gmsh element type (used by the v2 parser to
// classify elements that are not tagged with an entity hierarchy).
[[nodiscard]] int dim_from_gmsh_elem_type(int gmshCode) noexcept
{
    switch (gmshCode) {
    case 1:
        return 1;
    case 2:
    case 3:
        return 2;
    case 4:
    case 5:
    case 6:
    case 7:
        return 3;
    case 15:
        return 0;
    default:
        return -1;
    }
}

} // namespace

// -- main parser --------------------------------------------------------------

GmshReadResult parse_gmsh_msh_string(const std::string& text, std::string sourceHint)
{
    GmshReadResult r;
    r.mesh.sourceFormat = "gmsh_msh";
    r.mesh.sourcePath = sourceHint;

    LineSource src(text);
    std::string line;
    std::size_t lineNo = 0;

    auto fail = [&](const std::string& msg) {
        r.ok = false;
        r.error = "[" + sourceHint + ":" + std::to_string(lineNo) + "] " + msg;
        return r;
    };

    // -- $MeshFormat (must come first) ----------------------------------------
    bool isV2 = false;
    if (!src.next(line, lineNo) || line != "$MeshFormat") {
        return fail("expected $MeshFormat as first section");
    }
    if (!src.next(line, lineNo)) {
        return fail("unexpected EOF inside $MeshFormat");
    }
    {
        auto toks = split_ws(line);
        if (toks.size() < 3) {
            return fail("$MeshFormat must have 3 tokens (version file-type data-size)");
        }
        double version = 0.0;
        long long fileType = 0;
        if (!to_double(toks[0], version) || !to_int(toks[1], fileType)) {
            return fail("$MeshFormat: malformed version / file-type");
        }
        if (fileType != 0) {
            return fail("binary Gmsh .msh is not supported (file-type=" + std::to_string(fileType)
                        + "); please re-export as ASCII");
        }
        if (version >= 2.0 && version < 3.0) {
            isV2 = true;
        } else if (version >= 4.0 && version < 5.0) {
            isV2 = false;
        } else {
            return fail("unsupported Gmsh format version " + toks[0]
                        + " (supported: 2.x and 4.x ASCII)");
        }
    }
    if (!src.next(line, lineNo) || line != "$EndMeshFormat") {
        return fail("expected $EndMeshFormat");
    }

    // -- $PhysicalNames / $Entities are optional; collected if present --------
    std::vector<PhysicalName> physicals;
    std::vector<Entity> entities; // all dims pooled

    // -- node table -----------------------------------------------------------
    // Gmsh node tags are 1-based and may be sparse; we build a tag -> idx map.
    std::unordered_map<long long, NodeIdx> tagToIdx;
    UnstructuredZone zone;
    zone.name = "gmsh_zone";

    // Pre-create a working section list keyed by (entityDim, entityTag,
    // gmshType) so multiple element blocks with the same key merge.
    struct SecKey
    {
        int entityDim;
        long long entityTag;
        int gmshType;
        bool operator==(const SecKey& o) const noexcept
        {
            return entityDim == o.entityDim && entityTag == o.entityTag && gmshType == o.gmshType;
        }
    };
    struct SecKeyHash
    {
        std::size_t operator()(const SecKey& k) const noexcept
        {
            return std::hash<long long>{}((static_cast<long long>(k.entityDim) << 56)
                                          ^ (static_cast<long long>(k.gmshType) << 48)
                                          ^ k.entityTag);
        }
    };
    std::unordered_map<SecKey, std::size_t, SecKeyHash> secIndex;

    // (entityDim, entityTag) -> [section-index, ...] for boundary derivation.
    std::unordered_map<long long, std::vector<std::size_t>> entityKey_to_sections;
    auto pack_entity_key = [](int dim, long long tag) -> long long {
        return (static_cast<long long>(dim) << 56) ^ tag;
    };

    std::uint8_t maxDim = 0;

    // -- section-pump loop ----------------------------------------------------
    while (src.next(line, lineNo)) {
        if (line == "$PhysicalNames") {
            if (!src.next(line, lineNo))
                return fail("EOF in $PhysicalNames");
            long long n = 0;
            if (!to_int(line, n) || n < 0)
                return fail("$PhysicalNames count is not a non-negative integer");
            physicals.reserve(static_cast<std::size_t>(n));
            for (long long i = 0; i < n; ++i) {
                if (!src.next(line, lineNo))
                    return fail("EOF inside $PhysicalNames body");
                auto toks = split_ws(line);
                if (toks.size() < 3)
                    return fail("$PhysicalNames row needs dim tag \"name\"");
                PhysicalName pn;
                long long d = 0;
                if (!to_int(toks[0], d) || !to_int(toks[1], pn.tag))
                    return fail("$PhysicalNames row has malformed dim/tag");
                pn.dim = static_cast<int>(d);
                pn.name = strip_quotes(toks[2]);
                physicals.push_back(std::move(pn));
            }
            if (!src.next(line, lineNo) || line != "$EndPhysicalNames")
                return fail("expected $EndPhysicalNames");
            continue;
        }

        if (line == "$Entities") {
            if (!src.next(line, lineNo))
                return fail("EOF in $Entities");
            auto toks = split_ws(line);
            if (toks.size() < 4)
                return fail("$Entities header needs 4 counts");
            long long counts[4] = {0, 0, 0, 0};
            for (int d = 0; d < 4; ++d) {
                if (!to_int(toks[d], counts[d]) || counts[d] < 0)
                    return fail("$Entities count is not a non-negative integer");
            }
            for (int dim = 0; dim < 4; ++dim) {
                for (long long i = 0; i < counts[dim]; ++i) {
                    if (!src.next(line, lineNo))
                        return fail("EOF inside $Entities body");
                    auto erow = split_ws(line);
                    if (erow.empty())
                        return fail("$Entities row is empty");
                    Entity e;
                    e.dim = dim;
                    long long t = 0;
                    if (!to_int(erow[0], t))
                        return fail("$Entities row: malformed entity tag");
                    e.tag = t;
                    // Skip bbox: points have 3 coords, others have 6.
                    // Then numPhysicalTags + physTags, then for dim>0
                    // numBoundingEntities + bounding tags. We do not need
                    // the bounding tags so we parse only up to the
                    // physical-tag list.
                    std::size_t k = 1 + static_cast<std::size_t>(dim == 0 ? 3 : 6);
                    if (k >= erow.size())
                        return fail("$Entities row too short");
                    long long nPhys = 0;
                    if (!to_int(erow[k], nPhys) || nPhys < 0)
                        return fail("$Entities: malformed numPhysicalTags");
                    ++k;
                    for (long long p = 0; p < nPhys; ++p) {
                        if (k >= erow.size())
                            return fail("$Entities: physical-tag list truncated");
                        long long pt = 0;
                        if (!to_int(erow[k++], pt))
                            return fail("$Entities: malformed physical tag");
                        e.physicalTags.push_back(pt);
                    }
                    entities.push_back(std::move(e));
                }
            }
            if (!src.next(line, lineNo) || line != "$EndEntities")
                return fail("expected $EndEntities");
            continue;
        }

        if (line == "$Nodes") {
            if (!src.next(line, lineNo))
                return fail("EOF in $Nodes");
            if (isV2) {
                long long n = 0;
                if (!to_int(line, n) || n < 0)
                    return fail("$Nodes v2 count is not a non-negative integer");
                zone.x.reserve(static_cast<std::size_t>(n));
                zone.y.reserve(static_cast<std::size_t>(n));
                zone.z.reserve(static_cast<std::size_t>(n));
                tagToIdx.reserve(static_cast<std::size_t>(n));
                for (long long i = 0; i < n; ++i) {
                    if (!src.next(line, lineNo))
                        return fail("EOF inside $Nodes (v2 rows)");
                    auto row = split_ws(line);
                    if (row.size() < 4)
                        return fail("$Nodes v2 row needs tag x y z");
                    long long t = 0;
                    double x = 0.0, y = 0.0, z = 0.0;
                    if (!to_int(row[0], t) || !to_double(row[1], x) || !to_double(row[2], y)
                        || !to_double(row[3], z))
                        return fail("$Nodes v2 row malformed");
                    const NodeIdx idx = static_cast<NodeIdx>(zone.x.size());
                    zone.x.push_back(x);
                    zone.y.push_back(y);
                    zone.z.push_back(z);
                    tagToIdx.emplace(t, idx);
                }
                if (!src.next(line, lineNo) || line != "$EndNodes")
                    return fail("expected $EndNodes");
                continue;
            }
            auto toks = split_ws(line);
            if (toks.size() < 4)
                return fail("$Nodes header needs 4 ints");
            long long numBlocks = 0, numNodes = 0, minTag = 0, maxTag = 0;
            if (!to_int(toks[0], numBlocks) || !to_int(toks[1], numNodes)
                || !to_int(toks[2], minTag) || !to_int(toks[3], maxTag))
                return fail("$Nodes header has malformed integers");
            if (numBlocks < 0 || numNodes < 0)
                return fail("$Nodes header has negative counts");

            zone.x.reserve(static_cast<std::size_t>(numNodes));
            zone.y.reserve(static_cast<std::size_t>(numNodes));
            zone.z.reserve(static_cast<std::size_t>(numNodes));
            tagToIdx.reserve(static_cast<std::size_t>(numNodes));

            for (long long b = 0; b < numBlocks; ++b) {
                if (!src.next(line, lineNo))
                    return fail("EOF inside $Nodes body (block header)");
                auto bhdr = split_ws(line);
                if (bhdr.size() < 4)
                    return fail("$Nodes block header needs 4 ints");
                long long entityDim = 0;
                long long entityTag = 0;
                long long parametric = 0;
                long long nInBlock = 0;
                if (!to_int(bhdr[0], entityDim) || !to_int(bhdr[1], entityTag)
                    || !to_int(bhdr[2], parametric) || !to_int(bhdr[3], nInBlock))
                    return fail("$Nodes block header malformed");
                if (parametric != 0)
                    return fail("parametric node blocks are not supported");
                if (nInBlock < 0)
                    return fail("$Nodes block has negative count");

                std::vector<long long> tags;
                tags.reserve(static_cast<std::size_t>(nInBlock));
                for (long long i = 0; i < nInBlock; ++i) {
                    if (!src.next(line, lineNo))
                        return fail("EOF in $Nodes tag list");
                    long long t = 0;
                    auto trow = split_ws(line);
                    if (trow.empty() || !to_int(trow[0], t))
                        return fail("$Nodes tag row malformed");
                    tags.push_back(t);
                }
                for (long long i = 0; i < nInBlock; ++i) {
                    if (!src.next(line, lineNo))
                        return fail("EOF in $Nodes coord list");
                    auto crow = split_ws(line);
                    if (crow.size() < 3)
                        return fail("$Nodes coord row needs 3 values");
                    double x = 0.0, y = 0.0, z = 0.0;
                    if (!to_double(crow[0], x) || !to_double(crow[1], y) || !to_double(crow[2], z))
                        return fail("$Nodes coord row malformed");
                    const NodeIdx idx = static_cast<NodeIdx>(zone.x.size());
                    zone.x.push_back(x);
                    zone.y.push_back(y);
                    zone.z.push_back(z);
                    tagToIdx.emplace(tags[static_cast<std::size_t>(i)], idx);
                }
            }
            if (!src.next(line, lineNo) || line != "$EndNodes")
                return fail("expected $EndNodes");
            continue;
        }

        if (line == "$Elements") {
            if (!src.next(line, lineNo))
                return fail("EOF in $Elements");
            if (isV2) {
                long long n = 0;
                if (!to_int(line, n) || n < 0)
                    return fail("$Elements v2 count is not a non-negative integer");
                std::unordered_set<long long> seenEntKeys;
                for (long long i = 0; i < n; ++i) {
                    if (!src.next(line, lineNo))
                        return fail("EOF in $Elements (v2 rows)");
                    auto row = split_ws(line);
                    if (row.size() < 3)
                        return fail("$Elements v2 row needs elm-num elm-type num-tags ...");
                    long long elNum = 0, elType = 0, numTags = 0;
                    if (!to_int(row[0], elNum) || !to_int(row[1], elType)
                        || !to_int(row[2], numTags) || numTags < 0)
                        return fail("$Elements v2 row: malformed leading ints");
                    if (row.size() < static_cast<std::size_t>(3 + numTags))
                        return fail("$Elements v2 row: tag list truncated");
                    long long physTag = 0;
                    if (numTags >= 1) {
                        if (!to_int(row[3], physTag))
                            return fail("$Elements v2 row: malformed physical tag");
                    }
                    if (elType == 15) {
                        // 1-node point: silently skipped (no topology).
                        continue;
                    }
                    const ElementType simType = gmsh_element_type(static_cast<int>(elType));
                    if (simType == ElementType::Unknown)
                        return fail("unsupported Gmsh element type code " + std::to_string(elType));
                    const int dim = dim_from_gmsh_elem_type(static_cast<int>(elType));
                    if (dim < 0)
                        return fail("internal: dim<0 for elType " + std::to_string(elType));
                    const std::uint8_t vpe = vertices_per_element(simType);
                    if (vpe == 0)
                        return fail("internal: vpe==0");
                    const std::size_t expected = 3 + static_cast<std::size_t>(numTags) + vpe;
                    if (row.size() < expected)
                        return fail("$Elements v2 row has too few nodes for "
                                    + std::string(element_type_name(simType)));

                    // For v2 we substitute (dim, physTag) for the missing
                    // (entityDim, entityTag) pair.
                    SecKey key{dim, physTag, static_cast<int>(elType)};
                    std::size_t secIdx;
                    if (auto it = secIndex.find(key); it != secIndex.end()) {
                        secIdx = it->second;
                    } else {
                        ElementSection sec;
                        sec.type = simType;
                        sec.name = std::string(element_type_name(simType)) + "_dim"
                                   + std::to_string(dim) + "_phys" + std::to_string(physTag);
                        secIdx = zone.sections.size();
                        zone.sections.push_back(std::move(sec));
                        secIndex.emplace(key, secIdx);
                        entityKey_to_sections[pack_entity_key(dim, physTag)].push_back(secIdx);
                    }
                    // Synthesise one entity record per (dim, physTag) so
                    // the post-loop boundary derivation code path is
                    // shared with v4.
                    const long long entKey = pack_entity_key(dim, physTag);
                    if (seenEntKeys.insert(entKey).second) {
                        Entity e;
                        e.dim = dim;
                        e.tag = physTag;
                        e.physicalTags.push_back(physTag);
                        entities.push_back(std::move(e));
                    }
                    ElementSection& sec = zone.sections[secIdx];
                    for (std::uint8_t k = 0; k < vpe; ++k) {
                        long long t = 0;
                        if (!to_int(row[3 + numTags + k], t))
                            return fail("$Elements v2 row: node tag is not an int");
                        auto it = tagToIdx.find(t);
                        if (it == tagToIdx.end())
                            return fail("$Elements v2 row: unknown node tag " + std::to_string(t));
                        sec.nodes.push_back(it->second);
                    }
                    if (dim > maxDim)
                        maxDim = static_cast<std::uint8_t>(dim);
                }
                if (!src.next(line, lineNo) || line != "$EndElements")
                    return fail("expected $EndElements");
                continue;
            }
            auto toks = split_ws(line);
            if (toks.size() < 4)
                return fail("$Elements header needs 4 ints");
            long long numBlocks = 0, numEl = 0, minTag = 0, maxTag = 0;
            if (!to_int(toks[0], numBlocks) || !to_int(toks[1], numEl) || !to_int(toks[2], minTag)
                || !to_int(toks[3], maxTag))
                return fail("$Elements header has malformed integers");
            if (numBlocks < 0 || numEl < 0)
                return fail("$Elements header has negative counts");

            for (long long b = 0; b < numBlocks; ++b) {
                if (!src.next(line, lineNo))
                    return fail("EOF inside $Elements (block header)");
                auto bhdr = split_ws(line);
                if (bhdr.size() < 4)
                    return fail("$Elements block header needs 4 ints");
                long long entityDim = 0, entityTag = 0;
                long long elType = 0, nInBlock = 0;
                if (!to_int(bhdr[0], entityDim) || !to_int(bhdr[1], entityTag)
                    || !to_int(bhdr[2], elType) || !to_int(bhdr[3], nInBlock))
                    return fail("$Elements block header malformed");
                if (nInBlock < 0)
                    return fail("$Elements block has negative count");

                if (elType == 15) {
                    // 1-node point: skip rows.
                    for (long long i = 0; i < nInBlock; ++i) {
                        if (!src.next(line, lineNo))
                            return fail("EOF inside skipped $Elements rows");
                    }
                    continue;
                }
                const ElementType simType = gmsh_element_type(static_cast<int>(elType));
                if (simType == ElementType::Unknown) {
                    return fail("unsupported Gmsh element type code " + std::to_string(elType)
                                + " (high-order / partitioned elements are "
                                  "not yet implemented in Phase 23 Pass 1)");
                }
                const std::uint8_t vpe = vertices_per_element(simType);
                if (vpe == 0)
                    return fail("internal: vpe==0");

                // Locate-or-create the section for this (entityDim, entityTag,
                // elType) triple.
                SecKey key{static_cast<int>(entityDim), entityTag, static_cast<int>(elType)};
                std::size_t secIdx;
                if (auto it = secIndex.find(key); it != secIndex.end()) {
                    secIdx = it->second;
                } else {
                    ElementSection sec;
                    sec.type = simType;
                    sec.name = std::string(element_type_name(simType)) + "_dim"
                               + std::to_string(entityDim) + "_ent" + std::to_string(entityTag);
                    secIdx = zone.sections.size();
                    zone.sections.push_back(std::move(sec));
                    secIndex.emplace(key, secIdx);
                    entityKey_to_sections[pack_entity_key(static_cast<int>(entityDim), entityTag)]
                        .push_back(secIdx);
                }
                ElementSection& sec = zone.sections[secIdx];
                sec.nodes.reserve(sec.nodes.size() + static_cast<std::size_t>(nInBlock) * vpe);

                for (long long i = 0; i < nInBlock; ++i) {
                    if (!src.next(line, lineNo))
                        return fail("EOF in $Elements rows");
                    auto erow = split_ws(line);
                    if (erow.size() < static_cast<std::size_t>(1 + vpe))
                        return fail("$Elements row has too few nodes for "
                                    + std::string(element_type_name(simType)));
                    // Skip erow[0] (element tag).
                    for (std::uint8_t k = 0; k < vpe; ++k) {
                        long long t = 0;
                        if (!to_int(erow[1 + k], t))
                            return fail("$Elements row: node tag is not an int");
                        auto it = tagToIdx.find(t);
                        if (it == tagToIdx.end())
                            return fail("$Elements row: unknown node tag " + std::to_string(t));
                        sec.nodes.push_back(it->second);
                    }
                }

                if (entityDim > maxDim)
                    maxDim = static_cast<std::uint8_t>(entityDim);
            }
            if (!src.next(line, lineNo) || line != "$EndElements")
                return fail("expected $EndElements");
            continue;
        }

        // Unknown section: skip to its $EndXxx (Gmsh ignores forward-
        // compatible sections this way).
        if (line.size() > 1 && line.front() == '$') {
            std::string endTag = "$End" + line.substr(1);
            std::size_t guard = 0;
            while (src.next(line, lineNo)) {
                if (line == endTag)
                    break;
                if (++guard > 10'000'000)
                    return fail("runaway unknown section (no matching " + endTag + ")");
            }
            continue;
        }
        return fail("unexpected token outside section: '" + line + "'");
    }

    if (zone.x.empty())
        return fail("Gmsh file contained no nodes");
    if (zone.sections.empty())
        return fail("Gmsh file contained no supported elements");

    // -- derive boundary patches from physical groups -------------------------
    //
    // Surface-dim = mesh-dim - 1.  For a 3D mesh we group 2D physical
    // groups into BoundaryPatches; for a 2D mesh we group 1D physical
    // groups.  Each patch's faceElementIndices lists the section indices
    // (within zone.sections) that belong to entities carrying that
    // physical tag.
    if (maxDim >= 2) {
        const int surfDim = static_cast<int>(maxDim) - 1;

        // physTag -> patch
        std::unordered_map<long long, std::size_t> physTag_to_patch;
        auto get_patch = [&](long long pt) -> BoundaryPatch& {
            auto it = physTag_to_patch.find(pt);
            if (it != physTag_to_patch.end())
                return zone.boundaries[it->second];
            BoundaryPatch bp;
            bp.bcType = "wall";
            // Look up name; default to "physical_<id>".
            bp.name = "physical_" + std::to_string(pt);
            for (const auto& pn : physicals) {
                if (pn.dim == surfDim && pn.tag == pt) {
                    bp.name = pn.name.empty() ? bp.name : pn.name;
                    break;
                }
            }
            const std::size_t idx = zone.boundaries.size();
            zone.boundaries.push_back(std::move(bp));
            physTag_to_patch.emplace(pt, idx);
            return zone.boundaries[idx];
        };

        for (const auto& e : entities) {
            if (e.dim != surfDim || e.physicalTags.empty())
                continue;
            auto itSec = entityKey_to_sections.find(pack_entity_key(e.dim, e.tag));
            if (itSec == entityKey_to_sections.end())
                continue;
            for (long long pt : e.physicalTags) {
                BoundaryPatch& bp = get_patch(pt);
                for (std::size_t s : itSec->second) {
                    bp.faceElementIndices.push_back(static_cast<std::uint32_t>(s));
                }
            }
        }
        // Deterministic order: by name.
        std::sort(zone.boundaries.begin(),
                  zone.boundaries.end(),
                  [](const BoundaryPatch& a, const BoundaryPatch& b) { return a.name < b.name; });
    }

    r.mesh.zones.push_back(std::move(zone));
    r.dimension = maxDim ? maxDim : static_cast<std::uint8_t>(1);
    r.ok = true;
    return r;
}

GmshReadResult read_gmsh_msh(const std::string& path)
{
    GmshReadResult r;
    std::string text = read_whole_file(path);
    if (text.empty()) {
        r.ok = false;
        r.error = "could not open or read '" + path + "'";
        return r;
    }
    return parse_gmsh_msh_string(text, path);
}

} // namespace simall::io
