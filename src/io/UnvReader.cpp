// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/UnvReader.cpp
// Phase  : 23 Pass 2
//
// Two-pass parser for the SDRC/I-DEAS Universal mesh format:
//   Pass A : split the source into datasets (delimited by "-1" sentinels).
//   Pass B : per known dataset number, walk the body as a stream of
//            whitespace-separated tokens and emit ImportedMesh state.
//
// Token-stream parsing (rather than column-precise FORTRAN parsing) is
// robust against the small whitespace differences between UNV producers
// (HyperMesh pads to 13-char fields, IcemCFD to 10-char, gmsh to 5-char).
// The reader is tolerant of an arbitrary number of leading-zero pads and
// blank lines between records as long as token counts match the spec.
// =============================================================================
#include "io/UnvReader.hpp"

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

// -- public element-descriptor translator ------------------------------------

ElementType unv_element_type(int feDescriptor) noexcept
{
    switch (feDescriptor) {
    case 11:
    case 21:
    case 22:
        return ElementType::Bar2;
    case 41:
    case 91:
        return ElementType::Tri3;
    case 44:
    case 94:
        return ElementType::Quad4;
    case 111:
        return ElementType::Tetra4;
    case 112:
        return ElementType::Penta6;
    case 115:
        return ElementType::Hexa8;
    default:
        return ElementType::Unknown;
    }
}

// -- public Fortran-D float normaliser ---------------------------------------

std::string unv_normalize_float_token(std::string s)
{
    for (char& c : s) {
        if (c == 'D' || c == 'd')
            c = 'e';
    }
    return s;
}

namespace
{

[[nodiscard]] std::string read_whole_file(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return {};
    std::ostringstream oss;
    oss << f.rdbuf();
    return oss.str();
}

[[nodiscard]] std::vector<std::string> split_lines(const std::string& src)
{
    std::vector<std::string> out;
    std::string cur;
    out.reserve(256);
    for (char c : src) {
        if (c == '\n') {
            if (!cur.empty() && cur.back() == '\r')
                cur.pop_back();
            out.push_back(std::move(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) {
        if (cur.back() == '\r')
            cur.pop_back();
        out.push_back(std::move(cur));
    }
    return out;
}

[[nodiscard]] std::string trim_copy(std::string s)
{
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i])))
        ++i;
    std::size_t j = s.size();
    while (j > i && std::isspace(static_cast<unsigned char>(s[j - 1])))
        --j;
    return s.substr(i, j - i);
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
    const std::string norm = unv_normalize_float_token(s);
    char* end = nullptr;
    out = std::strtod(norm.c_str(), &end);
    return end != norm.c_str() && (end == nullptr || *end == '\0');
}

// One UNV dataset (header line + body lines collected between -1 sentinels).
struct Dataset
{
    int number = 0;
    std::size_t startLine = 0; ///< 1-based source line of header
    std::vector<std::string> bodyLines;
};

// Split source into dataset blocks. Returns false on structural error.
[[nodiscard]] bool split_datasets(const std::vector<std::string>& lines,
                                  std::vector<Dataset>& outDatasets,
                                  std::string& err)
{
    enum class State
    {
        Outside,
        ExpectNumber,
        InBody
    };
    State st = State::Outside;
    Dataset cur;
    std::size_t lineNo = 0;
    for (const auto& raw : lines) {
        ++lineNo;
        const std::string tr = trim_copy(raw);
        if (st == State::Outside) {
            if (tr.empty())
                continue;
            if (tr == "-1") {
                st = State::ExpectNumber;
                continue;
            }
            err = "[unv:" + std::to_string(lineNo) + "] expected '-1' sentinel, got '" + tr + "'";
            return false;
        }
        if (st == State::ExpectNumber) {
            if (tr.empty())
                continue;
            long long n = 0;
            if (!to_int(tr, n)) {
                err = "[unv:" + std::to_string(lineNo) + "] expected dataset number, got '" + tr
                      + "'";
                return false;
            }
            cur.number = static_cast<int>(n);
            cur.startLine = lineNo;
            cur.bodyLines.clear();
            st = State::InBody;
            continue;
        }
        // State::InBody
        if (tr == "-1") {
            outDatasets.push_back(std::move(cur));
            cur = Dataset{};
            st = State::Outside;
            continue;
        }
        cur.bodyLines.push_back(raw);
    }
    if (st != State::Outside) {
        err = "[unv:" + std::to_string(lineNo) + "] unterminated dataset (missing trailing '-1')";
        return false;
    }
    return true;
}

// Token stream over a vector<string> body. Whitespace-separated, comments
// are not part of the UNV grammar so we just walk character-by-character.
class TokenCursor
{
public:
    explicit TokenCursor(const std::vector<std::string>& body)
    {
        for (const auto& l : body) {
            tokenize_line(l);
        }
    }
    bool eof() const noexcept { return pos_ >= toks_.size(); }
    std::size_t remaining() const noexcept
    {
        return pos_ >= toks_.size() ? 0 : toks_.size() - pos_;
    }
    const std::string& next() { return toks_[pos_++]; }

private:
    void tokenize_line(const std::string& l)
    {
        std::string cur;
        for (char c : l) {
            if (std::isspace(static_cast<unsigned char>(c))) {
                if (!cur.empty()) {
                    toks_.push_back(std::move(cur));
                    cur.clear();
                }
            } else {
                cur.push_back(c);
            }
        }
        if (!cur.empty())
            toks_.push_back(std::move(cur));
    }
    std::vector<std::string> toks_;
    std::size_t pos_ = 0;
};

} // namespace

// -- main parser --------------------------------------------------------------

UnvReadResult parse_unv_string(const std::string& text, std::string sourceHint)
{
    UnvReadResult r;
    r.mesh.sourceFormat = "ideas_unv";
    r.mesh.sourcePath = sourceHint;

    const auto lines = split_lines(text);

    std::vector<Dataset> datasets;
    {
        std::string splitErr;
        if (!split_datasets(lines, datasets, splitErr)) {
            r.ok = false;
            r.error = splitErr;
            return r;
        }
    }

    // -- working state --------------------------------------------------------
    UnstructuredZone zone;
    zone.name = "unv_zone";

    // UNV node labels can be sparse; map label -> 0-based idx into zone.x/y/z.
    std::unordered_map<long long, NodeIdx> nodeIdx;

    // Map element label -> (section index in zone.sections).
    std::unordered_map<long long, std::size_t> elemSec;

    // Map FE-descriptor key -> section index (one section per descriptor).
    std::unordered_map<int, std::size_t> descSec;

    // Permanent groups collected from datasets 2467 / 2477.
    struct GroupRow
    {
        std::string name;
        std::vector<long long> elementLabels;
    };
    std::vector<GroupRow> groups;

    std::uint8_t maxDim = 0;

    auto fail = [&](const Dataset& ds, const std::string& msg) {
        r.ok = false;
        r.error = "[" + sourceHint + ":" + std::to_string(ds.startLine) + " dataset "
                  + std::to_string(ds.number) + "] " + msg;
        return r;
    };

    auto dim_of = [](ElementType t) -> std::uint8_t {
        switch (t) {
        case ElementType::Bar2:
            return 1;
        case ElementType::Tri3:
        case ElementType::Quad4:
            return 2;
        case ElementType::Tetra4:
        case ElementType::Pyra5:
        case ElementType::Penta6:
        case ElementType::Hexa8:
            return 3;
        default:
            return 0;
        }
    };

    // -- per-dataset dispatch -------------------------------------------------
    for (const auto& ds : datasets) {
        switch (ds.number) {
        case 2411: { // Nodes
            TokenCursor tc(ds.bodyLines);
            while (!tc.eof()) {
                if (tc.remaining() < 7) {
                    return fail(ds, "truncated node record");
                }
                long long label = 0;
                long long expCs = 0, dispCs = 0, color = 0;
                if (!to_int(tc.next(), label) || !to_int(tc.next(), expCs)
                    || !to_int(tc.next(), dispCs) || !to_int(tc.next(), color)) {
                    return fail(ds, "node header is not 4 ints");
                }
                double x = 0.0, y = 0.0, z = 0.0;
                if (!to_double(tc.next(), x) || !to_double(tc.next(), y)
                    || !to_double(tc.next(), z)) {
                    return fail(ds, "node coords are not 3 floats");
                }
                const NodeIdx idx = static_cast<NodeIdx>(zone.x.size());
                zone.x.push_back(x);
                zone.y.push_back(y);
                zone.z.push_back(z);
                nodeIdx.emplace(label, idx);
            }
            break;
        }

        case 2412: { // Elements
            TokenCursor tc(ds.bodyLines);
            while (!tc.eof()) {
                if (tc.remaining() < 6) {
                    return fail(ds, "truncated element header");
                }
                long long label = 0, fe = 0, pp = 0, mp = 0, col = 0, n = 0;
                if (!to_int(tc.next(), label) || !to_int(tc.next(), fe) || !to_int(tc.next(), pp)
                    || !to_int(tc.next(), mp) || !to_int(tc.next(), col) || !to_int(tc.next(), n)) {
                    return fail(ds, "element header is not 6 ints");
                }
                const int feCode = static_cast<int>(fe);
                const auto et = unv_element_type(feCode);
                if (et == ElementType::Unknown) {
                    return fail(ds, "unsupported FE descriptor " + std::to_string(feCode));
                }
                const std::uint8_t vpe = vertices_per_element(et);
                if (vpe == 0 || vpe != n) {
                    return fail(ds,
                                "FE descriptor " + std::to_string(feCode) + " declares num_nodes="
                                    + std::to_string(n) + " but expected " + std::to_string(vpe));
                }
                // Beam descriptors carry a 3-int orientation record.
                if (feCode == 11 || feCode == 21 || feCode == 22) {
                    if (tc.remaining() < 3) {
                        return fail(ds, "beam record missing orientation");
                    }
                    long long o1 = 0, o2 = 0, o3 = 0;
                    if (!to_int(tc.next(), o1) || !to_int(tc.next(), o2)
                        || !to_int(tc.next(), o3)) {
                        return fail(ds, "beam orientation is not 3 ints");
                    }
                }
                if (tc.remaining() < vpe) {
                    return fail(ds, "element connectivity truncated");
                }
                // Locate-or-create the per-descriptor section.
                std::size_t secIdx;
                if (auto it = descSec.find(feCode); it != descSec.end()) {
                    secIdx = it->second;
                } else {
                    ElementSection sec;
                    sec.type = et;
                    sec.name = std::string(element_type_name(et)) + "_fe" + std::to_string(feCode);
                    secIdx = zone.sections.size();
                    zone.sections.push_back(std::move(sec));
                    descSec.emplace(feCode, secIdx);
                }
                ElementSection& sec = zone.sections[secIdx];
                for (std::uint8_t k = 0; k < vpe; ++k) {
                    long long tag = 0;
                    if (!to_int(tc.next(), tag)) {
                        return fail(ds, "node label is not an int");
                    }
                    auto it = nodeIdx.find(tag);
                    if (it == nodeIdx.end()) {
                        return fail(ds, "unknown node label " + std::to_string(tag));
                    }
                    sec.nodes.push_back(it->second);
                }
                elemSec.emplace(label, secIdx);
                if (const auto d = dim_of(et); d > maxDim)
                    maxDim = d;
            }
            break;
        }

        case 2467:
        case 2477: { // Permanent Groups
            // Header record: 8 ints. Last is n_entities.
            // Name record: free-form (up to 40 chars), single line.
            // Entity records: 4 ints per entity, packed 2-per-line.
            // We iterate until the body is exhausted (multiple groups
            // can share one dataset block).
            std::size_t li = 0;
            while (li < ds.bodyLines.size()) {
                // Skip blank lines between groups.
                while (li < ds.bodyLines.size() && trim_copy(ds.bodyLines[li]).empty())
                    ++li;
                if (li >= ds.bodyLines.size())
                    break;

                // Header (8 ints, possibly spanning whitespace).
                TokenCursor headerTc({ds.bodyLines[li++]});
                if (headerTc.remaining() < 8) {
                    return fail(ds, "group header is not 8 ints");
                }
                long long hdr[8] = {0, 0, 0, 0, 0, 0, 0, 0};
                for (int k = 0; k < 8; ++k) {
                    if (!to_int(headerTc.next(), hdr[k])) {
                        return fail(ds, "group header has a non-int field");
                    }
                }
                const long long nEntities = hdr[7];
                if (nEntities < 0) {
                    return fail(ds, "group has negative entity count");
                }
                if (li >= ds.bodyLines.size()) {
                    return fail(ds, "group missing name line");
                }
                const std::string groupName = trim_copy(ds.bodyLines[li++]);
                GroupRow gr;
                gr.name = groupName.empty() ? "unnamed_group" : groupName;

                // Read ceil(nEntities/2) lines.
                const long long nLines = (nEntities + 1) / 2;
                long long remaining = nEntities;
                for (long long row = 0; row < nLines; ++row) {
                    if (li >= ds.bodyLines.size()) {
                        return fail(ds, "group entity rows truncated");
                    }
                    TokenCursor rowTc({ds.bodyLines[li++]});
                    const long long perRow = std::min<long long>(2, remaining);
                    for (long long e = 0; e < perRow; ++e) {
                        if (rowTc.remaining() < 4) {
                            return fail(ds, "group entity row missing fields");
                        }
                        long long entType = 0, entTag = 0, leaf = 0, comp = 0;
                        if (!to_int(rowTc.next(), entType) || !to_int(rowTc.next(), entTag)
                            || !to_int(rowTc.next(), leaf) || !to_int(rowTc.next(), comp)) {
                            return fail(ds, "group entity row has non-int");
                        }
                        // entType == 8 is "FE element". Other types
                        // (7 node, 11 dof set, ...) are dropped from
                        // this pass's boundary derivation.
                        if (entType == 8) {
                            gr.elementLabels.push_back(entTag);
                        }
                    }
                    remaining -= perRow;
                }
                groups.push_back(std::move(gr));
            }
            break;
        }

        default:
            // Forward-compat: silently skip unknown datasets
            // (151 header, 164 units, 2420 coord sys, etc.).
            break;
        }
        if (!r.ok && !r.error.empty())
            return r; // propagate early fail
    }

    if (zone.x.empty()) {
        r.ok = false;
        r.error = "UNV file contained no nodes (dataset 2411 missing)";
        return r;
    }
    if (zone.sections.empty()) {
        r.ok = false;
        r.error = "UNV file contained no supported elements (dataset 2412 missing)";
        return r;
    }

    // -- derive boundary patches from element-only groups at surfDim ----------
    if (maxDim >= 2) {
        const std::uint8_t surfDim = static_cast<std::uint8_t>(maxDim - 1);
        for (const auto& g : groups) {
            if (g.elementLabels.empty())
                continue;
            std::unordered_set<std::size_t> sectionsHit;
            bool allSurface = true;
            for (long long lab : g.elementLabels) {
                auto it = elemSec.find(lab);
                if (it == elemSec.end()) {
                    allSurface = false;
                    break;
                }
                const auto et = zone.sections[it->second].type;
                if (dim_of(et) != surfDim) {
                    allSurface = false;
                    break;
                }
                sectionsHit.insert(it->second);
            }
            if (!allSurface || sectionsHit.empty())
                continue;
            BoundaryPatch bp;
            bp.name = g.name;
            bp.bcType = "wall";
            std::vector<std::size_t> sorted(sectionsHit.begin(), sectionsHit.end());
            std::sort(sorted.begin(), sorted.end());
            for (std::size_t s : sorted) {
                bp.faceElementIndices.push_back(static_cast<std::uint32_t>(s));
            }
            zone.boundaries.push_back(std::move(bp));
        }
        std::sort(zone.boundaries.begin(),
                  zone.boundaries.end(),
                  [](const BoundaryPatch& a, const BoundaryPatch& b) { return a.name < b.name; });
    }

    r.mesh.zones.push_back(std::move(zone));
    r.dimension = maxDim ? maxDim : static_cast<std::uint8_t>(1);
    r.ok = true;
    return r;
}

UnvReadResult read_unv(const std::string& path)
{
    UnvReadResult r;
    std::string text = read_whole_file(path);
    if (text.empty()) {
        r.ok = false;
        r.error = "could not open or read '" + path + "'";
        return r;
    }
    return parse_unv_string(text, path);
}

} // namespace simall::io
