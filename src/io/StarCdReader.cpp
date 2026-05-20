// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/StarCdReader.cpp
// Phase  : 23 Pass 3
//
// Line-oriented parser for the Star-CD / Star-CCM+ ProSTAR ASCII mesh
// triple. The three files share nothing but the integer node IDs they
// reference, so we parse each independently and stitch via a label->idx
// map populated from the .vrt pass.
// =============================================================================
#include "io/StarCdReader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string_view>
#include <unordered_map>

namespace simall::io {

namespace {

[[nodiscard]] std::string read_whole_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return {};
    std::ostringstream oss;
    oss << f.rdbuf();
    return oss.str();
}

[[nodiscard]] std::vector<std::string> split_lines(const std::string& src) {
    std::vector<std::string> out;
    std::string              cur;
    out.reserve(256);
    for (char c : src) {
        if (c == '\n') {
            if (!cur.empty() && cur.back() == '\r') cur.pop_back();
            out.push_back(std::move(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) {
        if (cur.back() == '\r') cur.pop_back();
        out.push_back(std::move(cur));
    }
    return out;
}

[[nodiscard]] std::vector<std::string> tokenize(const std::string& line) {
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

[[nodiscard]] bool is_blank_or_comment(const std::string& tr) {
    return tr.empty() || tr[0] == '#' || tr[0] == '!' || tr[0] == '*';
}

[[nodiscard]] std::string trim_copy(std::string s) {
    std::size_t i = 0;
    while (i < s.size()
           && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    std::size_t j = s.size();
    while (j > i
           && std::isspace(static_cast<unsigned char>(s[j - 1]))) --j;
    return s.substr(i, j - i);
}

[[nodiscard]] bool to_int(const std::string& s, long long& out) {
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

void dedup_preserve_order(const std::vector<long long>& in,
                          std::vector<long long>&       out) {
    out.clear();
    out.reserve(in.size());
    for (long long v : in) {
        bool seen = false;
        for (long long u : out) if (u == v) { seen = true; break; }
        if (!seen) out.push_back(v);
    }
}

}  // namespace

// -- public classifiers ------------------------------------------------------

ElementType starcd_classify_cell(const std::vector<long long>& nodes8,
                                  std::vector<long long>&       unique) noexcept {
    if (nodes8.size() != 8) return ElementType::Unknown;
    dedup_preserve_order(nodes8, unique);
    switch (unique.size()) {
        case 4: return ElementType::Tetra4;
        case 5: return ElementType::Pyra5;
        case 6: return ElementType::Penta6;
        case 8: return ElementType::Hexa8;
        default: return ElementType::Unknown;
    }
}

ElementType starcd_classify_face(const std::vector<long long>& nodes4,
                                  std::vector<long long>&       unique) noexcept {
    if (nodes4.size() != 4) return ElementType::Unknown;
    dedup_preserve_order(nodes4, unique);
    switch (unique.size()) {
        case 3: return ElementType::Tri3;
        case 4: return ElementType::Quad4;
        default: return ElementType::Unknown;
    }
}

// -- main parser --------------------------------------------------------------

StarCdReadResult parse_starcd_strings(const std::string& vrtText,
                                       const std::string& celText,
                                       const std::string& bndText,
                                       std::string        sourceHint) {
    StarCdReadResult r;
    r.mesh.sourceFormat = "starcd_prostar";
    r.mesh.sourcePath   = sourceHint;

    UnstructuredZone zone;
    zone.name = "starcd_zone";

    // -- .vrt --------------------------------------------------------------
    std::unordered_map<long long, NodeIdx> nodeIdx;
    {
        const auto vlines = split_lines(vrtText);
        std::size_t ln = 0;
        for (const auto& raw : vlines) {
            ++ln;
            const std::string tr = trim_copy(raw);
            if (is_blank_or_comment(tr)) continue;
            auto toks = tokenize(raw);
            if (toks.size() < 4) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".vrt:" + std::to_string(ln)
                        + "] vertex line needs 4 tokens (id x y z)";
                return r;
            }
            long long label = 0;
            double    x = 0.0, y = 0.0, z = 0.0;
            if (!to_int(toks[0], label) || !to_double(toks[1], x)
                || !to_double(toks[2], y) || !to_double(toks[3], z)) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".vrt:" + std::to_string(ln)
                        + "] malformed vertex tokens";
                return r;
            }
            const NodeIdx idx = static_cast<NodeIdx>(zone.x.size());
            zone.x.push_back(x);
            zone.y.push_back(y);
            zone.z.push_back(z);
            if (!nodeIdx.emplace(label, idx).second) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".vrt:" + std::to_string(ln)
                        + "] duplicate vertex id " + std::to_string(label);
                return r;
            }
        }
    }

    if (zone.x.empty()) {
        r.ok    = false;
        r.error = "Star-CD .vrt file contained no vertices";
        return r;
    }

    // -- .cel --------------------------------------------------------------
    // Map ElementType -> section index within zone.sections.
    std::unordered_map<int, std::size_t> volSecIdx;  // key = static_cast<int>(type)
    auto get_vol_section = [&](ElementType t) -> std::size_t {
        const int key = static_cast<int>(t);
        if (auto it = volSecIdx.find(key); it != volSecIdx.end()) return it->second;
        ElementSection sec;
        sec.type = t;
        sec.name = std::string("vol_") + element_type_name(t);
        const std::size_t i = zone.sections.size();
        zone.sections.push_back(std::move(sec));
        volSecIdx.emplace(key, i);
        return i;
    };

    std::uint8_t maxDim = 0;
    auto dim_of = [](ElementType t) -> std::uint8_t {
        switch (t) {
            case ElementType::Bar2:   return 1;
            case ElementType::Tri3:
            case ElementType::Quad4:  return 2;
            case ElementType::Tetra4:
            case ElementType::Pyra5:
            case ElementType::Penta6:
            case ElementType::Hexa8:  return 3;
            default:                  return 0;
        }
    };

    {
        const auto clines = split_lines(celText);
        std::size_t ln = 0;
        for (const auto& raw : clines) {
            ++ln;
            const std::string tr = trim_copy(raw);
            if (is_blank_or_comment(tr)) continue;
            auto toks = tokenize(raw);
            // Minimum 9 tokens: cell_id + 8 nodes.  Optional trailing
            // cell_type and region_id are accepted and ignored (the
            // connectivity is authoritative for type detection).
            if (toks.size() < 9) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".cel:" + std::to_string(ln)
                        + "] cell line needs >= 9 tokens (id + 8 nodes)";
                return r;
            }
            long long             cellId = 0;
            std::vector<long long> n8(8, 0);
            if (!to_int(toks[0], cellId)) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".cel:" + std::to_string(ln)
                        + "] cell id is not an int";
                return r;
            }
            for (int k = 0; k < 8; ++k) {
                if (!to_int(toks[1 + k], n8[k])) {
                    r.ok    = false;
                    r.error = "[" + sourceHint + ".cel:" + std::to_string(ln)
                            + "] node tag is not an int";
                    return r;
                }
            }
            std::vector<long long> unique;
            const ElementType      et = starcd_classify_cell(n8, unique);
            if (et == ElementType::Unknown) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".cel:" + std::to_string(ln)
                        + "] cannot infer element type from "
                          "distinct-node count " + std::to_string(unique.size())
                        + " (supported: 4=tet, 5=pyramid, 6=prism, 8=hex)";
                return r;
            }
            const std::size_t si    = get_vol_section(et);
            ElementSection&   sec   = zone.sections[si];
            for (long long tag : unique) {
                auto it = nodeIdx.find(tag);
                if (it == nodeIdx.end()) {
                    r.ok    = false;
                    r.error = "[" + sourceHint + ".cel:" + std::to_string(ln)
                            + "] unknown vertex id " + std::to_string(tag);
                    return r;
                }
                sec.nodes.push_back(it->second);
            }
            if (const auto d = dim_of(et); d > maxDim) maxDim = d;
        }
    }

    if (volSecIdx.empty()) {
        r.ok    = false;
        r.error = "Star-CD .cel file contained no cells";
        return r;
    }

    // -- .bnd --------------------------------------------------------------
    // One ElementSection per unique patch name; BoundaryPatch references it.
    if (!bndText.empty()) {
        std::unordered_map<std::string, std::size_t> patchSecIdx;
        std::unordered_map<std::string, std::size_t> patchBndIdx;
        auto get_patch = [&](const std::string& name, ElementType t)
            -> std::pair<std::size_t, std::size_t> {
            auto itSec = patchSecIdx.find(name);
            std::size_t si;
            if (itSec != patchSecIdx.end()) {
                si = itSec->second;
            } else {
                ElementSection sec;
                sec.type = t;
                sec.name = "patch_" + name;
                si       = zone.sections.size();
                zone.sections.push_back(std::move(sec));
                patchSecIdx.emplace(name, si);

                BoundaryPatch bp;
                bp.name   = name;
                bp.bcType = "wall";
                bp.faceElementIndices.push_back(
                    static_cast<std::uint32_t>(si));
                const std::size_t bi = zone.boundaries.size();
                zone.boundaries.push_back(std::move(bp));
                patchBndIdx.emplace(name, bi);
            }
            return {si, patchBndIdx[name]};
        };

        const auto blines = split_lines(bndText);
        std::size_t ln = 0;
        for (const auto& raw : blines) {
            ++ln;
            const std::string tr = trim_copy(raw);
            if (is_blank_or_comment(tr)) continue;
            auto toks = tokenize(raw);
            if (toks.size() < 6) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".bnd:" + std::to_string(ln)
                        + "] boundary line needs >= 6 tokens (id + 4 "
                          "nodes + region_id)";
                return r;
            }
            long long bId = 0, regionId = 0;
            std::vector<long long> n4(4, 0);
            if (!to_int(toks[0], bId)) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".bnd:" + std::to_string(ln)
                        + "] boundary id is not an int";
                return r;
            }
            for (int k = 0; k < 4; ++k) {
                if (!to_int(toks[1 + k], n4[k])) {
                    r.ok    = false;
                    r.error = "[" + sourceHint + ".bnd:" + std::to_string(ln)
                            + "] node tag is not an int";
                    return r;
                }
            }
            if (!to_int(toks[5], regionId)) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".bnd:" + std::to_string(ln)
                        + "] region_id is not an int";
                return r;
            }
            // Optional trailing tokens: type_name (toks[6]) and
            // patch_name (toks[7..]). The patch identity is the
            // patch_name when present, otherwise "region_<region_id>".
            std::string patchName;
            if (toks.size() >= 8) {
                patchName = toks[7];
                for (std::size_t k = 8; k < toks.size(); ++k) {
                    patchName.push_back('_');
                    patchName.append(toks[k]);
                }
            } else {
                patchName = "region_" + std::to_string(regionId);
            }
            std::vector<long long> unique;
            const ElementType      ft = starcd_classify_face(n4, unique);
            if (ft == ElementType::Unknown) {
                r.ok    = false;
                r.error = "[" + sourceHint + ".bnd:" + std::to_string(ln)
                        + "] cannot infer face type from distinct-node "
                          "count " + std::to_string(unique.size())
                        + " (supported: 3=tri, 4=quad)";
                return r;
            }
            auto [si, bi]      = get_patch(patchName, ft);
            ElementSection& sec = zone.sections[si];
            if (sec.type != ft) {
                // Mixed-type patches not supported in this pass: report.
                r.ok    = false;
                r.error = "[" + sourceHint + ".bnd:" + std::to_string(ln)
                        + "] patch '" + patchName + "' mixes "
                        + element_type_name(sec.type) + " and "
                        + element_type_name(ft) + " faces (Phase 23 Pass 3 "
                          "requires single-type patches; split the patch "
                          "in your pre-processor)";
                return r;
            }
            for (long long tag : unique) {
                auto it = nodeIdx.find(tag);
                if (it == nodeIdx.end()) {
                    r.ok    = false;
                    r.error = "[" + sourceHint + ".bnd:" + std::to_string(ln)
                            + "] unknown vertex id " + std::to_string(tag);
                    return r;
                }
                sec.nodes.push_back(it->second);
            }
            (void)bi;
        }

        std::sort(zone.boundaries.begin(), zone.boundaries.end(),
                  [](const BoundaryPatch& a, const BoundaryPatch& b) {
                      return a.name < b.name;
                  });
    }

    r.mesh.zones.push_back(std::move(zone));
    r.dimension = maxDim ? maxDim : static_cast<std::uint8_t>(1);
    r.ok        = true;
    return r;
}

StarCdReadResult read_starcd(const std::string& vrtPath,
                              const std::string& celPath,
                              const std::string& bndPath) {
    StarCdReadResult r;
    const std::string vrt = read_whole_file(vrtPath);
    const std::string cel = read_whole_file(celPath);
    if (vrt.empty()) {
        r.ok    = false;
        r.error = "could not open .vrt file '" + vrtPath + "'";
        return r;
    }
    if (cel.empty()) {
        r.ok    = false;
        r.error = "could not open .cel file '" + celPath + "'";
        return r;
    }
    std::string bnd;
    if (!bndPath.empty()) {
        bnd = read_whole_file(bndPath);
        if (bnd.empty()) {
            r.ok    = false;
            r.error = "could not open .bnd file '" + bndPath + "'";
            return r;
        }
    }
    return parse_starcd_strings(vrt, cel, bnd, vrtPath);
}

}  // namespace simall::io
