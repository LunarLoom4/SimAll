// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/BlockMeshDict.cpp
// Phase  : 23 Pass 14 - blockMeshDict parser + MultiblockHex builder.
// =============================================================================
#include "meshing/BlockMeshDict.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace simall::meshing
{

namespace
{

// ---------------------------------------------------------------------------
// Tokeniser
// ---------------------------------------------------------------------------
struct Tokeniser
{
    std::string_view src;
    std::size_t pos = 0;
    std::size_t line = 1;

    void skip_ws_and_comments()
    {
        while (pos < src.size()) {
            const char c = src[pos];
            if (c == '\n') {
                ++line;
                ++pos;
                continue;
            }
            if (std::isspace(static_cast<unsigned char>(c))) {
                ++pos;
                continue;
            }
            if (c == '/' && pos + 1 < src.size() && src[pos + 1] == '/') {
                while (pos < src.size() && src[pos] != '\n')
                    ++pos;
                continue;
            }
            if (c == '/' && pos + 1 < src.size() && src[pos + 1] == '*') {
                pos += 2;
                while (pos + 1 < src.size() && !(src[pos] == '*' && src[pos + 1] == '/')) {
                    if (src[pos] == '\n')
                        ++line;
                    ++pos;
                }
                if (pos + 1 < src.size())
                    pos += 2;
                continue;
            }
            return;
        }
    }

    bool eof()
    {
        skip_ws_and_comments();
        return pos >= src.size();
    }

    char peek()
    {
        skip_ws_and_comments();
        return pos < src.size() ? src[pos] : '\0';
    }

    bool match(char c)
    {
        skip_ws_and_comments();
        if (pos < src.size() && src[pos] == c) {
            ++pos;
            return true;
        }
        return false;
    }

    void expect(char c)
    {
        if (!match(c)) {
            std::ostringstream os;
            os << "blockMeshDict line " << line << ": expected '" << c << "'";
            if (pos < src.size())
                os << " near '" << src[pos] << "'";
            throw std::runtime_error(os.str());
        }
    }

    // Word: alpha/digit/_ ; used for identifiers AND keywords like "hex".
    std::string word()
    {
        skip_ws_and_comments();
        const std::size_t start = pos;
        while (pos < src.size()) {
            const char c = src[pos];
            if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.')
                ++pos;
            else
                break;
        }
        if (start == pos) {
            std::ostringstream os;
            os << "blockMeshDict line " << line << ": expected identifier";
            throw std::runtime_error(os.str());
        }
        return std::string(src.substr(start, pos - start));
    }

    double number()
    {
        skip_ws_and_comments();
        const std::size_t start = pos;
        if (pos < src.size() && (src[pos] == '+' || src[pos] == '-'))
            ++pos;
        while (pos < src.size()) {
            const char c = src[pos];
            if (std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == 'e' || c == 'E'
                || c == '+' || c == '-')
                ++pos;
            else
                break;
        }
        if (start == pos) {
            std::ostringstream os;
            os << "blockMeshDict line " << line << ": expected number";
            throw std::runtime_error(os.str());
        }
        return std::strtod(std::string(src.substr(start, pos - start)).c_str(), nullptr);
    }

    std::uint32_t uint()
    {
        const double v = number();
        if (v < 0.0) {
            std::ostringstream os;
            os << "blockMeshDict line " << line << ": expected non-negative integer";
            throw std::runtime_error(os.str());
        }
        return static_cast<std::uint32_t>(v + 0.5);
    }
};

// ---------------------------------------------------------------------------
// Section parsers
// ---------------------------------------------------------------------------
void parse_vertices(Tokeniser& t, std::vector<util::Vec3d>& out)
{
    t.expect('(');
    while (t.peek() != ')') {
        t.expect('(');
        const double x = t.number();
        const double y = t.number();
        const double z = t.number();
        t.expect(')');
        out.push_back({x, y, z});
    }
    t.expect(')');
    if (t.peek() == ';')
        t.expect(';');
}

void parse_blocks(Tokeniser& t, std::vector<BlockMeshBlock>& out)
{
    t.expect('(');
    while (t.peek() != ')') {
        const std::string kw = t.word();
        if (kw != "hex") {
            std::ostringstream os;
            os << "blockMeshDict line " << t.line << ": unsupported block type '" << kw
               << "' (only 'hex' supported)";
            throw std::runtime_error(os.str());
        }
        BlockMeshBlock b;
        t.expect('(');
        for (int i = 0; i < 8; ++i)
            b.vertices[i] = t.uint();
        t.expect(')');

        // Optional block name between vertex list and division list.
        if (t.peek() != '(') {
            b.name = t.word();
        }

        t.expect('(');
        for (int i = 0; i < 3; ++i)
            b.divisions[i] = t.uint();
        t.expect(')');

        // Optional grading clause.  We accept simpleGrading (g g g);
        // edgeGrading (.. 12 values ..) is rejected.
        if (t.peek() != ')' && t.peek() != 'h') {
            // Peek a word to identify the grading style.
            const std::size_t savePos = t.pos;
            const std::size_t saveLine = t.line;
            t.skip_ws_and_comments();
            if (t.pos < t.src.size() && std::isalpha(static_cast<unsigned char>(t.src[t.pos]))) {
                const std::string g = t.word();
                if (g == "simpleGrading") {
                    t.expect('(');
                    for (int i = 0; i < 3; ++i)
                        b.grading[i] = t.number();
                    t.expect(')');
                } else if (g == "edgeGrading") {
                    std::ostringstream os;
                    os << "blockMeshDict line " << t.line
                       << ": edgeGrading not supported in this pass";
                    throw std::runtime_error(os.str());
                } else {
                    // Not a grading keyword -- roll back.
                    t.pos = savePos;
                    t.line = saveLine;
                }
            }
        }
        out.push_back(std::move(b));
    }
    t.expect(')');
    if (t.peek() == ';')
        t.expect(';');
}

void parse_patch_body(Tokeniser& t, BlockMeshPatch& p)
{
    t.expect('{');
    while (t.peek() != '}') {
        const std::string key = t.word();
        if (key == "type") {
            p.type = t.word();
            t.expect(';');
        } else if (key == "faces") {
            t.expect('(');
            while (t.peek() != ')') {
                t.expect('(');
                std::array<std::uint32_t, 4> face{};
                for (int i = 0; i < 4; ++i)
                    face[i] = t.uint();
                t.expect(')');
                p.faces.push_back(face);
            }
            t.expect(')');
            if (t.peek() == ';')
                t.expect(';');
        } else {
            // Skip unknown key  value(s) up to ';'.
            while (!t.eof() && t.peek() != ';' && t.peek() != '}') {
                if (t.peek() == '(') {
                    int depth = 0;
                    do {
                        const char c = t.peek();
                        if (c == '(')
                            ++depth;
                        else if (c == ')')
                            --depth;
                        ++t.pos;
                    } while (depth > 0 && t.pos < t.src.size());
                } else {
                    ++t.pos;
                }
            }
            if (t.peek() == ';')
                t.expect(';');
        }
    }
    t.expect('}');
}

void parse_boundary(Tokeniser& t, std::vector<BlockMeshPatch>& out)
{
    t.expect('(');
    while (t.peek() != ')') {
        BlockMeshPatch p;
        p.name = t.word();
        parse_patch_body(t, p);
        out.push_back(std::move(p));
    }
    t.expect(')');
    if (t.peek() == ';')
        t.expect(';');
}

} // namespace

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------
BlockMeshDict parse_block_mesh_dict(std::string_view text)
{
    Tokeniser t{text};
    BlockMeshDict d;
    while (!t.eof()) {
        const std::string key = t.word();
        if (key == "convertToMeters") {
            d.convertToMeters = t.number();
            t.expect(';');
        } else if (key == "vertices") {
            parse_vertices(t, d.vertices);
        } else if (key == "blocks") {
            parse_blocks(t, d.blocks);
        } else if (key == "boundary") {
            parse_boundary(t, d.patches);
        } else if (key == "edges" || key == "mergePatchPairs" || key == "FoamFile"
                   || key == "defaultPatch") {
            // Skip the trailing dictionary or list body.
            t.skip_ws_and_comments();
            if (t.peek() == '{') {
                int depth = 0;
                do {
                    const char c = t.peek();
                    if (c == '{')
                        ++depth;
                    else if (c == '}')
                        --depth;
                    ++t.pos;
                } while (depth > 0 && t.pos < t.src.size());
            } else if (t.peek() == '(') {
                int depth = 0;
                do {
                    const char c = t.peek();
                    if (c == '(')
                        ++depth;
                    else if (c == ')')
                        --depth;
                    ++t.pos;
                } while (depth > 0 && t.pos < t.src.size());
            }
            if (t.peek() == ';')
                t.expect(';');
        } else {
            std::ostringstream os;
            os << "blockMeshDict line " << t.line << ": unknown top-level key '" << key << "'";
            throw std::runtime_error(os.str());
        }
    }
    return d;
}

BlockMeshDict load_block_mesh_dict(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("blockMeshDict: cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    return parse_block_mesh_dict(ss.str());
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------
namespace
{

void emit_double(std::ostream& os, double v)
{
    // Round-trip-friendly: trim trailing zeros so 1.0 prints as "1" but keep
    // full precision for non-integer values.
    if (v == static_cast<long long>(v)) {
        os << static_cast<long long>(v);
    } else {
        std::ostringstream tmp;
        tmp.precision(17);
        tmp << v;
        os << tmp.str();
    }
}

} // namespace

std::string write_block_mesh_dict(const BlockMeshDict& dict)
{
    std::ostringstream os;
    os << "// blockMeshDict (generated by SimAll Beta)\n\n";
    os << "convertToMeters ";
    emit_double(os, dict.convertToMeters);
    os << ";\n\n";

    os << "vertices\n(\n";
    for (std::size_t i = 0; i < dict.vertices.size(); ++i) {
        const auto& v = dict.vertices[i];
        os << "    (";
        emit_double(os, v.x);
        os << ' ';
        emit_double(os, v.y);
        os << ' ';
        emit_double(os, v.z);
        os << ")    // " << i << '\n';
    }
    os << ");\n\n";

    os << "blocks\n(\n";
    for (const auto& b : dict.blocks) {
        os << "    hex (";
        for (int i = 0; i < 8; ++i) {
            if (i)
                os << ' ';
            os << b.vertices[i];
        }
        os << ')';
        if (!b.name.empty())
            os << ' ' << b.name;
        os << " (" << b.divisions[0] << ' ' << b.divisions[1] << ' ' << b.divisions[2]
           << ") simpleGrading (";
        emit_double(os, b.grading[0]);
        os << ' ';
        emit_double(os, b.grading[1]);
        os << ' ';
        emit_double(os, b.grading[2]);
        os << ")\n";
    }
    os << ");\n\n";

    if (!dict.patches.empty()) {
        os << "boundary\n(\n";
        for (const auto& p : dict.patches) {
            os << "    " << p.name << "\n    {\n";
            os << "        type " << (p.type.empty() ? "patch" : p.type) << ";\n";
            os << "        faces\n        (\n";
            for (const auto& f : p.faces) {
                os << "            (" << f[0] << ' ' << f[1] << ' ' << f[2] << ' ' << f[3] << ")\n";
            }
            os << "        );\n    }\n";
        }
        os << ");\n";
    }
    return os.str();
}

void save_block_mesh_dict(const BlockMeshDict& dict, const std::string& path)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("blockMeshDict: cannot write " + path);
    out << write_block_mesh_dict(dict);
}

// ---------------------------------------------------------------------------
// Multiblock builder
// ---------------------------------------------------------------------------
namespace
{

// HexBlock 6-face vertex indexing (matches MultiblockHex face zone convention).
//   0: imin  (-x)   1: imax  (+x)
//   2: jmin  (-y)   3: jmax  (+y)
//   4: kmin  (-z)   5: kmax  (+z)
constexpr int kHexFaceCorners[6][4] = {
    {0, 4, 7, 3}, // imin
    {1, 2, 6, 5}, // imax
    {0, 1, 5, 4}, // jmin
    {3, 7, 6, 2}, // jmax
    {0, 3, 2, 1}, // kmin
    {4, 5, 6, 7}  // kmax
};

// Sort a 4-tuple in ascending order and pack into a 64-bit key.  Vertex
// indices are 16-bit (max 65535 vertices -- adequate for any sane
// blockMeshDict).
std::uint64_t face_key(std::uint32_t a, std::uint32_t b, std::uint32_t c, std::uint32_t d) noexcept
{
    std::array<std::uint32_t, 4> v{a, b, c, d};
    std::sort(v.begin(), v.end());
    return (static_cast<std::uint64_t>(v[0])) | (static_cast<std::uint64_t>(v[1]) << 16)
           | (static_cast<std::uint64_t>(v[2]) << 32) | (static_cast<std::uint64_t>(v[3]) << 48);
}

} // namespace

BlockMeshBuildStats build_multiblock(const BlockMeshDict& dict, MultiblockHex& out)
{
    BlockMeshBuildStats stats;

    // Build face-key index: (sorted-4-vertex key) -> (blockIndex, faceSlot).
    struct BlockFaceLoc
    {
        std::uint32_t block;
        std::uint32_t face;
    };
    std::unordered_map<std::uint64_t, BlockFaceLoc> faceIndex;
    faceIndex.reserve(dict.blocks.size() * 6 * 2);

    std::vector<HexBlock> built;
    built.reserve(dict.blocks.size());
    for (std::size_t bi = 0; bi < dict.blocks.size(); ++bi) {
        const BlockMeshBlock& b = dict.blocks[bi];
        HexBlock hb;
        hb.name = b.name;
        for (int c = 0; c < 8; ++c) {
            const std::uint32_t vid = b.vertices[c];
            if (vid >= dict.vertices.size()) {
                std::ostringstream os;
                os << "blockMeshDict: block " << bi << " corner " << c << " references vertex "
                   << vid << " but only " << dict.vertices.size() << " defined";
                throw std::runtime_error(os.str());
            }
            const util::Vec3d& v = dict.vertices[vid];
            hb.corners[c] = {
                v.x * dict.convertToMeters, v.y * dict.convertToMeters, v.z * dict.convertToMeters};
        }
        hb.divisions = b.divisions;
        hb.grading = b.grading;
        hb.faceZones = {0, 0, 0, 0, 0, 0};

        for (int f = 0; f < 6; ++f) {
            const auto& fc = kHexFaceCorners[f];
            const std::uint64_t key = face_key(
                b.vertices[fc[0]], b.vertices[fc[1]], b.vertices[fc[2]], b.vertices[fc[3]]);
            faceIndex[key] =
                BlockFaceLoc{static_cast<std::uint32_t>(bi), static_cast<std::uint32_t>(f)};
        }
        built.push_back(std::move(hb));
    }

    // Assign patch zones.  Zone id = 1 + patchIndex.
    for (std::size_t pi = 0; pi < dict.patches.size(); ++pi) {
        const BlockMeshPatch& p = dict.patches[pi];
        const std::uint32_t zone = static_cast<std::uint32_t>(pi + 1);
        for (const auto& f : p.faces) {
            const std::uint64_t key = face_key(f[0], f[1], f[2], f[3]);
            const auto it = faceIndex.find(key);
            if (it == faceIndex.end()) {
                ++stats.orphanFaces;
                continue;
            }
            built[it->second.block].faceZones[it->second.face] = zone;
            ++stats.patchFaces;
        }
    }

    for (auto& hb : built)
        out.add_block(std::move(hb));
    stats.blocks = dict.blocks.size();
    stats.patches = dict.patches.size();
    return stats;
}

} // namespace simall::meshing
