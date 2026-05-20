// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/FaceBasedMeshReader.cpp
//
// Face-based polyhedral mesh directory ASCII reader.  Tokenizes each of
// the five files once into (kind,text) tokens with `(`, `)`, `{`, `}`,
// `;` as single-char delimiters, then walks the stream with per-file
// parsers.  Cell classification follows the face-signature rule
// documented in the header.
// =============================================================================
#include "io/FaceBasedMeshReader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string_view>
#include <unordered_map>

namespace simall::io
{

ElementType classify_face_based_cell(std::uint32_t nTri, std::uint32_t nQuad) noexcept
{
    if (nTri == 4 && nQuad == 0)
        return ElementType::Tetra4;
    if (nTri == 0 && nQuad == 6)
        return ElementType::Hexa8;
    if (nTri == 4 && nQuad == 1)
        return ElementType::Pyra5;
    if (nTri == 2 && nQuad == 3)
        return ElementType::Penta6;
    return ElementType::Unknown;
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

// -- tokenizer ---------------------------------------------------------------

enum class TokKind
{
    Word,
    LParen,
    RParen,
    LBrace,
    RBrace,
    Semi,
    EndOfFile
};

struct Tok
{
    TokKind kind = TokKind::EndOfFile;
    std::string text;
    std::size_t line = 0;
};

class MeshTokenLex
{
public:
    explicit MeshTokenLex(const std::string& src) { tokenize(src); }
    const Tok& peek() const { return at(pos_); }
    const Tok& next() { return at(pos_++); }
    bool eof() const { return pos_ >= toks_.size(); }
    std::size_t pos() const { return pos_; }
    void set_pos(std::size_t p) { pos_ = p; }

private:
    static const Tok& eofTok()
    {
        static const Tok t{TokKind::EndOfFile, {}, 0};
        return t;
    }
    const Tok& at(std::size_t p) const { return p < toks_.size() ? toks_[p] : eofTok(); }
    void tokenize(const std::string& s)
    {
        std::size_t i = 0;
        std::size_t ln = 1;
        const auto N = s.size();
        std::string cur;
        auto flush = [&](std::size_t lineWhere) {
            if (!cur.empty()) {
                toks_.push_back({TokKind::Word, std::move(cur), lineWhere});
                cur.clear();
            }
        };
        while (i < N) {
            const char c = s[i];
            if (c == '\n') {
                flush(ln);
                ++ln;
                ++i;
                continue;
            }
            if (std::isspace(static_cast<unsigned char>(c))) {
                flush(ln);
                ++i;
                continue;
            }
            if (c == '/' && i + 1 < N && s[i + 1] == '/') {
                flush(ln);
                while (i < N && s[i] != '\n')
                    ++i;
                continue;
            }
            if (c == '/' && i + 1 < N && s[i + 1] == '*') {
                flush(ln);
                i += 2;
                while (i + 1 < N && !(s[i] == '*' && s[i + 1] == '/')) {
                    if (s[i] == '\n')
                        ++ln;
                    ++i;
                }
                if (i + 1 < N)
                    i += 2;
                continue;
            }
            TokKind k = TokKind::Word;
            switch (c) {
            case '(':
                k = TokKind::LParen;
                break;
            case ')':
                k = TokKind::RParen;
                break;
            case '{':
                k = TokKind::LBrace;
                break;
            case '}':
                k = TokKind::RBrace;
                break;
            case ';':
                k = TokKind::Semi;
                break;
            default:
                break;
            }
            if (k != TokKind::Word) {
                flush(ln);
                toks_.push_back({k, std::string(1, c), ln});
                ++i;
                continue;
            }
            cur.push_back(c);
            ++i;
        }
        flush(ln);
    }
    std::vector<Tok> toks_;
    std::size_t pos_ = 0;
};

// Skip an optional `<HeaderWord> { ... }` dictionary header at the
// front of `lex` (if any).  If the leading token is not a bare word
// followed by `{`, the cursor is restored and parsing continues.
void skip_optional_header(MeshTokenLex& lex)
{
    const std::size_t mark = lex.pos();
    if (lex.eof())
        return;
    const Tok& t = lex.peek();
    if (t.kind != TokKind::Word)
        return;
    // A leading list size is a digit-only word -- do not consume it.
    if (std::all_of(t.text.begin(), t.text.end(), [](char c) {
            return std::isdigit(static_cast<unsigned char>(c));
        })) {
        return;
    }
    lex.next();
    if (lex.eof() || lex.peek().kind != TokKind::LBrace) {
        lex.set_pos(mark);
        return;
    }
    int depth = 0;
    while (!lex.eof()) {
        const Tok& u = lex.next();
        if (u.kind == TokKind::LBrace)
            ++depth;
        else if (u.kind == TokKind::RBrace) {
            --depth;
            if (depth == 0)
                return;
        }
    }
}

// Parse leading `N (`, returning N and leaving cursor positioned right
// after the opening `(`.
[[nodiscard]] bool parse_list_open(MeshTokenLex& lex, long long& count, std::string& err)
{
    if (lex.eof() || lex.peek().kind != TokKind::Word) {
        err = "expected list size before '('";
        return false;
    }
    if (!to_int(lex.next().text, count) || count < 0) {
        err = "list size is not a non-negative integer";
        return false;
    }
    if (lex.eof() || lex.peek().kind != TokKind::LParen) {
        err = "expected '(' after list size";
        return false;
    }
    lex.next();
    return true;
}

// -- per-file parsers ---------------------------------------------------------

struct PointXYZ
{
    double x, y, z;
};

[[nodiscard]] bool parse_points(const std::string& text,
                                std::vector<PointXYZ>& out,
                                std::string& err)
{
    MeshTokenLex lex(text);
    skip_optional_header(lex);
    long long n = 0;
    if (!parse_list_open(lex, n, err)) {
        err = "points: " + err;
        return false;
    }
    out.reserve(static_cast<std::size_t>(n));
    for (long long i = 0; i < n; ++i) {
        if (lex.eof() || lex.next().kind != TokKind::LParen) {
            err = "points: expected '(' before point " + std::to_string(i);
            return false;
        }
        double xyz[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k) {
            if (lex.eof() || lex.peek().kind != TokKind::Word
                || !to_double(lex.next().text, xyz[k])) {
                err = "points: malformed coordinate at point " + std::to_string(i);
                return false;
            }
        }
        if (lex.eof() || lex.next().kind != TokKind::RParen) {
            err = "points: expected ')' after point " + std::to_string(i);
            return false;
        }
        out.push_back({xyz[0], xyz[1], xyz[2]});
    }
    return true;
}

[[nodiscard]] bool parse_faces(const std::string& text,
                               std::vector<std::vector<long long>>& out,
                               std::string& err)
{
    MeshTokenLex lex(text);
    skip_optional_header(lex);
    long long n = 0;
    if (!parse_list_open(lex, n, err)) {
        err = "faces: " + err;
        return false;
    }
    out.reserve(static_cast<std::size_t>(n));
    for (long long i = 0; i < n; ++i) {
        if (lex.eof() || lex.peek().kind != TokKind::Word) {
            err = "faces: expected vertex count at face " + std::to_string(i);
            return false;
        }
        long long nv = 0;
        if (!to_int(lex.next().text, nv) || nv < 3) {
            err = "faces: face " + std::to_string(i) + " has invalid vertex count";
            return false;
        }
        if (lex.eof() || lex.next().kind != TokKind::LParen) {
            err = "faces: expected '(' at face " + std::to_string(i);
            return false;
        }
        std::vector<long long> v;
        v.reserve(static_cast<std::size_t>(nv));
        for (long long k = 0; k < nv; ++k) {
            if (lex.eof() || lex.peek().kind != TokKind::Word) {
                err = "faces: expected vertex at face " + std::to_string(i);
                return false;
            }
            long long vi = 0;
            if (!to_int(lex.next().text, vi) || vi < 0) {
                err = "faces: malformed vertex at face " + std::to_string(i);
                return false;
            }
            v.push_back(vi);
        }
        if (lex.eof() || lex.next().kind != TokKind::RParen) {
            err = "faces: expected ')' at face " + std::to_string(i);
            return false;
        }
        out.push_back(std::move(v));
    }
    return true;
}

[[nodiscard]] bool parse_int_list(const std::string& text,
                                  std::vector<long long>& out,
                                  std::string& err,
                                  const char* what)
{
    MeshTokenLex lex(text);
    skip_optional_header(lex);
    long long n = 0;
    if (!parse_list_open(lex, n, err)) {
        err = std::string(what) + ": " + err;
        return false;
    }
    out.reserve(static_cast<std::size_t>(n));
    for (long long i = 0; i < n; ++i) {
        if (lex.eof() || lex.peek().kind != TokKind::Word) {
            err = std::string(what) + ": expected integer at row " + std::to_string(i);
            return false;
        }
        long long v = 0;
        if (!to_int(lex.next().text, v)) {
            err = std::string(what) + ": malformed integer at row " + std::to_string(i);
            return false;
        }
        out.push_back(v);
    }
    return true;
}

struct PatchSpec
{
    std::string name;
    std::string type;
    long long nFaces = 0;
    long long startFace = 0;
};

[[nodiscard]] bool parse_boundary(const std::string& text,
                                  std::vector<PatchSpec>& out,
                                  std::string& err)
{
    MeshTokenLex lex(text);
    skip_optional_header(lex);
    long long n = 0;
    if (!parse_list_open(lex, n, err)) {
        err = "boundary: " + err;
        return false;
    }
    out.reserve(static_cast<std::size_t>(n));
    for (long long i = 0; i < n; ++i) {
        if (lex.eof() || lex.peek().kind != TokKind::Word) {
            err = "boundary: expected patch name at patch " + std::to_string(i);
            return false;
        }
        PatchSpec ps;
        ps.name = lex.next().text;
        if (lex.eof() || lex.next().kind != TokKind::LBrace) {
            err = "boundary: expected '{' after patch name '" + ps.name + "'";
            return false;
        }
        while (!lex.eof() && lex.peek().kind != TokKind::RBrace) {
            if (lex.peek().kind != TokKind::Word) {
                err = "boundary: expected key in patch '" + ps.name + "'";
                return false;
            }
            const std::string key = lex.next().text;
            std::vector<std::string> vals;
            while (!lex.eof() && lex.peek().kind != TokKind::Semi) {
                if (lex.peek().kind == TokKind::RBrace)
                    break;
                vals.push_back(lex.next().text);
            }
            if (!lex.eof() && lex.peek().kind == TokKind::Semi)
                lex.next();
            if (key == "type" && !vals.empty()) {
                ps.type = vals.front();
            } else if (key == "nFaces" && !vals.empty()) {
                long long v = 0;
                if (!to_int(vals.front(), v) || v < 0) {
                    err = "boundary: malformed nFaces in '" + ps.name + "'";
                    return false;
                }
                ps.nFaces = v;
            } else if (key == "startFace" && !vals.empty()) {
                long long v = 0;
                if (!to_int(vals.front(), v) || v < 0) {
                    err = "boundary: malformed startFace in '" + ps.name + "'";
                    return false;
                }
                ps.startFace = v;
            }
            // Other keys are silently accepted and discarded.
        }
        if (lex.eof() || lex.next().kind != TokKind::RBrace) {
            err = "boundary: missing '}' for patch '" + ps.name + "'";
            return false;
        }
        out.push_back(std::move(ps));
    }
    return true;
}

} // namespace

// -- main parser --------------------------------------------------------------

FaceBasedMeshReadResult parse_face_based_mesh_strings(const std::string& pointsText,
                                                      const std::string& facesText,
                                                      const std::string& ownerText,
                                                      const std::string& neighbourText,
                                                      const std::string& boundaryText,
                                                      std::string sourceHint)
{
    FaceBasedMeshReadResult r;
    r.mesh.sourceFormat = "face_based_polymesh_dir";
    r.mesh.sourcePath = sourceHint;

    std::string err;

    std::vector<PointXYZ> points;
    if (!parse_points(pointsText, points, err)) {
        r.ok = false;
        r.error = "[" + sourceHint + "] " + err;
        return r;
    }
    if (points.empty()) {
        r.ok = false;
        r.error = "[" + sourceHint + "] points file is empty";
        return r;
    }

    std::vector<std::vector<long long>> faces;
    if (!parse_faces(facesText, faces, err)) {
        r.ok = false;
        r.error = "[" + sourceHint + "] " + err;
        return r;
    }
    if (faces.empty()) {
        r.ok = false;
        r.error = "[" + sourceHint + "] faces file is empty";
        return r;
    }

    std::vector<long long> owner;
    if (!parse_int_list(ownerText, owner, err, "owner")) {
        r.ok = false;
        r.error = "[" + sourceHint + "] " + err;
        return r;
    }
    if (owner.size() != faces.size()) {
        r.ok = false;
        r.error = "[" + sourceHint + "] owner.size() (" + std::to_string(owner.size())
                  + ") != faces.size() (" + std::to_string(faces.size()) + ")";
        return r;
    }

    std::vector<long long> neighbour;
    if (!neighbourText.empty()) {
        if (!parse_int_list(neighbourText, neighbour, err, "neighbour")) {
            r.ok = false;
            r.error = "[" + sourceHint + "] " + err;
            return r;
        }
        if (neighbour.size() > faces.size()) {
            r.ok = false;
            r.error = "[" + sourceHint + "] neighbour.size() > faces.size()";
            return r;
        }
    }

    std::vector<PatchSpec> patches;
    if (!boundaryText.empty()) {
        if (!parse_boundary(boundaryText, patches, err)) {
            r.ok = false;
            r.error = "[" + sourceHint + "] " + err;
            return r;
        }
    }

    long long maxCell = -1;
    for (long long c : owner)
        if (c > maxCell)
            maxCell = c;
    for (long long c : neighbour)
        if (c > maxCell)
            maxCell = c;
    if (maxCell < 0) {
        r.ok = false;
        r.error = "[" + sourceHint + "] no cells in owner";
        return r;
    }
    const std::size_t numCells = static_cast<std::size_t>(maxCell) + 1;

    for (std::size_t f = 0; f < faces.size(); ++f) {
        for (long long v : faces[f]) {
            if (v < 0 || static_cast<std::size_t>(v) >= points.size()) {
                r.ok = false;
                r.error = "[" + sourceHint + "] face " + std::to_string(f) + " references vertex "
                          + std::to_string(v) + " but only " + std::to_string(points.size())
                          + " points are defined";
                return r;
            }
        }
    }

    std::vector<std::vector<std::size_t>> cellFaces(numCells);
    for (std::size_t f = 0; f < faces.size(); ++f) {
        const long long oc = owner[f];
        if (oc < 0 || static_cast<std::size_t>(oc) >= numCells) {
            r.ok = false;
            r.error = "[" + sourceHint + "] face " + std::to_string(f) + " owner out of range";
            return r;
        }
        cellFaces[static_cast<std::size_t>(oc)].push_back(f);
        if (f < neighbour.size()) {
            const long long nc = neighbour[f];
            if (nc < 0 || static_cast<std::size_t>(nc) >= numCells) {
                r.ok = false;
                r.error =
                    "[" + sourceHint + "] face " + std::to_string(f) + " neighbour out of range";
                return r;
            }
            cellFaces[static_cast<std::size_t>(nc)].push_back(f);
        }
    }

    UnstructuredZone zone;
    zone.name = "imported_zone";
    zone.x.reserve(points.size());
    zone.y.reserve(points.size());
    zone.z.reserve(points.size());
    for (const auto& p : points) {
        zone.x.push_back(p.x);
        zone.y.push_back(p.y);
        zone.z.push_back(p.z);
    }

    std::unordered_map<int, std::size_t> volSec;
    auto get_vol_section = [&](ElementType t) -> std::size_t {
        const int key = static_cast<int>(t);
        if (auto it = volSec.find(key); it != volSec.end())
            return it->second;
        ElementSection sec;
        sec.type = t;
        sec.name = std::string("vol_") + element_type_name(t) + "_unordered";
        const std::size_t i = zone.sections.size();
        zone.sections.push_back(std::move(sec));
        volSec.emplace(key, i);
        return i;
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

    std::uint8_t maxDim = 0;

    for (std::size_t c = 0; c < numCells; ++c) {
        std::uint32_t nTri = 0, nQuad = 0;
        for (auto f : cellFaces[c]) {
            const auto sz = faces[f].size();
            if (sz == 3)
                ++nTri;
            else if (sz == 4)
                ++nQuad;
            else {
                r.ok = false;
                r.error = "[" + sourceHint + "] cell " + std::to_string(c)
                          + " uses an unsupported face with " + std::to_string(sz)
                          + " vertices (this pass supports tri/quad only)";
                return r;
            }
        }
        const ElementType et = classify_face_based_cell(nTri, nQuad);
        if (et == ElementType::Unknown) {
            r.ok = false;
            r.error = "[" + sourceHint + "] cell " + std::to_string(c)
                      + " has unsupported face signature (" + std::to_string(nTri) + " tri / "
                      + std::to_string(nQuad)
                      + " quad); supported cells: "
                        "Tetra4(4t,0q), Hexa8(0t,6q), Pyra5(4t,1q), "
                        "Penta6(2t,3q)";
            return r;
        }
        const std::uint8_t vpe = vertices_per_element(et);
        std::vector<long long> uniq;
        uniq.reserve(vpe);
        for (auto f : cellFaces[c]) {
            for (long long v : faces[f]) {
                bool seen = false;
                for (long long u : uniq)
                    if (u == v) {
                        seen = true;
                        break;
                    }
                if (!seen)
                    uniq.push_back(v);
            }
        }
        if (uniq.size() != vpe) {
            r.ok = false;
            r.error = "[" + sourceHint + "] cell " + std::to_string(c) + " classified as "
                      + element_type_name(et) + " but face vertices produced "
                      + std::to_string(uniq.size()) + " distinct points (expected "
                      + std::to_string(vpe) + ")";
            return r;
        }
        ElementSection& sec = zone.sections[get_vol_section(et)];
        for (long long v : uniq) {
            sec.nodes.push_back(static_cast<NodeIdx>(v));
        }
        if (const auto d = dim_of(et); d > maxDim)
            maxDim = d;
    }

    for (const auto& ps : patches) {
        if (ps.nFaces == 0)
            continue;
        if (ps.startFace < 0 || ps.startFace + ps.nFaces > static_cast<long long>(faces.size())) {
            r.ok = false;
            r.error = "[" + sourceHint + "] patch '" + ps.name + "' face range ["
                      + std::to_string(ps.startFace) + ", "
                      + std::to_string(ps.startFace + ps.nFaces)
                      + ") is out of bounds (faces.size()=" + std::to_string(faces.size()) + ")";
            return r;
        }
        const auto firstSz = faces[static_cast<std::size_t>(ps.startFace)].size();
        ElementType et = (firstSz == 3)   ? ElementType::Tri3
                         : (firstSz == 4) ? ElementType::Quad4
                                          : ElementType::Unknown;
        if (et == ElementType::Unknown) {
            r.ok = false;
            r.error = "[" + sourceHint + "] patch '" + ps.name + "' first face has "
                      + std::to_string(firstSz) + " vertices (only 3/4 supported)";
            return r;
        }
        for (long long off = 0; off < ps.nFaces; ++off) {
            const auto& fv = faces[static_cast<std::size_t>(ps.startFace + off)];
            if (fv.size() != firstSz) {
                r.ok = false;
                r.error = "[" + sourceHint + "] patch '" + ps.name + "' is non-homogeneous (face "
                          + std::to_string(ps.startFace + off) + " has " + std::to_string(fv.size())
                          + " vertices, expected " + std::to_string(firstSz) + ")";
                return r;
            }
        }
        ElementSection sec;
        sec.type = et;
        sec.name = "patch_" + ps.name;
        sec.nodes.reserve(static_cast<std::size_t>(ps.nFaces * firstSz));
        for (long long off = 0; off < ps.nFaces; ++off) {
            for (long long v : faces[static_cast<std::size_t>(ps.startFace + off)]) {
                sec.nodes.push_back(static_cast<NodeIdx>(v));
            }
        }
        const std::size_t si = zone.sections.size();
        zone.sections.push_back(std::move(sec));

        BoundaryPatch bp;
        bp.name = ps.name;
        bp.bcType = ps.type.empty() ? std::string("patch") : ps.type;
        bp.faceElementIndices.push_back(static_cast<std::uint32_t>(si));
        zone.boundaries.push_back(std::move(bp));
    }
    std::sort(zone.boundaries.begin(),
              zone.boundaries.end(),
              [](const BoundaryPatch& a, const BoundaryPatch& b) { return a.name < b.name; });

    r.mesh.zones.push_back(std::move(zone));
    r.dimension = maxDim ? maxDim : static_cast<std::uint8_t>(1);
    r.ok = true;
    return r;
}

FaceBasedMeshReadResult read_face_based_mesh_dir(const std::string& dirPath)
{
    FaceBasedMeshReadResult r;
    auto join = [&](const char* name) {
        std::string p = dirPath;
        if (!p.empty() && p.back() != '/' && p.back() != '\\')
            p.push_back('/');
        p.append(name);
        return p;
    };
    const std::string pts = read_whole_file(join("points"));
    const std::string fcs = read_whole_file(join("faces"));
    const std::string own = read_whole_file(join("owner"));
    const std::string nei = read_whole_file(join("neighbour"));
    const std::string bnd = read_whole_file(join("boundary"));
    if (pts.empty() || fcs.empty() || own.empty()) {
        r.ok = false;
        r.error = "face-based mesh dir '" + dirPath + "' missing required points/faces/owner";
        return r;
    }
    return parse_face_based_mesh_strings(pts, fcs, own, nei, bnd, dirPath);
}

} // namespace simall::io
