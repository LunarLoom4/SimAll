// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/VtuReader.cpp
// Phase  : 23 Pass 9
//
// ASCII VTK Unstructured Grid (`.vtu`) reader.  Header documents the format
// and constraints accepted.  All XML scanning is implemented inline (no
// third-party XML dependency): the file is small, the schema is fixed.
// =============================================================================

#include "io/VtuReader.hpp"

#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace simall::io
{

// =============================================================================
// Minimal XML mini-DOM.  Sufficient for the VTU subset documented in the
// header; supports: <?xml ?> prolog, comments, CDATA, self-closing tags,
// double- and single-quoted attribute values, and the five standard XML
// entities (&amp; &lt; &gt; &quot; &apos;).  Namespaces are not interpreted
// (colon characters are treated as ordinary name characters).
// =============================================================================
namespace
{

struct XmlAttr
{
    std::string name;
    std::string value;
};

struct XmlNode
{
    std::string name;
    std::vector<XmlAttr> attrs;
    std::string text; // concatenated CDATA + char data
    std::vector<XmlNode> children;
};

inline bool is_ws(char c) noexcept
{
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}
inline bool is_name_start(char c) noexcept
{
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_';
}
inline bool is_name_char(char c) noexcept
{
    return is_name_start(c) || (c >= '0' && c <= '9') || c == ':' || c == '.' || c == '-';
}

void skip_ws(std::string_view s, std::size_t& p)
{
    while (p < s.size() && is_ws(s[p]))
        ++p;
}

bool starts_with(std::string_view s, std::size_t p, std::string_view what)
{
    return p + what.size() <= s.size() && std::memcmp(s.data() + p, what.data(), what.size()) == 0;
}

bool skip_misc(std::string_view s, std::size_t& p, std::string& err)
{
    while (true) {
        skip_ws(s, p);
        if (starts_with(s, p, "<!--")) {
            auto end = s.find("-->", p + 4);
            if (end == std::string_view::npos) {
                err = "unterminated XML comment";
                return false;
            }
            p = end + 3;
            continue;
        }
        if (starts_with(s, p, "<?")) {
            auto end = s.find("?>", p + 2);
            if (end == std::string_view::npos) {
                err = "unterminated XML processing instruction";
                return false;
            }
            p = end + 2;
            continue;
        }
        if (starts_with(s, p, "<!DOCTYPE")) {
            auto end = s.find('>', p + 9);
            if (end == std::string_view::npos) {
                err = "unterminated DOCTYPE";
                return false;
            }
            p = end + 1;
            continue;
        }
        break;
    }
    return true;
}

bool parse_name(std::string_view s, std::size_t& p, std::string& out, std::string& err)
{
    if (p >= s.size() || !is_name_start(s[p])) {
        err = "expected XML name";
        return false;
    }
    auto start = p;
    while (p < s.size() && is_name_char(s[p]))
        ++p;
    out.assign(s.data() + start, p - start);
    return true;
}

std::string xml_unescape(std::string_view in)
{
    std::string out;
    out.reserve(in.size());
    for (std::size_t i = 0; i < in.size();) {
        if (in[i] == '&') {
            auto sc = in.find(';', i + 1);
            if (sc != std::string_view::npos) {
                auto ent = in.substr(i + 1, sc - i - 1);
                if (ent == "amp") {
                    out += '&';
                    i = sc + 1;
                    continue;
                }
                if (ent == "lt") {
                    out += '<';
                    i = sc + 1;
                    continue;
                }
                if (ent == "gt") {
                    out += '>';
                    i = sc + 1;
                    continue;
                }
                if (ent == "quot") {
                    out += '"';
                    i = sc + 1;
                    continue;
                }
                if (ent == "apos") {
                    out += '\'';
                    i = sc + 1;
                    continue;
                }
                if (!ent.empty() && ent[0] == '#') {
                    // numeric character reference (rare in VTU but cheap)
                    int base = 10;
                    std::string_view digits = ent.substr(1);
                    if (!digits.empty() && (digits[0] == 'x' || digits[0] == 'X')) {
                        base = 16;
                        digits.remove_prefix(1);
                    }
                    int code = 0;
                    bool ok = !digits.empty();
                    for (char c : digits) {
                        int d = -1;
                        if (c >= '0' && c <= '9')
                            d = c - '0';
                        else if (base == 16 && c >= 'a' && c <= 'f')
                            d = 10 + c - 'a';
                        else if (base == 16 && c >= 'A' && c <= 'F')
                            d = 10 + c - 'A';
                        if (d < 0 || d >= base) {
                            ok = false;
                            break;
                        }
                        code = code * base + d;
                    }
                    if (ok && code > 0 && code < 128) {
                        out += static_cast<char>(code);
                        i = sc + 1;
                        continue;
                    }
                }
            }
        }
        out += in[i++];
    }
    return out;
}

bool parse_attrs(std::string_view s, std::size_t& p, std::vector<XmlAttr>& attrs, std::string& err)
{
    while (true) {
        skip_ws(s, p);
        if (p >= s.size()) {
            err = "unexpected EOF inside tag";
            return false;
        }
        if (s[p] == '>' || s[p] == '/')
            return true;
        XmlAttr a;
        if (!parse_name(s, p, a.name, err))
            return false;
        skip_ws(s, p);
        if (p >= s.size() || s[p] != '=') {
            err = "expected '=' after attribute name '" + a.name + "'";
            return false;
        }
        ++p;
        skip_ws(s, p);
        if (p >= s.size() || (s[p] != '"' && s[p] != '\'')) {
            err = "expected quoted value for attribute '" + a.name + "'";
            return false;
        }
        char quote = s[p++];
        auto vstart = p;
        while (p < s.size() && s[p] != quote)
            ++p;
        if (p >= s.size()) {
            err = "unterminated attribute value for '" + a.name + "'";
            return false;
        }
        a.value = xml_unescape(std::string_view(s.data() + vstart, p - vstart));
        ++p; // consume closing quote
        attrs.push_back(std::move(a));
    }
}

bool parse_element(std::string_view s, std::size_t& p, XmlNode& node, std::string& err)
{
    if (p >= s.size() || s[p] != '<') {
        err = "expected '<' to start element";
        return false;
    }
    ++p;
    if (!parse_name(s, p, node.name, err))
        return false;
    if (!parse_attrs(s, p, node.attrs, err))
        return false;
    skip_ws(s, p);
    if (p >= s.size()) {
        err = "unexpected EOF in opening tag <" + node.name + ">";
        return false;
    }
    if (s[p] == '/') {
        ++p;
        if (p >= s.size() || s[p] != '>') {
            err = "expected '>' after '/' in <" + node.name + "/>";
            return false;
        }
        ++p;
        return true;
    }
    if (s[p] != '>') {
        err = "expected '>' to close opening tag <" + node.name + ">";
        return false;
    }
    ++p;

    // Parse content
    while (true) {
        // Char data up to next '<'
        auto tstart = p;
        while (p < s.size() && s[p] != '<')
            ++p;
        if (p > tstart) {
            node.text.append(s.data() + tstart, p - tstart);
        }
        if (p >= s.size()) {
            err = "unterminated element <" + node.name + ">";
            return false;
        }
        if (starts_with(s, p, "</")) {
            p += 2;
            std::string close_name;
            if (!parse_name(s, p, close_name, err))
                return false;
            if (close_name != node.name) {
                err = "mismatched close tag </" + close_name + "> for <" + node.name + ">";
                return false;
            }
            skip_ws(s, p);
            if (p >= s.size() || s[p] != '>') {
                err = "expected '>' in close tag </" + close_name + ">";
                return false;
            }
            ++p;
            return true;
        }
        if (starts_with(s, p, "<!--")) {
            auto end = s.find("-->", p + 4);
            if (end == std::string_view::npos) {
                err = "unterminated comment inside <" + node.name + ">";
                return false;
            }
            p = end + 3;
            continue;
        }
        if (starts_with(s, p, "<![CDATA[")) {
            auto end = s.find("]]>", p + 9);
            if (end == std::string_view::npos) {
                err = "unterminated CDATA inside <" + node.name + ">";
                return false;
            }
            node.text.append(s.data() + p + 9, end - p - 9);
            p = end + 3;
            continue;
        }
        // child element
        XmlNode child;
        if (!parse_element(s, p, child, err))
            return false;
        node.children.push_back(std::move(child));
    }
}

bool parse_xml(const std::string& text, XmlNode& root, std::string& err)
{
    std::string_view s(text);
    std::size_t p = 0;
    if (!skip_misc(s, p, err))
        return false;
    if (p >= s.size()) {
        err = "empty XML document";
        return false;
    }
    if (!parse_element(s, p, root, err))
        return false;
    if (!skip_misc(s, p, err))
        return false;
    // Trailing data is permitted to be just whitespace/comments; ignore.
    return true;
}

const XmlNode* find_child(const XmlNode& parent, std::string_view name)
{
    for (auto const& c : parent.children) {
        if (c.name == name)
            return &c;
    }
    return nullptr;
}

const XmlAttr* find_attr(const XmlNode& n, std::string_view name)
{
    for (auto const& a : n.attrs) {
        if (a.name == name)
            return &a;
    }
    return nullptr;
}

bool attr_or_empty(const XmlNode& n, std::string_view name, std::string& out)
{
    if (auto* a = find_attr(n, name)) {
        out = a->value;
        return true;
    }
    out.clear();
    return false;
}

// =============================================================================
// Numeric token parsing (whitespace-separated ASCII).
// =============================================================================
bool parse_doubles(const std::string& text,
                   std::size_t expected,
                   std::vector<double>& out,
                   std::string& err)
{
    out.clear();
    out.reserve(expected);
    std::istringstream is(text);
    double v;
    while (is >> v)
        out.push_back(v);
    if (!is.eof()) {
        // Did not consume all input cleanly — find offending token.
        is.clear();
        std::string bad;
        is >> bad;
        if (!bad.empty()) {
            err = "non-numeric token '" + bad + "' in floating-point DataArray";
            return false;
        }
    }
    if (out.size() != expected) {
        err = "expected " + std::to_string(expected) + " floats in DataArray but parsed "
              + std::to_string(out.size());
        return false;
    }
    return true;
}

bool parse_longs(const std::string& text,
                 std::size_t expected,
                 std::vector<long long>& out,
                 std::string& err)
{
    out.clear();
    out.reserve(expected);
    std::istringstream is(text);
    long long v;
    while (is >> v)
        out.push_back(v);
    if (!is.eof()) {
        is.clear();
        std::string bad;
        is >> bad;
        if (!bad.empty()) {
            err = "non-integer token '" + bad + "' in integer DataArray";
            return false;
        }
    }
    if (out.size() != expected) {
        err = "expected " + std::to_string(expected) + " integers in DataArray but parsed "
              + std::to_string(out.size());
        return false;
    }
    return true;
}

// =============================================================================
// VTK cell-type → SimAll element-type table.
// =============================================================================
constexpr int kVtkVertex = 1;
constexpr int kVtkLine = 3;
constexpr int kVtkTriangle = 5;
constexpr int kVtkQuad = 9;
constexpr int kVtkTetra = 10;
constexpr int kVtkHexahedron = 12;
constexpr int kVtkWedge = 13;
constexpr int kVtkPyramid = 14;

} // anonymous namespace

ElementType vtu_element_type(int vtkCellType) noexcept
{
    switch (vtkCellType) {
    case kVtkLine:
        return ElementType::Bar2;
    case kVtkTriangle:
        return ElementType::Tri3;
    case kVtkQuad:
        return ElementType::Quad4;
    case kVtkTetra:
        return ElementType::Tetra4;
    case kVtkHexahedron:
        return ElementType::Hexa8;
    case kVtkWedge:
        return ElementType::Penta6;
    case kVtkPyramid:
        return ElementType::Pyra5;
    default:
        return ElementType::Unknown;
    }
}

// =============================================================================
// Main parse driver.
// =============================================================================
VtuReadResult parse_vtu_string(const std::string& text, std::string hint)
{
    VtuReadResult r;

    XmlNode root;
    std::string err;
    if (!parse_xml(text, root, err)) {
        r.error = "VTU '" + hint + "': XML parse error: " + err;
        return r;
    }

    if (root.name != "VTKFile") {
        r.error = "VTU '" + hint + "': root element is <" + root.name + ">, expected <VTKFile>";
        return r;
    }

    std::string vtkType;
    attr_or_empty(root, "type", vtkType);
    if (vtkType != "UnstructuredGrid") {
        r.error = "VTU '" + hint + "': <VTKFile type=\"" + vtkType
                  + "\"> not supported (only 'UnstructuredGrid')";
        return r;
    }

    std::string compressor;
    if (attr_or_empty(root, "compressor", compressor) && !compressor.empty()) {
        r.error = "VTU '" + hint + "': compressor='" + compressor
                  + "' is not supported (re-export uncompressed ASCII)";
        return r;
    }

    if (find_child(root, "AppendedData") != nullptr) {
        r.error = "VTU '" + hint
                  + "': <AppendedData> blocks are not "
                    "supported (re-export with format=\"ascii\")";
        return r;
    }

    const XmlNode* grid = find_child(root, "UnstructuredGrid");
    if (!grid) {
        r.error = "VTU '" + hint + "': missing <UnstructuredGrid> child";
        return r;
    }

    // Reject multi-piece datasets (parallel-output convenience).
    int pieceCount = 0;
    for (auto const& c : grid->children) {
        if (c.name == "Piece")
            ++pieceCount;
    }
    if (pieceCount == 0) {
        r.error = "VTU '" + hint + "': <UnstructuredGrid> contains no <Piece>";
        return r;
    }
    if (pieceCount > 1) {
        r.error = "VTU '" + hint + "': multi-piece datasets (" + std::to_string(pieceCount)
                  + " pieces) are not supported; merge pieces first";
        return r;
    }

    const XmlNode* piece = find_child(*grid, "Piece");

    std::string npStr, ncStr;
    attr_or_empty(*piece, "NumberOfPoints", npStr);
    attr_or_empty(*piece, "NumberOfCells", ncStr);
    if (npStr.empty() || ncStr.empty()) {
        r.error = "VTU '" + hint
                  + "': <Piece> missing NumberOfPoints "
                    "or NumberOfCells attribute";
        return r;
    }
    long long nPoints = 0, nCells = 0;
    try {
        nPoints = std::stoll(npStr);
        nCells = std::stoll(ncStr);
    } catch (...) {
        r.error = "VTU '" + hint
                  + "': <Piece> has non-integer "
                    "NumberOfPoints/NumberOfCells";
        return r;
    }
    if (nPoints < 0 || nCells < 0) {
        r.error = "VTU '" + hint + "': negative NumberOfPoints/NumberOfCells";
        return r;
    }

    // ----- <Points> ---------------------------------------------------------
    const XmlNode* points = find_child(*piece, "Points");
    if (!points) {
        r.error = "VTU '" + hint + "': <Piece> missing <Points>";
        return r;
    }
    const XmlNode* pointsDA = find_child(*points, "DataArray");
    if (!pointsDA) {
        r.error = "VTU '" + hint + "': <Points> missing <DataArray>";
        return r;
    }
    {
        std::string fmt;
        attr_or_empty(*pointsDA, "format", fmt);
        if (fmt != "ascii") {
            r.error = "VTU '" + hint + "': <Points>/<DataArray format=\"" + fmt
                      + "\"> not supported (must be \"ascii\")";
            return r;
        }
        std::string ncomp;
        attr_or_empty(*pointsDA, "NumberOfComponents", ncomp);
        if (!ncomp.empty() && ncomp != "3") {
            r.error = "VTU '" + hint
                      + "': <Points>/<DataArray "
                        "NumberOfComponents=\""
                      + ncomp + "\"> not supported (must be 3)";
            return r;
        }
    }

    std::vector<double> coords;
    if (!parse_doubles(pointsDA->text, static_cast<std::size_t>(nPoints) * 3, coords, err)) {
        r.error = "VTU '" + hint + "': <Points>: " + err;
        return r;
    }

    // ----- <Cells> ---------------------------------------------------------
    const XmlNode* cells = find_child(*piece, "Cells");
    if (!cells) {
        r.error = "VTU '" + hint + "': <Piece> missing <Cells>";
        return r;
    }
    const XmlNode* connDA = nullptr;
    const XmlNode* offsetDA = nullptr;
    const XmlNode* typesDA = nullptr;
    for (auto const& c : cells->children) {
        if (c.name != "DataArray")
            continue;
        std::string nm;
        attr_or_empty(c, "Name", nm);
        if (nm == "connectivity")
            connDA = &c;
        else if (nm == "offsets")
            offsetDA = &c;
        else if (nm == "types")
            typesDA = &c;
    }
    if (!connDA || !offsetDA || !typesDA) {
        r.error = "VTU '" + hint
                  + "': <Cells> requires DataArrays named "
                    "'connectivity', 'offsets', and 'types'";
        return r;
    }
    for (auto const* da : {connDA, offsetDA, typesDA}) {
        std::string fmt;
        attr_or_empty(*da, "format", fmt);
        if (fmt != "ascii") {
            std::string nm;
            attr_or_empty(*da, "Name", nm);
            r.error = "VTU '" + hint + "': <Cells>/<DataArray Name=\"" + nm + "\" format=\"" + fmt
                      + "\"> not supported (must be \"ascii\")";
            return r;
        }
    }

    std::vector<long long> offsets;
    if (!parse_longs(offsetDA->text, static_cast<std::size_t>(nCells), offsets, err)) {
        r.error = "VTU '" + hint + "': <Cells>/offsets: " + err;
        return r;
    }
    std::vector<long long> vtkTypes;
    if (!parse_longs(typesDA->text, static_cast<std::size_t>(nCells), vtkTypes, err)) {
        r.error = "VTU '" + hint + "': <Cells>/types: " + err;
        return r;
    }

    std::size_t connExpected = nCells > 0 ? static_cast<std::size_t>(offsets.back()) : 0;
    std::vector<long long> connectivity;
    if (!parse_longs(connDA->text, connExpected, connectivity, err)) {
        r.error = "VTU '" + hint + "': <Cells>/connectivity: " + err;
        return r;
    }

    // ----- Build UnstructuredZone ------------------------------------------
    UnstructuredZone zone;
    zone.name = "vtu_zone";
    zone.x.resize(static_cast<std::size_t>(nPoints));
    zone.y.resize(static_cast<std::size_t>(nPoints));
    zone.z.resize(static_cast<std::size_t>(nPoints));
    for (long long i = 0; i < nPoints; ++i) {
        zone.x[static_cast<std::size_t>(i)] = coords[3 * i + 0];
        zone.y[static_cast<std::size_t>(i)] = coords[3 * i + 1];
        zone.z[static_cast<std::size_t>(i)] = coords[3 * i + 2];
    }

    // Group cells by element type into ElementSections, in the order each
    // type first appears.  This keeps the section ordering stable and
    // independent of cell numbering.
    std::unordered_map<int, std::size_t> typeToSecIdx;
    std::uint8_t maxDim = 0;
    long long prevOffset = 0;
    for (long long ci = 0; ci < nCells; ++ci) {
        long long off = offsets[static_cast<std::size_t>(ci)];
        long long npe = off - prevOffset;
        int vtkT = static_cast<int>(vtkTypes[static_cast<std::size_t>(ci)]);

        if (off < prevOffset) {
            r.error = "VTU '" + hint + "': non-monotonic 'offsets' at cell " + std::to_string(ci);
            return r;
        }
        if (static_cast<std::size_t>(off) > connectivity.size()) {
            r.error = "VTU '" + hint + "': 'offsets' references position " + std::to_string(off)
                      + " past end of connectivity (" + std::to_string(connectivity.size()) + ")";
            return r;
        }

        if (vtkT == kVtkVertex) {
            prevOffset = off; // silently skip standalone vertices
            continue;
        }

        ElementType et = vtu_element_type(vtkT);
        if (et == ElementType::Unknown) {
            r.error = "VTU '" + hint + "': unsupported VTK cell type " + std::to_string(vtkT)
                      + " at cell " + std::to_string(ci)
                      + " (high-order / polyhedral cells are not yet implemented)";
            return r;
        }
        auto expectedNpe = static_cast<long long>(vertices_per_element(et));
        if (npe != expectedNpe) {
            r.error = "VTU '" + hint + "': cell " + std::to_string(ci) + " declared as VTK type "
                      + std::to_string(vtkT) + " (" + element_type_name(et) + ") expects "
                      + std::to_string(expectedNpe) + " nodes but has " + std::to_string(npe);
            return r;
        }

        auto it = typeToSecIdx.find(vtkT);
        std::size_t sidx;
        if (it == typeToSecIdx.end()) {
            ElementSection sec;
            sec.type = et;
            sec.name = std::string("vtu_") + element_type_name(et);
            zone.sections.push_back(std::move(sec));
            sidx = zone.sections.size() - 1;
            typeToSecIdx.emplace(vtkT, sidx);
        } else {
            sidx = it->second;
        }

        auto& sec = zone.sections[sidx];
        for (long long k = prevOffset; k < off; ++k) {
            long long ix = connectivity[static_cast<std::size_t>(k)];
            if (ix < 0 || ix >= nPoints) {
                r.error = "VTU '" + hint + "': cell " + std::to_string(ci)
                          + " references out-of-range node index " + std::to_string(ix)
                          + " (nPoints=" + std::to_string(nPoints) + ")";
                return r;
            }
            sec.nodes.push_back(static_cast<NodeIdx>(ix));
        }

        // track top dimension to expose via VtuReadResult::dimension
        std::uint8_t d = 0;
        switch (et) {
        case ElementType::Bar2:
            d = 1;
            break;
        case ElementType::Tri3:
        case ElementType::Quad4:
            d = 2;
            break;
        case ElementType::Tetra4:
        case ElementType::Pyra5:
        case ElementType::Penta6:
        case ElementType::Hexa8:
            d = 3;
            break;
        default:
            d = 0;
            break;
        }
        if (d > maxDim)
            maxDim = d;

        prevOffset = off;
    }

    if (maxDim == 0)
        maxDim = 3; // empty / vertex-only mesh - assume volume
    r.dimension = maxDim;

    r.mesh.sourceFormat = "vtu";
    r.mesh.sourcePath = hint;
    r.mesh.zones.push_back(std::move(zone));
    r.ok = true;
    return r;
}

VtuReadResult read_vtu(const std::string& path)
{
    VtuReadResult r;
    std::ifstream fs(path, std::ios::binary);
    if (!fs) {
        r.error = "VTU: cannot open '" + path + "' for reading";
        return r;
    }
    std::ostringstream ss;
    ss << fs.rdbuf();
    return parse_vtu_string(ss.str(), path);
}

} // namespace simall::io
