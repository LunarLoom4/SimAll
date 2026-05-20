// =============================================================================
// SimAll Beta — GUI Core (Qt-free)
// File   : src/gui_core/DockPerspective.cpp
// =============================================================================
#include "gui_core/DockPerspective.hpp"

#include <charconv>
#include <sstream>

namespace simall::gui_core
{

const DockPlacement* DockPerspective::find(std::string_view id) const noexcept
{
    for (auto& p : placements)
        if (p.panelId == id)
            return &p;
    return nullptr;
}
DockPlacement* DockPerspective::find(std::string_view id) noexcept
{
    for (auto& p : placements)
        if (p.panelId == id)
            return &p;
    return nullptr;
}
DockPlacement& DockPerspective::upsert(const std::string& id)
{
    if (auto* p = find(id))
        return *p;
    placements.push_back({id});
    return placements.back();
}

std::string area_to_string(DockArea a)
{
    switch (a) {
    case DockArea::Left:
        return "Left";
    case DockArea::Right:
        return "Right";
    case DockArea::Top:
        return "Top";
    case DockArea::Bottom:
        return "Bottom";
    case DockArea::Floating:
        return "Floating";
    case DockArea::Central:
        return "Central";
    case DockArea::Hidden:
        return "Hidden";
    }
    return "Hidden";
}

std::optional<DockArea> area_from_string(std::string_view s)
{
    if (s == "Left")
        return DockArea::Left;
    if (s == "Right")
        return DockArea::Right;
    if (s == "Top")
        return DockArea::Top;
    if (s == "Bottom")
        return DockArea::Bottom;
    if (s == "Floating")
        return DockArea::Floating;
    if (s == "Central")
        return DockArea::Central;
    if (s == "Hidden")
        return DockArea::Hidden;
    return std::nullopt;
}

std::string serialize(const DockPerspective& p)
{
    std::ostringstream os;
    os << "[perspective]\n"
       << "name=" << p.name << '\n'
       << "qt=" << p.qtBlob << '\n';
    for (const auto& pl : p.placements) {
        os << "[panel]\n"
           << "id=" << pl.panelId << '\n'
           << "area=" << area_to_string(pl.area) << '\n'
           << "visible=" << (pl.visible ? '1' : '0') << '\n'
           << "floating=" << (pl.floating ? '1' : '0') << '\n'
           << "rect=" << pl.x << ',' << pl.y << ',' << pl.w << ',' << pl.h << '\n'
           << "tab=" << pl.tabGroupId << ',' << pl.tabIndex << '\n';
    }
    return os.str();
}

namespace
{

inline std::string_view trim(std::string_view s)
{
    while (!s.empty()
           && (s.front() == ' ' || s.front() == '\t' || s.front() == '\r' || s.front() == '\n'))
        s.remove_prefix(1);
    while (!s.empty()
           && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r' || s.back() == '\n'))
        s.remove_suffix(1);
    return s;
}

bool parse_int(std::string_view s, int& out)
{
    auto r = std::from_chars(s.data(), s.data() + s.size(), out);
    return r.ec == std::errc{};
}

void split_csv(std::string_view s, std::vector<std::string>& out)
{
    out.clear();
    size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == ',') {
            out.emplace_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
}

} // namespace

std::optional<DockPerspective> deserialize(std::string_view txt)
{
    DockPerspective p;
    bool inPanel = false;
    DockPlacement cur;
    std::vector<std::string> csv;
    auto flush_panel = [&] {
        if (inPanel) {
            p.placements.push_back(cur);
            cur = {};
            inPanel = false;
        }
    };
    size_t lineStart = 0;
    for (size_t i = 0; i <= txt.size(); ++i) {
        if (i == txt.size() || txt[i] == '\n') {
            auto line = trim(txt.substr(lineStart, i - lineStart));
            lineStart = i + 1;
            if (line.empty())
                continue;
            if (line == "[perspective]") {
                flush_panel();
                continue;
            }
            if (line == "[panel]") {
                flush_panel();
                inPanel = true;
                continue;
            }
            auto eq = line.find('=');
            if (eq == std::string_view::npos)
                continue;
            auto key = line.substr(0, eq);
            auto val = line.substr(eq + 1);
            if (!inPanel) {
                if (key == "name")
                    p.name = std::string(val);
                else if (key == "qt")
                    p.qtBlob = std::string(val);
            } else {
                if (key == "id")
                    cur.panelId = std::string(val);
                else if (key == "area") {
                    auto a = area_from_string(val);
                    if (a)
                        cur.area = *a;
                } else if (key == "visible")
                    cur.visible = (val == "1");
                else if (key == "floating")
                    cur.floating = (val == "1");
                else if (key == "rect") {
                    split_csv(val, csv);
                    if (csv.size() == 4) {
                        parse_int(csv[0], cur.x);
                        parse_int(csv[1], cur.y);
                        parse_int(csv[2], cur.w);
                        parse_int(csv[3], cur.h);
                    }
                } else if (key == "tab") {
                    split_csv(val, csv);
                    if (csv.size() == 2) {
                        parse_int(csv[0], cur.tabGroupId);
                        parse_int(csv[1], cur.tabIndex);
                    }
                }
            }
        }
    }
    flush_panel();
    return p;
}

} // namespace simall::gui_core
