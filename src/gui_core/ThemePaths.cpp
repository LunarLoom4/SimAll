// =============================================================================
// SimAll Beta — GUI Core (Qt-free)
// File   : src/gui_core/ThemePaths.cpp
// =============================================================================
#include "gui_core/ThemePaths.hpp"

#include <cstdlib>
#include <string>

namespace simall::gui_core {

namespace fs = std::filesystem;

namespace {

std::optional<fs::path> env_path(const char* var) {
#ifdef _MSC_VER
    char*  raw   = nullptr;
    size_t sz    = 0;
    if (_dupenv_s(&raw, &sz, var) != 0 || !raw) return std::nullopt;
    fs::path p(raw);
    std::free(raw);
    return p;
#else
    if (const char* v = std::getenv(var)) return fs::path(v);
    return std::nullopt;
#endif
}

void push_candidates(std::vector<fs::path>& out, const fs::path& root) {
    if (root.empty()) return;
    out.push_back(root);
    out.push_back(root / "resources");
    out.push_back(root.parent_path() / "resources");
    out.push_back(root.parent_path().parent_path() / "resources");
}

}  // namespace

ThemeAssetLocations locate_theme_assets(const fs::path& executableDir,
                                        const fs::path& sourceRoot,
                                        std::string_view themeName) {
    ThemeAssetLocations out;
    std::vector<fs::path> roots;
    if (auto env = env_path("SIMALL_THEME_DIR")) roots.push_back(*env);
    push_candidates(roots, executableDir);
    push_candidates(roots, sourceRoot);
#ifndef _WIN32
    roots.emplace_back("/usr/share/simall");
    roots.emplace_back("/usr/local/share/simall");
#endif
    const std::string qssName = "simall_" + std::string(themeName) + ".qss";

    for (const auto& r : roots) {
        std::error_code ec;
        if (!fs::exists(r, ec)) { out.searchedDirectories.push_back(r); continue; }

        const fs::path stylesDir = r / "styles";
        const fs::path fontsDir  = r / "fonts" / "Inter";
        out.searchedDirectories.push_back(stylesDir);
        out.searchedDirectories.push_back(fontsDir);

        if (!out.stylesheet) {
            fs::path candidate = stylesDir / qssName;
            if (fs::exists(candidate, ec)) out.stylesheet = candidate;
        }
        if (fs::exists(fontsDir, ec)) {
            for (auto& e : fs::directory_iterator(fontsDir, ec)) {
                const auto& p = e.path();
                if (p.extension() == ".otf" || p.extension() == ".ttf")
                    out.fonts.push_back(p);
            }
        }
        if (out.stylesheet && !out.fonts.empty()) break;
    }
    return out;
}

}  // namespace simall::gui_core
