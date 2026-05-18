// =============================================================================
// SimAll Beta — GUI Core (Qt-free)
// File   : src/gui_core/ThemePaths.hpp
//
// File-system search for theme assets (QSS sheet and Inter font files).
// Lookup order (first hit wins):
//   1. $SIMALL_THEME_DIR (env override, if set)
//   2. <executable_dir>/resources/styles  +  /resources/fonts
//   3. <source_root>/resources/styles     +  /resources/fonts
//   4. /usr/share/simall/...                  (POSIX install prefix)
//
// Returning std::optional<std::filesystem::path> makes the Qt-side
// ThemeManager fail gracefully (it falls back to the built-in dark palette
// without QSS overrides) when the asset directory is not deployed.
// =============================================================================
#pragma once

#include <filesystem>
#include <optional>
#include <vector>

namespace simall::gui_core {

struct ThemeAssetLocations {
    std::optional<std::filesystem::path> stylesheet;          // simall_dark.qss
    std::vector<std::filesystem::path>   fonts;               // Inter-*.otf
    std::vector<std::filesystem::path>   searchedDirectories; // diagnostic
};

[[nodiscard]] ThemeAssetLocations locate_theme_assets(
    const std::filesystem::path& executableDir = {},
    const std::filesystem::path& sourceRoot    = {},
    std::string_view             themeName     = "dark");

}  // namespace simall::gui_core
