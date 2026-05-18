// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/ThemeManager.hpp
//
// Centralised Qt theme application:
//   * Loads the Inter font family from disk if available (otherwise the
//     editor uses the platform default sans).
//   * Installs the SimAll dark palette.
//   * Applies the QSS stylesheet from resources/styles/simall_dark.qss.
//   * Exposes accent colours so dialogs can pick up the same hues without
//     hard-coding hex strings.
// =============================================================================
#pragma once

#include <QColor>
#include <QString>
#include <vector>

namespace simall::gui {

enum class Theme { Dark, Light, HighContrast };

struct ThemeAccents {
    QColor primary    {0x00, 0x7A, 0xCC};
    QColor success    {0x16, 0x9F, 0x6B};
    QColor warning    {0xE0, 0xA6, 0x00};
    QColor danger     {0xE5, 0x2B, 0x50};
    QColor surface    {0x1E, 0x1E, 0x1E};
    QColor surface2   {0x25, 0x25, 0x26};
    QColor surface3   {0x33, 0x33, 0x34};
    QColor text       {0xE8, 0xE8, 0xE8};
    QColor textMuted  {0x9A, 0x9A, 0x9A};
};

class ThemeManager {
public:
    static ThemeManager& instance();

    // Apply the theme to QApplication: palette, font, stylesheet.
    // `sourceRoot` is passed to ThemePaths so dev builds find assets in-tree.
    void apply(Theme theme, const QString& sourceRoot = {});

    [[nodiscard]] Theme              current() const noexcept { return current_; }
    [[nodiscard]] const ThemeAccents& accents() const noexcept { return accents_; }
    [[nodiscard]] QString             font_family() const noexcept { return family_; }
    [[nodiscard]] QStringList         search_log() const noexcept { return searchLog_; }

private:
    ThemeManager() = default;
    void apply_dark_palette();
    void apply_light_palette();
    void apply_high_contrast_palette();
    void load_inter_font(const QString& sourceRoot);
    void load_stylesheet(Theme t, const QString& sourceRoot);

    Theme         current_ = Theme::Dark;
    ThemeAccents  accents_;
    QString       family_  = "Segoe UI";
    QStringList   searchLog_;
};

}  // namespace simall::gui
