// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/ThemeManager.cpp
// =============================================================================
#include "gui/ThemeManager.hpp"
#include "gui_core/ThemePaths.hpp"
#include "core/Logger.hpp"

#include <QApplication>
#include <QFile>
#include <QFont>
#include <QFontDatabase>
#include <QPalette>
#include <QTextStream>

namespace simall::gui {

ThemeManager& ThemeManager::instance() {
    static ThemeManager s;
    return s;
}

void ThemeManager::apply(Theme theme, const QString& sourceRoot) {
    current_ = theme;
    switch (theme) {
        case Theme::Dark:          apply_dark_palette();          break;
        case Theme::Light:         apply_light_palette();         break;
        case Theme::HighContrast:  apply_high_contrast_palette(); break;
    }
    load_inter_font(sourceRoot);
    load_stylesheet(theme, sourceRoot);
}

void ThemeManager::apply_dark_palette() {
    accents_ = ThemeAccents{};   // defaults already dark
    QPalette pal;
    pal.setColor(QPalette::Window,          accents_.surface);
    pal.setColor(QPalette::WindowText,      accents_.text);
    pal.setColor(QPalette::Base,            accents_.surface2);
    pal.setColor(QPalette::AlternateBase,   QColor(0x2D, 0x2D, 0x30));
    pal.setColor(QPalette::Text,            accents_.text);
    pal.setColor(QPalette::Button,          accents_.surface3);
    pal.setColor(QPalette::ButtonText,      accents_.text);
    pal.setColor(QPalette::ToolTipBase,     QColor(0x25, 0x25, 0x26));
    pal.setColor(QPalette::ToolTipText,     accents_.text);
    pal.setColor(QPalette::Highlight,       accents_.primary);
    pal.setColor(QPalette::HighlightedText, Qt::white);
    pal.setColor(QPalette::Link,            accents_.primary);
    pal.setColor(QPalette::Disabled, QPalette::Text,       accents_.textMuted);
    pal.setColor(QPalette::Disabled, QPalette::ButtonText, accents_.textMuted);
    qApp->setPalette(pal);
}

void ThemeManager::apply_light_palette() {
    accents_.surface   = QColor(0xFA, 0xFA, 0xFA);
    accents_.surface2  = QColor(0xFF, 0xFF, 0xFF);
    accents_.surface3  = QColor(0xEC, 0xEC, 0xEC);
    accents_.text      = QColor(0x1B, 0x1B, 0x1B);
    accents_.textMuted = QColor(0x6E, 0x6E, 0x6E);
    QPalette pal;
    pal.setColor(QPalette::Window,        accents_.surface);
    pal.setColor(QPalette::WindowText,    accents_.text);
    pal.setColor(QPalette::Base,          accents_.surface2);
    pal.setColor(QPalette::Text,          accents_.text);
    pal.setColor(QPalette::Button,        accents_.surface3);
    pal.setColor(QPalette::ButtonText,    accents_.text);
    pal.setColor(QPalette::Highlight,     accents_.primary);
    pal.setColor(QPalette::HighlightedText, Qt::white);
    qApp->setPalette(pal);
}

void ThemeManager::apply_high_contrast_palette() {
    accents_.surface   = Qt::black;
    accents_.surface2  = Qt::black;
    accents_.surface3  = QColor(0x20, 0x20, 0x20);
    accents_.text      = Qt::white;
    accents_.textMuted = QColor(0xCC, 0xCC, 0xCC);
    accents_.primary   = QColor(0xFF, 0xD7, 0x00);
    QPalette pal;
    pal.setColor(QPalette::Window,        accents_.surface);
    pal.setColor(QPalette::WindowText,    accents_.text);
    pal.setColor(QPalette::Base,          accents_.surface2);
    pal.setColor(QPalette::Text,          accents_.text);
    pal.setColor(QPalette::Button,        accents_.surface3);
    pal.setColor(QPalette::ButtonText,    accents_.text);
    pal.setColor(QPalette::Highlight,     accents_.primary);
    pal.setColor(QPalette::HighlightedText, Qt::black);
    qApp->setPalette(pal);
}

void ThemeManager::load_inter_font(const QString& sourceRoot) {
    const auto loc = gui_core::locate_theme_assets(
        QCoreApplication::applicationDirPath().toStdString(),
        sourceRoot.toStdString(),
        current_ == Theme::Dark ? "dark"
                                : current_ == Theme::Light ? "light" : "highcontrast");

    searchLog_.clear();
    for (auto& p : loc.searchedDirectories)
        searchLog_ << QString::fromStdString(p.string());

    bool loaded = false;
    for (auto& f : loc.fonts) {
        int id = QFontDatabase::addApplicationFont(QString::fromStdString(f.string()));
        if (id >= 0) {
            const QStringList fams = QFontDatabase::applicationFontFamilies(id);
            if (!fams.isEmpty()) { family_ = fams.first(); loaded = true; }
        }
    }
    if (!loaded) {
        // Inter not deployed — fall back to platform-native UI sans.
        family_ =
#ifdef Q_OS_WIN
            "Segoe UI"
#elif defined(Q_OS_MAC)
            ".AppleSystemUIFont"
#else
            "Noto Sans"
#endif
            ;
    }
    QFont f(family_, 9);
    qApp->setFont(f);
    SIMALL_LOG_INFO("Theme", "Applied font family: ", family_.toStdString(),
                    " (Inter ", (loaded ? "loaded" : "missing — fallback"), ")");
}

void ThemeManager::load_stylesheet(Theme t, const QString& sourceRoot) {
    const char* name = t == Theme::Dark         ? "dark"
                     : t == Theme::Light        ? "light"
                                                : "highcontrast";
    const auto loc = gui_core::locate_theme_assets(
        QCoreApplication::applicationDirPath().toStdString(),
        sourceRoot.toStdString(), name);
    if (!loc.stylesheet) {
        SIMALL_LOG_INFO("Theme", "Stylesheet '", name,
                        "' not found — palette-only mode");
        qApp->setStyleSheet({});
        return;
    }
    QFile f(QString::fromStdString(loc.stylesheet->string()));
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
        SIMALL_LOG_INFO("Theme", "Failed to open stylesheet ",
                        loc.stylesheet->string());
        return;
    }
    QTextStream in(&f);
    qApp->setStyleSheet(in.readAll());
    SIMALL_LOG_INFO("Theme", "Loaded stylesheet: ", loc.stylesheet->string());
}

}  // namespace simall::gui
