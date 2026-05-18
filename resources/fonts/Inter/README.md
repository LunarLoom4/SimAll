# Inter font

SimAll's GUI ships with the **Inter** typeface as the preferred application
font.  Drop the official OTF/TTF files into this directory and they will be
auto-loaded at startup by `simall::gui::ThemeManager::load_inter_font`.

Expected filenames (any subset is fine — the GUI loads everything matching
`*.otf` or `*.ttf` from this folder):

    Inter-Thin.otf
    Inter-ExtraLight.otf
    Inter-Light.otf
    Inter-Regular.otf
    Inter-Medium.otf
    Inter-SemiBold.otf
    Inter-Bold.otf
    Inter-ExtraBold.otf
    Inter-Black.otf

Download: https://rsms.me/inter/   (SIL OFL 1.1 — redistribution permitted).

If the directory is empty, the GUI falls back to the platform-native UI
sans (Segoe UI on Windows, .AppleSystemUIFont on macOS, Noto Sans on Linux).
