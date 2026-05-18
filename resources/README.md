# Resource Bundle

| Path | Purpose |
|---|---|
| `simall_resources.qrc` | Qt resource manifest (compiled into the app) |
| `themes/dark.qss`      | dark-mode stylesheet (spec §3.6) |
| `icons/`               | UI icons (SVG, scale-clean) |
| `splash/splash.svg`    | startup splash (1024 × 640) |
| `fonts/`               | bundled fonts (see below) |

## Fonts

The default UI font is **Inter** with **Segoe UI** as fallback
(spec §3.7). Inter is OFL-licensed and may be vendored under
`fonts/inter/`. To avoid blowing up the repo size, the fonts are
fetched on first build by `scripts/bootstrap_fonts.ps1` from
the official `rsms/inter` release tarball.

If Inter is unavailable at runtime, Qt automatically falls back to
the next entry in the font stack defined in `dark.qss`.

## Icon style

All SVGs use the dark-theme palette and a 1.2–1.4 px stroke at
16 × 16 baseline. App icon (`icons/app.svg`) is exported to multi-
size `.ico` (`scripts/build_icons.ps1`) for Windows packaging.

## Compiling resources

`simall_resources.qrc` is referenced from `src/gui/CMakeLists.txt`
via `qt_add_resources(simall_gui simall_resources PREFIX "/" FILES …)`.
