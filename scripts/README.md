# Helper Scripts

| Script | Purpose |
|---|---|
| `build.ps1` / `build.sh`      | one-shot configure + build wrapper |
| `format.ps1`                  | clang-format the entire tree (or `-Check` for CI) |
| `run_regression.py`           | ctest driver + JUnit converter |
| `bootstrap_third_party.ps1`   | `git submodule update --init --recursive --depth 1` |
| `bootstrap_fonts.ps1`         | fetch the OFL Inter font release |
| `regen_protos.ps1`            | regenerate `*.pb.h/.pb.cc` from `src/io/proto/` |
| `verify_data.ps1`             | check `data/manifest.json` SHA-256 hashes (TODO) |
| `check_dep_graph.py`          | enforce module dependency order (TODO) |
| `check_includes.py`           | enforce forbidden-includes rules (TODO) |
