# SimAll Beta 1.0.0-rc1 — Release Notes

**Release date:** 2026-05-18
**Status:** Release candidate (gating: CI green on Windows/Linux × Debug/Release × MSVC/GCC/Clang).

## Highlights

* First release candidate marking the end of the 20-week industrial-CFD
  build-out.  All subsystems listed in the master architecture blueprint
  are present, individually unit-tested, and exercised by the regression
  pyramid (regression / verification / golden / performance).
* Public plugin SDK is now frozen at ABI version 1.  Three reference
  plugins ship alongside the runtime as worked examples.
* Six new verification driver applications (`simall_naca`,
  `simall_shocktube`, `simall_rb`, `simall_flameD`, `simall_cylinder`,
  `simall_pipe`) join the existing LDC, BFS, and TGV CLIs.
* Cross-platform CPack packaging refined: Windows NSIS installer plus ZIP,
  macOS DragNDrop plus TGZ, Linux DEB plus TGZ — all sliced into
  `runtime`, `devel`, and `docs` components.

## Supported toolchains

| OS              | Compiler   | Build types                |
|-----------------|------------|----------------------------|
| Windows 11 22H2 | MSVC 19.36 | Debug, Release             |
| Ubuntu 22.04    | GCC 11+    | Release, RelWithDebInfo    |
| Ubuntu 22.04    | Clang 14+  | Release, ASan+UBSan        |
| macOS 13+       | AppleClang | Release (community-tested) |

## Upgrade notes

* Project version is now `1.0.0-rc1`.  Downstream consumers using
  `find_package(SimAll CONFIG)` must allow major version `1`.
* Plugin ABI is locked.  Existing out-of-tree plugins recompiled against
  the v1 SDK header remain compatible until the next ABI bump.

## Known issues

* MPI back-end requires OpenMPI ≥ 4.1 on Linux; older releases will
  configure but trigger runtime deprecation warnings.
* CUDA back-end is opt-in via `-DSIMALL_ENABLE_CUDA=ON`; the default CPU
  fall-back is used by the CI matrix.

## Acknowledgements

Reference data used in the verification suite is reproduced from
peer-reviewed sources; the full attribution list lives in
`docs/validation/report.md`.
