# Changelog

All notable changes to SimAll Beta are documented here.  The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project adheres
to [Semantic Versioning](https://semver.org/).

## [1.0.0-rc1] — 2026-05-18

### Added — Week 20 (Plugins, packaging, applications, docs)
- Header-only plugin SDK (`plugins/SimAllPluginSdk.hpp`) with `PluginBase`,
  `CategorisedPlugin`, and the `SIMALL_DECLARE_PLUGIN(Class)` macro.
- Three reference plugins: `simall_curve_postproc`, `simall_csv_writer`,
  `simall_hello_solver`, installed under `lib/simall/plugins`.
- Six new verification driver applications: `simall_naca`, `simall_shocktube`,
  `simall_rb`, `simall_flameD`, `simall_cylinder`, `simall_pipe`.
- LICENSE (Apache-2.0), CHANGELOG, and RELEASE_NOTES_1.0-rc1.md.
- User manual (`docs/user/manual.md`), theory guide
  (`docs/theory/governing_equations.md`), validation report
  (`docs/validation/report.md`), and third-party notices.
- CPack `LICENSE`, `RESOURCES`, and component-level descriptions wired
  through the existing `simall_setup_cpack` helper.

### Added — Week 19 (Validation, regression, perf, CI)
- `tests/regression/` with the 10 canonical CFD benchmark headers (LDC,
  BFS, channel, TGV, NACA, shocktube, RB, FlameD, cylinder, pipe).
- `tests/verification/` MMS suite with FD5 Laplacian convergence proof.
- `tests/performance/` micro-benchmark harness (SAXPY / dot / matvec).
- `tests/golden/` FNV-1a hash gates against quantised reference vectors.
- CI matrix extended with ASan+UBSan job and ctest label-filtered passes.

### Added — Week 18 (Long-tail modules)
- Adjoint (discrete-GMRES on Aᵀ + continuous shape sensitivity).
- Optimization: MMA, CMA-ES, full-factorial / LHS / Sobol DOE, Kriging.
- Morphing: trivariate Bernstein-Bezier FFD lattice.
- Dynamics: 6-DOF coupler with explicit / loose-implicit / under-relaxed modes.
- Rotating: sliding-mesh driver with inverse-distance donor weighting.
- AMR: rank-based mark+enforce orchestrator.
- IBM: Mohd-Yusof / Uhlmann direct forcing with Roma 3-pt δ-kernel.
- EMag: linear-tet A-V element matrices + CSR assembler.
- Acoustics: Curle dipole integrator (complements existing FW-H).
- ROM: operator inference + ECSW NNLS hyper-reduction.

### Added — Weeks 1–17
See `docs/CHANGELOG_PREVIOUS.md` for the per-week itemisation.

### Changed
- Project version bumped from 0.1.0 to 1.0.0-rc1.

### Known limitations
- GPU back-end requires CUDA 12.x; `SIMALL_ENABLE_CUDA=OFF` falls back to
  the CPU SpMV path.
- MPI is optional; serial fall-back paths are exercised in CI.
- VTK 9.2+ is required for the GUI; headless builds drop the
  `SimAll::Visualization` target.
