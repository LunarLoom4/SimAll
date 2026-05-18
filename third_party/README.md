# Third-Party Dependencies

SimAll Beta links a fixed set of third-party libraries. Each is
vendored either as a `git submodule` under this directory or
fetched on-demand by `FetchContent`. **Never** copy source into the
tree by hand — always vendor through one of these two channels so
licence audit and version bumps remain mechanical.

## Pinned versions

| Library  | Version | Channel        | License | Purpose |
|---|---|---|---|---|
| Eigen           | 3.4.0    | FetchContent | MPL 2.0     | linear algebra primitives (required) |
| Catch2          | v3.5.4   | FetchContent | BSL 1.0     | unit-test framework |
| fmt             | 10.2.1   | FetchContent | MIT         | formatting / logging |
| spdlog          | 1.13.0   | FetchContent | MIT         | structured logging |
| pybind11        | 2.12.0   | FetchContent | BSD 3-clause | Python bindings |
| METIS           | 5.1.0    | FetchContent | Apache 2.0  | mesh partitioner |
| nlohmann/json   | 3.11.3   | FetchContent | MIT         | project file I/O |
| VTK             | 9.3.x    | system / superbuild | BSD 3-clause | visualization |
| Qt              | 6.5+     | system       | LGPL v3 / commercial | GUI |
| OpenCASCADE     | 7.7.x    | system / superbuild | LGPL v2.1 + exception | CAD kernel |
| MPI             | OpenMPI 4 / MS-MPI | system | as upstream | distributed memory |
| TBB             | 2021.x   | system       | Apache 2.0  | threading |
| CUDA Toolkit    | 12.x     | system       | NVIDIA EULA | GPU acceleration |
| PETSc           | 3.20.x   | system       | BSD 2-clause | optional linear solvers |

## How to vendor

1. Add a new entry to `submodules.cmake` (FetchContent) **or** add a
   git submodule under `third_party/<name>` and link it via
   `add_subdirectory(third_party/<name>)`.
2. Update the **License** column above and add the upstream `LICENSE`
   text under `third_party/<name>/LICENSE`.
3. Bump the pinned version on a single line in `submodules.cmake`.
4. Run `scripts/bootstrap_third_party.ps1` to sync.

## Licence audit

A nightly CI job runs `scripts/license_audit.py` (TODO) which walks
each vendored module's `LICENSE` file and verifies it matches the
table above.
