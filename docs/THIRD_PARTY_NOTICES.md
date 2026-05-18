# Third-party notices

SimAll Beta links against and bundles several third-party components.
This file enumerates the runtime / build-time dependencies along with
their upstream licenses.

| Component | Version | License | Role |
|-----------|---------|---------|------|
| Eigen3            | ≥ 3.4    | MPL-2.0  | Dense linear algebra |
| Intel TBB         | ≥ 2021   | Apache-2.0 | Task scheduler (optional) |
| OpenMP            | ≥ 4.5    | (compiler) | Shared-memory parallelism (optional) |
| OpenMPI / MS-MPI  | ≥ 4.1    | BSD-3-Clause / MS-EULA | Distributed-memory parallelism (optional) |
| CUDA Toolkit      | ≥ 12.0   | NVIDIA EULA | GPU back-end (optional) |
| VTK               | ≥ 9.2    | BSD-3-Clause | Rendering / mesh data structures |
| OpenCASCADE       | ≥ 7.6    | LGPL-2.1 | CAD kernel |
| Qt6               | ≥ 6.5    | LGPL-3.0 | GUI shell |
| spdlog            | ≥ 1.12   | MIT | Logger backend |
| pybind11          | ≥ 2.11   | BSD-3-Clause | Python bindings (optional) |
| Catch2            | 3.5.4    | BSL-1.0 | Unit-test framework |
| HDF5              | ≥ 1.12   | BSD-style | Result store back-end (optional) |
| CGNS              | ≥ 4.3    | zlib-style | Native CFD interchange (optional) |

Reference data shipped under `tests/regression/cases/` is reproduced from
peer-reviewed publications; full citations live in
`docs/validation/report.md`.
