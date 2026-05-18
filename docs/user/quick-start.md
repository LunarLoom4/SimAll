# Quick Start

## Prerequisites

| Platform | Compiler | CMake | Qt | VTK | OpenCASCADE |
|---|---|---|---|---|---|
| Windows 10/11 x64 | Visual Studio 2022 (17.6+) | 3.22+ | 6.5+ | 9.2+ | 7.6+ |
| Ubuntu 22.04 x64  | gcc 11 / clang 14          | 3.22+ | 6.5+ | 9.2+ | 7.6+ |
| Red Hat 9 x64     | gcc 12 / icpx               | 3.22+ | 6.5+ | 9.2+ | 7.6+ |

Install optional accelerators as needed: MPI, OpenMP, TBB, CUDA 12+,
PETSc 3.20+.

## Build

```powershell
git clone https://example.com/SimAll_Beta.git
cd SimAll_Beta
./scripts/bootstrap_third_party.ps1   # fetch Eigen, Catch2, fmt, spdlog, ...
cmake --preset default
cmake --build build --config Release -j
ctest --preset default --output-on-failure
```

On Linux:

```bash
./scripts/build.sh linux-release
ctest --preset linux-release
```

## Run a tutorial

```powershell
./build/bin/simall_cavity   # lid-driven cavity Re=1000
./build/bin/simall_bfs      # backward-facing step
./build/bin/simall_taylor_green  # 3-D Taylor–Green vortex Re=1600
./build/bin/simall_beta     # full GUI
```

## Open the GUI

`simall_beta` launches the Qt GUI. Use **File → Open Sample** to
load `data/stl/unit_cube.stl` and explore the ribbon.
