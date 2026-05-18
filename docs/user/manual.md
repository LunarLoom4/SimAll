# SimAll Beta — User Manual

> Version 1.0.0-rc1.  This manual targets end users of the SimAll Beta
> industrial CFD/CAE platform.  Developers should consult the SDK guide
> under `docs/sdk/` and the architecture overview under `docs/architecture/`.

## 1. Installation

### Windows
1. Download `SimAll-Beta-1.0.0-rc1-win64.exe` (NSIS installer) or
   `SimAll-Beta-1.0.0-rc1-win64.zip` (portable archive).
2. Run the installer; choose `Runtime`, `SDK / Development`, and
   `Documentation` components as needed.
3. The installer adds `bin/` to the PATH on request.

### Linux (DEB)
```bash
sudo dpkg -i simall-beta_1.0.0-rc1_amd64.deb
sudo apt-get install -f      # resolve runtime depends
```

### macOS (DragNDrop)
Drag `SimAll Beta.app` from the mounted disk image into `/Applications`.

## 2. Quick start

Run the lid-driven cavity acceptance driver:

```bash
simall_cavity --re 100 --nx 64 --ny 64 --out cavity_Re100.vtk
```

Open `cavity_Re100.vtk` in the bundled GUI (`SimAll`) or ParaView to
inspect the velocity / pressure field.

## 3. Verification driver suite

The following CLIs ship in the runtime component and emit the canonical
reference values used to gate the solver:

| Tool                | Reference benchmark           |
|---------------------|-------------------------------|
| `simall_cavity`     | Ghia 1982 lid-driven cavity   |
| `simall_taylor_green` | TGV analytic decay          |
| `simall_bfs`        | Armaly 1983 BFS               |
| `simall_naca`       | Thin-airfoil-theory NACA-0012 |
| `simall_shocktube`  | Sod exact Riemann solution    |
| `simall_rb`         | Rayleigh-Bénard onset         |
| `simall_flameD`     | Sandia Flame D centreline     |
| `simall_cylinder`   | Roshko + Henderson cylinder   |
| `simall_pipe`       | Hagen-Poiseuille pipe profile |

Each tool prints its acceptance tolerance so operators can scope-check
the build before launching a long production run.

## 4. GUI workflow

1. Launch `SimAll`.
2. Use the **Project** ribbon to create a new project (`.simall` file).
3. Import geometry via **CAD → Import** (STEP/IGES/STL).
4. Mesh via **Meshing → Generate** (Cartesian / hex-sweep / polyhedral).
5. Configure physics under **Solver**, **Turbulence**, **Multiphase**,
   etc.  The 13-tab ribbon exposes one panel per subsystem.
6. Run via **Solve → Start** and monitor residuals in the
   **Solver Monitor** dock.
7. Post-process with **Visualization → Filters** (contour, streamline,
   isosurface, LIC, volume raycast).
8. Save results to `.simall` (native) or export to CGNS / Tecplot / VTK.

## 5. Scripting

```python
import simall as sa
prj = sa.Project.open("cavity.simall")
prj.solver.simple.run(steps=500)
prj.save()
```

The Python REPL is also available via **View → Python Console** inside
the GUI.

## 6. Plugins

Drop a third-party `.dll` / `.so` / `.dylib` into
`<install>/lib/simall/plugins/` and restart SimAll.  Plugin ABI is
v1 — plugins compiled against the v1 SDK remain compatible until the
next major release.

## 7. Support

* Documentation portal: <https://example.org/simall/docs>
* Issue tracker:        <https://example.org/simall/issues>
* Mailing list:         `simall-users@example.org`
