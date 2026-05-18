# Layer Dependency Graph

The dependency graph below is enforced by [CMakeLists.txt](../../CMakeLists.txt)
ordering. Each subsystem may depend ONLY on subsystems declared earlier.

```text
utilities → core → io → parallel → gpu → materials → cad → meshing
        → solver → turbulence
                → heat_transfer
                → radiation
                → combustion
                → multiphase
                → particles
                → morphing → ibm → porous → rotating → dynamics
                → acoustics → emag → rom → zones → amr
                → scripting → optimization → adjoint
        → visualization → gui
```

A static check (`scripts/check_dep_graph.py`) parses
`CMakeLists.txt` and asserts no edge goes backwards.

## Hot-path forbidden inclusions

The following includes are FORBIDDEN inside any solver, meshing, or
turbulence translation unit:

- `<QtWidgets>`, `<QtGui>`, `<QtCore>` — GUI primitives never leak.
- `<vtk*.h>`  except where visualization is the explicit purpose.
- `<TopoDS.hxx>` and other OCCT headers — physics modules must not
  see CAD types; meshing wraps them through `cad/CadKernel`.

`scripts/check_includes.py` enforces this at CI time.
