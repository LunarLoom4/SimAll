# Scripting

SimAll Beta embeds Python 3.10+ via pybind11. The interpreter is
exposed in the GUI's bottom dock and from the command line:

```powershell
simall_beta --script my_workflow.py
```

## Module map

```python
import simall

project   = simall.Project.open("cavity.simallproj")
mesh      = project.mesh
physics   = project.physics
solver    = project.solver
monitors  = project.monitors
```

## Example — automated parametric sweep

```python
import simall
for Re in (100, 400, 1000, 3200):
    proj = simall.Project.open("cavity.simallproj")
    proj.materials["air"].viscosity = 1.0 / Re
    proj.solver.max_iterations = 5000
    proj.solver.run()
    proj.export_cgns(f"cavity_Re{Re}.cgns")
```

Bindings live under [src/scripting/](../../src/scripting); add a new
binding by deriving from `scripting::IBinding`.
