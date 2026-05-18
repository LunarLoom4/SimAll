"""lid_driven_cavity.py — Canonical 2-D lid-driven cavity benchmark.

Sets up a 1 m x 1 m square domain with a moving top lid at U=1 m/s,
steady incompressible SIMPLE, Re = 1000 (kinematic viscosity 1e-3).
Drives the run through CommandBus commands so the entire setup is
captured by any active macro recorder.
"""
import simall

ctx = simall.ctx()
bus = simall.CommandBus.instance()

def cmd(name, **kw):
    c = simall.Command()
    c.name = name
    c.args = [(k, str(v)) for k, v in kw.items()]
    return c

# --- Domain & mesh ----------------------------------------------------------
bus.dispatch(cmd("Geometry.CreateBox",
                 xmin=0, ymin=0, zmin=0,
                 xmax=1, ymax=1, zmax=0.01))
bus.dispatch(cmd("Mesh.Surface", size=0.02, growth=1.1))
bus.dispatch(cmd("Mesh.Volume",  algorithm="hex_dominant", quality=0.4))

# --- Material ---------------------------------------------------------------
bus.dispatch(cmd("Material.SetFluid",
                 density=1.0, viscosity=1.0e-3))

# --- Boundary conditions ----------------------------------------------------
bus.dispatch(cmd("BC.Set", zone="top",    kind="wall", u=1.0, v=0.0, w=0.0))
bus.dispatch(cmd("BC.Set", zone="bottom", kind="wall"))
bus.dispatch(cmd("BC.Set", zone="left",   kind="wall"))
bus.dispatch(cmd("BC.Set", zone="right",  kind="wall"))
bus.dispatch(cmd("BC.Set", zone="front",  kind="symmetry"))
bus.dispatch(cmd("BC.Set", zone="back",   kind="symmetry"))

# --- Solver -----------------------------------------------------------------
bus.dispatch(cmd("Solver.Setup",
                 algorithm="SIMPLE", relax_p=0.3, relax_u=0.7,
                 tolerance=1e-6, max_iter=2000))

ctx.set("case.name", "lid_driven_cavity_Re1000")
bus.dispatch(cmd("Run.Start"))
