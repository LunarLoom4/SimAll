"""batch_run.py — Parametric Reynolds sweep over the lid-driven cavity.

Demonstrates programmatic case generation: for each Reynolds number, set
the viscosity, run the solver, export results to a numbered folder.
"""
import os
import simall

bus  = simall.CommandBus.instance()
ctx  = simall.ctx()

def cmd(name, **kw):
    c = simall.Command()
    c.name = name
    c.args = [(k, str(v)) for k, v in kw.items()]
    return c

REYNOLDS_LIST = [100, 400, 1000, 3200]
OUT_ROOT      = os.path.join(ctx.get("project.path") or ".", "studies", "Re_sweep")

for Re in REYNOLDS_LIST:
    nu  = 1.0 / Re                       # L=1, U=1 -> Re = 1/nu
    out = os.path.join(OUT_ROOT, "Re_%d" % Re)
    simall.log("=== Re = %d (nu = %.3e) ===" % (Re, nu))

    bus.dispatch(cmd("Material.SetFluid", density=1.0, viscosity=nu))
    bus.dispatch(cmd("Solver.Setup",      algorithm="SIMPLE", tolerance=1e-6, max_iter=4000))
    bus.dispatch(cmd("Run.Start"))
    bus.dispatch(cmd("Results.ExportEnsight", path=out))

simall.log("Re sweep complete: %d cases." % len(REYNOLDS_LIST))
