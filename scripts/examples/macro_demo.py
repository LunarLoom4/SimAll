"""macro_demo.py — Synthetic macro replay.

Constructs a small list of `simall.Command` objects and dispatches them
through the CommandBus.  This mirrors what the GUI MacroRecorder writes
out when the user clicks "Save macro".
"""
import simall

bus = simall.CommandBus.instance()

def cmd(name, **kwargs):
    c = simall.Command()
    c.name = name
    c.args = [(k, str(v)) for k, v in kwargs.items()]
    return c

macro = [
    cmd("File.Open",           path="cases/cavity.simall"),
    cmd("Mesh.Surface",        size=0.005, growth=1.2),
    cmd("Mesh.Volume",         algorithm="tetra", quality=0.3),
    cmd("Solver.Setup",        algorithm="SIMPLE", tolerance=1e-5),
    cmd("Run.Start"),
]
for c in macro:
    if bus.has_handler(c.name):
        bus.dispatch(c)
    else:
        simall.log("(no handler) " + c.to_python_call())
