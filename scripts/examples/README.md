# SimAll example scripts

These scripts run inside the SimAll embedded Python interpreter (Python
3.9+ via pybind11).  They are also valid standalone Python files for
linting purposes — they import the `simall` module which is only present
inside the SimAll process.

| File | Demonstrates |
|------|-------------|
| `hello.py`                  | Basic logging and context-bag access. |
| `udf_velocity_profile.py`   | Registering a UDF that returns a per-point velocity vector. |
| `macro_demo.py`             | Replaying a recorded macro to set up a job programmatically. |
| `lid_driven_cavity.py`      | End-to-end script for the canonical 2-D lid-driven cavity benchmark. |
| `batch_run.py`              | Loop over a parametric design study (Re sweep). |

Run from inside SimAll with:

```python
exec(open(r"scripts/examples/lid_driven_cavity.py").read())
```

or from the command line via the headless harness:

```
simall-cli --script scripts/examples/lid_driven_cavity.py
```
