# SimAll Beta — Documentation

Welcome. The documentation is organised as:

- **Architecture** — strict layer graph, threading model, plugin model.
- **Pipelines** — CAD → Mesh → Solver → Post-processing.
- **Theory** — finite-volume formulation and each physics module.
- **User Guide** — quick-start, workflow tree, BCs, monitors, scripting.
- **Plugin SDK** — authoring, ABI, examples.
- **Validation** — canonical test cases with tolerances.
- **API reference** — Doxygen-generated browsable API.

Build locally:

```powershell
pip install mkdocs-material mkdoxy
doxygen docs/Doxyfile
mkdocs serve -f docs/mkdocs.yml
```

