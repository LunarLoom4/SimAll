# Plugin Examples

Reference plugins live under [plugins/](../../plugins). Each
demonstrates a single extension surface.

| Folder | Demonstrates |
|---|---|
| `plugins/example_turbulence`  | custom one-equation turbulence model |
| `plugins/example_material`    | temperature-dependent viscosity material |
| `plugins/example_exporter`    | export to a custom XML format |
| `plugins/example_filter`      | post-processing filter (Q-criterion) |
| `plugins/example_binding`     | Python binding for a custom field |

To build the examples, set `-DSIMALL_BUILD_PLUGINS=ON` (default).
