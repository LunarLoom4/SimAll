# Plugin ABI Contract

The plugin ABI is defined by `SIMALL_PLUGIN_ABI` (in
[plugins/IPlugin.hpp](../../plugins/IPlugin.hpp)) and the layouts of
all interfaces listed below.

## ABI rules

- ABI is **C-only at the export boundary**. The `CreatePlugin` symbol
  must use `extern "C"` and return `IPlugin*`.
- Polymorphic interfaces are stable inside one major version; their
  vtables MUST NOT have methods inserted or reordered.
- All cross-DLL strings are `const char*` UTF-8 with `'\0'` terminator.
- Containers (`std::vector`, etc.) MUST NOT cross the ABI boundary.
  Use span-style `(T const*, size_t)` pairs.

## Frozen interfaces (ABI v1)

| Interface | Header |
|---|---|
| `IPlugin`                 | `plugins/IPlugin.hpp` |
| `turbulence::ITurbulenceModel` | `src/turbulence/ITurbulenceModel.hpp` |
| `combustion::ICombustionModel` | `src/combustion/ICombustionModel.hpp` |
| `materials::IMaterial`         | `src/materials/IMaterial.hpp` |
| `io::IExporter`                | `src/io/IExporter.hpp` |
| `visualization::IFilter`       | `src/visualization/IFilter.hpp` |
| `scripting::IBinding`          | `src/scripting/IBinding.hpp` |
