# Plugin Model

A plugin is a shared library (`*.dll` / `*.so` / `*.dylib`) that
exports a single C symbol:

```cpp
extern "C" SIMALL_PLUGIN_EXPORT simall::IPlugin* CreatePlugin();
```

See [plugins/IPlugin.hpp](../../plugins/IPlugin.hpp) for the contract.

## Plugin categories

| Category | Base interface | Registered with |
|---|---|---|
| Turbulence model | `turbulence::ITurbulenceModel` | `turbulence::TurbulenceRegistry` |
| Combustion model | `combustion::ICombustionModel` | `combustion::CombustionRegistry` |
| Material        | `materials::IMaterial`         | `materials::MaterialDatabase` |
| Exporter        | `io::IExporter`                | `io::ExporterRegistry` |
| Postprocess filter | `visualization::IFilter`     | `visualization::FilterRegistry` |
| Script binding  | `scripting::IBinding`          | `scripting::Repl` |

## Lifecycle

1. `core::PluginLoader::scan(<directory>)` enumerates `*.simallplugin`.
2. ABI version is checked against `SIMALL_PLUGIN_ABI`.
3. `CreatePlugin()` is invoked; the returned `IPlugin*` is stored in
   `core::PluginRegistry`.
4. On shutdown the registry calls `IPlugin::shutdown()` and unloads
   the library.

## ABI stability

Bumping the ABI requires changing `SIMALL_PLUGIN_ABI` and rebuilding
all plugins. Major releases are allowed to break ABI; minor releases
must remain compatible.
