# Plugin SDK — Authoring a Plugin

A plugin is built as a `SHARED` CMake target that links against
`simall::core`, `simall::solver`, and any subsystems it extends.

## Minimal CMake

```cmake
add_library(myplugin SHARED MyTurbModel.cpp)
target_link_libraries(myplugin PRIVATE simall::core simall::turbulence)
set_target_properties(myplugin PROPERTIES
  PREFIX ""
  SUFFIX ".simallplugin")
install(TARGETS myplugin DESTINATION plugins)
```

## Skeleton

```cpp
#include <simall/plugins/IPlugin.hpp>
#include <simall/turbulence/ITurbulenceModel.hpp>

class MyTurb : public simall::turbulence::ITurbulenceModel {
    // ... initialise, advance, eddyViscosity ...
};

class MyPlugin : public simall::IPlugin {
public:
    const char* name()    const override { return "MyTurb"; }
    const char* version() const override { return "1.0.0"; }
    int         abi()     const override { return SIMALL_PLUGIN_ABI; }
    bool initialize(simall::core::IServiceLocator& sl) override {
        sl.get<simall::turbulence::TurbulenceRegistry>()
          .registerModel<MyTurb>("my-turb");
        return true;
    }
    void shutdown() override {}
};

extern "C" SIMALL_PLUGIN_EXPORT simall::IPlugin* CreatePlugin() {
    return new MyPlugin();
}
```

## Loading

Drop `myplugin.simallplugin` into `<install>/plugins/`. The
`PluginLoader` auto-discovers and registers it on startup.
