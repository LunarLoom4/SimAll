simall_add_module(core
    SOURCES
        Application.cpp
        ThreadPool.cpp
        TaskGraph.cpp
        PluginLoader.cpp
        PluginRegistry.cpp
        ConfigStore.cpp
        CrashHandler.cpp
        ResourceMonitor.cpp
        UnitSystem.cpp
        RotatingFileSink.cpp
        BoundaryConditionRegistry.cpp
    HEADERS
        Application.hpp
        Command.hpp
        ConfigStore.hpp
        CrashHandler.hpp
        EventBus.hpp
        IProject.hpp
        Logger.hpp
        PluginLoader.hpp
        PluginRegistry.hpp
        ResourceMonitor.hpp
        RotatingFileSink.hpp
        ServiceLocator.hpp
        TaskGraph.hpp
        ThreadPool.hpp
        UnitSystem.hpp
        BoundaryConditionRegistry.hpp
    PUBLIC_DEPS
        SimAll::Utilities)

# Platform-specific link requirements for runtime introspection libs.
if(WIN32)
    target_link_libraries(simall_core PRIVATE dbghelp Psapi)
else()
    find_package(Threads REQUIRED)
    target_link_libraries(simall_core PUBLIC Threads::Threads ${CMAKE_DL_LIBS})
endif()

# `core/PluginLoader.hpp` and `core/PluginRegistry.hpp` include
# "plugins/IPlugin.hpp" / "plugins/SimAllPluginSdk.hpp", which live at the
# repository root (outside src/).  Expose the root as a public include
# directory so any consumer that pulls in core's plugin headers can resolve
# them without having to know the plugins folder lives one level up.
target_include_directories(simall_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/simall>)
