set(_simall_scripting_sources
    ScriptEngine.cpp
    ScriptContext.cpp
    Repl.cpp
    CommandBus.cpp
    MacroRecorder.cpp
    UdfHost.cpp
    PyBindings.cpp
    PythonEngine.cpp)

set(_simall_scripting_headers
    ScriptEngine.hpp
    ScriptContext.hpp
    Repl.hpp
    CommandBus.hpp
    MacroRecorder.hpp
    UdfHost.hpp
    PyBindings.hpp
    PythonEngine.hpp)

simall_add_module(scripting
    SOURCES     ${_simall_scripting_sources}
    HEADERS     ${_simall_scripting_headers}
    PUBLIC_DEPS SimAll::Core SimAll::Utilities)

# Embed Python only when both Python3 + pybind11 resolved cleanly.  The
# stub branch in PyBindings.cpp remains compiled otherwise so the library
# is always link-complete.
if(SIMALL_HAVE_PYTHON AND TARGET pybind11::embed)
    target_compile_definitions(simall_scripting PUBLIC SIMALL_HAVE_PYTHON=1)
    target_link_libraries(simall_scripting PUBLIC pybind11::embed)
endif()
