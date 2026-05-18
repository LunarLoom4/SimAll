simall_add_module(zones
    SOURCES FanHexZone.cpp SelectionManager.cpp
    HEADERS FanHexZone.hpp NamedSelection.hpp SelectionManager.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Solver)
