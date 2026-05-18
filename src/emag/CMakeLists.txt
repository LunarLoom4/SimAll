simall_add_module(emag
    SOURCES Mhd.cpp AvFormulation.cpp
    HEADERS Mhd.hpp AvFormulation.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Solver)
