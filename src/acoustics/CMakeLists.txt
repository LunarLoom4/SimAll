simall_add_module(acoustics
    SOURCES FwhSurface.cpp CurleAcoustics.cpp
    HEADERS FwhSurface.hpp CurleAcoustics.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Solver)
