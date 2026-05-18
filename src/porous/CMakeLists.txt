simall_add_module(porous
    SOURCES PorousMedia.cpp
    HEADERS PorousMedia.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Solver)
