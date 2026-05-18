simall_add_module(ibm
    SOURCES ImmersedBoundary.cpp DirectForcingIbm.cpp
    HEADERS ImmersedBoundary.hpp DirectForcingIbm.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Solver)
