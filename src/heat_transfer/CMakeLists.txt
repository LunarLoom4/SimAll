simall_add_module(heat_transfer
    SOURCES HeatTransfer.cpp ConjugateHeatTransfer.cpp Solidification.cpp
    HEADERS HeatTransfer.hpp ConjugateHeatTransfer.hpp Solidification.hpp
    PUBLIC_DEPS SimAll::Solver SimAll::Meshing)
