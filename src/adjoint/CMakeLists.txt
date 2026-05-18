simall_add_module(adjoint
    SOURCES Adjoint.cpp DiscreteAdjoint.cpp ContinuousAdjoint.cpp
    HEADERS Adjoint.hpp DiscreteAdjoint.hpp ContinuousAdjoint.hpp
    PUBLIC_DEPS SimAll::Solver)
