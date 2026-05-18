simall_add_module(rotating
    SOURCES MultipleReferenceFrame.cpp SlidingMeshDriver.cpp
    HEADERS MultipleReferenceFrame.hpp SlidingMeshDriver.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Solver)
