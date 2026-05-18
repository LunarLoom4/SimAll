simall_add_module(optimization
    SOURCES Optimization.cpp SimpTopology.cpp
            MmaOptimizer.cpp CmaEs.cpp Doe.cpp Kriging.cpp
    HEADERS Optimization.hpp SimpTopology.hpp
            MmaOptimizer.hpp CmaEs.hpp Doe.hpp Kriging.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Meshing SimAll::Solver)
