simall_add_module(amr
    SOURCES AdaptiveRefinement.cpp RefinementApplier.cpp AmrApplier.cpp
    HEADERS AdaptiveRefinement.hpp RefinementApplier.hpp AmrApplier.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Solver)
