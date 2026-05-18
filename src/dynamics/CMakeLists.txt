simall_add_module(dynamics
    SOURCES SixDof.cpp SixDofCoupler.cpp
    HEADERS SixDof.hpp SixDofCoupler.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities)
