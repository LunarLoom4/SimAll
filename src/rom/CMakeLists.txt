simall_add_module(rom
    SOURCES PodDmd.cpp OperatorInference.cpp EcswHyperReduction.cpp
    HEADERS PodDmd.hpp OperatorInference.hpp EcswHyperReduction.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities)
