simall_add_module(morphing
    SOURCES RbfMorpher.cpp FfdBox.cpp
    HEADERS RbfMorpher.hpp FfdBox.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing)
