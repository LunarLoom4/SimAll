simall_add_module(materials
    SOURCES
        MaterialDatabase.cpp
        Janaf.cpp
        SolidsCatalog.cpp
        FluidsCatalog.cpp
    HEADERS
        Material.hpp
        RealGasEos.hpp
        TransportProperties.hpp
        PolynomialFit.hpp
        TabulatedProperty.hpp
        Janaf.hpp
        AnisotropicConductivity.hpp
        SolidsCatalog.hpp
        FluidsCatalog.hpp
        MaterialAssignment.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities)
