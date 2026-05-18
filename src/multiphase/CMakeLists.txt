simall_add_module(multiphase
    SOURCES Multiphase.cpp VolumeOfFluid.cpp PlicReconstruction.cpp CavitationSchnerrSauer.cpp PopulationBalance.cpp SurfaceTensionCsf.cpp LeeMassTransfer.cpp CavitationKunz.cpp CavitationZwart.cpp MixtureModel.cpp DriftFlux.cpp EulerianMultifluid.cpp IsoAdvector.cpp RpiWallBoiling.cpp
    HEADERS Multiphase.hpp VolumeOfFluid.hpp PlicReconstruction.hpp CavitationSchnerrSauer.hpp PopulationBalance.hpp SurfaceTensionCsf.hpp LeeMassTransfer.hpp CavitationKunz.hpp CavitationZwart.hpp MixtureModel.hpp DriftFlux.hpp EulerianMultifluid.hpp IsoAdvector.hpp RpiWallBoiling.hpp
    PUBLIC_DEPS SimAll::Solver)
