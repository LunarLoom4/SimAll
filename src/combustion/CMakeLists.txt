simall_add_module(combustion
    SOURCES Combustion.cpp LaminarFiniteRate.cpp EddyDissipationConcept.cpp
            ChemkinParser.cpp CanteraParser.cpp ThermoNasaParser.cpp FgmFlamelet.cpp SlfmFlamelet.cpp
            TransportedPdf.cpp GEquation.cpp TfcModel.cpp
            NoxThermalPromptFuel.cpp SootMossBrookes.cpp SootMomMethod.cpp
    HEADERS Combustion.hpp LaminarFiniteRate.hpp EddyDissipationConcept.hpp
            ChemkinParser.hpp CanteraParser.hpp ThermoNasaParser.hpp FgmFlamelet.hpp SlfmFlamelet.hpp
            TransportedPdf.hpp GEquation.hpp TfcModel.hpp
            NoxThermalPromptFuel.hpp SootMossBrookes.hpp SootMomMethod.hpp
    PUBLIC_DEPS SimAll::Solver)
