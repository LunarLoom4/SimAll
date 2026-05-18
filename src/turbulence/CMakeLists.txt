simall_add_module(turbulence
    SOURCES Models.cpp KOmegaSST.cpp KEpsilonStandard.cpp SpalartAllmaras.cpp Smagorinsky.cpp Wale.cpp GammaReThetaTransition.cpp SaDdes.cpp WallFunctions.cpp KEpsilonRng.cpp KEpsilonRealizable.cpp KOmegaStandard.cpp DynamicSmagorinsky.cpp IDDES.cpp SasSst.cpp LRR_Reynolds_Stress.cpp KKLOmegaTransition.cpp
    HEADERS ITurbulenceModel.hpp KOmegaSST.hpp KEpsilonStandard.hpp SpalartAllmaras.hpp Smagorinsky.hpp Wale.hpp GammaReThetaTransition.hpp SaDdes.hpp WallFunctions.hpp KEpsilonRng.hpp KEpsilonRealizable.hpp KOmegaStandard.hpp DynamicSmagorinsky.hpp IDDES.hpp SasSst.hpp LRR_Reynolds_Stress.hpp KKLOmegaTransition.hpp
    PUBLIC_DEPS SimAll::Solver SimAll::Meshing)
