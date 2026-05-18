// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/Models.cpp
// Phase  : 8.1–8.5
//
// Concrete model skeletons (one TU per family). Each model declares its
// transport equations + closures; the actual matrix assembly is performed in
// dedicated implementation files added during Phase 8 hand-off.
// =============================================================================
#include "turbulence/ITurbulenceModel.hpp"
#include "turbulence/KOmegaSST.hpp"
#include "turbulence/KEpsilonStandard.hpp"
#include "turbulence/SpalartAllmaras.hpp"
#include "turbulence/Smagorinsky.hpp"
#include "turbulence/Wale.hpp"
#include "turbulence/KEpsilonRng.hpp"
#include "turbulence/KEpsilonRealizable.hpp"
#include "turbulence/KOmegaStandard.hpp"
#include "turbulence/DynamicSmagorinsky.hpp"
#include "turbulence/IDDES.hpp"
#include "turbulence/SasSst.hpp"
#include "turbulence/LRR_Reynolds_Stress.hpp"
#include "turbulence/KKLOmegaTransition.hpp"
#include "core/Logger.hpp"

namespace simall::turbulence {

namespace {

struct Registrar {
    Registrar() {
        auto& r = TurbulenceRegistry::instance();
        r.register_model("kEpsilonStandard", []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<KEpsilonStandard_Full>(); });
        r.register_model("kOmegaSST",        []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<KOmegaSST_Full>(); });
        r.register_model("SpalartAllmaras",  []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<SpalartAllmaras_Full>(); });
        r.register_model("smagorinsky",      []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<Smagorinsky_LES>(); });
        r.register_model("WALE",             []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<WALE_LES>(); });
        r.register_model("kEpsilonRNG",        []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<KEpsilonRng_Full>(); });
        r.register_model("kEpsilonRealizable", []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<KEpsilonRealizable_Full>(); });
        r.register_model("kOmegaStandard",     []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<KOmegaStandard_Full>(); });
        r.register_model("dynamicSmagorinsky", []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<DynamicSmagorinsky_LES>(); });
        r.register_model("SA-IDDES",           []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<IDDES_Full>(); });
        r.register_model("SAS-SST",            []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<SasSst_Full>(); });
        r.register_model("LRR-RSM",            []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<LRR_Reynolds_Stress_Full>(); });
        r.register_model("kKLOmegaTransition", []() -> std::unique_ptr<ITurbulenceModel> {
            return std::make_unique<KKLOmegaTransition_Full>(); });
        // γ-Reθ_t (Langtry-Menter 2009) is a *coupled* model that augments
        // k-ω SST via the gamma_eff scalar field — see
        // GammaReThetaTransition.hpp. It is not a stand-alone closure and is
        // therefore activated by selecting "kOmegaSST+gammaReTheta" once the
        // host enables transition coupling, rather than via factory creation.
        SIMALL_LOG_INFO("Turbulence",
            "Registered models: k-ε (Std/RNG/Realizable), k-ω (Std/SST), SA, SA-DDES, SA-IDDES, "
            "Smagorinsky, dynSmag, WALE, SAS-SST, LRR-RSM, k-kL-ω, γ-Reθ (coupled to SST)");
    }
};
static Registrar s_registrar;

}  // namespace

}  // namespace simall::turbulence
