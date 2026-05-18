// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/SolidsCatalog.cpp
// =============================================================================
#include "materials/SolidsCatalog.hpp"

namespace simall::materials {

SolidsCatalog::SolidsCatalog() {
    entries_ = {
      //  name                          rho       cp     k      E         nu     alpha       eps
      { "Aluminum-6061",                2700.0,  896.0, 167.0, 6.89e10,  0.33,  23.6e-6,    0.09 },
      { "Aluminum-pure",                2702.0,  903.0, 237.0, 7.0e10,   0.33,  23.1e-6,    0.04 },
      { "Steel-AISI1010",               7832.0,  434.0,  63.9, 2.0e11,   0.29,  12.0e-6,    0.20 },
      { "Steel-Stainless-304",          7900.0,  477.0,  14.9, 1.93e11,  0.27,  17.3e-6,    0.30 },
      { "Steel-Stainless-316L",         8000.0,  500.0,  16.3, 1.93e11,  0.30,  16.0e-6,    0.28 },
      { "Copper",                       8933.0,  385.0, 401.0, 1.10e11,  0.34,  16.5e-6,    0.04 },
      { "Brass-CartridgeBrass",         8530.0,  380.0, 110.0, 1.10e11,  0.34,  19.9e-6,    0.04 },
      { "Titanium-Ti6Al4V",             4430.0,  526.0,   6.7, 1.14e11,  0.34,   8.6e-6,    0.30 },
      { "Inconel-718",                  8190.0,  435.0,  11.4, 2.0e11,   0.29,  13.0e-6,    0.28 },
      { "Nickel",                       8900.0,  444.0,  90.7, 2.0e11,   0.31,  13.4e-6,    0.06 },
      { "Silicon",                      2329.0,  712.0, 148.0, 1.30e11,  0.28,   2.6e-6,    0.30 },
      { "Silicon-Carbide",              3210.0,  675.0,  120.0,4.1e11,   0.14,   4.0e-6,    0.85 },
      { "Alumina-99",                   3960.0,  880.0,  35.0, 3.7e11,   0.22,   8.1e-6,    0.85 },
      { "Glass-BorosilicatePyrex",      2225.0,  835.0,   1.4, 6.4e10,   0.20,   3.3e-6,    0.90 },
      { "Concrete-Structural",          2300.0,  880.0,   1.4, 3.0e10,   0.20,  10.0e-6,    0.94 },
      { "BrickRed",                     1920.0,  835.0,   0.72,2.0e10,   0.20,   6.0e-6,    0.93 },
      { "Wood-Oak",                      545.0, 2385.0,   0.17,1.1e10,   0.37,   5.0e-6,    0.90 },
      { "Rubber-Silicone",              1100.0, 1500.0,   0.20,5.0e6,    0.49, 270.0e-6,    0.94 },
      { "Polymer-PTFE",                 2200.0, 1050.0,   0.25,5.0e8,    0.46, 135.0e-6,    0.92 },
      { "Polymer-Epoxy",                1180.0, 1700.0,   0.20,3.5e9,    0.35,  60.0e-6,    0.85 },
      { "Ice-273K",                      920.0, 2040.0,   2.20,9.3e9,    0.33,  51.0e-6,    0.97 },
      { "Graphite-PolyXfine",           1950.0,  710.0, 120.0, 1.0e10,   0.20,   7.9e-6,    0.85 },
      { "Foam-PUExpanded",                40.0, 1300.0,   0.024,5.0e6,   0.30,  60.0e-6,    0.90 },
    };
}

const SolidProperties* SolidsCatalog::find(const std::string& name) const noexcept {
    for (auto const& e : entries_)
        if (e.name == name) return &e;
    return nullptr;
}

}  // namespace simall::materials
