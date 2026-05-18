// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/FluidsCatalog.hpp
// Phase  : 15.6 — Reference fluid library (Newtonian, single-species).
//
// Reference-quality properties @ 300 K, 1 atm.  Where Sutherland viscosity
// data are available the Sutherland parameters are included so the catalog
// entry can directly populate a Material with a PropertyFunction.
//
// Pulled from CRC Handbook, NIST WebBook, White (Viscous Fluid Flow).
// =============================================================================
#pragma once

#include "materials/TransportProperties.hpp"

#include <string>
#include <vector>

namespace simall::materials {

struct FluidProperties {
    std::string name;
    double      density;              // kg/m^3
    double      viscosity;            // Pa·s
    double      conductivity;         // W/(m·K)
    double      specificHeat;         // J/(kg·K)
    double      molecularWeight;      // kg/mol
    double      gamma;                // c_p / c_v   (gases only; -1 for liquids)
    double      surfaceTension;       // N/m  (liquid–air; -1 for gases)
    double      vaporPressure300K;    // Pa   (liquids; 0 for gases)
    bool        isGas;
    SutherlandViscosity    sutherlandVisc;     // valid if isGas
    SutherlandConductivity sutherlandCond;     // valid if isGas
};

class FluidsCatalog {
public:
    FluidsCatalog();

    const FluidProperties* find(const std::string& name) const noexcept;
    const std::vector<FluidProperties>& all() const noexcept { return entries_; }

private:
    std::vector<FluidProperties> entries_;
};

}  // namespace simall::materials
