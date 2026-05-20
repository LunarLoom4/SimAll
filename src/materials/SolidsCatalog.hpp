// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/SolidsCatalog.hpp
// Phase  : 15.5 — Engineering solids catalog (structural, conductive).
//
// Curated catalog of common engineering solids with reference-quality
// thermophysical properties at 300 K (representative ranges per
// CRC Handbook / Incropera & DeWitt / Touloukian).  Used by
// MaterialAssignment.setSolid() to populate FieldRegistry::scalar("k") etc.
// =============================================================================
#pragma once

#include <string>
#include <vector>

namespace simall::materials
{

struct SolidProperties
{
    std::string name;
    double density;          // kg/m^3
    double specificHeat;     // J/(kg·K)
    double conductivity;     // W/(m·K)  (isotropic; anisotropic via AnisotropicConductivity)
    double youngsModulus;    // Pa
    double poissonsRatio;    // —
    double thermalExpansion; // 1/K
    double emissivity;       // —      (polished surface assumption)
};

/// Built-in solids catalog. All values @ 300 K, representative.
class SolidsCatalog
{
public:
    SolidsCatalog();

    const SolidProperties* find(const std::string& name) const noexcept;
    const std::vector<SolidProperties>& all() const noexcept { return entries_; }

private:
    std::vector<SolidProperties> entries_;
};

} // namespace simall::materials
