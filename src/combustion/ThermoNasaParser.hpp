// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/ThermoNasaParser.hpp
// Phase  : 11.11 — NASA polynomial thermodynamics parser.
//
// Supports two record formats:
//
//  - NASA-7 (Gordon-McBride 1971, two temperature ranges, 15 coefficients):
//      cp/R   = a1 + a2 T + a3 T² + a4 T³ + a5 T⁴
//      h/(RT) = a1 + a2 T/2 + a3 T²/3 + a4 T³/4 + a5 T⁴/5 + a6/T
//      s/R    = a1 ln T + a2 T + a3 T²/2 + a4 T³/3 + a5 T⁴/4 + a7
//
//  - NASA-9 (Burcat-Ruscic 2005, two or more ranges, 7 polynomial + 2 integ.):
//      cp/R   = a1/T² + a2/T + a3 + a4 T + a5 T² + a6 T³ + a7 T⁴
//      h/(RT) = -a1/T² + a2 ln(T)/T + a3 + a4 T/2 + a5 T²/3 + a6 T³/4
//               + a7 T⁴/5 + b1/T
//      s/R    = -a1/(2 T²) - a2/T + a3 ln(T) + a4 T + a5 T²/2
//               + a6 T³/3 + a7 T⁴/4 + b2
//
// Both formats follow strict 80-column FORTRAN layout; reader is whitespace-
// tolerant within columns.  Temperature ranges are inclusive of T_low/T_high
// per range; selection picks the smallest range containing T.
//
// All cp / h / s returned in SI (J/(kg·K), J/kg, J/(kg·K)).
// =============================================================================
#pragma once

#include <map>
#include <string>
#include <vector>

namespace simall::combustion
{

struct NasaTempRange
{
    double Tlow = 0.0;
    double Thigh = 0.0;
    double a[9]{}; // NASA-7 uses [0..6], NASA-9 uses [0..6]+b1=a[7]+b2=a[8]
};

struct NasaSpeciesThermo
{
    std::string name;
    std::map<std::string, int> composition; // element symbol → atom count
    double molarMass = 0.0;                 // g/mol (computed)
    bool isNasa9 = false;
    std::vector<NasaTempRange> ranges; // 2 for NASA-7, ≥2 for NASA-9
};

class ThermoNasaParser
{
public:
    /// Parse a Chemkin-style THERMO block or a NASA-9 burcat-style file.
    /// Auto-detects format by record-line structure.
    std::vector<NasaSpeciesThermo> parse_file(const std::string& path) const;
    std::vector<NasaSpeciesThermo> parse_string(const std::string& text) const;

    /// Specific heat at constant pressure (J/(kg·K)) using the NASA polynomial.
    static double cp(const NasaSpeciesThermo& sp, double T);
    /// Enthalpy h (J/kg) (includes formation enthalpy reference T = 298.15 K).
    static double h(const NasaSpeciesThermo& sp, double T);
    /// Entropy s (J/(kg·K)) at T (and reference pressure 1 atm).
    static double s(const NasaSpeciesThermo& sp, double T);

    /// Universal gas constant (SI).
    static constexpr double Ru = 8.314462618; // J/(mol·K)
};

} // namespace simall::combustion
