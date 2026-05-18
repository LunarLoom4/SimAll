// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/ChemkinParser.hpp
// Phase  : 11.10 — Chemkin-II mechanism file parser.
//
// Reads a Chemkin-II text file (.inp / .mech) with four canonical blocks:
//
//   ELEMENTS   <symbols...>   END
//   SPECIES    <names...>     END
//   THERMO     <NASA-7 records>   END     (optional, often in a sister file)
//   REACTIONS  [units]
//      <reaction equation>   A   beta   Ea
//      [DUPLICATE]
//      [REV / A_r beta_r Ea_r]
//      [LOW  / A0 beta0 Ea0]   (Lindemann/Troe fall-off)
//      [TROE / a T*** T* [T**]]
//      [THIRD-BODY enhancement +H2/2.0/+H2O/16/...]
//   END
//
// Output:
//   - Vector<ChemkinElement>  (symbol + atomic weight)
//   - Vector<ChemkinSpecies>  (name, composition map, optional thermo ref)
//   - Vector<ChemkinReaction> (equation, stoich coeffs, Arrhenius, REV/LOW/TROE)
//
// Unit handling: parses leading "REACTIONS [CAL/MOLE|KCAL/MOLE|JOULES/MOLE|
// KELVINS] [MOLES|MOLECULES]" header and converts internally to SI (J/mol, mol).
// =============================================================================
#pragma once

#include <map>
#include <string>
#include <vector>

namespace simall::combustion {

struct ChemkinElement {
    std::string symbol;
    double      atomicWeight = 0.0;     // [g/mol]
};

struct ChemkinSpecies {
    std::string                 name;
    std::map<std::string, int>  composition;   // element symbol → atom count
    bool                        hasThermo = false;
    int                         thermoIndex = -1;
};

enum class ReactionDirection { Forward, Reversible };

struct ChemkinArrhenius {
    double A    = 0.0;
    double beta = 0.0;
    double Ea   = 0.0;       // [J/mol]
};

struct TroeFallOff {
    bool   enabled = false;
    double a = 0.0, T3 = 0.0, T1 = 0.0, T2 = 0.0;   // a + (1-a)e^{-T/T***} + a e^{-T/T*} + e^{-T**/T}
};

struct ChemkinReaction {
    std::string                       equation;       // raw "A + B = C + D"
    std::map<std::string, double>     reactants;      // species → ν
    std::map<std::string, double>     products;
    ReactionDirection                 direction = ReactionDirection::Forward;
    ChemkinArrhenius                  fwd;
    ChemkinArrhenius                  rev;            // when REV/ given
    bool                              hasReverse = false;
    bool                              isLindemann = false;
    ChemkinArrhenius                  low;            // fall-off low-pressure
    TroeFallOff                       troe;
    bool                              isDuplicate = false;
    std::map<std::string, double>     thirdBodyEff;   // enhancement factors
};

enum class EnergyUnit { CalPerMole, KCalPerMole, JoulesPerMole, KJoulesPerMole, Kelvins };
enum class QuantityUnit { Moles, Molecules };

struct ChemkinMechanism {
    std::vector<ChemkinElement>  elements;
    std::vector<ChemkinSpecies>  species;
    std::vector<ChemkinReaction> reactions;
    EnergyUnit                   energyUnit   = EnergyUnit::CalPerMole;
    QuantityUnit                 quantityUnit = QuantityUnit::Moles;
};

class ChemkinParser {
public:
    /// Throws std::runtime_error on I/O failure or fundamentally malformed
    /// input.  Returns mechanism with SI-normalised Ea (J/mol) on success.
    ChemkinMechanism parse_file(const std::string& path) const;

    /// Parse a complete mechanism from an already-loaded string.
    ChemkinMechanism parse_string(const std::string& text) const;
};

}  // namespace simall::combustion
