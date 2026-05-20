// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/Janaf.hpp
// Phase  : 15.3 — JANAF / NASA-7 polynomial thermodynamics.
//
// Stores the standard CHEMKIN/NASA-7 two-range fit:
//
//   c_p^0(T)/R     = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
//   h^0(T)/(R T)   = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T
//   s^0(T)/R       = a1 ln T + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4 + a7
//
// Two coefficient sets:  low range  [Tmin, Tmid],  high range  [Tmid, Tmax].
// =============================================================================
#pragma once

#include <array>
#include <cmath>
#include <iosfwd>
#include <istream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::materials
{

inline constexpr double kRgasUniversalJanaf = 8.314462618; // J/(mol·K)

struct JanafSpecies
{
    std::string name;
    double molecularWeight = 0.0; // kg/mol
    double Tmin = 200.0;
    double Tmid = 1000.0;
    double Tmax = 6000.0;
    std::array<double, 7> lowCoeffs{};  // [Tmin, Tmid]
    std::array<double, 7> highCoeffs{}; // [Tmid, Tmax]

    const std::array<double, 7>& pick(double T) const noexcept
    {
        return (T < Tmid) ? lowCoeffs : highCoeffs;
    }

    /// Specific heat at constant pressure  c_p  [J/(kg·K)].
    double cp(double T) const noexcept
    {
        const auto& a = pick(T);
        const double cpoverR =
            a[0] + a[1] * T + a[2] * T * T + a[3] * T * T * T + a[4] * T * T * T * T;
        return cpoverR * kRgasUniversalJanaf / molecularWeight;
    }

    /// Total enthalpy  h  [J/kg]  (reference @ 298.15 K folded into a6).
    double enthalpy(double T) const noexcept
    {
        const auto& a = pick(T);
        const double hoverRT = a[0] + 0.5 * a[1] * T + (1.0 / 3.0) * a[2] * T * T
                               + 0.25 * a[3] * T * T * T + 0.20 * a[4] * T * T * T * T + a[5] / T;
        return hoverRT * kRgasUniversalJanaf * T / molecularWeight;
    }

    /// Standard entropy  s  [J/(kg·K)]  (no pressure correction here).
    double entropy(double T) const noexcept
    {
        const auto& a = pick(T);
        const double soverR = a[0] * std::log(T) + a[1] * T + 0.5 * a[2] * T * T
                              + (1.0 / 3.0) * a[3] * T * T * T + 0.25 * a[4] * T * T * T * T + a[6];
        return soverR * kRgasUniversalJanaf / molecularWeight;
    }
};

/// Parser for CHEMKIN-format THERMO blocks (NASA-7).  Accepts an input
/// stream positioned at the first species line (after the optional
/// "THERMO ALL" header and global T-range line).  Stops at "END".
///
/// Returns species keyed by uppercased name.
class JanafThermoParser
{
public:
    std::unordered_map<std::string, JanafSpecies> parse(std::istream& in);

    /// Same, but reads a single species record (4 fixed-format lines).
    /// Returns true if the record parsed cleanly, false on EOF/END.
    bool parseRecord(std::istream& in, JanafSpecies& out);

    const std::vector<std::string>& errors() const noexcept { return errors_; }

private:
    std::vector<std::string> errors_;
};

} // namespace simall::materials
