// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/UnitSystem.hpp
// Phase  : 1.3 (APPLICATION CORE → Unit System)
//
// Tagged-quantity unit conversions. All solver and CAD math runs in SI
// internally; this layer converts between SI, USCS (foot-pound-second),
// and CGS for I/O and the GUI. Per-quantity registries map a string unit
// label ("mm", "psi", "BTU/(h·ft·°F)", …) to a linear `value * scale +
// offset` transform from SI.
//
// Usage:
//   auto v_si = units::convert(Quantity::Pressure, 14.7, "psi", UnitSystem::SI);
//   auto v_us = units::convert(Quantity::Length,  1.0,  UnitSystem::SI, "in");
//
// Temperature uses an affine transform (offset != 0 for Celsius/Fahrenheit).
// =============================================================================
#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace simall::core::units {

enum class UnitSystem : int { SI, USCS, CGS };

enum class Quantity : int {
    Length,
    Mass,
    Time,
    Temperature,
    Velocity,
    Pressure,
    Force,
    Energy,
    Power,
    Density,
    DynamicViscosity,
    KinematicViscosity,
    Conductivity,
    SpecificHeat,
    HeatFlux,
    MassFlow,
    VolumeFlow,
    _Count
};

/// Pretty name (e.g. "Pressure").
std::string_view quantityName(Quantity q);

/// SI unit symbol for a quantity (e.g. "Pa" for Pressure).
std::string_view siSymbol(Quantity q);

/// USCS / CGS conventional symbol for a quantity.
std::string_view conventionalSymbol(Quantity q, UnitSystem sys);

/// Convert `value` of `q` expressed in unit `fromLabel` into SI base.
/// Returns nullopt if the label is not registered for that quantity.
std::optional<double> toSI(Quantity q, double value, std::string_view fromLabel);

/// Convert SI `value` of `q` into the unit `toLabel`.
std::optional<double> fromSI(Quantity q, double siValue, std::string_view toLabel);

/// Convert between any two labels: a → SI → b.
std::optional<double> convert(Quantity q, double value,
                              std::string_view fromLabel,
                              std::string_view toLabel);

/// Convert SI value into the system-default unit (e.g. m → ft for USCS).
double fromSI(Quantity q, double siValue, UnitSystem toSys);
double toSI(Quantity q, double value, UnitSystem fromSys);

/// Register a custom unit label for a quantity: si = value * scale + offset.
void registerUnit(Quantity q, std::string label, double scale, double offset = 0.0);

/// Enumerate every label currently registered for a quantity.
std::vector<std::string> labelsFor(Quantity q);

}  // namespace simall::core::units
