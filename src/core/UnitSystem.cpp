// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/UnitSystem.cpp
// =============================================================================
#include "core/UnitSystem.hpp"

#include <array>
#include <cstddef>
#include <mutex>
#include <unordered_map>

namespace simall::core::units
{

namespace
{

struct UnitDef
{
    double scale = 1.0;
    double offset = 0.0;
};
using LabelMap = std::unordered_map<std::string, UnitDef>;

constexpr std::size_t kQ = static_cast<std::size_t>(Quantity::_Count);

struct Registry
{
    std::mutex mtx;
    std::array<LabelMap, kQ> maps;
};

Registry& reg()
{
    static Registry r;
    static std::once_flag once;
    std::call_once(once, [] {
        auto add = [&](Quantity q, const char* lbl, double s, double o = 0.0) {
            r.maps[static_cast<std::size_t>(q)][lbl] = {s, o};
        };
        // Identity SI entries.
        add(Quantity::Length, "m", 1.0);
        add(Quantity::Length, "cm", 0.01);
        add(Quantity::Length, "mm", 1e-3);
        add(Quantity::Length, "um", 1e-6);
        add(Quantity::Length, "km", 1000.0);
        add(Quantity::Length, "in", 0.0254);
        add(Quantity::Length, "ft", 0.3048);
        add(Quantity::Length, "yd", 0.9144);
        add(Quantity::Length, "mi", 1609.344);

        add(Quantity::Mass, "kg", 1.0);
        add(Quantity::Mass, "g", 1e-3);
        add(Quantity::Mass, "mg", 1e-6);
        add(Quantity::Mass, "lb", 0.45359237);
        add(Quantity::Mass, "slug", 14.5939029);

        add(Quantity::Time, "s", 1.0);
        add(Quantity::Time, "ms", 1e-3);
        add(Quantity::Time, "us", 1e-6);
        add(Quantity::Time, "min", 60.0);
        add(Quantity::Time, "h", 3600.0);

        add(Quantity::Temperature, "K", 1.0, 0.0);
        add(Quantity::Temperature, "C", 1.0, 273.15);
        add(Quantity::Temperature, "F", 5.0 / 9.0, 459.67 * 5.0 / 9.0);
        add(Quantity::Temperature, "R", 5.0 / 9.0, 0.0);

        add(Quantity::Velocity, "m/s", 1.0);
        add(Quantity::Velocity, "km/h", 1000.0 / 3600.0);
        add(Quantity::Velocity, "mph", 0.44704);
        add(Quantity::Velocity, "ft/s", 0.3048);
        add(Quantity::Velocity, "knot", 0.514444);

        add(Quantity::Pressure, "Pa", 1.0);
        add(Quantity::Pressure, "kPa", 1000.0);
        add(Quantity::Pressure, "MPa", 1e6);
        add(Quantity::Pressure, "bar", 1e5);
        add(Quantity::Pressure, "atm", 101325.0);
        add(Quantity::Pressure, "psi", 6894.757293168);
        add(Quantity::Pressure, "torr", 133.322387415);

        add(Quantity::Force, "N", 1.0);
        add(Quantity::Force, "kN", 1000.0);
        add(Quantity::Force, "lbf", 4.4482216152605);
        add(Quantity::Force, "dyn", 1e-5);

        add(Quantity::Energy, "J", 1.0);
        add(Quantity::Energy, "kJ", 1000.0);
        add(Quantity::Energy, "MJ", 1e6);
        add(Quantity::Energy, "cal", 4.184);
        add(Quantity::Energy, "kcal", 4184.0);
        add(Quantity::Energy, "BTU", 1055.05585262);
        add(Quantity::Energy, "erg", 1e-7);

        add(Quantity::Power, "W", 1.0);
        add(Quantity::Power, "kW", 1000.0);
        add(Quantity::Power, "MW", 1e6);
        add(Quantity::Power, "hp", 745.6998715822702);
        add(Quantity::Power, "BTU/h", 0.29307107017222);

        add(Quantity::Density, "kg/m^3", 1.0);
        add(Quantity::Density, "g/cm^3", 1000.0);
        add(Quantity::Density, "lb/ft^3", 16.018463373960140);

        add(Quantity::DynamicViscosity, "Pa.s", 1.0);
        add(Quantity::DynamicViscosity, "cP", 1e-3);
        add(Quantity::DynamicViscosity, "P", 0.1);
        add(Quantity::DynamicViscosity, "lb/(ft.s)", 1.488163943569553);

        add(Quantity::KinematicViscosity, "m^2/s", 1.0);
        add(Quantity::KinematicViscosity, "cSt", 1e-6);
        add(Quantity::KinematicViscosity, "St", 1e-4);
        add(Quantity::KinematicViscosity, "ft^2/s", 0.09290304);

        add(Quantity::Conductivity, "W/(m.K)", 1.0);
        add(Quantity::Conductivity, "BTU/(h.ft.F)", 1.730734666372968);
        add(Quantity::Conductivity, "cal/(s.cm.K)", 418.4);

        add(Quantity::SpecificHeat, "J/(kg.K)", 1.0);
        add(Quantity::SpecificHeat, "kJ/(kg.K)", 1000.0);
        add(Quantity::SpecificHeat, "BTU/(lb.F)", 4186.8);

        add(Quantity::HeatFlux, "W/m^2", 1.0);
        add(Quantity::HeatFlux, "BTU/(h.ft^2)", 3.15459074506);

        add(Quantity::MassFlow, "kg/s", 1.0);
        add(Quantity::MassFlow, "g/s", 1e-3);
        add(Quantity::MassFlow, "lb/s", 0.45359237);
        add(Quantity::MassFlow, "lb/h", 0.45359237 / 3600.0);

        add(Quantity::VolumeFlow, "m^3/s", 1.0);
        add(Quantity::VolumeFlow, "L/s", 1e-3);
        add(Quantity::VolumeFlow, "L/min", 1e-3 / 60.0);
        add(Quantity::VolumeFlow, "gpm", 6.30901964e-5); // US gallon per minute
        add(Quantity::VolumeFlow, "cfm", 4.71947443e-4); // cubic foot per minute
    });
    return r;
}

const char* kNames[kQ] = {"Length",
                          "Mass",
                          "Time",
                          "Temperature",
                          "Velocity",
                          "Pressure",
                          "Force",
                          "Energy",
                          "Power",
                          "Density",
                          "DynamicViscosity",
                          "KinematicViscosity",
                          "Conductivity",
                          "SpecificHeat",
                          "HeatFlux",
                          "MassFlow",
                          "VolumeFlow"};
const char* kSi[kQ] = {"m",
                       "kg",
                       "s",
                       "K",
                       "m/s",
                       "Pa",
                       "N",
                       "J",
                       "W",
                       "kg/m^3",
                       "Pa.s",
                       "m^2/s",
                       "W/(m.K)",
                       "J/(kg.K)",
                       "W/m^2",
                       "kg/s",
                       "m^3/s"};
const char* kUscs[kQ] = {"ft",
                         "lb",
                         "s",
                         "F",
                         "ft/s",
                         "psi",
                         "lbf",
                         "BTU",
                         "hp",
                         "lb/ft^3",
                         "lb/(ft.s)",
                         "ft^2/s",
                         "BTU/(h.ft.F)",
                         "BTU/(lb.F)",
                         "BTU/(h.ft^2)",
                         "lb/s",
                         "cfm"};
const char* kCgs[kQ] = {"cm",
                        "g",
                        "s",
                        "K",
                        "cm/s",
                        "Pa",
                        "dyn",
                        "erg",
                        "W",
                        "g/cm^3",
                        "P",
                        "St",
                        "cal/(s.cm.K)",
                        "J/(kg.K)",
                        "W/m^2",
                        "g/s",
                        "m^3/s"};

const char* defaultLabel(Quantity q, UnitSystem sys)
{
    auto i = static_cast<std::size_t>(q);
    switch (sys) {
    case UnitSystem::SI:
        return kSi[i];
    case UnitSystem::USCS:
        return kUscs[i];
    case UnitSystem::CGS:
        return kCgs[i];
    }
    return kSi[i];
}

} // namespace

std::string_view quantityName(Quantity q)
{
    return kNames[static_cast<std::size_t>(q)];
}
std::string_view siSymbol(Quantity q)
{
    return kSi[static_cast<std::size_t>(q)];
}
std::string_view conventionalSymbol(Quantity q, UnitSystem sys)
{
    return defaultLabel(q, sys);
}

std::optional<double> toSI(Quantity q, double value, std::string_view fromLabel)
{
    auto& r = reg();
    std::lock_guard lk(r.mtx);
    auto& m = r.maps[static_cast<std::size_t>(q)];
    auto it = m.find(std::string(fromLabel));
    if (it == m.end())
        return std::nullopt;
    return value * it->second.scale + it->second.offset;
}

std::optional<double> fromSI(Quantity q, double siValue, std::string_view toLabel)
{
    auto& r = reg();
    std::lock_guard lk(r.mtx);
    auto& m = r.maps[static_cast<std::size_t>(q)];
    auto it = m.find(std::string(toLabel));
    if (it == m.end())
        return std::nullopt;
    return (siValue - it->second.offset) / it->second.scale;
}

std::optional<double> convert(Quantity q,
                              double value,
                              std::string_view fromLabel,
                              std::string_view toLabel)
{
    auto si = toSI(q, value, fromLabel);
    if (!si)
        return std::nullopt;
    return fromSI(q, *si, toLabel);
}

double fromSI(Quantity q, double siValue, UnitSystem toSys)
{
    auto v = fromSI(q, siValue, defaultLabel(q, toSys));
    return v ? *v : siValue;
}
double toSI(Quantity q, double value, UnitSystem fromSys)
{
    auto v = toSI(q, value, defaultLabel(q, fromSys));
    return v ? *v : value;
}

void registerUnit(Quantity q, std::string label, double scale, double offset)
{
    auto& r = reg();
    std::lock_guard lk(r.mtx);
    r.maps[static_cast<std::size_t>(q)][std::move(label)] = {scale, offset};
}

std::vector<std::string> labelsFor(Quantity q)
{
    auto& r = reg();
    std::lock_guard lk(r.mtx);
    auto const& m = r.maps[static_cast<std::size_t>(q)];
    std::vector<std::string> out;
    out.reserve(m.size());
    for (auto const& [k, _] : m)
        out.push_back(k);
    return out;
}

} // namespace simall::core::units
