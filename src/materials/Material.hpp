// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/Material.hpp
// Phase  : 15 (MATERIAL DATABASE)
//
// Material model with temperature-dependent properties (polynomial fits,
// Sutherland law, tabulated interpolation per spec section 15).
// =============================================================================
#pragma once

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::materials
{

enum class PropertyKind
{
    Density,
    Viscosity,
    Conductivity,
    SpecificHeat,
    MolecularWeight,
    SurfaceTension,
    ElectricalConductivity,
    Emissivity
};

class PropertyFunction
{
public:
    using Fn = std::function<double(double T /*K*/, double p /*Pa*/)>;
    PropertyFunction() = default;
    PropertyFunction(double constant) : fn_([constant](double, double) { return constant; }) {}
    PropertyFunction(Fn fn) : fn_(std::move(fn)) {}
    double operator()(double T, double p) const { return fn_ ? fn_(T, p) : 0.0; }

private:
    Fn fn_;
};

struct SpeciesEntry
{
    std::string name;
    double molecularWeight;
};

class Material
{
public:
    explicit Material(std::string name) : name_(std::move(name)) {}

    const std::string& name() const noexcept { return name_; }
    void set(PropertyKind k, PropertyFunction f) { props_[k] = std::move(f); }
    double evaluate(PropertyKind k, double T, double p) const
    {
        auto it = props_.find(k);
        return it == props_.end() ? 0.0 : it->second(T, p);
    }

    std::vector<SpeciesEntry>& species() noexcept { return species_; }

private:
    std::string name_;
    std::unordered_map<PropertyKind, PropertyFunction> props_;
    std::vector<SpeciesEntry> species_;
};

class MaterialDatabase
{
public:
    Material& add(std::string name);
    Material* find(const std::string& name);
    const auto& all() const noexcept { return store_; }

private:
    std::unordered_map<std::string, Material> store_;
};

} // namespace simall::materials
