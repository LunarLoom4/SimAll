// =============================================================================
// SimAll Beta — GUI Core (Qt-free)
// File   : src/gui_core/PropertyDescriptor.cpp
// =============================================================================
#include "gui_core/PropertyDescriptor.hpp"

#include <cmath>
#include <filesystem>
#include <sstream>

namespace simall::gui_core {

namespace {

ValidationResult validate_numeric(const PropertyDescriptor& d,
                                  double v) {
    if (!std::isfinite(v))
        return ValidationResult::failure(d.display() + " must be finite");
    if (v < d.minimum)
        return ValidationResult::failure(d.display() + " below minimum " +
                                         std::to_string(d.minimum));
    if (v > d.maximum)
        return ValidationResult::failure(d.display() + " above maximum " +
                                         std::to_string(d.maximum));
    return ValidationResult::success();
}

ValidationResult validate_enum(const PropertyDescriptor& d, int idx) {
    if (idx < 0 || static_cast<size_t>(idx) >= d.enumOptions.size())
        return ValidationResult::failure(d.display() +
                                         " enum index out of range");
    return ValidationResult::success();
}

ValidationResult validate_vec3(const PropertyDescriptor& d,
                               const std::array<double,3>& v) {
    for (int i = 0; i < 3; ++i) {
        auto r = validate_numeric(d, v[i]);
        if (!r.ok) return r;
    }
    if (d.type == PropertyType::Range) {
        if (v[0] > v[1])
            return ValidationResult::failure(d.display() +
                                             " range min > max");
    }
    if (d.type == PropertyType::Color) {
        for (int i = 0; i < 3; ++i)
            if (v[i] < 0.0 || v[i] > 1.0)
                return ValidationResult::failure(d.display() +
                                                  " color component outside [0,1]");
    }
    return ValidationResult::success();
}

ValidationResult validate_path(const PropertyDescriptor& d,
                               const std::string& p) {
    if (p.empty())
        return ValidationResult::failure(d.display() + " path is empty");
    if (d.pathMustExist && !std::filesystem::exists(p))
        return ValidationResult::failure(d.display() +
                                         " path does not exist: " + p);
    return ValidationResult::success();
}

}  // namespace

ValidationResult PropertyDescriptor::validate(const Variant& c) const {
    // built-in checks first
    switch (type) {
        case PropertyType::Double:
        case PropertyType::Int: {
            if (!variant_holds_numeric(c))
                return ValidationResult::failure(display() + " expected numeric");
            const double v = variant_as_double(c);
            auto r = validate_numeric(*this, v);
            if (!r.ok) return r;
            break;
        }
        case PropertyType::Bool: {
            if (!std::holds_alternative<bool>(c))
                return ValidationResult::failure(display() + " expected bool");
            break;
        }
        case PropertyType::Enum: {
            if (!std::holds_alternative<int>(c))
                return ValidationResult::failure(display() + " expected enum index (int)");
            auto r = validate_enum(*this, std::get<int>(c));
            if (!r.ok) return r;
            break;
        }
        case PropertyType::String: {
            if (!std::holds_alternative<std::string>(c))
                return ValidationResult::failure(display() + " expected string");
            break;
        }
        case PropertyType::Vec3:
        case PropertyType::Color:
        case PropertyType::Range: {
            if (!std::holds_alternative<std::array<double,3>>(c))
                return ValidationResult::failure(display() + " expected vec3-like");
            auto r = validate_vec3(*this, std::get<std::array<double,3>>(c));
            if (!r.ok) return r;
            break;
        }
        case PropertyType::FilePath: {
            if (!std::holds_alternative<std::string>(c))
                return ValidationResult::failure(display() + " expected path string");
            auto r = validate_path(*this, std::get<std::string>(c));
            if (!r.ok) return r;
            break;
        }
    }
    if (customValidator) {
        auto r = customValidator(c);
        if (!r.ok) return r;
    }
    return ValidationResult::success();
}

ValidationResult PropertyDescriptor::commit(Variant candidate) {
    auto r = validate(candidate);
    if (!r.ok) return r;
    currentValue = std::move(candidate);
    if (onChanged) onChanged(currentValue);
    return ValidationResult::success();
}

PropertyDescriptor* find(std::vector<PropertyDescriptor>& bag,
                          std::string_view name) noexcept {
    for (auto& d : bag) if (d.propertyName == name) return &d;
    return nullptr;
}

const PropertyDescriptor* find(const std::vector<PropertyDescriptor>& bag,
                                std::string_view name) noexcept {
    for (auto& d : bag) if (d.propertyName == name) return &d;
    return nullptr;
}

std::string to_textual_snapshot(const std::vector<PropertyDescriptor>& bag) {
    std::ostringstream ss;
    for (const auto& d : bag) {
        ss << d.propertyName << " = ";
        std::visit([&](const auto& v) {
            using T = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<T, double>)
                ss << v;
            else if constexpr (std::is_same_v<T, int>)
                ss << v;
            else if constexpr (std::is_same_v<T, bool>)
                ss << (v ? "true" : "false");
            else if constexpr (std::is_same_v<T, std::string>)
                ss << '"' << v << '"';
            else if constexpr (std::is_same_v<T, std::array<double,3>>)
                ss << v[0] << ',' << v[1] << ',' << v[2];
        }, d.currentValue);
        ss << '\n';
    }
    return ss.str();
}

}  // namespace simall::gui_core
