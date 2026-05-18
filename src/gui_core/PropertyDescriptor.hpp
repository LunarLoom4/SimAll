// =============================================================================
// SimAll Beta — GUI Core (Qt-free)
// File   : src/gui_core/PropertyDescriptor.hpp
//
// Metadata-only descriptor for a single editable property.  Pure C++, no Qt,
// no resource dependencies — it is consumed both by the desktop widget builder
// (`gui::PropertyEditorV2`) and by the headless scripting / CLI batch layer
// (W16 PyBindings + W20 plugins).
//
// Design intent:
//   * Strongly-typed via std::variant.  No string-based property protocol.
//   * Validation is descriptor-driven (`validate()`); the widget never needs
//     to know solver/units semantics.
//   * Visibility and enablement are *functions of the current value set*,
//     not of widget state — so headless validators behave identically.
//   * Edit events fan out through `onChanged` + a typed `Bus` so the
//     workflow tree can refresh derived state.
// =============================================================================
#pragma once

#include <array>
#include <functional>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace simall::gui_core {

// -- type tag -----------------------------------------------------------------
enum class PropertyType : uint8_t {
    Double = 0,
    Int,
    Bool,
    Enum,
    String,
    Vec3,
    FilePath,   // requires `fileFilter` and `pathMustExist`
    Color,      // stored as Vec3 (RGB, linear, 0..1)
    Range       // pair-of-doubles; stored as Vec3 {min, max, 0}
};

// -- value -------------------------------------------------------------------
using Variant = std::variant<double, int, bool, std::string,
                             std::array<double, 3>>;

inline bool variant_holds_numeric(const Variant& v) {
    return std::holds_alternative<double>(v) || std::holds_alternative<int>(v);
}

inline double variant_as_double(const Variant& v) {
    if (auto* d = std::get_if<double>(&v)) return *d;
    if (auto* i = std::get_if<int>(&v))    return static_cast<double>(*i);
    return 0.0;
}

// -- validation result -------------------------------------------------------
struct ValidationResult {
    bool        ok = true;
    std::string message;   // empty when ok
    static ValidationResult success()                       { return {true, {}}; }
    static ValidationResult failure(std::string msg)        { return {false, std::move(msg)}; }
};

// -- descriptor --------------------------------------------------------------
//
// All fields are public-by-design: GUI code reads them straight; scripting
// code constructs them by aggregate init.
//
struct PropertyDescriptor {
    // identity
    std::string  propertyName;
    std::string  displayName;       // optional; falls back to propertyName
    std::string  category;          // optional; controls grouping in v2 editor

    // type + default
    PropertyType type             = PropertyType::Double;
    Variant      defaultValue     = double{0};
    Variant      currentValue     = double{0};

    // numeric bounds (apply to Double / Int / Range)
    double       minimum          = -1e300;
    double       maximum          =  1e300;
    double       step             = 0.0;     // 0 → editor default

    // enum
    std::vector<std::string> enumOptions;

    // file-path
    std::string  fileFilter;        // "STEP (*.step *.stp);;All (*)"
    bool         pathMustExist     = false;

    // metadata
    std::string  tooltip;
    std::string  units;
    std::string  helpUrl;
    bool         requiresRestart   = false;
    bool         advanced          = false;  // hidden unless Show Advanced

    // dynamic predicates — pure functions of the *property bag* they live in.
    // The editor passes `this` as the lookup context; descriptors closure
    // over whichever PropertyDescriptor* references they need.
    std::function<bool()>                 visibilityCondition;
    std::function<bool()>                 enabledCondition;

    // custom validator runs after built-in min/max/enum checks.
    std::function<ValidationResult(const Variant&)> customValidator;

    // notification — invoked AFTER successful commit.
    std::function<void(const Variant&)>   onChanged;

    // -- API ---------------------------------------------------------------
    [[nodiscard]] const std::string& display() const noexcept {
        return displayName.empty() ? propertyName : displayName;
    }
    [[nodiscard]] bool is_visible() const {
        return !visibilityCondition || visibilityCondition();
    }
    [[nodiscard]] bool is_enabled() const {
        return !enabledCondition || enabledCondition();
    }

    // validate a candidate value against the descriptor's constraints.
    // Used by the GUI BEFORE commit and by scripting code BEFORE assignment.
    [[nodiscard]] ValidationResult validate(const Variant& candidate) const;

    // commit `candidate` if valid.  Returns the same ValidationResult.
    ValidationResult commit(Variant candidate);
};

// -- bag helpers --------------------------------------------------------------
//
// Find a descriptor by `propertyName` in a vector.  Returns nullptr if missing.
// Used by visibilityCondition closures, e.g.:
//     d.visibilityCondition = [&bag]{
//         auto* p = simall::gui_core::find(bag, "turbulence.model");
//         return p && std::get<int>(p->currentValue) != 0;
//     };
//
[[nodiscard]] PropertyDescriptor*       find(std::vector<PropertyDescriptor>& bag,
                                              std::string_view name) noexcept;
[[nodiscard]] const PropertyDescriptor* find(const std::vector<PropertyDescriptor>& bag,
                                              std::string_view name) noexcept;

// Snapshot the bag's name → currentValue pairs in textual form (for project
// serialization, undo history, and diff display).  Numeric values are
// printed with std::to_string; vectors are "x,y,z"; booleans are "true"/"false".
[[nodiscard]] std::string to_textual_snapshot(const std::vector<PropertyDescriptor>& bag);

}  // namespace simall::gui_core
