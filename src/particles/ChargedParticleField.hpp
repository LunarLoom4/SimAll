// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ChargedParticleField.hpp
// Phase  : 13.17 — Coulomb (electrostatic) body force on Lagrangian parcels.
//
// Each parcel carries a fixed charge q_p (C).  The acceleration imparted
// by an externally-prescribed electric field E (V/m) is
//
//   a = q_p E / m_p
//
// The field is either:
//   * uniform vector (props.E_uniform), or
//   * three cell-centred fields "E_x", "E_y", "E_z" looked up in the
//     FieldRegistry (set by an external Poisson solver or BC).
//
// Optional magnetic Lorentz term (B uniform) is supported via
//   F_mag = q_p · (v × B)
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "particles/LagrangianTracker.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

namespace simall::particles {

struct ChargedParticleProps {
    double      charge_per_parcel = 1.0e-12;  // [C]
    util::Vec3d E_uniform{0,0,0};             // V/m, used when useFieldRegistry==false
    util::Vec3d B_uniform{0,0,0};             // T,    Lorentz magnetic
    bool        useFieldRegistry = false;     // if true, read "E_x/y/z"
    bool        includeMagnetic  = false;
};

class ChargedParticleField {
public:
    void initialize(const meshing::Mesh& mesh, ChargedParticleProps props);

    /// Integrates ∆v = (q/m)(E + v × B) · dt for every active parcel.
    /// Returns the maximum |Δv| applied this step.
    double apply(double dt, LagrangianTracker& tracker,
                 solver::FieldRegistry& fields);

    const ChargedParticleProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    ChargedParticleProps p_{};
};

}  // namespace simall::particles
