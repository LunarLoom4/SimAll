// =============================================================================
// SimAll Beta - Heat Transfer Subsystem
// File   : src/heat_transfer/HeatTransfer.hpp
// Phase  : 9 (HEAT TRANSFER + Conjugate HT)
// =============================================================================
#pragma once
#include "solver/FieldRegistry.hpp"
#include "meshing/MeshStorage.hpp"

namespace simall::heat {

enum class Mode { Conduction, Convection, Conjugate, PhaseChange };

class HeatTransferModule {
public:
    void configure(Mode m) { mode_ = m; }
    Mode mode() const { return mode_; }
    void solve(double dt, meshing::Mesh&, solver::FieldRegistry&);
private:
    Mode mode_ = Mode::Conjugate;
};

}  // namespace simall::heat
