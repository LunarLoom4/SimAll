#include "heat_transfer/HeatTransfer.hpp"
namespace simall::heat {
void HeatTransferModule::solve(double, meshing::Mesh&, solver::FieldRegistry&) {
    // CHT interface flux continuity is enforced in Phase 9 implementation.
}
}
