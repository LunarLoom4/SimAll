// Phase 12 (MULTIPHASE).
#pragma once
#include "solver/FieldRegistry.hpp"
namespace simall::multiphase {
enum class Model { VOF, Mixture, Eulerian, DPM, Cavitation, Spray };
class MultiphaseModule {
public:
    void configure(Model m) { model_ = m; }
    void solve(double dt, solver::FieldRegistry&);
private:
    Model model_ = Model::VOF;
};
}
