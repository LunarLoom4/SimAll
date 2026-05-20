// Phase 10 (RADIATION) — model interfaces.
#pragma once
#include "solver/FieldRegistry.hpp"
namespace simall::radiation
{
enum class Model
{
    P1,
    DiscreteOrdinates,
    Rosseland,
    MonteCarlo,
    SurfaceToSurface
};
class RadiationModule
{
public:
    void configure(Model m) { model_ = m; }
    void solve(double dt, solver::FieldRegistry&);

private:
    Model model_ = Model::P1;
};
} // namespace simall::radiation
