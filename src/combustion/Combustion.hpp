// Phase 11 (COMBUSTION).
#pragma once
#include "solver/FieldRegistry.hpp"

#include <string>
#include <vector>
namespace simall::combustion
{

struct ReactionStep
{
    std::string equation;
    double Ea, A, beta;
};
struct Mechanism
{
    std::vector<ReactionStep> steps;
};

enum class Model
{
    SpeciesTransport,
    FiniteRate,
    EddyDissipation,
    Flamelet,
    PDF,
    Premixed,
    NonPremixed,
    NOx,
    Soot
};

class CombustionModule
{
public:
    void load_chemkin(const std::string& path);
    void configure(Model m) { model_ = m; }
    void solve(double dt, solver::FieldRegistry&);

private:
    Model model_ = Model::EddyDissipation;
    Mechanism mech_;
};
} // namespace simall::combustion
