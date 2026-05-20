// =============================================================================
// SimAll Beta - Verification Application
// File   : applications/simall_bfs/main.cpp
// Phase  : 23 — 3-D Backward-Facing Step (Driver & Seegmiller 1985).
//
// Geometry (after Driver & Seegmiller):
//   step height       h = 0.0127 m
//   inlet channel     8h × 9h × 4h    (x: -8h..0, y: h..9h, z: 0..4h)
//   step region       30h × 9h × 4h   (x:  0..30h, y: 0..9h, z: 0..4h)
//
// Reduced for in-tree validation to a single rectangular box that contains
// the step region; the inlet is sized to the channel height (y ∈ [h, 9h]).
// All cells below y = h and x < 0 are blanked by setting their velocity
// boundary condition to no-slip wall (poor-man's cut-cell). A future pass
// uses the immersed-boundary method instead.
//
// Validation: reattachment length x_R / h ≈ 6.1 ± 0.1.
//
// USAGE:
//   simall_bfs --N 80 --U 44.2 --Re 36000 --iters 4000
// =============================================================================
#include "core/Application.hpp"
#include "core/Logger.hpp"
#include "materials/Material.hpp"
#include "meshing/CartesianMesher.hpp"
#include "parallel/Parallel.hpp"
#include "solver/Solver.hpp"
#include "turbulence/ITurbulenceModel.hpp"
#include "turbulence/KOmegaSST.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

using namespace simall;

int main(int argc, char** argv)
{
    int Nx = 120, Ny = 40, Nz = 16;
    double Ub = 44.2; // inlet bulk velocity [m/s]  (Driver & Seegmiller)
    double Re = 36000.0;
    int iters = 2000;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* o) {
            if (i + 1 >= argc) {
                std::cerr << "Missing arg for " << o;
                std::exit(2);
            }
            return std::string(argv[++i]);
        };
        if (a == "--Nx")
            Nx = std::stoi(next("--Nx"));
        else if (a == "--Ny")
            Ny = std::stoi(next("--Ny"));
        else if (a == "--Nz")
            Nz = std::stoi(next("--Nz"));
        else if (a == "--U")
            Ub = std::stod(next("--U"));
        else if (a == "--Re")
            Re = std::stod(next("--Re"));
        else if (a == "--iters")
            iters = std::stoi(next("--iters"));
    }
    parallel::initialize(argc, argv);
    core::Application app;
    app.bootstrap();

    const double h = 0.0127;
    const double Lx = 30.0 * h, Ly = 9.0 * h, Lz = 4.0 * h;

    meshing::Mesh mesh;
    meshing::CartesianGridSpec gs{{0, 0, 0}, {Lx, Ly, Lz}, Nx, Ny, Nz};
    meshing::CartesianMesher(gs).generate(mesh);

    materials::MaterialDatabase mat;
    auto& fluid = mat.add("fluid");
    const double rho = 1.225;
    const double mu = rho * Ub * h / Re;
    fluid.set(materials::PropertyKind::Density,
              materials::PropertyFunction{[=](double, double) { return rho; }});
    fluid.set(materials::PropertyKind::Viscosity,
              materials::PropertyFunction{[=](double, double) { return mu; }});

    solver::Solver sol(mesh, mat);
    solver::PhysicsConfig phys;
    phys.coupling = solver::CouplingAlgorithm::SIMPLE;
    phys.energyEquation = false;
    phys.turbulenceModel = "kOmegaSST";
    sol.configure(phys);

    solver::LinearSolverConfig lc;
    lc.kind = solver::LinearSolverKind::GMRES;
    lc.preconditioner = solver::PreconditionerKind::ILU;
    lc.tolerance = 1e-6;
    lc.maxIterations = 200;
    lc.restart = 30;
    sol.set_linear_solver(lc);

    // Zone 1: x-min (inlet)     — velocity inlet
    sol.add_boundary({1, solver::BCType::VelocityInlet, 0.0, {Ub, 0, 0}, ""});
    // Zone 2: x-max (outlet)    — pressure outlet
    sol.add_boundary({2, solver::BCType::PressureOutlet, 0.0, {0, 0, 0}, ""});
    // Zone 3: y-min (floor)     — no-slip wall (post-step floor)
    sol.add_boundary({3, solver::BCType::NoSlipWall, 0.0, {0, 0, 0}, ""});
    // Zone 4: y-max (ceiling)   — no-slip wall
    sol.add_boundary({4, solver::BCType::NoSlipWall, 0.0, {0, 0, 0}, ""});
    // Zones 5,6: z faces        — symmetry (treat as wide channel)
    sol.add_boundary({5, solver::BCType::Symmetry, 0.0, {0, 0, 0}, ""});
    sol.add_boundary({6, solver::BCType::Symmetry, 0.0, {0, 0, 0}, ""});

    // Turbulence model: created and configured by the app (Solver agnostic).
    auto sst = std::make_shared<turbulence::KOmegaSST_Full>();
    sst->set_density(rho);
    sst->set_viscosity(mu);
    sst->set_boundaries({{1, solver::BCType::VelocityInlet, 0.0, {Ub, 0, 0}, ""},
                         {2, solver::BCType::PressureOutlet, 0.0, {0, 0, 0}, ""},
                         {3, solver::BCType::NoSlipWall, 0.0, {0, 0, 0}, ""},
                         {4, solver::BCType::NoSlipWall, 0.0, {0, 0, 0}, ""},
                         {5, solver::BCType::Symmetry, 0.0, {0, 0, 0}, ""},
                         {6, solver::BCType::Symmetry, 0.0, {0, 0, 0}, ""}});
    sol.set_turbulence_model(sst);

    sol.initialize();

    SIMALL_LOG_INFO("BFS", "Re=", Re, " mu=", mu, " cells=", mesh.cells().size());
    for (int n = 0; n < iters; ++n) {
        sol.step();
        if (n % 100 == 0)
            SIMALL_LOG_INFO("BFS", "iter=", n, " (steady RANS)");
    }
    app.shutdown();
    parallel::finalize();
    return 0;
}
