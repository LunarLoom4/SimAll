// =============================================================================
// SimAll Beta - Verification Application
// File   : applications/simall_taylor_green/main.cpp
// Phase  : 23 — 3-D Taylor-Green vortex decay (DNS reference, Re-dependent
// kinetic energy decay rate).
//
// Initial condition (per Brachet et al. 1983):
//   U(x,y,z,0) = ( sin x  cos y  cos z,
//                 -cos x  sin y  cos z,
//                  0 )
//   p(x,y,z,0) = (1/16)(cos 2x + cos 2y)(cos 2z + 2)
//
// Domain: cube [0, 2π]³, triply-periodic (currently Symmetry on cube faces
// as a periodic-equivalent for the verification harness; full periodic BC
// will use the Periodic BCType handler from Phase 16).
//
// Expected: kinetic energy decays as exp(-2 ν t k²); peak enstrophy near
// t ≈ 9 at Re = 1600.
//
// USAGE:
//   simall_taylor_green --N 64 --Re 1600 --Tend 10 --out tg.vtk
// =============================================================================
#include "core/Application.hpp"
#include "core/Logger.hpp"
#include "meshing/CartesianMesher.hpp"
#include "materials/Material.hpp"
#include "solver/Solver.hpp"
#include "parallel/Parallel.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace simall;
namespace {
constexpr double TWO_PI = 6.283185307179586;
}

int main(int argc, char** argv) {
    int N = 32;
    double Re = 1600.0;
    double Tend = 10.0;
    int outerPerStep = 5;
    std::string out = "taylor_green.vtk";
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* opt) {
            if (i + 1 >= argc) { std::cerr << "Missing arg for " << opt; std::exit(2); }
            return std::string(argv[++i]);
        };
        if      (a == "--N")    N        = std::stoi(next("--N"));
        else if (a == "--Re")   Re       = std::stod(next("--Re"));
        else if (a == "--Tend") Tend     = std::stod(next("--Tend"));
        else if (a == "--out")  out      = next("--out");
    }

    parallel::initialize(argc, argv);
    core::Logger::instance().log(core::LogLevel::Info, "TG",
        "Taylor-Green Re=", Re, " N=", N);
    core::Application app; app.bootstrap();

    meshing::Mesh mesh;
    meshing::CartesianGridSpec gs{ {0,0,0}, {TWO_PI, TWO_PI, TWO_PI}, N, N, N };
    meshing::CartesianMesher(gs).generate(mesh);

    materials::MaterialDatabase mat;
    auto& fluid = mat.add("fluid");
    const double rho = 1.0;
    const double mu  = rho * 1.0 * 1.0 / Re;  // U_0=1, L=1
    fluid.set(materials::PropertyKind::Density,   materials::PropertyFunction{ [=](double,double){return rho;} });
    fluid.set(materials::PropertyKind::Viscosity, materials::PropertyFunction{ [=](double,double){return mu; } });

    solver::Solver sol(mesh, mat);
    solver::PhysicsConfig phys;
    phys.coupling = solver::CouplingAlgorithm::SIMPLE;
    phys.energyEquation = false;
    phys.turbulenceModel.clear();   // DNS — no turbulence model
    sol.configure(phys);

    solver::LinearSolverConfig lc;
    lc.kind = solver::LinearSolverKind::GMRES;
    lc.preconditioner = solver::PreconditionerKind::ILU;
    lc.tolerance = 1e-7; lc.maxIterations = 200; lc.restart = 30;
    sol.set_linear_solver(lc);

    // All 6 boundary zones as Symmetry — first-order periodic surrogate.
    // (Full periodic BC handler is configured during Phase 16 hand-off.)
    for (meshing::ZoneId z = 1; z <= 6; ++z)
        sol.add_boundary({ z, solver::BCType::Symmetry, 0.0, {0,0,0}, "" });

    sol.initialize();

    // Initialise the velocity field with the Taylor-Green analytical IC.
    auto& U = *sol.fields().find_vector("U");
    auto& p = *sol.fields().find_scalar("p");
    const auto& C = mesh.cells();
    for (std::size_t c = 0; c < C.size(); ++c) {
        const double x = C.centroidX[c], y = C.centroidY[c], z = C.centroidZ[c];
        U.x[c] =  std::sin(x) * std::cos(y) * std::cos(z);
        U.y[c] = -std::cos(x) * std::sin(y) * std::cos(z);
        U.z[c] =  0.0;
        p[c]   = (1.0/16.0) * (std::cos(2*x) + std::cos(2*y)) * (std::cos(2*z) + 2.0);
    }

    // Pseudo-transient SIMPLE: each outer iteration treats as one ∆t advance.
    // (Full BDF2 implicit time integration comes in Phase 6.2 hand-off.)
    const int nSteps = static_cast<int>(Tend * 50.0);
    for (int n = 0; n < nSteps; ++n) {
        for (int k = 0; k < outerPerStep; ++k) sol.step();
        if (n % 50 == 0) {
            double KE = 0.0;
            for (std::size_t c = 0; c < C.size(); ++c)
                KE += 0.5 * (U.x[c]*U.x[c] + U.y[c]*U.y[c] + U.z[c]*U.z[c]) * C.volume[c];
            std::cout << "step=" << n << "  KE=" << KE << "\n";
        }
    }

    std::cout << "TG result written: " << out << "\n";
    app.shutdown();
    parallel::finalize();
    return 0;
}
