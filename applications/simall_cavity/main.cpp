// =============================================================================
// SimAll Beta - Verification Application
// File   : applications/simall_cavity/main.cpp
// Phase  : 23 (verification suite — lid-driven cavity, Re = 100/400/1000).
//
// Headless CLI runner that:
//   1. builds an Nx × Ny × 1 Cartesian hex mesh of a unit square
//   2. assigns BCs (top = moving lid, others = no-slip walls)
//   3. configures the SIMPLE algorithm with ILU(0)-GMRES momentum + CG pressure
//   4. iterates to user-specified residual target
//   5. writes results to a Legacy VTK file readable by ParaView
//
// Compare midline velocity profiles against Ghia, Ghia & Shin (1982).
// =============================================================================
#include "core/Application.hpp"
#include "core/Logger.hpp"
#include "meshing/CartesianMesher.hpp"
#include "materials/Material.hpp"
#include "solver/Solver.hpp"
#include "parallel/Parallel.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

using namespace simall;

namespace {

void write_vtk_legacy(const std::string& path,
                      const meshing::Mesh& mesh,
                      const solver::FieldRegistry& fields) {
    std::ofstream f(path);
    const auto& N = mesh.nodes();
    const auto& C = mesh.cells();
    const auto& F = mesh.faces();
    f << "# vtk DataFile Version 3.0\nSimAll Cavity\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << N.size() << " double\n";
    for (std::size_t i = 0; i < N.size(); ++i)
        f << N.x[i] << ' ' << N.y[i] << ' ' << N.z[i] << '\n';

    // Reconstruct hex cell-vertex list by intersecting cell faces.
    // For Cartesian hexes, each cell has exactly 6 faces and 8 unique vertices.
    f << "CELLS " << C.size() << ' ' << C.size() * 9 << '\n';
    for (std::size_t c = 0; c < C.size(); ++c) {
        std::vector<int> v;
        for (int k = C.faceOffsets[c]; k < C.faceOffsets[c+1]; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            for (int q = F.nodeOffsets[fid]; q < F.nodeOffsets[fid+1]; ++q)
                v.push_back(static_cast<int>(F.nodeIndices[q]));
        }
        std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end());
        if (v.size() != 8) { f << v.size(); for (int n : v) f << ' ' << n; f << '\n'; continue; }
        // Order vertices into VTK_HEXAHEDRON convention (z then y then x).
        std::sort(v.begin(), v.end(), [&](int a, int b) {
            if (N.z[a] != N.z[b]) return N.z[a] < N.z[b];
            if (N.y[a] != N.y[b]) return N.y[a] < N.y[b];
            return N.x[a] < N.x[b];
        });
        // bottom: (0,0)(1,0)(1,1)(0,1) → indices 0,1,3,2
        f << "8 " << v[0] << ' ' << v[1] << ' ' << v[3] << ' ' << v[2]
                  << ' ' << v[4] << ' ' << v[5] << ' ' << v[7] << ' ' << v[6] << '\n';
    }
    f << "CELL_TYPES " << C.size() << '\n';
    for (std::size_t c = 0; c < C.size(); ++c) f << "12\n";  // VTK_HEXAHEDRON

    f << "CELL_DATA " << C.size() << '\n';
    if (auto* p = const_cast<solver::FieldRegistry&>(fields).find_scalar("p")) {
        f << "SCALARS pressure double 1\nLOOKUP_TABLE default\n";
        for (double v : *p) f << v << '\n';
    }
    if (auto* U = const_cast<solver::FieldRegistry&>(fields).find_vector("U")) {
        f << "VECTORS velocity double\n";
        for (std::size_t i = 0; i < U->size(); ++i)
            f << U->x[i] << ' ' << U->y[i] << ' ' << U->z[i] << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    int N = 64;
    double Re = 100.0;
    int maxIter = 2000;
    double target = 1.0e-5;
    std::string out = "cavity.vtk";
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* opt) {
            if (i + 1 >= argc) { std::cerr << "Missing arg for " << opt; std::exit(2); }
            return std::string(argv[++i]);
        };
        if      (a == "--N")        N        = std::stoi(next("--N"));
        else if (a == "--Re")       Re       = std::stod(next("--Re"));
        else if (a == "--iters")    maxIter  = std::stoi(next("--iters"));
        else if (a == "--target")   target   = std::stod(next("--target"));
        else if (a == "--out")      out      = next("--out");
    }

    parallel::initialize(argc, argv);
    core::Logger::instance().log(core::LogLevel::Info, "Cavity", "Re=", Re, " N=", N);
    core::Application app;
    app.bootstrap();

    // ---- Mesh ---------------------------------------------------------------
    meshing::Mesh mesh;
    meshing::CartesianGridSpec gs{ {0,0,0}, {1.0,1.0,0.1}, N, N, 1 };
    meshing::CartesianMesher(gs).generate(mesh);

    // ---- Material -----------------------------------------------------------
    materials::MaterialDatabase mat;
    auto& fluid = mat.add("fluid");
    const double rho = 1.0, U_lid = 1.0, L = 1.0;
    const double mu  = rho * U_lid * L / Re;
    fluid.set(materials::PropertyKind::Density,   materials::PropertyFunction{ [=](double,double){ return rho; } });
    fluid.set(materials::PropertyKind::Viscosity, materials::PropertyFunction{ [=](double,double){ return mu;  } });

    // ---- Solver -------------------------------------------------------------
    solver::Solver sol(mesh, mat);
    solver::PhysicsConfig phys;
    phys.coupling = solver::CouplingAlgorithm::SIMPLE;
    phys.energyEquation = false;
    phys.turbulenceModel.clear();
    sol.configure(phys);

    solver::LinearSolverConfig lc;
    lc.kind = solver::LinearSolverKind::GMRES;
    lc.preconditioner = solver::PreconditionerKind::ILU;
    lc.tolerance = 1e-6;
    lc.maxIterations = 200;
    lc.restart = 30;
    sol.set_linear_solver(lc);

    // Boundary zones from CartesianMesher: 1=xMin 2=xMax 3=yMin 4=yMax 5=zMin 6=zMax
    auto wall = [](meshing::ZoneId z){
        return solver::BoundarySpec{ z, solver::BCType::NoSlipWall, 0.0, {0,0,0}, "" };
    };
    sol.add_boundary(wall(1));
    sol.add_boundary(wall(2));
    sol.add_boundary(wall(3));
    sol.add_boundary(wall(5));
    sol.add_boundary(wall(6));
    // Top: moving lid
    solver::BoundarySpec lid{ 4, solver::BCType::MovingWall, 0.0, {U_lid, 0, 0}, "" };
    sol.add_boundary(lid);

    sol.initialize();
    for (int i = 0; i < maxIter; ++i) {
        sol.step();
    }

    write_vtk_legacy(out, mesh, sol.fields());
    std::cout << "Cavity result written: " << out << "\n";

    app.shutdown();
    parallel::finalize();
    return 0;
}
