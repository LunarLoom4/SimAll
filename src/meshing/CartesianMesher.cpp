// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/CartesianMesher.cpp
// =============================================================================
#include "meshing/CartesianMesher.hpp"

#include "core/Logger.hpp"
#include "meshing/Connectivity.hpp"

namespace simall::meshing
{

void CartesianMesher::generate(Mesh& mesh)
{
    const int Nx = spec_.Nx, Ny = spec_.Ny, Nz = spec_.Nz;
    const double dx = spec_.extent.x / Nx;
    const double dy = spec_.extent.y / Ny;
    const double dz = spec_.extent.z / Nz;

    NodeStorage nodes;
    const int Px = Nx + 1, Py = Ny + 1, Pz = Nz + 1;
    nodes.reserve(std::size_t(Px) * Py * Pz);
    auto pid = [&](int i, int j, int k) { return NodeId(((k * Py) + j) * Px + i); };
    for (int k = 0; k < Pz; ++k)
        for (int j = 0; j < Py; ++j)
            for (int i = 0; i < Px; ++i) {
                nodes.x.push_back(spec_.origin.x + i * dx);
                nodes.y.push_back(spec_.origin.y + j * dy);
                nodes.z.push_back(spec_.origin.z + k * dz);
            }

    std::vector<CellDescriptor> cells;
    cells.reserve(std::size_t(Nx) * Ny * Nz);

    // VTK_HEXAHEDRON vertex ordering, outward-normal faces.
    // n0 = (i,j,k),     n1 = (i+1,j,k),   n2 = (i+1,j+1,k),  n3 = (i,j+1,k)
    // n4 = (i,j,k+1),   n5 = (i+1,j,k+1), n6 = (i+1,j+1,k+1),n7 = (i,j+1,k+1)
    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                NodeId n0 = pid(i, j, k);
                NodeId n1 = pid(i + 1, j, k);
                NodeId n2 = pid(i + 1, j + 1, k);
                NodeId n3 = pid(i, j + 1, k);
                NodeId n4 = pid(i, j, k + 1);
                NodeId n5 = pid(i + 1, j, k + 1);
                NodeId n6 = pid(i + 1, j + 1, k + 1);
                NodeId n7 = pid(i, j + 1, k + 1);

                CellDescriptor cd;
                cd.faces = {
                    {n0, n3, n2, n1}, // -Z (back)
                    {n4, n5, n6, n7}, // +Z (front)
                    {n0, n1, n5, n4}, // -Y (bottom)
                    {n3, n7, n6, n2}, // +Y (top)
                    {n0, n4, n7, n3}, // -X (left)
                    {n1, n2, n6, n5}, // +X (right)
                };
                cells.push_back(std::move(cd));
            }

    ConnectivityBuilder::build(mesh, nodes, cells);

    // Boundary-zone classification using face centroid coordinates.
    const auto& F = mesh.faces();
    const ZoneId zLeft = mesh.add_zone("xMin", true);
    const ZoneId zRight = mesh.add_zone("xMax", true);
    const ZoneId zBottom = mesh.add_zone("yMin", true);
    const ZoneId zTop = mesh.add_zone("yMax", true);
    const ZoneId zBack = mesh.add_zone("zMin", true);
    const ZoneId zFront = mesh.add_zone("zMax", true);
    const double tol = 1e-9;
    for (std::size_t f = 0; f < F.size(); ++f) {
        if (F.neighbor[f] != kBoundaryCell)
            continue;
        const double cx = F.centroidX[f], cy = F.centroidY[f], cz = F.centroidZ[f];
        if (std::abs(cx - spec_.origin.x) < tol)
            mesh.faces().boundaryZone[f] = zLeft;
        else if (std::abs(cx - (spec_.origin.x + spec_.extent.x)) < tol)
            mesh.faces().boundaryZone[f] = zRight;
        else if (std::abs(cy - spec_.origin.y) < tol)
            mesh.faces().boundaryZone[f] = zBottom;
        else if (std::abs(cy - (spec_.origin.y + spec_.extent.y)) < tol)
            mesh.faces().boundaryZone[f] = zTop;
        else if (std::abs(cz - spec_.origin.z) < tol)
            mesh.faces().boundaryZone[f] = zBack;
        else if (std::abs(cz - (spec_.origin.z + spec_.extent.z)) < tol)
            mesh.faces().boundaryZone[f] = zFront;
    }

    SIMALL_LOG_INFO("Mesh",
                    "Cartesian grid generated: ",
                    Nx,
                    "x",
                    Ny,
                    "x",
                    Nz,
                    " (",
                    mesh.cells().size(),
                    " cells)");
}

} // namespace simall::meshing
