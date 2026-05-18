// =============================================================================
// SimAll Beta - Lumped Zone Models Subsystem
// File   : src/zones/FanHexZone.cpp
// =============================================================================
#include "zones/FanHexZone.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::zones {

namespace {
util::Vec3d normalize(util::Vec3d v) {
    const double n = v.norm();
    return (n > 1e-30) ? util::Vec3d{v.x/n, v.y/n, v.z/n} : util::Vec3d{0,0,0};
}
}

// -----------------------------------------------------------------------------
//   FanZone
// -----------------------------------------------------------------------------
void FanZone::initialize(const meshing::Mesh& m, std::vector<FanZoneSpec> specs) {
    mesh_  = &m;
    specs_ = std::move(specs);
    const std::size_t nC = m.cells().size();
    cellZoneId_.assign(nC, 0);
    // Mesh storage does not carry per-cell zone IDs natively; the project
    // attaches cell-zone tags as a scalar registry field "cellZone" in
    // [io/Project]. We honour that convention if present, otherwise zone 0
    // (no membership) is assumed and the user is expected to drive the fan
    // through the throughFaceZone alone.
    throughFaces_.assign(specs_.size(), {});
    const auto& Ff = m.faces();
    for (std::size_t i = 0; i < specs_.size(); ++i) {
        for (std::size_t f = 0; f < Ff.size(); ++f) {
            if (Ff.neighbor[f] == meshing::kBoundaryCell
             && Ff.boundaryZone[f] == specs_[i].throughFaceZone) {
                throughFaces_[i].push_back(f);
            }
        }
    }
}

void FanZone::apply(solver::FieldRegistry& F) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    auto& Sfan = F.vector("S_FanMom", nC);
    std::fill(Sfan.x.begin(), Sfan.x.end(), 0.0);
    std::fill(Sfan.y.begin(), Sfan.y.end(), 0.0);
    std::fill(Sfan.z.begin(), Sfan.z.end(), 0.0);

    const auto* U = F.find_vector("U");
    if (!U) return;
    const auto& Ff = mesh_->faces();
    const auto& Cs = mesh_->cells();
    const auto* cellZone = F.find_scalar("cellZone");   // optional integer tags
    for (std::size_t i = 0; i < specs_.size(); ++i) {
        const auto& sp = specs_[i];
        const util::Vec3d axis = normalize(sp.axis);
        // Compute volumetric flow through the through-face zone (signed
        // along the fan axis).
        double Q = 0.0;
        for (meshing::FaceId f : throughFaces_[i]) {
            const meshing::CellId c = Ff.owner[f];
            const double un = U->x[c]*axis.x + U->y[c]*axis.y + U->z[c]*axis.z;
            const double A  = std::sqrt(Ff.areaX[f]*Ff.areaX[f]
                                       + Ff.areaY[f]*Ff.areaY[f]
                                       + Ff.areaZ[f]*Ff.areaZ[f]);
            Q += un * A;
        }
        const double dp = sp.curve.c0 + sp.curve.c1*Q + sp.curve.c2*Q*Q;
        // Distribute Δp as body force across cell zone.
        double zoneVol = 0.0;
        for (std::size_t c = 0; c < nC; ++c) {
            const auto z = cellZone ? static_cast<meshing::ZoneId>((*cellZone)[c]) : 0;
            if (z == sp.cellZone) zoneVol += Cs.volume[c];
        }
        if (zoneVol < 1e-30) continue;
        const double fmag = dp / zoneVol;       // N/m³
        for (std::size_t c = 0; c < nC; ++c) {
            const auto z = cellZone ? static_cast<meshing::ZoneId>((*cellZone)[c]) : 0;
            if (z != sp.cellZone) continue;
            Sfan.x[c] += fmag * axis.x;
            Sfan.y[c] += fmag * axis.y;
            Sfan.z[c] += fmag * axis.z;
        }
    }
}

// -----------------------------------------------------------------------------
//   HexZone
// -----------------------------------------------------------------------------
void HexZone::initialize(const meshing::Mesh& m, std::vector<HexZoneSpec> specs) {
    mesh_ = &m; specs_ = std::move(specs);
    cellZoneId_.assign(m.cells().size(), 0);
}

void HexZone::apply(solver::FieldRegistry& F) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    auto& Shx = F.scalar("S_HexEn",  nC);
    auto& Shm = F.vector("S_HexMom", nC);
    std::fill(Shx.begin(), Shx.end(), 0.0);
    std::fill(Shm.x.begin(), Shm.x.end(), 0.0);
    std::fill(Shm.y.begin(), Shm.y.end(), 0.0);
    std::fill(Shm.z.begin(), Shm.z.end(), 0.0);

    const auto* U = F.find_vector("U");
    const auto* T = F.find_scalar("T");
    const auto* cellZone = F.find_scalar("cellZone");
    if (!U || !T) return;
    for (const auto& sp : specs_) {
        for (std::size_t c = 0; c < nC; ++c) {
            const auto z = cellZone ? static_cast<meshing::ZoneId>((*cellZone)[c]) : 0;
            if (z != sp.cellZone) continue;
            const double Umag = std::sqrt(U->x[c]*U->x[c]+U->y[c]*U->y[c]+U->z[c]*U->z[c]);
            const double h = sp.h_a + sp.h_b * std::pow(std::max(Umag, 1e-12), sp.h_c);
            Shx[c] += h * (sp.T_ambient - (*T)[c]);
            const double F_darcy = sp.darcy + sp.forchheimer * Umag;
            Shm.x[c] -= F_darcy * U->x[c];
            Shm.y[c] -= F_darcy * U->y[c];
            Shm.z[c] -= F_darcy * U->z[c];
        }
    }
}

}  // namespace simall::zones
