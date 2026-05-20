// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/AmrAdaptor.cpp
// =============================================================================
#include "meshing/AmrAdaptor.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace simall::meshing
{

void AmrAdaptor::initialize(Mesh& mesh, AmrProps props)
{
    mesh_ = &mesh;
    p_ = props;
    level_.clear();
    parent_.clear();
    SIMALL_LOG_INFO("Meshing",
                    "AmrAdaptor init: maxLevel=",
                    p_.maxLevel,
                    " refineThr=",
                    p_.refineThreshold,
                    " coarsenThr=",
                    p_.coarsenThreshold);
}

int AmrAdaptor::cell_level(CellId c) const
{
    const auto it = level_.find(c);
    return it == level_.end() ? 0 : it->second;
}

void AmrAdaptor::enforce_balance(std::vector<int>& flag, AmrReport& rep)
{
    if (!p_.balance21 || !mesh_)
        return;
    const auto& Ff = mesh_->faces();
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::size_t f = 0; f < Ff.size(); ++f) {
            const auto o = Ff.owner[f];
            const auto n = Ff.neighbor[f];
            if (n == kBoundaryCell)
                continue;
            const int lo = cell_level(o) + (flag[o] > 0 ? 1 : 0);
            const int ln = cell_level(n) + (flag[n] > 0 ? 1 : 0);
            if (lo - ln >= 2 && flag[n] <= 0) {
                flag[n] = 1;
                ++rep.cellsBalanced;
                changed = true;
            } else if (ln - lo >= 2 && flag[o] <= 0) {
                flag[o] = 1;
                ++rep.cellsBalanced;
                changed = true;
            }
        }
    }
}

AmrReport AmrAdaptor::adapt(solver::FieldRegistry& fields, const std::string& tagFieldName)
{
    AmrReport rep{};
    if (!mesh_)
        return rep;
    auto* tagPtr = fields.find_scalar(tagFieldName);
    if (!tagPtr) {
        SIMALL_LOG_INFO("Meshing", "AmrAdaptor: tag field '", tagFieldName, "' not found");
        return rep;
    }
    auto& C = mesh_->cells();
    const std::size_t nC = C.size();
    const auto& tag = *tagPtr;
    if (tag.size() != nC) {
        SIMALL_LOG_INFO("Meshing", "AmrAdaptor: tag size mismatch ", tag.size(), " vs cells ", nC);
        return rep;
    }
    // 1 = refine, -1 = coarsen, 0 = keep.
    std::vector<int> flag(nC, 0);
    for (std::size_t c = 0; c < nC; ++c) {
        const double v = std::abs(tag[c]);
        if (v >= p_.refineThreshold && cell_level(c) < p_.maxLevel)
            flag[c] = 1;
        else if (v <= p_.coarsenThreshold && cell_level(c) > 0)
            flag[c] = -1;
    }
    enforce_balance(flag, rep);

    // Apply refinement: append 8 child cells per marked parent.  The
    // parent cell is retained with zero volume to keep CellId stable for
    // downstream tools; tags carrying field data are migrated to children
    // by uniform copy in subsequent solver iterations.
    const std::size_t parentBase = C.volume.size();
    for (std::size_t c = 0; c < parentBase; ++c) {
        if (flag[c] != 1)
            continue;
        const double vol = C.volume[c];
        const double childVol = vol * 0.125;
        const double cx = C.centroidX[c];
        const double cy = C.centroidY[c];
        const double cz = C.centroidZ[c];
        const double h = std::cbrt(vol) * 0.25;
        for (int k = 0; k < 8; ++k) {
            const double sx = (k & 1) ? +h : -h;
            const double sy = (k & 2) ? +h : -h;
            const double sz = (k & 4) ? +h : -h;
            const CellId child = static_cast<CellId>(C.volume.size());
            C.volume.push_back(childVol);
            C.centroidX.push_back(cx + sx);
            C.centroidY.push_back(cy + sy);
            C.centroidZ.push_back(cz + sz);
            // Empty per-child face CSR — face topology owner re-stitch
            // belongs to a downstream refinement remesher.
            C.faceOffsets.push_back(C.faceOffsets.back());
            parent_[child] = static_cast<CellId>(c);
            level_[child] = cell_level(static_cast<CellId>(c)) + 1;
        }
        C.volume[c] = 0.0; // mark parent inactive
        ++rep.cellsRefined;
    }

    // Apply coarsening: when all 8 siblings carry -1, restore parent
    // volume to Σ child vol and zero out siblings.
    std::unordered_map<CellId, std::vector<CellId>> kids;
    for (const auto& [child, par] : parent_)
        kids[par].push_back(child);
    for (auto& [par, siblings] : kids) {
        if (siblings.size() < 8)
            continue;
        bool allFlagged = true;
        double sumV = 0.0;
        for (auto s : siblings) {
            if (s >= flag.size() || flag[s] != -1) {
                allFlagged = false;
                break;
            }
            sumV += C.volume[s];
        }
        if (!allFlagged)
            continue;
        for (auto s : siblings)
            C.volume[s] = 0.0;
        C.volume[par] = sumV;
        ++rep.cellsCoarsened;
    }
    SIMALL_LOG_INFO("Meshing",
                    "AmrAdapt: refined=",
                    rep.cellsRefined,
                    " coarsened=",
                    rep.cellsCoarsened,
                    " balanced=",
                    rep.cellsBalanced,
                    " newCellCount=",
                    C.volume.size());
    return rep;
}

} // namespace simall::meshing
