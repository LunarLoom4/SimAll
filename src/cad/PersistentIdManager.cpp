// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/PersistentIdManager.cpp
//
// Matching algorithm
// ------------------
// 1. snapshot(): walks every node in the topology graph that has a non-null
//    occHandle. For each, computes:
//       - centroid via BRepGProp on a face/edge/solid;
//       - bbox via BRepBndLib;
//       - measure (area for faces, length for edges, volume for solids,
//         0 for vertices).
//
// 2. remapTo(): for each entity in the new topology, computes the same
//    fingerprint, then searches snapshot prints of equal `type` for one
//    whose centroid lies within `positionTol` AND whose measure is within
//    `measureRelTol` relative deviation. Ambiguous matches (more than one
//    candidate satisfying the test) leave the new id fresh and the old id
//    becomes "unmatched" so the caller (named-selection manager) can warn.
//
// Complexity: O(N²) acceptable for typical part counts (≤ 10⁵ faces).
// =============================================================================
#include "cad/PersistentIdManager.hpp"
#include "cad/ShapeHandleInternal.hpp"
#include "core/Logger.hpp"

#include <BRep_Tool.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Pnt.hxx>

#include <cmath>
#include <limits>

namespace simall::cad {

namespace {

util::BoundingBox occBox(const TopoDS_Shape& s) {
    Bnd_Box bb;
    BRepBndLib::Add(s, bb);
    util::BoundingBox out;
    if (bb.IsVoid()) return out;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    bb.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    out.expand({xmin, ymin, zmin});
    out.expand({xmax, ymax, zmax});
    return out;
}

PersistentFingerprint fingerprintFor(const TopoDS_Shape& s,
                                     TopologyType        type,
                                     util::PersistentId  id) {
    PersistentFingerprint p;
    p.id   = id;
    p.type = type;
    p.bbox = occBox(s);

    GProp_GProps props;
    switch (type) {
        case TopologyType::Vertex: {
            gp_Pnt v = BRep_Tool::Pnt(TopoDS::Vertex(s));
            p.centroid = {v.X(), v.Y(), v.Z()};
            p.measure  = 0.0;
            break;
        }
        case TopologyType::Edge:
            BRepGProp::LinearProperties(s, props);
            p.measure  = props.Mass();
            p.centroid = {props.CentreOfMass().X(), props.CentreOfMass().Y(), props.CentreOfMass().Z()};
            break;
        case TopologyType::Face:
            BRepGProp::SurfaceProperties(s, props);
            p.measure  = props.Mass();
            p.centroid = {props.CentreOfMass().X(), props.CentreOfMass().Y(), props.CentreOfMass().Z()};
            break;
        case TopologyType::Solid:
            BRepGProp::VolumeProperties(s, props);
            p.measure  = props.Mass();
            p.centroid = {props.CentreOfMass().X(), props.CentreOfMass().Y(), props.CentreOfMass().Z()};
            break;
        default:
            p.centroid = p.bbox.center();
            p.measure  = 0.0;
            break;
    }
    return p;
}

}  // namespace

void PersistentIdManager::snapshot(const ShapeHandle& shape) {
    prints_.clear();
    if (!shape.valid()) return;
    const auto& g = shape.topology();
    prints_.reserve(g.size());
    for (auto const& [id, n] : g.nodes()) {
        if (!n.occHandle) continue;
        const TopoDS_Shape& s = *static_cast<const TopoDS_Shape*>(n.occHandle);
        prints_.push_back(fingerprintFor(s, n.type, id));
    }
}

RemapResult PersistentIdManager::remapTo(ShapeHandle& newShape,
                                         double positionTol,
                                         double measureRelTol) {
    RemapResult out;
    if (!newShape.valid()) return out;

    auto& g = newShape.topology();

    // Build fingerprints for the new topology BEFORE any id mutation.
    struct NewEntry {
        util::PersistentId id;   // current id in new graph (will be replaced)
        PersistentFingerprint fp;
    };
    std::vector<NewEntry> news;
    news.reserve(g.size());
    for (auto const& [id, n] : g.nodes()) {
        if (!n.occHandle) continue;
        const TopoDS_Shape& s = *static_cast<const TopoDS_Shape*>(n.occHandle);
        news.push_back({id, fingerprintFor(s, n.type, id)});
    }

    // Track which old prints are still available.
    std::vector<bool> oldClaimed(prints_.size(), false);

    // Build the per-new-id reassignment plan.
    std::unordered_map<util::PersistentId, util::PersistentId> newToOld;
    newToOld.reserve(news.size());

    for (auto const& e : news) {
        std::size_t best     = std::numeric_limits<std::size_t>::max();
        double      bestDist = positionTol;
        int         tied     = 0;

        for (std::size_t i = 0; i < prints_.size(); ++i) {
            if (oldClaimed[i]) continue;
            auto const& op = prints_[i];
            if (op.type != e.fp.type) continue;

            // Measure check: 0-measure (vertices) skip ratio test.
            if (op.measure > 0.0 && e.fp.measure > 0.0) {
                double rel = std::abs(op.measure - e.fp.measure) /
                             std::max(op.measure, e.fp.measure);
                if (rel > measureRelTol) continue;
            }
            double dx = op.centroid.x - e.fp.centroid.x;
            double dy = op.centroid.y - e.fp.centroid.y;
            double dz = op.centroid.z - e.fp.centroid.z;
            double d  = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (d > positionTol) continue;

            if (d < bestDist - 1e-12) {
                best = i; bestDist = d; tied = 1;
            } else if (std::abs(d - bestDist) <= 1e-12) {
                ++tied;
            }
        }

        if (best != std::numeric_limits<std::size_t>::max() && tied == 1) {
            oldClaimed[best] = true;
            newToOld[e.id]   = prints_[best].id;
            out.oldToNew[prints_[best].id] = prints_[best].id;  // id preserved
        }
    }

    // Apply rewrite: replace each new node id with the matched old id where
    // available. TopologyGraph stores nodes_ keyed by id, so we rebuild the
    // underlying map. Uses the mutable_nodes() accessor exposed for this
    // purpose by Week 3.
    auto& nodes = g.mutable_nodes();
    std::unordered_map<util::PersistentId, TopologyNode> rebuilt;
    rebuilt.reserve(nodes.size());
    for (auto& [id, n] : nodes) {
        util::PersistentId newId = id;
        auto it = newToOld.find(id);
        if (it != newToOld.end()) {
            newId = it->second;
            n.id  = newId;
        } else {
            out.freshlyAllocated.push_back(id);
        }
        rebuilt.emplace(newId, std::move(n));
    }
    nodes.swap(rebuilt);

    for (std::size_t i = 0; i < prints_.size(); ++i) {
        if (!oldClaimed[i]) out.unmatchedOld.push_back(prints_[i].id);
    }

    SIMALL_LOG_INFO("CAD/IdMgr",
                    "remap: ", out.oldToNew.size(), " preserved, ",
                    out.unmatchedOld.size(), " orphaned, ",
                    out.freshlyAllocated.size(), " fresh");
    return out;
}

}  // namespace simall::cad
