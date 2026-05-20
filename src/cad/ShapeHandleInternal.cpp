// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/ShapeHandleInternal.cpp
// =============================================================================
#include "cad/ShapeHandleInternal.hpp"

#include <TopExp_Explorer.hxx>

namespace simall::cad
{

TopologyType occKindToTopologyType(TopAbs_ShapeEnum e) noexcept
{
    switch (e) {
    case TopAbs_VERTEX:
        return TopologyType::Vertex;
    case TopAbs_EDGE:
        return TopologyType::Edge;
    case TopAbs_WIRE:
        return TopologyType::Wire;
    case TopAbs_FACE:
        return TopologyType::Face;
    case TopAbs_SHELL:
        return TopologyType::Shell;
    case TopAbs_SOLID:
        return TopologyType::Solid;
    default:
        return TopologyType::Compound;
    }
}

void rebuildTopologyGraph(const TopoDS_Shape& shape, TopologyGraph& g)
{
    g = {};
    if (shape.IsNull())
        return;
    auto traverse = [&](TopAbs_ShapeEnum kind) {
        for (TopExp_Explorer it(shape, kind); it.More(); it.Next()) {
            const TopoDS_Shape& s = it.Current();
            g.add(occKindToTopologyType(kind), const_cast<TopoDS_Shape*>(&s));
        }
    };
    traverse(TopAbs_SOLID);
    traverse(TopAbs_SHELL);
    traverse(TopAbs_FACE);
    traverse(TopAbs_EDGE);
    traverse(TopAbs_VERTEX);
}

ShapeHandle makeHandle(TopoDS_Shape shape)
{
    ShapeHandle h;
    h.impl_->shape = std::move(shape);
    rebuildTopologyGraph(h.impl_->shape, h.impl_->graph);
    return h;
}

} // namespace simall::cad
