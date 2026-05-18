// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/Booleans.cpp
// =============================================================================
#include "cad/Booleans.hpp"
#include "cad/ShapeHandleInternal.hpp"
#include "core/Logger.hpp"

#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_BuilderAlgo.hxx>
#include <BOPAlgo_GlueEnum.hxx>
#include <TopTools_ListOfShape.hxx>
#include <stdexcept>

namespace simall::cad {

namespace {

template <class BopApi>
ShapeHandle runOp(const ShapeHandle& a, const ShapeHandle& b,
                  const BooleanOptions& o, const char* tag) {
    if (!a.valid() || !b.valid())
        throw std::runtime_error(std::string("Booleans/") + tag + ": invalid input");

    BopApi op;
    TopTools_ListOfShape la, lb;
    la.Append(ShapeHandleAccess::shape(a));
    lb.Append(ShapeHandleAccess::shape(b));
    op.SetArguments(la);
    op.SetTools(lb);
    if (o.fuzzyValue > 0.0) op.SetFuzzyValue(o.fuzzyValue);
    op.SetRunParallel(o.runParallel ? Standard_True : Standard_False);
    op.SetNonDestructive(o.nonDestructive ? Standard_True : Standard_False);
    op.SetGlue(o.gluePartial ? BOPAlgo_GlueShift : BOPAlgo_GlueOff);
    op.SetCheckInverted(o.checkInverted ? Standard_True : Standard_False);
    op.Build();
    if (op.HasErrors())
        throw std::runtime_error(std::string("Booleans/") + tag + ": OCC reported errors");
    TopoDS_Shape result = op.Shape();
    if (result.IsNull())
        throw std::runtime_error(std::string("Booleans/") + tag + ": null result");

    SIMALL_LOG_INFO("CAD/Bool", tag, " OK");
    return makeHandle(std::move(result));
}

}  // namespace

ShapeHandle Booleans::fuse  (const ShapeHandle& a, const ShapeHandle& b, const BooleanOptions& o)
{ return runOp<BRepAlgoAPI_Fuse>  (a, b, o, "fuse"); }
ShapeHandle Booleans::cut   (const ShapeHandle& a, const ShapeHandle& b, const BooleanOptions& o)
{ return runOp<BRepAlgoAPI_Cut>   (a, b, o, "cut"); }
ShapeHandle Booleans::common(const ShapeHandle& a, const ShapeHandle& b, const BooleanOptions& o)
{ return runOp<BRepAlgoAPI_Common>(a, b, o, "common"); }

ShapeHandle Booleans::fuseMany(const std::vector<const ShapeHandle*>& inputs,
                               const BooleanOptions& o) {
    if (inputs.empty()) throw std::runtime_error("Booleans/fuseMany: empty input set");
    if (inputs.size() == 1) {
        TopoDS_Shape s = ShapeHandleAccess::shape(*inputs.front());
        return makeHandle(std::move(s));
    }

    BRepAlgoAPI_BuilderAlgo op;
    TopTools_ListOfShape args;
    for (auto* h : inputs) {
        if (!h || !h->valid())
            throw std::runtime_error("Booleans/fuseMany: invalid input shape");
        args.Append(ShapeHandleAccess::shape(*h));
    }
    op.SetArguments(args);
    if (o.fuzzyValue > 0.0) op.SetFuzzyValue(o.fuzzyValue);
    op.SetRunParallel(o.runParallel ? Standard_True : Standard_False);
    op.SetNonDestructive(o.nonDestructive ? Standard_True : Standard_False);
    op.Build();
    if (op.HasErrors())
        throw std::runtime_error("Booleans/fuseMany: OCC reported errors");
    TopoDS_Shape result = op.Shape();
    if (result.IsNull())
        throw std::runtime_error("Booleans/fuseMany: null result");

    SIMALL_LOG_INFO("CAD/Bool", "fuseMany OK (", inputs.size(), " inputs)");
    return makeHandle(std::move(result));
}

}  // namespace simall::cad
