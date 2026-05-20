// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/CadKernel.cpp
// Phase  : 4 — public facade.
//
// CadKernel is now a thin orchestration layer: it dispatches imports to the
// dedicated reader modules (StepReader / IgesReader / StlSurfaceReader /
// BrepIo) and delegates healing / tessellation to their respective modules.
// Format extension detection lives here so callers retain the single-entry
// `kernel.import("foo.step")` ergonomics.
// =============================================================================
#include "cad/CadKernel.hpp"

#include "cad/BrepIo.hpp"
#include "cad/Healing.hpp"
#include "cad/IgesReader.hpp"
#include "cad/ObjReader.hpp"
#include "cad/ShapeHandleInternal.hpp"
#include "cad/StepReader.hpp"
#include "cad/StlSurfaceReader.hpp"
#include "cad/Tessellator.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>

#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

namespace simall::cad
{

// ---- ShapeHandle PIMPL plumbing --------------------------------------------
ShapeHandle::ShapeHandle() : impl_(std::make_unique<Impl>()) {}
ShapeHandle::~ShapeHandle() = default;
ShapeHandle::ShapeHandle(ShapeHandle&&) noexcept = default;
ShapeHandle& ShapeHandle::operator=(ShapeHandle&&) noexcept = default;

bool ShapeHandle::valid() const noexcept
{
    return impl_ && !impl_->shape.IsNull();
}
const TopologyGraph& ShapeHandle::topology() const noexcept
{
    return impl_->graph;
}
TopologyGraph& ShapeHandle::topology() noexcept
{
    return impl_->graph;
}

util::BoundingBox ShapeHandle::bounds() const
{
    util::BoundingBox bb;
    if (!valid())
        return bb;
    Bnd_Box b;
    BRepBndLib::Add(impl_->shape, b);
    if (b.IsVoid())
        return bb;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    b.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    bb.expand({xmin, ymin, zmin});
    bb.expand({xmax, ymax, zmax});
    return bb;
}

// ---- CadKernel orchestrator ------------------------------------------------
CadKernel::CadKernel() = default;
CadKernel::~CadKernel() = default;

namespace
{
std::string lowerExt(const std::string& path)
{
    auto ext = std::filesystem::path(path).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext;
}
} // namespace

ShapeHandle CadKernel::import(const std::string& path)
{
    const auto ext = lowerExt(path);
    ShapeHandle h;
    if (ext == ".step" || ext == ".stp")
        h = StepReader{}.read(path);
    else if (ext == ".iges" || ext == ".igs")
        h = IgesReader{}.read(path);
    else if (ext == ".stl")
        h = StlSurfaceReader{}.read(path);
    else if (ext == ".obj")
        h = ObjReader{}.read(path);
    else if (ext == ".brep" || ext == ".brp")
        h = BrepIo{}.read(path);
    else
        throw std::runtime_error("Unsupported CAD format: " + ext);

    SIMALL_LOG_INFO("CAD", "Imported ", path, " (", h.topology().size(), " topo nodes)");
    return h;
}

void CadKernel::heal(ShapeHandle& shape, const HealingOptions& opts)
{
    Healing{}.repair(shape, opts);
}

TriangleMesh CadKernel::tessellate(const ShapeHandle& shape, const TessellationParams& p)
{
    return Tessellator{}.tessellate(shape, p);
}

TriangleMesh CadKernel::adaptive_tessellate(const ShapeHandle& shape,
                                            double maxEdge,
                                            double epsilon)
{
    return Tessellator{}.adaptive(shape, maxEdge, epsilon);
}

} // namespace simall::cad
