// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/BrepIo.cpp
// =============================================================================
#include "cad/BrepIo.hpp"

#include "cad/ShapeHandleInternal.hpp"
#include "core/Logger.hpp"

#include <filesystem>
#include <stdexcept>

#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <TopoDS_Shape.hxx>

namespace simall::cad
{

ShapeHandle BrepIo::read(const std::string& path)
{
    if (!std::filesystem::exists(path))
        throw std::runtime_error("BrepIo: file not found: " + path);

    TopoDS_Shape shape;
    BRep_Builder builder;
    if (!BRepTools::Read(shape, path.c_str(), builder))
        throw std::runtime_error("BrepIo: BRepTools::Read failed for " + path);
    if (shape.IsNull())
        throw std::runtime_error("BrepIo: read produced a null shape");

    auto h = makeHandle(std::move(shape));
    SIMALL_LOG_INFO("CAD/BREP", "Loaded ", path, " (", h.topology().size(), " topo nodes)");
    return h;
}

bool BrepIo::write(const ShapeHandle& shape, const std::string& path) const
{
    if (!shape.valid())
        return false;
    try {
        return BRepTools::Write(ShapeHandleAccess::shape(shape), path.c_str()) != Standard_False;
    } catch (...) {
        return false;
    }
}

} // namespace simall::cad
