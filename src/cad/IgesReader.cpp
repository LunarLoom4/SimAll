// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/IgesReader.cpp
// =============================================================================
#include "cad/IgesReader.hpp"
#include "cad/ShapeHandleInternal.hpp"
#include "cad/Healing.hpp"
#include "core/Logger.hpp"

#include <IGESControl_Reader.hxx>
#include <Interface_Static.hxx>
#include <IFSelect_ReturnStatus.hxx>

#include <filesystem>
#include <stdexcept>

namespace simall::cad {

ShapeHandle IgesReader::read(const std::string& path, const IgesReadOptions& opts) {
    report_ = {};

    if (!std::filesystem::exists(path))
        throw std::runtime_error("IgesReader: file not found: " + path);

    IGESControl_Reader reader;
    Interface_Static::SetIVal("read.iges.bspline.continuity", opts.fixContinuity ? 1 : 0);
    Interface_Static::SetIVal("read.iges.faulty.entities",    1);

    IFSelect_ReturnStatus status = reader.ReadFile(path.c_str());
    if (status != IFSelect_RetDone)
        throw std::runtime_error("IgesReader: read failed for " + path +
                                 " (status=" + std::to_string(int(status)) + ")");

    if (opts.readVisible) reader.SetReadVisible(Standard_True);

    Standard_Integer nRoots = reader.NbRootsForTransfer();
    report_.entityCount = std::size_t(nRoots);
    Standard_Integer transferred = reader.TransferRoots();
    report_.transferred = std::size_t(transferred);
    if (transferred == 0)
        throw std::runtime_error("IgesReader: nothing transferred from " + path);

    TopoDS_Shape shape = reader.OneShape();
    if (shape.IsNull())
        throw std::runtime_error("IgesReader: shape is null after transfer");

    auto h = makeHandle(std::move(shape));
    if (opts.heal) Healing{}.repair(h);

    SIMALL_LOG_INFO("CAD/IGES", "Loaded ", path, " — ", transferred, "/", nRoots,
                    " roots, ", h.topology().size(), " topo nodes");
    return h;
}

}  // namespace simall::cad
