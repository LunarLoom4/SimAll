// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/IgesReader.hpp
// Phase  : 4.1 (CAD readers — IGES legacy)
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"

#include <string>
#include <vector>

namespace simall::cad
{

struct IgesReadOptions
{
    bool readNames = true;
    bool readBSpline = true;
    bool readVisible = true;
    bool heal = true; // IGES files almost always need it
    bool fixContinuity = true;
};

struct IgesReadReport
{
    std::size_t entityCount = 0;
    std::size_t transferred = 0;
    std::vector<std::string> warnings;
};

class IgesReader
{
public:
    ShapeHandle read(const std::string& path, const IgesReadOptions& opts = {});
    IgesReadReport lastReport() const noexcept { return report_; }

private:
    IgesReadReport report_;
};

} // namespace simall::cad
