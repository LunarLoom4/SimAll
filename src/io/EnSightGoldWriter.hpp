// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/EnSightGoldWriter.hpp
// Phase  : 4.5 — EnSight Gold (CEI) "C Binary" exporter.
//
// EnSight Gold is the de-facto industrial CFD post-processing format
// supported by ParaView, EnSight, FieldView, Tecplot 360, and most
// downstream visualisation tools.
//
// Layout produced (matches CEI EnSight Gold v6.x specification):
//
//   <case>.case               — text case file (sections FORMAT,
//                                GEOMETRY, VARIABLE, TIME)
//   <case>.geo                — binary geometry (header, part, coordinates,
//                                element block)
//   <case>.<var>.<step>       — binary per-variable file
//
// Supports per-cell scalar and vector fields; arbitrary polyhedral
// connectivity emitted via the "nfaced" element block.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <string>
#include <vector>

namespace simall::io
{

class EnSightGoldWriter
{
public:
    /// Initialise the case directory. `caseDir` is created if it does not
    /// exist. `caseName` becomes the .case file base name.
    bool open(const std::string& caseDir, const std::string& caseName);

    /// Append a time step. Writes:
    ///   geometry on first call,
    ///   one variable file per scalar / vector field,
    ///   updates the .case file's TIME section.
    bool write_step(const meshing::Mesh& mesh,
                    const solver::FieldRegistry& fields,
                    double simulationTime);

    /// Flush and close the case file.
    bool close();

private:
    bool write_geometry(const meshing::Mesh& mesh);
    bool write_variable_scalar(const std::string& name, const solver::ScalarField& f, int step);
    bool write_variable_vector(const std::string& name, const solver::VectorField& f, int step);
    bool flush_case_file();

    std::string caseDir_;
    std::string caseName_;
    std::vector<double> times_;
    std::vector<std::string> scalarVars_;
    std::vector<std::string> vectorVars_;
    bool geometryWritten_ = false;
    bool isOpen_ = false;
};

} // namespace simall::io
