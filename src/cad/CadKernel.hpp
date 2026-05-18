// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/CadKernel.hpp
// Phase  : 4 (CAD KERNEL FAÇADE)
//
// Public surface of the OpenCASCADE wrapper. The rest of the codebase MUST
// only include this header — never <STEPControl_Reader.hxx> directly.
// =============================================================================
#pragma once

#include "TopologyGraph.hpp"
#include "utilities/MathTypes.hpp"

#include <memory>
#include <string>
#include <vector>

namespace simall::cad {

struct TessellationParams {
    double deflection = 1.0e-3;   // chordal tolerance (m)
    double angle      = 0.35;     // rad (~20 deg)
    bool   relative   = true;
    int    minTrianglesPerFace = 4;
};

struct HealingOptions {
    double tolerance        = 1.0e-6;
    bool   stitchGaps       = true;
    bool   removeSlivers    = true;
    bool   collapseTinyEdges = true;
    bool   harmonizeTolerances = true;
};

struct TriangleMesh {
    std::vector<util::Vec3d> points;
    std::vector<std::array<std::uint32_t, 3>> triangles;
    std::vector<util::PersistentId>            triangleFaceId;   // ←  picking key
};

class ShapeHandle {
public:
    ShapeHandle();
    ~ShapeHandle();
    ShapeHandle(ShapeHandle&&) noexcept;
    ShapeHandle& operator=(ShapeHandle&&) noexcept;
    ShapeHandle(const ShapeHandle&)            = delete;
    ShapeHandle& operator=(const ShapeHandle&) = delete;

    bool valid() const noexcept;
    const TopologyGraph& topology() const noexcept;
    TopologyGraph&       topology()       noexcept;

    util::BoundingBox bounds() const;

private:
    friend class  CadKernel;
    friend struct ShapeHandleAccess;       // implementation-side helper
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class CadKernel {
public:
    CadKernel();
    ~CadKernel();

    /// STEP/IGES/STL importer dispatch by extension.
    ShapeHandle import(const std::string& path);

    /// Repair operations (Phase 4.3).
    void heal(ShapeHandle& shape, const HealingOptions& opts = {});

    /// Convert exact geometry → triangles for visualization (Phase 4.4).
    /// Triangle ↔ face mapping is preserved in TriangleMesh::triangleFaceId.
    TriangleMesh tessellate(const ShapeHandle& shape, const TessellationParams& p = {});

    /// Curvature-adaptive tessellation: h = min(hmax, sqrt(2εR)). Section 5.6.
    TriangleMesh adaptive_tessellate(const ShapeHandle& shape,
                                     double maxEdgeLength,
                                     double chordalEpsilon);
};

}  // namespace simall::cad
