// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/Booleans.hpp
// Phase  : 4.2 (Boolean operations)
//
// Boolean set operations on solids and shells. SimAll exposes the three
// canonical ops (fuse / cut / common) plus a multi-input convenience for
// large assemblies. Underlying implementation is OCC's BOPAlgo, requested
// in parallel mode so it scales with ThreadPool worker count.
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"

#include <vector>

namespace simall::cad {

struct BooleanOptions {
    double fuzzyValue       = 0.0;   // 0 → OCC default tolerance
    bool   runParallel      = true;
    bool   gluePartial      = false; // BOPAlgo glue mode for non-overlapping touching faces
    bool   checkInverted    = true;
    bool   nonDestructive   = true;  // copy inputs so source shapes are preserved
};

class Booleans {
public:
    ShapeHandle fuse  (const ShapeHandle& a, const ShapeHandle& b, const BooleanOptions& o = {});
    ShapeHandle cut   (const ShapeHandle& a, const ShapeHandle& b, const BooleanOptions& o = {});
    ShapeHandle common(const ShapeHandle& a, const ShapeHandle& b, const BooleanOptions& o = {});

    /// Fuse many shapes at once (more efficient than pairwise loop).
    ShapeHandle fuseMany(const std::vector<const ShapeHandle*>& inputs,
                         const BooleanOptions& o = {});
};

}  // namespace simall::cad
