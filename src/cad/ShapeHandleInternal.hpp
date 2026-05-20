// =============================================================================
// SimAll Beta - CAD Subsystem (INTERNAL)
// File   : src/cad/ShapeHandleInternal.hpp
// Phase  : 4 (CAD KERNEL — internal sharing header)
//
// This header is the ONLY place in the project that exposes the layout of
// `ShapeHandle::Impl`. It must only be included from translation units in
// src/cad/ — never from public headers and never from other subsystems.
// All OpenCASCADE-using CAD modules (readers, healing, booleans, …) include
// this header so they can manipulate the wrapped TopoDS_Shape without
// duplicating PIMPL boilerplate.
//
// Subsystem boundary rule (revised from the original CadKernel.cpp comment):
//   - The CAD subsystem (src/cad/*.cpp) MAY include OpenCASCADE headers.
//   - Every other subsystem accesses geometry strictly via the public
//     facades declared in cad/CadKernel.hpp / cad/Tessellator.hpp etc.
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"
#include "cad/TopologyGraph.hpp"

// --- OpenCASCADE core types ------------------------------------------------
#include <TopAbs.hxx>
#include <TopoDS_Shape.hxx>

namespace simall::cad
{

struct ShapeHandle::Impl
{
    TopoDS_Shape shape;
    TopologyGraph graph;
};

/// Internal-only accessor — friended by `ShapeHandle`. Lets every CAD
/// translation unit reach the underlying TopoDS_Shape without exposing
/// OpenCASCADE types in any public header.
struct ShapeHandleAccess
{
    static ShapeHandle::Impl& impl(ShapeHandle& h) noexcept { return *h.impl_; }
    static const ShapeHandle::Impl& impl(const ShapeHandle& h) noexcept { return *h.impl_; }
    static TopoDS_Shape& shape(ShapeHandle& h) noexcept { return h.impl_->shape; }
    static const TopoDS_Shape& shape(const ShapeHandle& h) noexcept { return h.impl_->shape; }
};

/// Translate an OCC topology kind into the SimAll enum used by TopologyGraph.
TopologyType occKindToTopologyType(TopAbs_ShapeEnum e) noexcept;

/// Walk every (vertex, edge, wire, face, shell, solid) of `shape` and push
/// it into `g` as a node, replacing any previous content.
void rebuildTopologyGraph(const TopoDS_Shape& shape, TopologyGraph& g);

/// Construct a fresh ShapeHandle that owns `shape`; its topology graph is
/// populated immediately. Use this in every reader / boolean / healing op
/// to avoid re-implementing the PIMPL plumbing.
ShapeHandle makeHandle(TopoDS_Shape shape);

} // namespace simall::cad
