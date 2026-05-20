// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ActorRegistry.hpp
// Phase  : 3 / Week 14 (scene-management trio: ActorRegistry)
//
// Owns the *kind* and *display state* of every actor that may end up on the
// scene.  Geometry payload is stored polymorphically (variant-of-views) so a
// downstream renderer can dispatch with a single switch.  Picking and
// recorders read directly from this registry; the SceneGraph only stores
// ActorIds.
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace simall::visualization
{

enum class ActorKind : std::uint8_t
{
    Unknown = 0,
    Surface, // SurfaceMesh payload
    Lines,   // LineSet payload (streamlines, contours, edges)
    Glyphs,  // GlyphSet payload (vector glyphs)
    Volume,  // VolumeMesh payload (raycast or interior wireframe)
    Image2D, // Image payload (LIC, scalar bar, screenshot preview)
    Overlay  // 2D HUD primitive owned by the back-end
};

using ActorPayload =
    std::variant<std::monostate, SurfaceMesh, LineSet, GlyphSet, VolumeMesh, Image>;

struct ActorDisplay
{
    bool visible = true;
    bool picking = true;
    bool showEdges = false;
    double opacity = 1.0;
    double lineWidth = 1.0;
    double pointSize = 1.0;
    Color4 solidColor{0.85f, 0.85f, 0.90f, 1.0f};
    bool useScalarColor = false; // true → consult TF
    TransferFunction tf{};
    double scalarMin = 0.0;
    double scalarMax = 1.0;
};

struct ActorRecord
{
    ActorId id = kInvalidActorId;
    std::string name;
    ActorKind kind = ActorKind::Unknown;
    ActorPayload payload;
    ActorDisplay display;
    util::BoundingBox bbox; // updated when payload is replaced
};

class ActorRegistry
{
public:
    ActorRegistry();

    /// Create an empty actor of the given kind; caller follows up with
    /// set_payload() and set_display().  Names need NOT be unique but
    /// must be non-empty.
    ActorId create(ActorKind kind, std::string name);

    /// Replace the geometry/image payload; bounding box is recomputed.  The
    /// payload kind must match the actor's declared kind, otherwise the
    /// call is rejected (returns false).
    bool set_payload(ActorId id, ActorPayload payload);

    void set_display(ActorId id, const ActorDisplay& d);

    /// Permanently remove the actor.  Detach from any SceneGraph node
    /// separately — registry has no back-pointer.
    bool destroy(ActorId id);

    /// Read-only access; copies are returned to keep callers thread-safe.
    std::optional<ActorRecord> get(ActorId id) const;
    std::vector<ActorId> list() const;
    std::vector<ActorId> list(ActorKind kind) const;
    std::size_t size() const noexcept;

    /// Bulk visibility / opacity edits used by the layer panel in the GUI.
    void set_visible_all(ActorKind kind, bool visible);

private:
    static util::BoundingBox compute_bbox(const ActorPayload& p);

    mutable std::mutex mu_;
    std::unordered_map<ActorId, ActorRecord> actors_;
    ActorId nextId_ = 1;
};

} // namespace simall::visualization
