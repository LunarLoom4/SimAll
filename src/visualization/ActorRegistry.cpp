// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ActorRegistry.cpp
// =============================================================================
#include "visualization/ActorRegistry.hpp"

namespace simall::visualization
{

ActorRegistry::ActorRegistry() = default;

ActorId ActorRegistry::create(ActorKind kind, std::string name)
{
    std::lock_guard lock(mu_);
    ActorRecord r;
    r.id = nextId_++;
    r.kind = kind;
    r.name = std::move(name);
    actors_.emplace(r.id, std::move(r));
    return nextId_ - 1;
}

bool ActorRegistry::set_payload(ActorId id, ActorPayload payload)
{
    std::lock_guard lock(mu_);
    auto it = actors_.find(id);
    if (it == actors_.end())
        return false;

    // Kind must match the variant alternative.
    const auto expected = it->second.kind;
    const bool ok =
        (expected == ActorKind::Surface && std::holds_alternative<SurfaceMesh>(payload))
        || (expected == ActorKind::Lines && std::holds_alternative<LineSet>(payload))
        || (expected == ActorKind::Glyphs && std::holds_alternative<GlyphSet>(payload))
        || (expected == ActorKind::Volume && std::holds_alternative<VolumeMesh>(payload))
        || (expected == ActorKind::Image2D && std::holds_alternative<Image>(payload))
        || (expected == ActorKind::Overlay && std::holds_alternative<std::monostate>(payload));
    if (!ok)
        return false;

    it->second.bbox = compute_bbox(payload);
    it->second.payload = std::move(payload);
    return true;
}

void ActorRegistry::set_display(ActorId id, const ActorDisplay& d)
{
    std::lock_guard lock(mu_);
    auto it = actors_.find(id);
    if (it != actors_.end())
        it->second.display = d;
}

bool ActorRegistry::destroy(ActorId id)
{
    std::lock_guard lock(mu_);
    return actors_.erase(id) != 0;
}

std::optional<ActorRecord> ActorRegistry::get(ActorId id) const
{
    std::lock_guard lock(mu_);
    auto it = actors_.find(id);
    if (it == actors_.end())
        return std::nullopt;
    return it->second;
}

std::vector<ActorId> ActorRegistry::list() const
{
    std::lock_guard lock(mu_);
    std::vector<ActorId> out;
    out.reserve(actors_.size());
    for (auto& [id, _] : actors_)
        out.push_back(id);
    return out;
}

std::vector<ActorId> ActorRegistry::list(ActorKind kind) const
{
    std::lock_guard lock(mu_);
    std::vector<ActorId> out;
    for (auto& [id, rec] : actors_) {
        if (rec.kind == kind)
            out.push_back(id);
    }
    return out;
}

std::size_t ActorRegistry::size() const noexcept
{
    std::lock_guard lock(mu_);
    return actors_.size();
}

void ActorRegistry::set_visible_all(ActorKind kind, bool visible)
{
    std::lock_guard lock(mu_);
    for (auto& [_, rec] : actors_) {
        if (rec.kind == kind)
            rec.display.visible = visible;
    }
}

util::BoundingBox ActorRegistry::compute_bbox(const ActorPayload& p)
{
    util::BoundingBox box;
    std::visit(
        [&](auto const& payload) {
            using T = std::decay_t<decltype(payload)>;
            if constexpr (std::is_same_v<T, SurfaceMesh> || std::is_same_v<T, VolumeMesh>
                          || std::is_same_v<T, LineSet>) {
                for (auto const& p_ : payload.points)
                    box.expand(p_);
            } else if constexpr (std::is_same_v<T, GlyphSet>) {
                for (auto const& p_ : payload.anchor)
                    box.expand(p_);
            }
        },
        p);
    return box;
}

} // namespace simall::visualization
