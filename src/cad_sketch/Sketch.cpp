// =============================================================================
// SimAll Beta -- Parametric CAD / Sketcher
// File   : src/cad_sketch/Sketch.cpp
// Phase  : 23 Pass 23.1
// =============================================================================
#include "cad_sketch/Sketch.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

namespace simall::cad::sketch {

// --- parameters --------------------------------------------------------------

ParameterId Sketch::add_parameter(double initial_value, bool fixed, std::string name) {
    Parameter p;
    p.id    = next_pid_++;
    p.value = initial_value;
    p.fixed = fixed;
    p.name  = std::move(name);
    pid_index_.emplace(p.id, params_.size());
    params_.push_back(std::move(p));
    return params_.back().id;
}

void Sketch::set_parameter(ParameterId pid, double v) { param_ref(pid).value = v; }
double Sketch::get_parameter(ParameterId pid) const   { return param_ref(pid).value; }
void Sketch::fix_parameter(ParameterId pid, bool f)   { param_ref(pid).fixed = f; }

const Parameter& Sketch::parameter(ParameterId pid) const { return param_ref(pid); }

Parameter& Sketch::param_ref(ParameterId pid) {
    auto it = pid_index_.find(pid);
    if (it == pid_index_.end()) throw std::out_of_range("Sketch: unknown ParameterId");
    return params_[it->second];
}
const Parameter& Sketch::param_ref(ParameterId pid) const {
    auto it = pid_index_.find(pid);
    if (it == pid_index_.end()) throw std::out_of_range("Sketch: unknown ParameterId");
    return params_[it->second];
}

// --- entities ----------------------------------------------------------------

EntityId Sketch::push_entity(EntityKind k, std::vector<ParameterId> pids) {
    SketchEntity e;
    e.id     = next_eid_++;
    e.kind   = k;
    e.params = std::move(pids);
    eid_index_.emplace(e.id, entities_.size());
    entities_.push_back(std::move(e));
    return entities_.back().id;
}

EntityId Sketch::add_point(Vec2 p) {
    return push_entity(EntityKind::Point,
        { add_parameter(p.x), add_parameter(p.y) });
}

EntityId Sketch::add_line(Vec2 a, Vec2 b) {
    return push_entity(EntityKind::Line,
        { add_parameter(a.x), add_parameter(a.y),
          add_parameter(b.x), add_parameter(b.y) });
}

EntityId Sketch::add_circle(Vec2 c, double r) {
    return push_entity(EntityKind::Circle,
        { add_parameter(c.x), add_parameter(c.y), add_parameter(r) });
}

EntityId Sketch::add_arc(Vec2 c, double r, double t0, double t1) {
    return push_entity(EntityKind::Arc,
        { add_parameter(c.x), add_parameter(c.y), add_parameter(r),
          add_parameter(t0),  add_parameter(t1) });
}

EntityId Sketch::add_spline(const std::vector<Vec2>& cps) {
    std::vector<ParameterId> pids;
    pids.reserve(cps.size() * 2);
    for (const auto& p : cps) {
        pids.push_back(add_parameter(p.x));
        pids.push_back(add_parameter(p.y));
    }
    return push_entity(EntityKind::Spline, std::move(pids));
}

const SketchEntity& Sketch::entity(EntityId eid) const { return entity_ref(eid); }

const SketchEntity& Sketch::entity_ref(EntityId eid) const {
    auto it = eid_index_.find(eid);
    if (it == eid_index_.end()) throw std::out_of_range("Sketch: unknown EntityId");
    return entities_[it->second];
}

Vec2 Sketch::point_value(EntityId eid) const {
    const auto& e = entity_ref(eid);
    assert(e.kind == EntityKind::Point && e.params.size() == 2);
    return { get_parameter(e.params[0]), get_parameter(e.params[1]) };
}

Vec2 Sketch::line_endpoint(EntityId eid, int which) const {
    const auto& e = entity_ref(eid);
    assert(e.kind == EntityKind::Line && e.params.size() == 4);
    const std::size_t off = (which == 0) ? 0 : 2;
    return { get_parameter(e.params[off]), get_parameter(e.params[off + 1]) };
}

Vec2 Sketch::circle_center(EntityId eid) const {
    const auto& e = entity_ref(eid);
    assert((e.kind == EntityKind::Circle || e.kind == EntityKind::Arc) && e.params.size() >= 3);
    return { get_parameter(e.params[0]), get_parameter(e.params[1]) };
}

double Sketch::circle_radius(EntityId eid) const {
    const auto& e = entity_ref(eid);
    assert((e.kind == EntityKind::Circle || e.kind == EntityKind::Arc) && e.params.size() >= 3);
    return get_parameter(e.params[2]);
}

// --- constraints -------------------------------------------------------------

ConstraintId Sketch::push_constraint(ConstraintKind k, EntityId a, EntityId b,
                                     std::optional<double> v) {
    SketchConstraint c;
    c.id    = next_cid_++;
    c.kind  = k;
    c.a     = a;
    c.b     = b;
    c.value = v;
    cid_index_.emplace(c.id, constraints_.size());
    constraints_.push_back(c);
    return c.id;
}

ConstraintId Sketch::add_coincident   (EntityId p1, EntityId p2)              { return push_constraint(ConstraintKind::Coincident,    p1, p2, std::nullopt); }
ConstraintId Sketch::add_horizontal   (EntityId l)                            { return push_constraint(ConstraintKind::Horizontal,    l,  kInvalidId, std::nullopt); }
ConstraintId Sketch::add_vertical     (EntityId l)                            { return push_constraint(ConstraintKind::Vertical,      l,  kInvalidId, std::nullopt); }
ConstraintId Sketch::add_parallel     (EntityId l1, EntityId l2)              { return push_constraint(ConstraintKind::Parallel,      l1, l2, std::nullopt); }
ConstraintId Sketch::add_perpendicular(EntityId l1, EntityId l2)              { return push_constraint(ConstraintKind::Perpendicular, l1, l2, std::nullopt); }
ConstraintId Sketch::add_distance     (EntityId p1, EntityId p2, double d)    { return push_constraint(ConstraintKind::Distance,      p1, p2, d); }
ConstraintId Sketch::add_angle        (EntityId l1, EntityId l2, double rad)  { return push_constraint(ConstraintKind::Angle,         l1, l2, rad); }
ConstraintId Sketch::add_radius       (EntityId c,  double r)                 { return push_constraint(ConstraintKind::Radius,        c,  kInvalidId, r); }
ConstraintId Sketch::add_point_on_line  (EntityId pt, EntityId line)          { return push_constraint(ConstraintKind::PointOnLine,   pt, line, std::nullopt); }
ConstraintId Sketch::add_point_on_circle(EntityId pt, EntityId circle)        { return push_constraint(ConstraintKind::PointOnCircle, pt, circle, std::nullopt); }
ConstraintId Sketch::add_equal_length (EntityId l1, EntityId l2)              { return push_constraint(ConstraintKind::EqualLength,   l1, l2, std::nullopt); }
ConstraintId Sketch::add_fix          (EntityId pt) {
    // Fix is encoded by pinning the underlying parameters; the constraint
    // record itself exists purely for diagnostic/undo bookkeeping and emits
    // no residual at solve time (see emit_constraint_residuals).
    const auto& ent = entity_ref(pt);
    for (ParameterId pid : ent.params) param_ref(pid).fixed = true;
    return push_constraint(ConstraintKind::Fix, pt, kInvalidId, std::nullopt);
}

const SketchConstraint& Sketch::constraint(ConstraintId cid) const {
    auto it = cid_index_.find(cid);
    if (it == cid_index_.end()) throw std::out_of_range("Sketch: unknown ConstraintId");
    return constraints_[it->second];
}

// --- mutation ---------------------------------------------------------------

bool Sketch::remove_entity(EntityId eid) {
    auto it = eid_index_.find(eid);
    if (it == eid_index_.end()) return false;

    // Drop dependent constraints first.
    constraints_.erase(std::remove_if(constraints_.begin(), constraints_.end(),
        [eid](const SketchConstraint& c) { return c.a == eid || c.b == eid; }),
        constraints_.end());

    entities_.erase(entities_.begin() + static_cast<std::ptrdiff_t>(it->second));

    // Rebuild indices (cheap: sketches are small).
    eid_index_.clear();
    for (std::size_t i = 0; i < entities_.size(); ++i) eid_index_.emplace(entities_[i].id, i);
    cid_index_.clear();
    for (std::size_t i = 0; i < constraints_.size(); ++i) cid_index_.emplace(constraints_[i].id, i);
    return true;
}

bool Sketch::remove_constraint(ConstraintId cid) {
    auto it = cid_index_.find(cid);
    if (it == cid_index_.end()) return false;
    constraints_.erase(constraints_.begin() + static_cast<std::ptrdiff_t>(it->second));
    cid_index_.clear();
    for (std::size_t i = 0; i < constraints_.size(); ++i) cid_index_.emplace(constraints_[i].id, i);
    return true;
}

void Sketch::clear() {
    params_.clear();
    entities_.clear();
    constraints_.clear();
    pid_index_.clear();
    eid_index_.clear();
    cid_index_.clear();
    next_pid_ = next_eid_ = next_cid_ = 1;
}

}  // namespace simall::cad::sketch
