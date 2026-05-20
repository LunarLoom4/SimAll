// =============================================================================
// SimAll Beta -- Parametric CAD / Sketcher
// File   : src/cad_sketch/Sketch.hpp
// Phase  : 23 Pass 23.1
//
// Core data model for a parametric 2D sketch:
//   * Parameter -- a single free scalar degree of freedom.
//   * SketchEntity -- line / arc / circle / spline / point, expressed as
//     a small bundle of ParameterId references into the sketch's parameter
//     pool (so every geometric coordinate the solver can move is a Parameter).
//   * SketchConstraint -- coincident / horizontal / vertical / parallel /
//     perpendicular / distance / angle / radius / point-on-line, expressed
//     as a small bundle of entity / parameter references plus an optional
//     scalar target value.
//   * Sketch -- the container that owns parameters, entities, constraints,
//     and assigns persistent monotonically-increasing ids.
//
// The model deliberately stores everything by integer id (not pointer) so
// the entire sketch is trivially copyable, serialisable, and snapshot-safe
// for the undo stack.  The constraint solver lives in ConstraintSolver.hpp
// and operates purely on these structs.
// =============================================================================
#pragma once

#include "cad_sketch/Vec2.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::cad::sketch
{

// -----------------------------------------------------------------------------
// Strongly-typed ids.  Plain uint32_t with named aliases keeps the data POD
// while still giving readable function signatures.
// -----------------------------------------------------------------------------
using ParameterId = std::uint32_t;
using EntityId = std::uint32_t;
using ConstraintId = std::uint32_t;
inline constexpr std::uint32_t kInvalidId = 0;

// -----------------------------------------------------------------------------
// Parameter -- one scalar unknown.  `fixed = true` means the solver treats
// the parameter as constant (it will not appear in the unknown vector x).
// -----------------------------------------------------------------------------
struct Parameter
{
    ParameterId id{kInvalidId};
    double value{0.0};
    bool fixed{false};
    std::string name{}; // optional human-readable label
};

// -----------------------------------------------------------------------------
// Entity model.  Each kind reserves a fixed number of ParameterId slots in
// `params`.  Splines reserve 2*N (one (x,y) pair per control point).
// -----------------------------------------------------------------------------
enum class EntityKind : std::uint8_t
{
    Point = 0,  // params = {px, py}
    Line = 1,   // params = {ax, ay, bx, by}
    Circle = 2, // params = {cx, cy, r}
    Arc = 3,    // params = {cx, cy, r, theta_start, theta_end}
    Spline = 4  // params = {x0, y0, x1, y1, ..., xN, yN}  (N+1 control pts)
};

struct SketchEntity
{
    EntityId id{kInvalidId};
    EntityKind kind{EntityKind::Point};
    std::vector<ParameterId> params{}; // small-vector in practice
    bool construction{false};          // construction geom (not exported)
};

// -----------------------------------------------------------------------------
// Constraint model.  Different kinds reference different entity/parameter
// slots; rather than a discriminated union we use a flat struct with optional
// fields and let the solver dispatch on `kind`.  This keeps memory layout
// trivial and serialisation a one-liner.
// -----------------------------------------------------------------------------
enum class ConstraintKind : std::uint8_t
{
    Coincident = 0,    // two Points coincide                          -> 2 residuals
    Horizontal = 1,    // Line is horizontal (ay == by)                -> 1 residual
    Vertical = 2,      // Line is vertical   (ax == bx)                -> 1 residual
    Parallel = 3,      // two Lines parallel                           -> 1 residual
    Perpendicular = 4, // two Lines perpendicular                      -> 1 residual
    Distance = 5,      // distance between two Points equals `value`   -> 1 residual
    Angle = 6,         // angle between two Lines equals `value` (rad) -> 1 residual
    Radius = 7,        // Circle or Arc radius equals `value`          -> 1 residual
    PointOnLine = 8,   // Point lies on Line                           -> 1 residual
    PointOnCircle = 9, // Point lies on Circle                         -> 1 residual
    EqualLength = 10,  // two Lines have equal length                  -> 1 residual
    Fix = 11           // Pin one Point in place at its current loc    -> 2 residuals
};

struct SketchConstraint
{
    ConstraintId id{kInvalidId};
    ConstraintKind kind{ConstraintKind::Coincident};

    // Up to two entity references; meaning depends on `kind`.
    EntityId a{kInvalidId};
    EntityId b{kInvalidId};

    // Optional scalar target (distance / angle / radius / equal-length-pair).
    std::optional<double> value{};
};

// -----------------------------------------------------------------------------
// Sketch container.
// -----------------------------------------------------------------------------
class Sketch
{
public:
    // ---- parameter management -------------------------------------------------
    ParameterId add_parameter(double initial_value, bool fixed = false, std::string name = {});
    void set_parameter(ParameterId pid, double v);
    double get_parameter(ParameterId pid) const;
    void fix_parameter(ParameterId pid, bool fixed = true);
    [[nodiscard]] const Parameter& parameter(ParameterId pid) const;
    [[nodiscard]] std::size_t parameter_count() const noexcept { return params_.size(); }
    [[nodiscard]] const std::vector<Parameter>& parameters() const noexcept { return params_; }

    // ---- entity factories -----------------------------------------------------
    // All factories internally allocate the parameter pool entries.
    EntityId add_point(Vec2 p);
    EntityId add_line(Vec2 a, Vec2 b);
    EntityId add_circle(Vec2 center, double radius);
    EntityId add_arc(Vec2 center, double radius, double theta_start, double theta_end);
    EntityId add_spline(const std::vector<Vec2>& control_points);

    // ---- entity inspection ----------------------------------------------------
    [[nodiscard]] const SketchEntity& entity(EntityId eid) const;
    [[nodiscard]] const std::vector<SketchEntity>& entities() const noexcept { return entities_; }

    // Convenience: extract endpoint / centre / radius of an entity by reading
    // the underlying parameters.
    [[nodiscard]] Vec2 point_value(EntityId eid) const;
    [[nodiscard]] Vec2 line_endpoint(EntityId eid, int which) const; // which = 0 or 1
    [[nodiscard]] Vec2 circle_center(EntityId eid) const;
    [[nodiscard]] double circle_radius(EntityId eid) const;

    // ---- constraint factories -------------------------------------------------
    ConstraintId add_coincident(EntityId p1, EntityId p2);
    ConstraintId add_horizontal(EntityId line);
    ConstraintId add_vertical(EntityId line);
    ConstraintId add_parallel(EntityId l1, EntityId l2);
    ConstraintId add_perpendicular(EntityId l1, EntityId l2);
    ConstraintId add_distance(EntityId p1, EntityId p2, double distance);
    ConstraintId add_angle(EntityId l1, EntityId l2, double radians);
    ConstraintId add_radius(EntityId circle, double radius);
    ConstraintId add_point_on_line(EntityId pt, EntityId line);
    ConstraintId add_point_on_circle(EntityId pt, EntityId circle);
    ConstraintId add_equal_length(EntityId l1, EntityId l2);
    ConstraintId add_fix(EntityId pt);

    // ---- constraint inspection -----------------------------------------------
    [[nodiscard]] const SketchConstraint& constraint(ConstraintId cid) const;
    [[nodiscard]] const std::vector<SketchConstraint>& constraints() const noexcept
    {
        return constraints_;
    }

    // ---- mutation -------------------------------------------------------------
    bool remove_entity(EntityId eid); // also drops dependent constraints
    bool remove_constraint(ConstraintId cid);

    void clear();

private:
    ParameterId next_pid_{1};
    EntityId next_eid_{1};
    ConstraintId next_cid_{1};
    std::vector<Parameter> params_{};
    std::vector<SketchEntity> entities_{};
    std::vector<SketchConstraint> constraints_{};
    std::unordered_map<ParameterId, std::size_t> pid_index_{};
    std::unordered_map<EntityId, std::size_t> eid_index_{};
    std::unordered_map<ConstraintId, std::size_t> cid_index_{};

    // Internal helpers.
    [[nodiscard]] Parameter& param_ref(ParameterId pid);
    [[nodiscard]] const Parameter& param_ref(ParameterId pid) const;
    [[nodiscard]] const SketchEntity& entity_ref(EntityId eid) const;
    EntityId push_entity(EntityKind k, std::vector<ParameterId> pids);
    ConstraintId push_constraint(ConstraintKind k, EntityId a, EntityId b, std::optional<double> v);
};

} // namespace simall::cad::sketch
