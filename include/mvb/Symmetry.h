#pragma once

#include "mvb/NamedShape.h"

#include <array>
#include <vector>

namespace mvb {

// A reflective coordinate-plane symmetry of an assembled shape.
enum class SymmetryPlane : int {
    None = 0,
    XY = 1,   // z = 0, normal = ±Z
    YZ = 2,   // x = 0, normal = ±X
    XZ = 3,   // y = 0, normal = ±Y
};

// Which side of the plane to keep when cutting.
enum class SymmetryHalf : int {
    Positive = 0,   // normal · r ≥ 0
    Negative = 1,
};

// Result of analyzing one candidate plane.
struct SymmetryCandidate {
    SymmetryPlane plane = SymmetryPlane::None;
    double volume_plus  = 0.0;    // ∫ over positive half [m³]
    double volume_minus = 0.0;    // ∫ over negative half [m³]
    double volume_asymmetry = 1.0; // |V+ - V-| / V_total
    bool   valid = false;          // asymmetry < tolerance
};

struct SymmetryResult {
    std::vector<SymmetryCandidate> candidates;
    std::vector<SymmetryPlane>     valid_planes;
};

// Aggregate shape bounding box — written to (xmin..zmax). Convenience for
// callers of cut_half_space.
struct ShapeBBox {
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
};

ShapeBBox aggregate_bbox(const std::vector<NamedShape>& shapes);

// Geometric symmetry test — for each of {XY, YZ, XZ}, intersects every
// shape with a half-space box on both sides and compares volumes. A plane
// is valid when |V_+ - V_-| / V_total < volume_tolerance.
//
// Pure geometry: has no knowledge of FEM concepts (winding axis, ports,
// etc.). Consumers apply domain-specific filters on the returned planes.
SymmetryResult analyze_symmetry(const std::vector<NamedShape>& shapes,
                                double volume_tolerance = 1e-3);

// Cut a single named shape to one half-space of the given plane. Returns a
// vector because a boolean Common against a half-space can yield multiple
// disconnected solids (typical for a core that wraps around the plane);
// each sub-solid carries the original name, with a "_0", "_1", ... suffix
// when the cut produces more than one.
//
// The bbox passed in sizes the half-space box; it must cover all shapes
// that will be cut (call `aggregate_bbox` on the full set first). The
// cut plane itself is always at the world origin (x=0, y=0, or z=0).
std::vector<NamedShape> cut_half_space(const NamedShape& in,
                                       SymmetryPlane plane,
                                       SymmetryHalf  half,
                                       const ShapeBBox& bbox);

// Overload: cut every shape in a vector, flattening the per-shape results.
std::vector<NamedShape> cut_half_space(const std::vector<NamedShape>& in,
                                       SymmetryPlane plane,
                                       SymmetryHalf  half,
                                       const ShapeBBox& bbox);

// Compose multiple half-space cuts in sequence (e.g. YZ then XY → quarter
// domain, or YZ then XY then XZ → octant). Shapes surviving each cut with
// non-empty volume are forwarded to the next cut; names accumulate
// suffixes only when genuine splits occur.
std::vector<NamedShape> cut_to_region(
    const std::vector<NamedShape>& in,
    const std::vector<std::pair<SymmetryPlane, SymmetryHalf>>& cuts,
    const ShapeBBox& bbox);

const char* to_string(SymmetryPlane p);

// Convenience: index of the normal and the two in-plane axes (0=X, 1=Y, 2=Z).
// For SymmetryPlane::YZ the normal is X (0) and in-plane axes are Y, Z.
void plane_axes(SymmetryPlane p, int& i_normal, int& i_in_a, int& i_in_b);

} // namespace mvb
