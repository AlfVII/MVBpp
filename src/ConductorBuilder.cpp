#include "mvb/ConductorBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"
#include <array>
#include <functional>
#include <set>
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <Geom_BezierCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <cstdio>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepOffsetAPI_MakePipeShell.hxx>
#include <BRepLib.hxx>
#include <GCE2d_MakeSegment.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <gp_Vec.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Circ.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Dir.hxx>
#include <gp_XY.hxx>
#include <gp_XYZ.hxx>
#include <BRepAlgoAPI_Check.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BOPAlgo_GlueEnum.hxx>
#include <TopTools_ListOfShape.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Solid.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Wire.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <iostream>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace mvb {

namespace {

// OPT-IN tight sweep tolerance (metres) for the single-body MakePipeShell sweep, via
// MVB_SWEEP_TOL (0 / unset = OCCT's default 1e-4, the long-validated behaviour).
//
// WHY IT EXISTS: MakePipeShell's default 3D approximation tolerance is 1e-4 m -- 100 um -- while a
// densely wound coil's inter-turn copper clearance is ~70 um (PQ5050: 38 turns of r = 0.315 mm at
// outer-diameter pitch; ETD34 likewise). At the default tolerance adjacent sheets of the ONE swept
// lateral surface interpenetrate IN THE CAD; BRepCheck does not catch it, and gmsh then fails with
// "PLC Error: a segment and a facet intersect" at coarse mesh targets and "overlapping facets on
// surface N surface N" (the surface against ITSELF) at fine ones -- measured on both designs down
// to a 0.115 mm conductor target where chord error (~5 um) is far below the gap. MVB_SWEEP_TOL=1e-6
// eliminates that self-intersection (ETD34's finest-rung failure signature disappeared and it
// meshed for 95 minutes of real work instead of dying in 40 s).
//
// WHY IT IS NOT THE DEFAULT: the tight fit is not uniformly safe -- on the EP oblong stadium spine
// it produced a body that passed every internal check (volume, BRepCheck) while a section had
// collapsed off a crossing. Making it default needs a discriminator that accepts the tight body
// only when it is genuinely sound (crossing-containment via the solid classifier was measured too
// entangled with the fillet construction to bolt on quickly); that integration is registered work,
// not a flag flip. Until then: set MVB_SWEEP_TOL=1e-6 for dense round-column windings.
static const double kSweepTol3d =
    std::getenv("MVB_SWEEP_TOL") ? std::atof(std::getenv("MVB_SWEEP_TOL")) : 0.0;

constexpr double kPi = std::numbers::pi;
constexpr double kTwoPi = 2.0 * std::numbers::pi;
// Curved primitives are collision-checked as sampled polylines; the sampling step is
// chosen so each chord sags inward by at most this fraction of the wire radius. This is
// a DISCRETIZATION-DENSITY parameter of the measurement (like polygonSegments), not a
// clearance policy.
constexpr double kMaxSagFraction = 0.02;
// The measurement's own error bound: a sampled polyline lies within this sag of the true
// curve, so a measured polyline-polyline distance can under-read the true centreline
// distance by at most the two prims' sag bounds. Derived from the sampling rule actually
// applied in samplePrim/curveSampleCount -- never a policy constant.
inline double samplingSag(double wireRadius) {
    return std::max(1e-6, kMaxSagFraction * wireRadius);
}
// THE one collision criterion (Alf, 2026-08-07: no allowances, no invented tolerances --
// MKF/MAS geometry is the source of truth, and when it is wrong we THROW, we never absorb
// it): two wires' bare copper may TOUCH but never interpenetrate, so the minimum
// centreline separation is EXACTLY the sum of the BARE (conducting) radii. Every layout
// decision (vertical-fan slots, drift windows, lead corridors) demands precisely this.
// At the CHECK, the measured distance is credited with the sampling sag bounds above, so
// the gate throws exactly when copper overlap is certain beyond measurement resolution.
constexpr double gateMinSeparation(double bareA, double bareB) {
    return bareA + bareB;
}
// The connection plane ("connections radial, 90 degrees from the turns, in the YZ
// plane"). MKF's 2D winding-window cross-section maps into 3D on the -Z side of the
// column: 2D (x = radial, y = axial) -> 3D (0, y, -(x + zoff)), where zoff = 0 for round
// columns and columnDepth - columnWidth for rectangular/oblong ones (the same clearance
// past the column lands the wire deeper on the Z faces than on the X faces). Every
// crossing sits on this plane; every drawn ConnectionReservedSpace rectangle (the
// pink/blue boxes of the Painter SVG) is replayed here verbatim — nothing is invented.
constexpr double kPlaneAz = kPi / 2.0;

// Azimuth convention (OCCT right-handed rotation about +Y, the column axis), generalised
// to a vertical axis through the horizontal point (cx, cz):
//   pos(r, y, az) = (cx + r cos az, y, cz - r sin az)
gp_Pnt azPointC(double cx, double cz, double r, double y, double az) {
    return gp_Pnt(cx + r * std::cos(az), y, cz - r * std::sin(az));
}

// Rodrigues rotation of v about the unit axis by ang.
gp_XYZ rotateXYZ(const gp_XYZ& v, const gp_XYZ& axisUnit, double ang) {
    return v * std::cos(ang) + axisUnit.Crossed(v) * std::sin(ang) +
           axisUnit * (axisUnit.Dot(v) * (1.0 - std::cos(ang)));
}

// ---------------------------------------------------------------------------------------
// Path model.
//   SEG    — straight segment (leads, racetrack faces, toroidal tubes and chords).
//   ARC3   — exact circular arc about an arbitrary axis: P(t) = c + rot(axis, t)·v0 for
//            t in [0, sweep] (racetrack corners, oblong caps, toroidal elbows).
//   SPIRAL — spiral about a VERTICAL axis through (cx, cz): azimuth linear in the
//            parameter, radius/height linear (round-column wraps) or cosine-blended
//            (blend=true, oblong cap ramps) so the end tangents are purely azimuthal.
//   BLEND  — S-curve between two points whose tangent is parallel to a common direction
//            at BOTH ends (the rectangular-column transition ramp): the lateral offset
//            follows a cosine blend, so the chain stays tangent-continuous through the
//            crossings while passing EXACTLY through them. (A straight ramp leaves
//            ~15-degree kinks at its junctions; OCCT's pipe-shell corner transitions
//            flare unboundedly on such kinks — observed: 9 mm spikes.)
struct Seg {
    gp_Pnt a, b;
};
struct Arc3 {
    gp_Pnt c;            // arc centre
    gp_XYZ axis{0, 1, 0};  // unit rotation axis (right-handed sweep)
    gp_XYZ v0{0, 0, 0};    // centre -> start vector, |v0| = bend radius
    double sweep = 0;
};
struct Spiral {
    double cx = 0, cz = 0;   // vertical axis position
    double r0 = 0, y0 = 0, az0 = 0;
    double r1 = 0, y1 = 0, az1 = 0;
    bool blend = false;      // cosine-blend r/y (end tangents purely azimuthal)
};
struct Blend {
    gp_Pnt a, b;
    gp_XYZ u{1, 0, 0};       // unit tangent direction at both ends
};
struct Primitive {
    enum Kind { SEG, ARC3, SPIRAL, BLEND } kind = SEG;
    Seg seg{};
    Arc3 arc{};
    Spiral spiral{};
    Blend blendc{};
    std::string label;
    // Electrical turn ordinal this primitive belongs to (entrance lead = first turn's,
    // exit lead = last turn's). A continuous wire legitimately contacts itself only
    // between CONSECUTIVE turns; the gate exempts same-conductor pairs with
    // |ordinal diff| <= 1 and checks everything farther apart.
    size_t turnOrdinal = 0;
    // Lead primitives are the replayed MKF connection routes. Emission splits runs at
    // lead<->wrap boundaries: those junctions are the only 90-degree corners of the path
    // (wrap chains are tangent-continuous), and OCCT's pipe-shell mitre flares unboundedly
    // across them, while lead chains sweep robustly as exact cylinders + sphere elbows.
    bool isLead = false;
    // An inter-layer link (MKF's blue ConnectionReservedSpace box): the straight radial
    // step that carries the wire from one concentric layer out to the next. Like a lead it
    // meets its neighbouring wraps at 90-degree corners, so it is isolated into its own run
    // and swept piecewise (cylinder + sphere elbows) rather than through the pipe-shell.
    bool isConnection = false;
};
struct ConductorPath {
    std::string name;
    double wireRadius = 0.0;
    // Bare-COPPER (conducting) envelope, independent of paintCoating. The collision gate uses THIS, not
    // the (possibly insulated) wireRadius: two turns whose enamel touches on a tight winding are
    // physically valid and mesh fine at the conducting radius, but overlapping COPPER is a real short /
    // an un-meshable solid and is always rejected. For a rect wire these are the copper half-extents.
    double condRadius = 0.0, condWidth = 0.0, condHeight = 0.0;
    // Seam rotation this path was rigidly rotated by about the column axis (windings/parallels are
    // staggered so leads don't collide) — deferred Z end-run planning attaches at this azimuth.
    double seamRot = 0.0;
    std::vector<Primitive> prims;
    // Whether this winding may be emitted as ONE swept solid (see emitConductor). True only for
    // ROUND and OBLONG concentric columns, whose centerlines MakePipeShell sweeps with the section
    // centered on the spine, honoring the MKF crossings. Excluded:
    //   - RECTANGULAR columns: the long straight racetrack sides flip OCCT's pipe frame, producing
    //     a valid-but-displaced pipe sitting a radius or two off the crossings.
    //   - TOROIDAL windings: the high-torsion hole-threading spine only closes under a torsion-
    //     stable frame, which shifts the section ~1-2% off the spine at the lead-adjacent crossings.
    // Both keep the EXACT per-run compound, which places every crossing exactly.
    bool singleBodyCapable = false;
    // Rectangular / planar wire: the swept profile is an oriented RECTANGLE (width = radial,
    // height = axial, matching TurnBuilder::build_rect_profile). On a ROUND column it sweeps as one
    // body with a fixed binormal (singleBodyCapable). On a RECT/OBLONG column the wide flat section
    // cannot be swept around the racetrack corners (its inner edge would collapse), so instead each
    // primitive is built as its own rect solid -- straights as prisms, corners as REVOLVED rect
    // (annular wedges, clean at any corner radius) -- and the solids fused (useRectSolids).
    bool isRectangular = false;
    bool useRectSolids = false;
    // ROUND wire on a toroid / rect column also goes through the per-primitive analytic path (SEG ->
    // cylinder, ARC3 -> torus segment), because unioning the swept BSpline pipes of the per-run
    // compound leaves seam slivers. Analytic solids fuse into ONE clean FEM-ready body. roundProfile
    // selects the circle section in that path; round/oblong columns keep the clean single-body sweep.
    bool roundProfile = false;
    bool toroidal = false;    // rect wire on a toroid: the section's "axial" axis is the local
                              // AZIMUTHAL direction (tangent to the big ring), not the column Y.
    double wireWidth = 0.0;   // radial extent [m]
    double wireHeight = 0.0;  // axial extent [m]
    bool femReady = false;    // true -> pay for the one-piece/conformal FEM geometry; false -> fast
                              // per-run compound for drawing (see ConductorBuilder::Options::femReady)
};

// --- capsule distance helpers ----------------------------------------------------------
double segSegDistance(const gp_Pnt& p1, const gp_Pnt& q1, const gp_Pnt& p2, const gp_Pnt& q2) {
    // Segment-segment minimum distance (Ericson, Real-Time Collision Detection, 5.1.9).
    gp_XYZ d1 = q1.XYZ() - p1.XYZ();
    gp_XYZ d2 = q2.XYZ() - p2.XYZ();
    gp_XYZ r = p1.XYZ() - p2.XYZ();
    double a = d1.SquareModulus(), e = d2.SquareModulus(), f = d2.Dot(r);
    double s = 0.0, t = 0.0;
    if (a <= 1e-30 && e <= 1e-30) return p1.Distance(p2);
    if (a <= 1e-30) {
        t = std::clamp(f / e, 0.0, 1.0);
    } else {
        double c = d1.Dot(r);
        if (e <= 1e-30) {
            s = std::clamp(-c / a, 0.0, 1.0);
        } else {
            double b = d1.Dot(d2);
            double denom = a * e - b * b;
            if (denom > 1e-30) s = std::clamp((b * f - c * e) / denom, 0.0, 1.0);
            t = (b * s + f) / e;
            if (t < 0.0)      { t = 0.0; s = std::clamp(-c / a, 0.0, 1.0); }
            else if (t > 1.0) { t = 1.0; s = std::clamp((b - c) / a, 0.0, 1.0); }
        }
    }
    gp_XYZ c1 = p1.XYZ() + d1 * s;
    gp_XYZ c2 = p2.XYZ() + d2 * t;
    return std::sqrt((c1 - c2).SquareModulus());
}

int curveSampleCount(double radius, double azSpan, double wireRadius) {
    double maxSag = samplingSag(wireRadius);
    double stepAz = (radius > 1e-12)
                        ? 2.0 * std::acos(std::clamp(1.0 - maxSag / radius, 0.0, 1.0))
                        : std::max(azSpan, 1e-3);
    stepAz = std::clamp(stepAz, 1e-3, 0.2);
    return std::max(2, static_cast<int>(std::ceil(azSpan / stepAz)) + 1);
}

std::vector<gp_Pnt> samplePrim(const Primitive& p, double wireRadius) {
    if (p.kind == Primitive::SEG) return {p.seg.a, p.seg.b};
    if (p.kind == Primitive::ARC3) {
        const Arc3& a = p.arc;
        double radius = a.v0.Modulus();
        int n = curveSampleCount(radius, std::abs(a.sweep), wireRadius);
        std::vector<gp_Pnt> pts;
        pts.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            double t = a.sweep * i / (n - 1);
            pts.push_back(gp_Pnt(a.c.XYZ() + rotateXYZ(a.v0, a.axis, t)));
        }
        return pts;
    }
    if (p.kind == Primitive::BLEND) {
        const Blend& bl = p.blendc;
        gp_XYZ d = bl.b.XYZ() - bl.a.XYZ();
        double L = d.Dot(bl.u);
        gp_XYZ perp = d - bl.u * L;
        double perpM = perp.Modulus();
        if (perpM < 1e-12) return {bl.a, bl.b};
        // Peak curvature of the cosine blend: kappa = perp * pi^2 / (2 L^2); sample with
        // the same sag rule as the arcs, on the tightest osculating radius.
        double rMin = 2.0 * L * L / (kPi * kPi * perpM);
        double maxSag = samplingSag(wireRadius);
        double step = 2.0 * std::sqrt(std::max(1e-18, 2.0 * rMin * maxSag));
        double arcLen = std::sqrt(L * L + perpM * perpM) * 1.1;
        int n = std::clamp(static_cast<int>(std::ceil(arcLen / step)) + 1, 8, 512);
        std::vector<gp_Pnt> pts;
        pts.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            double t = static_cast<double>(i) / (n - 1);
            pts.push_back(gp_Pnt(bl.a.XYZ() + bl.u * (L * t) +
                                 perp * (0.5 * (1.0 - std::cos(kPi * t)))));
        }
        return pts;
    }
    const Spiral& sp = p.spiral;
    int n = curveSampleCount(std::max(sp.r0, sp.r1), std::abs(sp.az1 - sp.az0), wireRadius);
    std::vector<gp_Pnt> pts;
    pts.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        double t = static_cast<double>(i) / (n - 1);
        double f = sp.blend ? 0.5 * (1.0 - std::cos(kPi * t)) : t;
        pts.push_back(azPointC(sp.cx, sp.cz, sp.r0 + (sp.r1 - sp.r0) * f,
                               sp.y0 + (sp.y1 - sp.y0) * f, sp.az0 + (sp.az1 - sp.az0) * t));
    }
    return pts;
}

// A correct ROUND sweep keeps every centreline point at the centre of its circular section, so the
// nearest solid wall is ~wireRadius away. A framing that drifts the section off the spine (the toroid
// high-torsion failure mode of MakePipeShell's corrected Frenet) leaves an interior centreline point
// closer to one wall. Verify the sampled interior points all clear the floor; if not, the single body
// is geometrically wrong even though it is watertight and the right volume -- reject it.
static bool centrelineStaysCentred(const TopoDS_Shape& solid,
                                   const std::vector<gp_Pnt>& interiorPts, double wireRadius) {
    const double floor = 0.85 * wireRadius;  // backstop vs gross drift; a good MakePipe holds >0.98
    for (const auto& pt : interiorPts) {
        TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(pt);
        BRepExtrema_DistShapeShape d(v, solid);
        if (!d.IsDone() || d.Value() < floor) return false;
    }
    return true;
}

double polyPolyDistance(const std::vector<gp_Pnt>& polyA, const std::vector<gp_Pnt>& polyB) {
    double best = std::numeric_limits<double>::max();
    for (size_t i = 0; i + 1 < polyA.size(); ++i)
        for (size_t j = 0; j + 1 < polyB.size(); ++j)
            best = std::min(best, segSegDistance(polyA[i], polyA[i + 1], polyB[j], polyB[j + 1]));
    return best;
}

std::pair<gp_Pnt, gp_Pnt> primEndpoints(const Primitive& p) {
    if (p.kind == Primitive::ARC3) {
        return {gp_Pnt(p.arc.c.XYZ() + p.arc.v0),
                gp_Pnt(p.arc.c.XYZ() + rotateXYZ(p.arc.v0, p.arc.axis, p.arc.sweep))};
    }
    if (p.kind == Primitive::SPIRAL) {
        return {azPointC(p.spiral.cx, p.spiral.cz, p.spiral.r0, p.spiral.y0, p.spiral.az0),
                azPointC(p.spiral.cx, p.spiral.cz, p.spiral.r1, p.spiral.y1, p.spiral.az1)};
    }
    if (p.kind == Primitive::BLEND) return {p.blendc.a, p.blendc.b};
    return {p.seg.a, p.seg.b};
}

bool shareEndpoint(const Primitive& a, const Primitive& b) {
    auto [a0, a1] = primEndpoints(a);
    auto [b0, b1] = primEndpoints(b);
    constexpr double tol = 1e-9;
    return a0.Distance(b0) < tol || a0.Distance(b1) < tol ||
           a1.Distance(b0) < tol || a1.Distance(b1) < tol;
}

// The user-mandated hard gate: any wire-envelope overlap anywhere throws. Nothing is
// ever moved to resolve a collision — the winding data (or MKF's blocking) must change.
// Rotate an entire conductor path about the column axis (+Y through the origin) by `theta`. Every
// concentric winding is drawn with its seam/leads on the SAME reference plane (kPlaneAz), so different
// windings' terminal leads land in the same plane and collide. Rotating a winding rigidly about Y just
// moves its seam to a different angle — the turns stay the same rings at the same radius/height, the
// magnetics are unchanged — so staggering each winding's angle sends their leads out at distinct
// azimuths and they no longer overlap. Points rotate in XZ (Y fixed); the spiral is azimuth-based
// (centre on the axis, cx=cz=0) so its start/end azimuths simply shift by theta.
inline gp_Pnt rotYpt(const gp_Pnt& p, double c, double s) {
    return gp_Pnt(p.X() * c + p.Z() * s, p.Y(), -p.X() * s + p.Z() * c);
}
inline gp_XYZ rotYvec(const gp_XYZ& v, double c, double s) {
    return gp_XYZ(v.X() * c + v.Z() * s, v.Y(), -v.X() * s + v.Z() * c);
}
void rotatePathAboutY(ConductorPath& path, double theta) {
    if (theta == 0.0) return;
    const double c = std::cos(theta), s = std::sin(theta);
    for (auto& pr : path.prims) {
        switch (pr.kind) {
            case Primitive::SEG:
                pr.seg.a = rotYpt(pr.seg.a, c, s); pr.seg.b = rotYpt(pr.seg.b, c, s); break;
            case Primitive::ARC3:
                pr.arc.c = rotYpt(pr.arc.c, c, s);
                pr.arc.axis = rotYvec(pr.arc.axis, c, s); pr.arc.v0 = rotYvec(pr.arc.v0, c, s); break;
            case Primitive::SPIRAL: {
                const double ncx = pr.spiral.cx * c + pr.spiral.cz * s;
                const double ncz = -pr.spiral.cx * s + pr.spiral.cz * c;
                pr.spiral.cx = ncx; pr.spiral.cz = ncz;
                pr.spiral.az0 += theta; pr.spiral.az1 += theta; break;
            }
            case Primitive::BLEND:
                pr.blendc.a = rotYpt(pr.blendc.a, c, s); pr.blendc.b = rotYpt(pr.blendc.b, c, s);
                pr.blendc.u = rotYvec(pr.blendc.u, c, s); break;
        }
    }
}

// Excise REDUNDANT EXCURSIONS from a conductor path: stretches that leave a point and come back
// to it, contributing no net progress. Two shapes of this were observed, and one rule covers both:
//   * a RETRACE pair -- A->B immediately followed by B->A;
//   * a zig-zag that traverses the SAME segment twice in the SAME direction with a closed loop in
//     between. Measured on the litz round-column fixture, the connection router emitted
//     ... -> y=-4.862, across, up to y=-0.940, across, back DOWN to y=-4.862, and across AGAIN on
//     the identical segment: prims 21 and 25 had identical endpoints and direction.
// Either way the wire would lay DOUBLE material along the repeated line, which no real conductor
// does, and topologically the excursion closes a loop. A spine with loops is not a simple open path,
// so BRepAdaptor_CompCurve cannot traverse it and both MakePipeShell and MakePipe abort with the
// deeply unhelpful "BRepAdaptor_Curve::No geometry" -- the single-body export then silently degrades
// to a per-turn compound. (That cost a long detour: every spine edge had a valid 3D curve and a sane
// parameter range; the fault was topology, not geometry.)
//
// Removal is LOSSLESS: the excursion starts and ends at the same point, so dropping it leaves the
// neighbours joined exactly as before -- verified on the fixture, where the primitive after the
// excursion started precisely where the one before it ended.
//
// Restricted to LEAD / CONNECTION primitives. A WRAP must never be touched: its turns are the
// physics, and a legitimately closed ring would look like an excursion to this rule.
std::size_t dropRedundantExcursions(ConductorPath& path) {
    constexpr double kTol = 1e-9;   // observed bit-identical
    std::size_t dropped = 0;
    bool again = true;
    while (again) {
        again = false;
        const std::size_t n = path.prims.size();
        for (std::size_t i = 0; i < n && !again; ++i) {
            if (!(path.prims[i].isLead || path.prims[i].isConnection)) continue;
            const gp_Pnt start_i = primEndpoints(path.prims[i]).first;
            for (std::size_t j = i; j < n; ++j) {
                if (!(path.prims[j].isLead || path.prims[j].isConnection)) break;
                if (primEndpoints(path.prims[j]).second.Distance(start_i) > kTol) continue;
                // prims [i..j] leave start_i and return to it: a closed excursion.
                path.prims.erase(path.prims.begin() + i, path.prims.begin() + j + 1);
                dropped += (j - i + 1);
                again = true;
                break;
            }
        }
    }
    return dropped;
}

void checkCollisions(const std::vector<ConductorPath>& paths) {
    // MVB_LEAD_NO_VALIDATE: diagnostic export mode (see emitToroLead) -- the emitted copper is
    // EXPECTED to interpenetrate; skipping the gate lets the STEP reach disk for inspection.
    if (std::getenv("MVB_LEAD_NO_VALIDATE")) {
        std::cerr << "[ConductorBuilder] MVB_LEAD_NO_VALIDATE set: collision gate SKIPPED -- "
                     "diagnostic geometry, not for FEM\n";
        return;
    }
    // Pre-sample every primitive once (rect/toroidal wraps are 9-10 primitives per turn;
    // re-sampling per pair would dominate the gate).
    std::vector<std::vector<std::vector<gp_Pnt>>> polys(paths.size());
    for (size_t ci = 0; ci < paths.size(); ++ci) {
        polys[ci].reserve(paths[ci].prims.size());
        for (const auto& pr : paths[ci].prims)
            polys[ci].push_back(samplePrim(pr, paths[ci].wireRadius));
    }
    for (size_t ci = 0; ci < paths.size(); ++ci) {
        for (size_t cj = ci; cj < paths.size(); ++cj) {
            const auto& A = paths[ci];
            const auto& B = paths[cj];
            const bool rectPair = A.isRectangular || B.isRectangular;
            // Gate on the CONDUCTING (copper) envelope: overlapping copper is the real fault (short /
            // un-meshable solid); enamel touching on a tight winding is fine. Insulated (paintCoating)
            // radii would false-reject every tightly-wound layer.
            double minGap = gateMinSeparation(A.condRadius, B.condRadius) -
                            samplingSag(A.wireRadius) - samplingSag(B.wireRadius);
            for (size_t i = 0; i < A.prims.size(); ++i) {
                size_t jStart = (ci == cj) ? i + 1 : 0;
                for (size_t j = jStart; j < B.prims.size(); ++j) {
                    const auto& pa = A.prims[i];
                    const auto& pb = B.prims[j];
                    if (ci == cj &&
                        (pa.turnOrdinal > pb.turnOrdinal ? pa.turnOrdinal - pb.turnOrdinal
                                                         : pb.turnOrdinal - pa.turnOrdinal) <= 1) {
                        continue;
                    }
                    if (ci == cj && shareEndpoint(pa, pb)) continue;
                    double d = polyPolyDistance(polys[ci][i], polys[cj][j]);
                    if (std::getenv("MVB_GATE_DIAG") && d < minGap) {
                        // Closest sampled pair, so a reported overlap can be located in space.
                        double best = 1e30; gp_Pnt ba, bb;
                        for (const auto& pa2 : polys[ci][i])
                            for (const auto& pb2 : polys[cj][j]) {
                                const double dd = pa2.SquareDistance(pb2);
                                if (dd < best) { best = dd; ba = pa2; bb = pb2; }
                            }
                        std::fprintf(stderr,
                            "[gate] '%s' vs '%s' d=%.4fmm need=%.4fmm closest A=(%.3f,%.3f,%.3f) "
                            "rA=%.3f  B=(%.3f,%.3f,%.3f) rB=%.3f\n",
                            pa.label.c_str(), pb.label.c_str(), d*1e3, minGap*1e3,
                            ba.X()*1e3, ba.Y()*1e3, ba.Z()*1e3, std::hypot(ba.X(), ba.Z())*1e3,
                            bb.X()*1e3, bb.Y()*1e3, bb.Z()*1e3, std::hypot(bb.X(), bb.Z())*1e3);
                    }
                    // Rectangular wire: the round capsule is wrong -- the section is width(radial) x
                    // height(axial). Decompose the closest-approach vector into its AXIAL and in-
                    // plane (radial/tangential) parts and require the section boxes to overlap in
                    // BOTH. The axial axis is the column Y for concentric windings, or the toroid's
                    // AZIMUTHAL tangent (the direction adjacent turns are spaced) at the contact.
                    // (Approximate: the in-plane part lumps radial with tangential; round pairs keep
                    // the exact capsule test above unchanged.)
                    if (rectPair && d < 0.5 * std::hypot(A.condWidth + A.condHeight,
                                                         B.condWidth + B.condHeight)) {
                        gp_Pnt ca, cb;
                        double best = std::numeric_limits<double>::max();
                        for (const auto& va : polys[ci][i])
                            for (const auto& vb : polys[cj][j]) {
                                double dd = va.SquareDistance(vb);
                                if (dd < best) { best = dd; ca = va; cb = vb; }
                            }
                        gp_XYZ sep = cb.XYZ() - ca.XYZ();
                        gp_XYZ axialDir(0, 1, 0);
                        if (A.toroidal || B.toroidal) {
                            gp_XYZ mid = (ca.XYZ() + cb.XYZ()) * 0.5;
                            gp_XYZ az(-mid.Z(), 0.0, mid.X());  // azimuthal tangent (ring in XZ)
                            if (az.Modulus() > 1e-9) axialDir = az / az.Modulus();
                        }
                        double axialSep = std::abs(sep.Dot(axialDir));
                        double inPlaneSep =
                            std::sqrt(std::max(0.0, sep.SquareModulus() - axialSep * axialSep));
                        double axA = A.isRectangular ? A.condHeight / 2.0 : A.condRadius;
                        double axB = B.isRectangular ? B.condHeight / 2.0 : B.condRadius;
                        double ipA = A.isRectangular ? A.condWidth / 2.0 : A.condRadius;
                        double ipB = B.isRectangular ? B.condWidth / 2.0 : B.condRadius;
                        const double sagAB =
                            samplingSag(A.wireRadius) + samplingSag(B.wireRadius);
                        bool overlap = axialSep < (axA + axB - sagAB) &&
                                       inPlaneSep < (ipA + ipB - sagAB);
                        if (!overlap) continue;
                        d = std::min(axialSep, inPlaneSep);  // report the tighter gap
                    } else if (rectPair) {
                        continue;  // beyond the loosest rectangle reach -> no overlap
                    } else if (d >= minGap) {
                        continue;
                    }
                    {
                        auto [aa, ab] = primEndpoints(pa);
                        auto [ba, bb] = primEndpoints(pb);
                        auto pt = [](const gp_Pnt& p) {
                            std::ostringstream o;
                            o << "(" << p.X() << "," << p.Y() << "," << p.Z() << ")";
                            return o.str();
                        };
                        std::ostringstream s;
                        s << "ConductorBuilder: collision between " << A.name << " ["
                          << pa.label << " " << pt(aa) << "->" << pt(ab) << "] and " << B.name
                          << " [" << pb.label << " " << pt(ba) << "->" << pt(bb)
                          << "]: centreline distance " << d << " m < "
                          << (A.condRadius + B.condRadius)
                          << " m (copper envelopes overlap). Turn positions are never moved; "
                             "fix the winding data or the MKF blocking (see MKF ABT #187).";
                        throw std::runtime_error(s.str());
                    }
                }
            }
        }
    }
}

// --- geometry emission ------------------------------------------------------------------
TopoDS_Wire wireProfileWire(const gp_Pnt& center, const gp_Dir& normal, double radius,
                            int segments) {
    // DETERMINISTIC section orientation. gp_Ax2(center, normal) lets OCC choose the X
    // direction, so two primitives meeting with the SAME tangent could still get n-gons rotated
    // relative to each other: their abutting faces then only partially coincide, and at a mitre
    // junction the pieces can miss entirely (measured: a 56-primitive rect-column conductor fell
    // into 7 components with no endpoint gaps and nothing dropped). Deriving X from a fixed world
    // reference makes the orientation a pure function of the tangent, so tangent neighbours share
    // an exactly coincident section -- which is what conformal abutment requires.
    gp_XYZ ref(0, 1, 0);
    if (std::abs(normal.XYZ().Dot(ref)) > 0.9) ref = gp_XYZ(1, 0, 0);
    gp_XYZ xdir = ref - normal.XYZ() * normal.XYZ().Dot(ref);
    gp_Ax2 plane(center, normal, gp_Dir(xdir));
    if (segments <= 0) {
        gp_Circ circ(plane, radius);
        return BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire();
    }
    BRepBuilderAPI_MakePolygon poly;
    gp_Dir dx = plane.XDirection();
    gp_Dir dy = plane.YDirection();
    double offset = kPi / segments;
    // INSCRIBED n-gon (vertices on the nominal circle): the faceted wire never protrudes past
    // the wire's real envelope, which is what the zero-tolerance pairwise-overlap contracts rely
    // on. Its flats consequently sit at the apothem, r*cos(pi/n) -- 0.981*r at n=16 -- so any
    // clearance or coverage check against a faceted section must use the apothem, not r.
    for (int i = 0; i < segments; ++i) {
        double ang = kTwoPi * i / segments + offset;
        gp_XYZ p = center.XYZ() + dx.XYZ() * (radius * std::cos(ang)) +
                   dy.XYZ() * (radius * std::sin(ang));
        poly.Add(gp_Pnt(p));
    }
    poly.Close();
    return poly.Wire();
}
TopoDS_Face wireProfile(const gp_Pnt& center, const gp_Dir& normal, double radius,
                        int segments) {
    return BRepBuilderAPI_MakeFace(wireProfileWire(center, normal, radius, segments)).Face();
}

// The same round profile, but as TWO half-circle edges instead of one closed circle.
// Sweeping a CLOSED circular edge produces a single PERIODIC B-spline surface carrying a seam
// curve, and gmsh cannot order the boundary loop of such a face -- it aborts the 2D stage with
// "the 1D mesh seems not to be forming a closed loop" (the offending face reported simply as
// "BSpline surface"). That is what made the real-winding conductor unmeshable: every other
// primitive is an analytic cylinder or torus, which mesh fine; only the swept axial transitions
// carried a seam. Split in two and the sweep yields ordinary patches with real boundaries -- the
// identical solid, now meshable. (ABT #332)
TopoDS_Wire wireProfileWireSplit(const gp_Pnt& center, const gp_Dir& normal, double radius,
                                 int pieces = 2) {
    gp_Circ circ(gp_Ax2(center, normal), radius);
    BRepBuilderAPI_MakeWire w;
    for (int i = 0; i < pieces; ++i)
        w.Add(BRepBuilderAPI_MakeEdge(circ, kTwoPi * i / pieces, kTwoPi * (i + 1) / pieces).Edge());
    return w.Wire();
}

// Oriented RECTANGLE profile for rectangular/planar wire, in the plane perpendicular to the wire's
// travel (tangent). The AXIAL in-plane axis is the column axis projected perpendicular to the
// tangent; the RADIAL axis is tangent x axial. width spans radial, height spans axial -- matching
// TurnBuilder::build_rect_profile so the real-winding section is consistent with the per-turn one.
TopoDS_Wire rectProfileWire(const gp_Pnt& center, const gp_Dir& tangent, const gp_Dir& axialAxis,
                            double width, double height) {
    gp_XYZ t = tangent.XYZ();
    gp_XYZ a = axialAxis.XYZ() - t * (axialAxis.XYZ().Dot(t));  // axial, ⊥ tangent
    if (a.Modulus() < 1e-9) {                                   // tangent ∥ axial: pick any ⊥ axis
        gp_XYZ ref = std::abs(t.X()) < 0.9 ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 0, 1);
        a = ref - t * (ref.Dot(t));
    }
    a.Normalize();
    gp_XYZ w = t.Crossed(a);                                    // radial
    w.Normalize();
    double hw = width / 2.0, hh = height / 2.0;
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(center.XYZ() - w * hw - a * hh));
    poly.Add(gp_Pnt(center.XYZ() + w * hw - a * hh));
    poly.Add(gp_Pnt(center.XYZ() + w * hw + a * hh));
    poly.Add(gp_Pnt(center.XYZ() - w * hw + a * hh));
    poly.Close();
    return poly.Wire();
}

// One analytic edge per primitive, so the whole conductor is ONE path:
//  - SEG: straight edge.
//  - ARC3: exact circle edge trimmed to [0, sweep].
//  - SPIRAL at constant radius: TRUE HELIX — a 2D line in the (U = azimuth, V = height)
//    parametric space of a cylinder about its vertical axis (the canonical OCCT
//    construction: Geom2d line on Geom_CylindricalSurface + BRepLib::BuildCurves3d).
//  - SPIRAL with varying radius (conical wraps, oblong cap ramps): BSpline through the
//    sampled centreline.
TopoDS_Edge primEdge(const Primitive& pr, double wireRadius) {
    if (pr.kind == Primitive::SEG) {
        if (pr.seg.a.Distance(pr.seg.b) < 1e-12) return TopoDS_Edge();
        return BRepBuilderAPI_MakeEdge(pr.seg.a, pr.seg.b).Edge();
    }
    if (pr.kind == Primitive::ARC3) {
        double radius = pr.arc.v0.Modulus();
        if (radius < 1e-12 || std::abs(pr.arc.sweep) < 1e-12) return TopoDS_Edge();
        try {
            // gp_Circ parametrisation P(u) = c + r (cos u XDir + sin u YDir) with
            // YDir = axis x XDir matches rotateXYZ (axis is perpendicular to v0).
            gp_Ax2 ax2(pr.arc.c, gp_Dir(pr.arc.axis), gp_Dir(pr.arc.v0));
            gp_Circ circ(ax2, radius);
            return BRepBuilderAPI_MakeEdge(circ, 0.0, pr.arc.sweep).Edge();
        } catch (const Standard_Failure&) {
            return TopoDS_Edge();
        }
    }
    if (pr.kind == Primitive::SPIRAL) {
        // Every spiral lives on an ANALYTIC surface of revolution about its vertical
        // axis: a cylinder when the radius is constant, a cone otherwise (radius and
        // height share the blend factor, so the (r, y) trace is a straight meridian).
        // Building the edge as a 2D curve in the surface's (U = azimuth, V) parameter
        // space (+ BRepLib::BuildCurves3d) is the one construction OCCT's pipe-shell
        // sweeps reliably — free-form interpolated 3D splines, though positionally
        // tight, produce unbounded surface flares when swept.
        const Spiral& sp = pr.spiral;
        double dr = sp.r1 - sp.r0;
        double dy = sp.y1 - sp.y0;
        bool constantRadius = std::abs(dr) < 1e-12;
        if (constantRadius || std::abs(dy) > 1e-12) {
            try {
                // Frame with U measured from +X toward -Z — exactly the azPointC()
                // convention: P(U,V) matches (cx + r cos U, y, cz - r sin U).
                gp_Ax3 axis(gp_Pnt(sp.cx, sp.y0, sp.cz), gp_Dir(0, 1, 0),
                            gp_Dir(1, 0, 0));
                Handle(Geom_Surface) surf;
                double v1;   // V at the end (V = 0 at the start)
                if (constantRadius) {
                    axis.SetLocation(gp_Pnt(sp.cx, 0, sp.cz));
                    surf = new Geom_CylindricalSurface(axis, sp.r0);
                    // Cylinder V is plain height.
                    v1 = dy;
                } else {
                    double alpha = std::atan(dr / dy);
                    surf = new Geom_ConicalSurface(axis, alpha, sp.r0);
                    v1 = dy / std::cos(alpha);
                }
                double v0 = constantRadius ? sp.y0 : 0.0;
                TopoDS_Edge e;
                if (!sp.blend) {
                    Handle(Geom2d_TrimmedCurve) seg2d =
                        GCE2d_MakeSegment(gp_Pnt2d(sp.az0, v0),
                                          gp_Pnt2d(sp.az1, v0 + v1)).Value();
                    e = BRepBuilderAPI_MakeEdge(seg2d, surf).Edge();
                } else {
                    // Cosine blend of V over U — end tangents purely azimuthal, so the
                    // junctions into the neighbouring on-station geometry are exact.
                    int n = curveSampleCount(std::max(sp.r0, sp.r1),
                                             std::abs(sp.az1 - sp.az0), wireRadius);
                    Handle(TColgp_HArray1OfPnt2d) arr =
                        new TColgp_HArray1OfPnt2d(1, n);
                    for (int i = 0; i < n; ++i) {
                        double t = static_cast<double>(i) / (n - 1);
                        double f = 0.5 * (1.0 - std::cos(kPi * t));
                        arr->SetValue(i + 1, gp_Pnt2d(sp.az0 + (sp.az1 - sp.az0) * t,
                                                      v0 + v1 * f));
                    }
                    Geom2dAPI_Interpolate interp(arr, Standard_False, 1e-12);
                    interp.Load(gp_Vec2d(1, 0), gp_Vec2d(1, 0), Standard_True);
                    interp.Perform();
                    if (!interp.IsDone() || interp.Curve().IsNull())
                        return TopoDS_Edge();
                    e = BRepBuilderAPI_MakeEdge(interp.Curve(), surf).Edge();
                }
                BRepLib::BuildCurves3d(e);
                return e;
            } catch (const Standard_Failure&) {
                return TopoDS_Edge();
            }
        }
        if (sp.blend) {
            // Flat blended spiral (radius varies at constant height): a 2D spiral in
            // the horizontal plane's parameter space, interpolated with exact azimuthal
            // end tangents.
            try {
                gp_Ax3 frame(gp_Pnt(sp.cx, sp.y0, sp.cz), gp_Dir(0, 1, 0),
                             gp_Dir(1, 0, 0));
                Handle(Geom_Plane) plane = new Geom_Plane(frame);
                // Plane P(U,V) = O + U XDir + V YDir with YDir = -Z, so the azPointC
                // trace maps to (U, V) = (r cos az, r sin az).
                int n = curveSampleCount(std::max(sp.r0, sp.r1),
                                         std::abs(sp.az1 - sp.az0), wireRadius);
                Handle(TColgp_HArray1OfPnt2d) arr = new TColgp_HArray1OfPnt2d(1, n);
                for (int i = 0; i < n; ++i) {
                    double t = static_cast<double>(i) / (n - 1);
                    double f = 0.5 * (1.0 - std::cos(kPi * t));
                    double az = sp.az0 + (sp.az1 - sp.az0) * t;
                    double r = sp.r0 + dr * f;
                    arr->SetValue(i + 1, gp_Pnt2d(r * std::cos(az), r * std::sin(az)));
                }
                gp_Vec2d ta(-std::sin(sp.az0), std::cos(sp.az0));
                gp_Vec2d tb(-std::sin(sp.az1), std::cos(sp.az1));
                Geom2dAPI_Interpolate interp(arr, Standard_False, 1e-12);
                interp.Load(ta, tb, Standard_True);
                interp.Perform();
                if (!interp.IsDone() || interp.Curve().IsNull()) return TopoDS_Edge();
                TopoDS_Edge e = BRepBuilderAPI_MakeEdge(interp.Curve(), plane).Edge();
                BRepLib::BuildCurves3d(e);
                return e;
            } catch (const Standard_Failure&) {
                return TopoDS_Edge();
            }
        }
        // Flat LINEAR spiral: falls through to the sampled BSpline below.
    }
    if (pr.kind == Primitive::BLEND) {
        // The S-blend is PLANAR (the longitudinal direction and the offset span one
        // plane), so it too builds as a 2D cosine in an analytic surface's parameter
        // space — the pipe-reliable construction (see the spiral comment above).
        const Blend& bl = pr.blendc;
        gp_XYZ d = bl.b.XYZ() - bl.a.XYZ();
        double L = d.Dot(bl.u);
        gp_XYZ perp = d - bl.u * L;
        double perpM = perp.Modulus();
        if (L < 1e-12) return TopoDS_Edge();
        if (perpM < 1e-12) return BRepBuilderAPI_MakeEdge(bl.a, bl.b).Edge();
        try {
            gp_Ax3 frame(bl.a, gp_Dir(bl.u.Crossed(perp)), gp_Dir(bl.u));
            Handle(Geom_Plane) plane = new Geom_Plane(frame);
            // Plane P(U,V) = O + U*XDir + V*YDir with XDir = u, YDir = axis x XDir =
            // unit(perp).
            double stepL = 2.0 * std::sqrt(std::max(
                1e-18, 4.0 * L * L / (kPi * kPi * perpM) *
                           samplingSag(wireRadius)));
            int n = std::clamp(static_cast<int>(std::ceil(L / stepL)) + 1, 8, 512);
            Handle(TColgp_HArray1OfPnt2d) arr = new TColgp_HArray1OfPnt2d(1, n);
            for (int i = 0; i < n; ++i) {
                double t = static_cast<double>(i) / (n - 1);
                arr->SetValue(i + 1, gp_Pnt2d(L * t,
                                              perpM * 0.5 * (1.0 - std::cos(kPi * t))));
            }
            Geom2dAPI_Interpolate interp(arr, Standard_False, 1e-12);
            interp.Load(gp_Vec2d(1, 0), gp_Vec2d(1, 0), Standard_True);
            interp.Perform();
            if (!interp.IsDone() || interp.Curve().IsNull()) return TopoDS_Edge();
            TopoDS_Edge e = BRepBuilderAPI_MakeEdge(interp.Curve(), plane).Edge();
            BRepLib::BuildCurves3d(e);
            return e;
        } catch (const Standard_Failure&) {
            return TopoDS_Edge();
        }
    }
    // Flat varying-radius spirals at exactly constant height: approximated BSpline
    // through the sampled centreline.
    auto pts = samplePrim(pr, wireRadius);
    if (pts.size() < 2) return TopoDS_Edge();
    try {
        TColgp_Array1OfPnt arr(1, static_cast<Standard_Integer>(pts.size()));
        for (size_t i = 0; i < pts.size(); ++i)
            arr.SetValue(static_cast<Standard_Integer>(i + 1), pts[i]);
        Handle(Geom_BSplineCurve) bs = GeomAPI_PointsToBSpline(arr).Curve();
        if (bs.IsNull()) return TopoDS_Edge();
        return BRepBuilderAPI_MakeEdge(bs).Edge();
    } catch (const Standard_Failure&) {
        return TopoDS_Edge();
    }
}

// Sweep a contiguous span of primitives as ONE PipeShell (single spine wire). The span
// must not revisit any point — OCCT's pipe-shell corrupts the heap in curve-on-surface
// projection on self-touching spines (observed in 7.9: free() abort inside
// ShapeConstruct_ProjectCurveOnSurface), so callers split at those points.
// Build the conductor centerline as a G1 (tangent-continuous) wire for the single-body sweep.
// Every straight lead/link and every helix stays EXACT (analytic edge, so the pipe-shell does
// not flare the way it would on a free-form spline). At each SHARP corner the two segments are
// trimmed back by a hair and bridged with a short tangent-matched cubic -- rounding the wire
// bend the way a real wound wire physically bends. With no C0 corners left, MakePipeShell has
// nothing to "transition", so it produces a watertight solid instead of leaving gaps.
// Only round windings (straight segments + cylindrical helices) are handled here; any other
// primitive returns a null wire and the caller keeps the exact per-run compound.
TopoDS_Wire buildFilletedWire(const Primitive* const* prims, size_t count, double wireRadius) {
    if (count == 0) return TopoDS_Wire();
    const bool dbg = std::getenv("MVB_DIAG") != nullptr;
    auto kindName = [](int k) {
        return k == Primitive::SEG ? "SEG" : k == Primitive::ARC3 ? "ARC3"
             : k == Primitive::SPIRAL ? "SPIRAL" : "BLEND";
    };
    auto entryDir = [&](const Primitive& p) {
        auto pts = samplePrim(p, wireRadius);
        gp_Vec v(pts[1].XYZ() - pts.front().XYZ());
        return v.Magnitude() > 1e-12 ? gp_Dir(v) : gp_Dir(1, 0, 0);
    };
    auto exitDir = [&](const Primitive& p) {
        auto pts = samplePrim(p, wireRadius);
        gp_Vec v(pts.back().XYZ() - pts[pts.size() - 2].XYZ());
        return v.Magnitude() > 1e-12 ? gp_Dir(v) : gp_Dir(1, 0, 0);
    };
    auto arcLen = [&](const Primitive& p) {
        auto pts = samplePrim(p, wireRadius);
        double L = 0;
        for (size_t k = 1; k < pts.size(); ++k) L += pts[k].Distance(pts[k - 1]);
        return L;
    };
    // Trim distance at each end of each primitive: a fillet is needed only where consecutive
    // segments meet at an angle (a helix chains tangent-continuously to the next helix).
    std::vector<double> startTrim(count, 0.0), endTrim(count, 0.0);
    std::vector<double> lens(count, 0.0);
    for (size_t i = 0; i < count; ++i) lens[i] = arcLen(*prims[i]);
    // Fillet trim distance. For a 90-degree corner the inscribed bend radius equals the trim, so
    // this must EXCEED the wire radius: a tube of radius r swept through a bend of radius < r
    // self-intersects, and MakePipeShell/MakePipe then fail (this file already enforces the same
    // rule elsewhere as minBend = wireRadius * 1.02). At the old 0.8*wireRadius every 90-degree
    // lead/link corner asked for a bend TIGHTER than the wire, which is why the litz round-column
    // spine -- valid edge by edge, simple, loop-free -- still could not be swept.
    const double target = wireRadius * 1.25;
    // Only SEG / ARC3 / non-blend SPIRAL can be shrunk by arc length. A BLEND or a
    // cosine-blended SPIRAL has a non-linear arc-length map and cannot be trimmed, but it is
    // built to enter/leave its neighbours tangentially, so it rarely lands at a real corner.
    // When it does, route the WHOLE fillet trim onto the trimmable neighbour (an asymmetric
    // one-sided fillet is still G1); only if BOTH sides are untrimmable do we give up.
    auto trimmable = [](const Primitive& p) {
        return p.kind == Primitive::SEG || p.kind == Primitive::ARC3 ||
               (p.kind == Primitive::SPIRAL && !p.spiral.blend);
    };
    for (size_t i = 0; i + 1 < count; ++i) {
        if (exitDir(*prims[i]).Angle(entryDir(*prims[i + 1])) < 0.05) continue;
        bool ti = trimmable(*prims[i]), tj = trimmable(*prims[i + 1]);
        if (!ti && !tj) {
            if (dbg) std::cerr << "[spine] bail: both sides untrimmable at corner i=" << i
                               << " (" << kindName(prims[i]->kind) << " -> "
                               << kindName(prims[i + 1]->kind) << ")\n";
            return TopoDS_Wire();
        }
        // At a lead/link <-> wrap junction the WRAP endpoint is an MKF turn/crossing position
        // that must stay EXACTLY on the copper centerline; fillet only the lead/link side so the
        // wire still passes through the crossing (the rounded bend lives on the lead's approach).
        bool iLead = prims[i]->isLead || prims[i]->isConnection;
        bool jLead = prims[i + 1]->isLead || prims[i + 1]->isConnection;
        bool mixed = iLead != jLead;
        if (ti && (!mixed || iLead))
            endTrim[i] = std::max(endTrim[i], std::min(target, 0.4 * lens[i]));
        if (tj && (!mixed || jLead))
            startTrim[i + 1] = std::max(startTrim[i + 1], std::min(target, 0.4 * lens[i + 1]));
    }
    auto trimPrim = [](Primitive p, double sTrim, double eTrim, double totalLen) {
        if (p.kind == Primitive::SEG) {
            gp_XYZ dir = p.seg.b.XYZ() - p.seg.a.XYZ();
            double L = dir.Modulus();
            if (L > 1e-12) {
                gp_XYZ u = dir / L;
                p.seg.a = gp_Pnt(p.seg.a.XYZ() + u * sTrim);
                p.seg.b = gp_Pnt(p.seg.b.XYZ() - u * eTrim);
            }
        } else if (p.kind == Primitive::ARC3) {
            // Arc length is exactly radius * angle, so trimming by arc length is a pure
            // rotation of the start vector and a shrink of the sweep angle.
            double radius = p.arc.v0.Modulus();
            if (radius > 1e-12) {
                double sgn = p.arc.sweep > 0 ? 1.0 : -1.0;
                double dS = sgn * sTrim / radius;   // advance the start angle
                double dE = sgn * eTrim / radius;   // retract the end angle
                p.arc.v0 = rotateXYZ(p.arc.v0, p.arc.axis, dS);
                p.arc.sweep -= (dS + dE);
            }
        } else {  // non-blend SPIRAL (cylindrical helix, conical wrap, or flat linear spiral)
            double azSpan = p.spiral.az1 - p.spiral.az0;
            double speed = totalLen / std::max(1e-12, std::abs(azSpan));
            double sgn = azSpan > 0 ? 1.0 : -1.0;
            double naz0 = p.spiral.az0 + sgn * (sTrim / speed);
            double naz1 = p.spiral.az1 - sgn * (eTrim / speed);
            double f0 = (naz0 - p.spiral.az0) / azSpan;
            double f1 = (naz1 - p.spiral.az0) / azSpan;
            // Non-blend spiral: radius and height vary linearly with the azimuth fraction
            // (straight meridian), so interpolate BOTH ends off the ORIGINAL endpoints.
            double y0o = p.spiral.y0, dy = p.spiral.y1 - p.spiral.y0;
            double r0o = p.spiral.r0, dr = p.spiral.r1 - p.spiral.r0;
            p.spiral.y0 = y0o + dy * f0;
            p.spiral.y1 = y0o + dy * f1;
            p.spiral.r0 = r0o + dr * f0;
            p.spiral.r1 = r0o + dr * f1;
            p.spiral.az0 = naz0;
            p.spiral.az1 = naz1;
        }
        return p;
    };
    // Build each trimmed edge and capture its EXACT endpoint tangents (analytic D1, not the
    // sampled finite difference) so the bridging cubics are truly tangent to the neighbours --
    // any residual kink makes MakePipeShell fail on the mixed-curve spine.
    std::vector<TopoDS_Edge> edges(count);
    std::vector<gp_Pnt> ep0(count), ep1(count);
    std::vector<gp_Dir> et0(count), et1(count);
    for (size_t i = 0; i < count; ++i) {
        Primitive tp = trimPrim(*prims[i], startTrim[i], endTrim[i], lens[i]);
        TopoDS_Edge e = primEdge(tp, wireRadius);
        if (e.IsNull()) {
            if (dbg) std::cerr << "[spine] bail: primEdge NULL i=" << i << " kind="
                               << kindName(prims[i]->kind) << "\n";
            return TopoDS_Wire();
        }
        BRepAdaptor_Curve c(e);
        gp_Pnt a, b;
        gp_Vec da, db;
        c.D1(c.FirstParameter(), a, da);
        c.D1(c.LastParameter(), b, db);
        if (da.Magnitude() < 1e-12 || db.Magnitude() < 1e-12) {
            if (dbg) std::cerr << "[spine] bail: zero D1 tangent i=" << i << " kind="
                               << kindName(prims[i]->kind) << "\n";
            return TopoDS_Wire();
        }
        if (e.Orientation() == TopAbs_REVERSED) {
            std::swap(a, b);
            gp_Vec t = da; da = db.Reversed(); db = t.Reversed();
        }
        edges[i] = e;
        ep0[i] = a; et0[i] = gp_Dir(da);
        ep1[i] = b; et1[i] = gp_Dir(db);
    }
    BRepBuilderAPI_MakeWire wireMaker;
    for (size_t i = 0; i < count; ++i) {
        // A corner was filleted here (either side trimmed, possibly asymmetric) -> the two
        // trimmed ends leave a gap; bridge it with a tangent-matched cubic.
        if (i > 0 && (startTrim[i] > 0.0 || endTrim[i - 1] > 0.0)) {
            gp_Pnt A = ep1[i - 1], B = ep0[i];
            double gap = A.Distance(B);
            if (gap > 1e-9) {
                double s = gap / 3.0;
                TColgp_Array1OfPnt poles(1, 4);
                poles(1) = A;
                poles(2) = gp_Pnt(A.XYZ() + et1[i - 1].XYZ() * s);
                poles(3) = gp_Pnt(B.XYZ() - et0[i].XYZ() * s);
                poles(4) = B;
                Handle(Geom_BezierCurve) bez = new Geom_BezierCurve(poles);
                wireMaker.Add(BRepBuilderAPI_MakeEdge(bez).Edge());
                if (!wireMaker.IsDone()) {
                    if (dbg) std::cerr << "[spine] bail: wireMaker fail on fillet bridge i="
                                       << i << " gap=" << gap << "\n";
                    return TopoDS_Wire();
                }
            }
        }
        wireMaker.Add(edges[i]);
        if (!wireMaker.IsDone()) {
            if (dbg) std::cerr << "[spine] bail: wireMaker fail adding edge i=" << i << " kind="
                               << kindName(prims[i]->kind) << "\n";
            return TopoDS_Wire();
        }
    }
    if (!wireMaker.IsDone()) return TopoDS_Wire();
    TopoDS_Wire spine = wireMaker.Wire();
    // WELD near-coincident endpoints. Consecutive primitives are built independently and their
    // shared endpoints agree only to sub-micron precision, which is BELOW the vertex tolerance
    // BRepBuilderAPI_MakeWire welds at -- so it happily returns IsDone() with every edge still
    // carrying its own pair of vertices. The result looks fine edge by edge and even passes
    // BRepCheck_Analyzer, but it is a BAG OF EDGES, not a wire: measured on the litz round-column
    // spine, edges=65 vertices=130 (a connected open wire must have edges+1). Both MakePipe and
    // MakePipeShell drive a BRepAdaptor_CompCurve over the spine, which cannot traverse that and
    // aborts with the thoroughly misleading "BRepAdaptor_Curve::No geometry". FixConnected welds
    // vertices within tolerance and makes it a real wire. (1e-6 m is the same misalignment scale the
    // conductor fuse needed.)
    {
        ShapeFix_Wire wf;
        wf.Load(spine);
        wf.SetPrecision(1e-6);
        wf.SetMaxTolerance(1e-5);
        wf.ClearModes();
        wf.FixReorder();                 // order the edges head-to-tail first
        wf.FixConnected(1e-5);           // then weld endpoints within tolerance
        if (!wf.Wire().IsNull()) spine = wf.Wire();
    }
    // NEVER hand back a spine carrying an edge with no 3D curve. MakePipeShell adapts every spine
    // edge, so one geometry-less edge aborts all three framing modes with
    // "BRepAdaptor_Curve::No geometry" -- which surfaces as a bare "sweep NULL" and reads as a
    // framing problem when it is actually a malformed spine.
    {
        // Wire-LEVEL health, not just per-edge: MakePipe/MakePipeShell both drive a
        // BRepAdaptor_CompCurve over the spine, which needs an ORDERED, CONNECTED wire. A wire that
        // is valid edge-by-edge but disconnected still aborts with "BRepAdaptor_Curve::No geometry".
        // For a connected OPEN wire, vertices == edges + 1.
        // A sweepable spine must be a SIMPLE OPEN PATH. Both MakePipe and MakePipeShell drive a
        // BRepAdaptor_CompCurve along it, which can only traverse one. Count UNIQUE vertices (an
        // explorer revisits each shared vertex once per edge, so it reports 2*edges even when the
        // wire is perfectly connected): a simple open path has vertices == edges + 1, and any
        // shortfall is that many closed LOOPS. Loops arise when BRepBuilderAPI_MakeWire
        // proximity-welds two vertices that are geometrically coincident but NOT topologically
        // consecutive -- i.e. the conductor path revisits a position. OCCT reports this only as
        // "BRepAdaptor_Curve::No geometry" from inside the sweep, which sends you looking at curves
        // and framing instead of topology (it cost a long detour: every edge had a valid 3D curve
        // and a sane parameter range). Diagnose it here instead.
        TopTools_IndexedMapOfShape em, vm;
        TopExp::MapShapes(spine, TopAbs_EDGE, em);
        TopExp::MapShapes(spine, TopAbs_VERTEX, vm);
        const int ne = em.Extent(), nv = vm.Extent();
        const int loops = ne - nv + 1;   // one connected component
        if (dbg)
            std::cerr << "[spine] wire: edges=" << ne << " unique vertices=" << nv
                      << " loops=" << loops
                      << " closed=" << (BRep_Tool::IsClosed(spine) ? 1 : 0)
                      << " valid=" << (BRepCheck_Analyzer(spine).IsValid() ? 1 : 0) << "\n";
        if (loops > 0 && dbg) {
            // The primitives, so a loop junction can be tied to the pair that produced it.
            for (size_t q = 0; q < count; ++q) {
                auto [qa, qb] = primEndpoints(*prims[q]);
                std::cerr << "[spine]   prim " << q << " " << kindName(prims[q]->kind)
                          << (prims[q]->isLead ? " LEAD" : "")
                          << (prims[q]->isConnection ? " CONN" : "")
                          << " (" << qa.X()*1e3 << "," << qa.Y()*1e3 << "," << qa.Z()*1e3
                          << ") -> (" << qb.X()*1e3 << "," << qb.Y()*1e3 << "," << qb.Z()*1e3
                          << ") mm\n";
            }
            // Which vertices are the loop junctions? A simple open path has exactly two vertices of
            // valence 1 (the free ends) and the rest valence 2. Anything higher is where MakeWire
            // welded non-consecutive edges together, and its coordinates say whether the PATH truly
            // revisits that point or two merely-nearby points were fused.
            TopTools_IndexedDataMapOfShapeListOfShape v2e;
            TopExp::MapShapesAndAncestors(spine, TopAbs_VERTEX, TopAbs_EDGE, v2e);
            for (int k = 1; k <= v2e.Extent(); ++k) {
                const int val = v2e(k).Extent();
                if (val == 2) continue;
                gp_Pnt P = BRep_Tool::Pnt(TopoDS::Vertex(v2e.FindKey(k)));
                std::cerr << "[spine]   vertex valence=" << val << " at ("
                          << P.X()*1e3 << ", " << P.Y()*1e3 << ", " << P.Z()*1e3 << ") mm\n";
            }
        }
        if (loops > 0) {
            if (dbg) std::cerr << "[spine] bail: spine is not a simple open path (" << loops
                               << " loop(s)) -- the conductor path revisits a position, so "
                               << "MakeWire welded non-consecutive vertices. No sweep can "
                               << "traverse it; caller falls back to the per-run compound.\n";
            return TopoDS_Wire();
        }
        int idx = 0;
        for (TopExp_Explorer e(spine, TopAbs_EDGE); e.More(); e.Next(), ++idx) {
            double f = 0.0, l = 0.0;
            Handle(Geom_Curve) c = BRep_Tool::Curve(TopoDS::Edge(e.Current()), f, l);
            if (dbg) {
                GProp_GProps lp; BRepGProp::LinearProperties(e.Current(), lp);
                std::cerr << "[spine] edge " << idx << " len=" << lp.Mass()*1e3 << "mm"
                          << " curve=" << (c.IsNull() ? "NULL" : c->DynamicType()->Name())
                          << " range=[" << f << "," << l << "]\n";
            }
            if (c.IsNull() || BRep_Tool::Degenerated(TopoDS::Edge(e.Current()))) {
                if (dbg) std::cerr << "[spine] bail: spine edge " << idx
                                   << " has no 3D curve (degenerated="
                                   << (BRep_Tool::Degenerated(TopoDS::Edge(e.Current())) ? 1 : 0)
                                   << ")\n";
                return TopoDS_Wire();
            }
        }
    }
    return spine;
}

// Sweep an already-built (G1) spine wire into a solid.
//  - Round wire: an exact circle profile. Frenet framing keeps the round section centered exactly
//    on the spine but cannot close a spine with inflection points; corrected Frenet parallel-
//    transports the frame and closes those, discrete is the last resort. (Round-section sweeps are
//    rotation-invariant, so any frame that closes keeps the section on the spine.)
//  - Rectangular wire: an ORIENTED rectangle profile, swept with a FIXED BINORMAL along the column
//    axis so the section's flat faces stay axial/radial the whole way round (round columns only, so
//    no primitive travels axially to make that frame flip). The section is NOT rotation-invariant,
//    so exactly one frame is correct -- no fallback ladder.
TopoDS_Shape sweepWire(const TopoDS_Wire& spine, const gp_Pnt& p0, const gp_Dir& t0,
                       double wireRadius, int profileSegments, bool rectangular = false,
                       double rectWidth = 0.0, double rectHeight = 0.0,
                       const gp_Dir& axialAxis = gp_Dir(0, 1, 0),
                       bool preferSimplePipe = false, bool tryFixedBinormal = false) {
    if (rectangular) {
        try {
            TopoDS_Wire prof = rectProfileWire(p0, t0, axialAxis, rectWidth, rectHeight);
            BRepOffsetAPI_MakePipeShell ps(spine);
            ps.SetMode(axialAxis);   // fixed binormal = column axis: section stays axial/radial
            ps.Add(prof);
            ps.Build();
            if (!ps.IsDone() || !ps.MakeSolid()) {
                if (std::getenv("MVB_DIAG"))
                    std::cerr << "[sweepWire] rect fixed-binormal status="
                              << (int)ps.GetStatus() << " isDone=" << ps.IsDone() << "\n";
                return TopoDS_Shape();
            }
            TopoDS_Shape s = ps.Shape();
            if (BRepCheck_Analyzer(s).IsValid()) return s;
            if (std::getenv("MVB_DIAG"))
                std::cerr << "[sweepWire] rect fixed-binormal: swept but INVALID\n";
        } catch (const Standard_Failure& f) {
            if (std::getenv("MVB_DIAG"))
                std::cerr << "[sweepWire] rect fixed-binormal threw: "
                          << (f.GetMessageString() ? f.GetMessageString() : "(null)") << "\n";
        }
        return TopoDS_Shape();
    }
    // A TOROID's hole-threading spine has high torsion: MakePipeShell's corrected Frenet closes a
    // watertight, right-volume body but drifts the round section a percent or two off the spine at
    // the lead-adjacent crossing. The simple MakePipe with a FACE profile keeps the section centred
    // there, so for toroids we use it EXCLUSIVELY -- never the corrected-Frenet fallback, which would
    // silently ship an off-centre crossing. If MakePipe can't close (complex spread spine), return
    // null and let the caller drop to the exact per-run compound.
    if (preferSimplePipe) {
        try {
            TopoDS_Face prof = wireProfile(p0, t0, wireRadius, profileSegments);
            BRepOffsetAPI_MakePipe mp(spine, prof);
            mp.Build();
            if (mp.IsDone() && !mp.Shape().IsNull() && BRepCheck_Analyzer(mp.Shape()).IsValid())
                return mp.Shape();
        } catch (const Standard_Failure&) {
        }
        return TopoDS_Shape();
    }
    const bool sweepDiag = std::getenv("MVB_DIAG") != nullptr;
    static const char* kModeName[3] = {"Frenet", "correctedFrenet", "discrete"};
    // Attribute a "No geometry" throw to the SPINE or the PROFILE: both are adapted inside the
    // sweep, so the exception alone does not say which.
    if (sweepDiag) {
        try {
            TopoDS_Wire dp = wireProfileWire(p0, t0, wireRadius, profileSegments);
            int npe = 0, nbad = 0;
            for (TopExp_Explorer e(dp, TopAbs_EDGE); e.More(); e.Next()) {
                ++npe; double f, l;
                if (BRep_Tool::Curve(TopoDS::Edge(e.Current()), f, l).IsNull()) ++nbad;
            }
            std::cerr << "[sweepWire] profile: edges=" << npe << " without-curve=" << nbad
                      << " r=" << wireRadius << " segs=" << profileSegments << "\n";
        } catch (const Standard_Failure& e) {
            std::cerr << "[sweepWire] profile CONSTRUCTION threw (" << e.GetMessageString() << ")\n";
        }
        int nse = 0, nsbad = 0;
        for (TopExp_Explorer e(spine, TopAbs_EDGE); e.More(); e.Next()) {
            ++nse; double f, l;
            if (BRep_Tool::Curve(TopoDS::Edge(e.Current()), f, l).IsNull()) ++nsbad;
        }
        std::cerr << "[sweepWire] spine: edges=" << nse << " without-curve=" << nsbad << "\n";
    }
    // FIXED-BINORMAL frame first when requested (rect-column whole-path sweep): on a racetrack
    // spine every wrap turns about the column axis, so the binormal IS the axis -- an exact,
    // drift-free frame. The Frenet family reframes along the spine and was measured to displace
    // the section 18 um on a 56-prim spine: enough to close a 19.5 um inter-wrap clearance to
    // 1.5 um and make an otherwise sound body unmeshable. Rotation is immaterial for the round
    // section; only the POSITIONING stability matters. Falls through to the Frenet ladder if the
    // fixed frame cannot build (e.g. a spine edge tangent parallel to the axis).
    if (tryFixedBinormal) {
        try {
            TopoDS_Wire prof = wireProfileWire(p0, t0, wireRadius, profileSegments);
            BRepOffsetAPI_MakePipeShell ps(spine);
            if (kSweepTol3d > 0.0) ps.SetTolerance(kSweepTol3d, kSweepTol3d, 1e-2);
            ps.SetMode(axialAxis);
            ps.Add(prof);
            ps.Build();
            if (ps.IsDone() && ps.MakeSolid()) {
                TopoDS_Shape s = ps.Shape();
                if (BRepCheck_Analyzer(s).IsValid()) {
                    if (sweepDiag) std::cerr << "[sweepWire] fixed-binormal: ACCEPTED\n";
                    return s;
                }
            }
            if (sweepDiag) std::cerr << "[sweepWire] fixed-binormal: failed, falling to Frenet\n";
        } catch (const Standard_Failure& e) {
            if (sweepDiag) std::cerr << "[sweepWire] fixed-binormal threw ("
                                     << e.GetMessageString() << ")\n";
        }
    }
    for (int mode = 0; mode < 3; ++mode) {
        try {
            TopoDS_Wire prof = wireProfileWire(p0, t0, wireRadius, profileSegments);
            BRepOffsetAPI_MakePipeShell ps(spine);
            if (kSweepTol3d > 0.0) ps.SetTolerance(kSweepTol3d, kSweepTol3d, 1e-2);
            if (mode == 0) ps.SetMode(Standard_True);         // Frenet
            else if (mode == 1) ps.SetMode(Standard_False);   // corrected Frenet (torsion-stable)
            else ps.SetDiscreteMode();                        // discrete
            ps.Add(prof);
            ps.Build();
            // Report WHICH stage each mode failed at. "sweep NULL" alone cannot distinguish a
            // MakePipeShell that never built from one that built an invalid solid, and those need
            // opposite fixes (spine vs profile/framing).
            if (!ps.IsDone()) {
                if (sweepDiag) std::cerr << "[sweepWire] mode " << kModeName[mode]
                                         << ": Build() not done\n";
                continue;
            }
            if (!ps.MakeSolid()) {
                if (sweepDiag) std::cerr << "[sweepWire] mode " << kModeName[mode]
                                         << ": MakeSolid() failed (shell did not close)\n";
                continue;
            }
            TopoDS_Shape s = ps.Shape();
            if (BRepCheck_Analyzer(s).IsValid()) return s;
            if (sweepDiag) {
                int nsol = 0;
                for (TopExp_Explorer e(s, TopAbs_SOLID); e.More(); e.Next()) ++nsol;
                GProp_GProps gp; BRepGProp::VolumeProperties(s, gp);
                std::cerr << "[sweepWire] mode " << kModeName[mode] << ": built but INVALID"
                          << " (solids=" << nsol << " vol=" << gp.Mass()*1e9 << "mm3)\n";
            }
        } catch (const Standard_Failure& e) {
            if (sweepDiag) std::cerr << "[sweepWire] mode " << kModeName[mode] << ": threw ("
                                     << e.GetMessageString() << ")\n";
            continue;
        }
    }
    // LAST RESORT for a round section: the SIMPLE pipe. MakePipeShell can fail on a long spine of
    // many short, alternating curve types -- a litz round-column winding lands 65 edges that
    // alternate BSpline wraps with sub-millimetre Bezier fillets and Line leads, and all three
    // framing modes abort with "BRepAdaptor_Curve::No geometry" even though every spine edge has a
    // valid 3D curve and a sane parameter range (verified edge by edge). MakePipe is a different
    // algorithm with a different failure envelope, and a round section is rotation-invariant so its
    // framing is immaterial. Caller still volume-checks the result, so a wrong body cannot slip in.
    try {
        TopoDS_Face prof = wireProfile(p0, t0, wireRadius, profileSegments);
        BRepOffsetAPI_MakePipe mp(spine, prof);
        mp.Build();
        if (mp.IsDone() && !mp.Shape().IsNull() && BRepCheck_Analyzer(mp.Shape()).IsValid()) {
            if (sweepDiag) std::cerr << "[sweepWire] simple MakePipe succeeded\n";
            return mp.Shape();
        }
        if (sweepDiag) std::cerr << "[sweepWire] simple MakePipe: done=" << mp.IsDone()
                                 << " null=" << mp.Shape().IsNull() << "\n";
    } catch (const Standard_Failure& e) {
        if (sweepDiag) std::cerr << "[sweepWire] simple MakePipe threw (" << e.GetMessageString()
                                 << ")\n";
    }
    return TopoDS_Shape();
}

// Forward: the run's junction tangents decide the sweep's transition mode (see sweepRun).
static gp_Dir primFwdStart(const Primitive& p, double r);
static gp_Dir primFwdEnd(const Primitive& p, double r);

TopoDS_Shape sweepRun(const Primitive* const* prims, size_t count, double wireRadius,
                      int wirePolygonSegments) {
    if (count == 0) return TopoDS_Shape();
    try {
        BRepBuilderAPI_MakeWire wireMaker;
        for (size_t i = 0; i < count; ++i) {
            TopoDS_Edge e = primEdge(*prims[i], wireRadius);
            if (e.IsNull()) continue;
            wireMaker.Add(e);
            if (!wireMaker.IsDone()) return TopoDS_Shape();
        }
        if (!wireMaker.IsDone()) return TopoDS_Shape();
        TopoDS_Wire spine = wireMaker.Wire();

        auto firstPts = samplePrim(*prims[0], wireRadius);
        gp_Vec t0(firstPts[1].XYZ() - firstPts[0].XYZ());
        if (t0.Magnitude() < 1e-12) return TopoDS_Shape();
        // The profile ALWAYS follows the requested faceting (Alf: one setting, no exceptions;
        // segments <= 0 is the explicit exact-round request). The old hardcoded exact circle
        // sped the sweep up but produced PERIODIC surfaces gmsh refuses to mesh, and mixed
        // round turns into faceted conductors whenever different paths built them.
        TopoDS_Wire prof =
            wireProfileWire(firstPts.front(), gp_Dir(t0), wireRadius, wirePolygonSegments);

        BRepOffsetAPI_MakePipeShell ps(spine);
        // TRANSITION MODE. Corner machinery is only correct where there IS a corner, and both of
        // OCCT 7.9's corner modes crash on spines that have none:
        //   RoundCorner  -- trim machinery needs edges longer than its rounding radius and
        //                   segfaults on the short MKF link segments
        //                   (BRepFill_TrimShellCorner::Perform).
        //   RightCorner  -- at a TANGENT junction the two section edges it compares are parallel;
        //                   for non line/circle sections (our swept BSplines) OCCT's
        //                   Extrema_ExtCC::PrepareParallelResult appends to mySqDist WITHOUT
        //                   appending to mypoints, so NbExt() reports 1 while mypoints is empty
        //                   and Points(1) reads out of bounds -- a SIGSEGV in a release build,
        //                   through BRepFill_TrimShellCorner::ChooseSection. Reproduced on
        //                   06_llc_xfmr_eq4128: two consecutive 'over dragback' spirals meet at
        //                   9e-7 deg across the station plane.
        // A run whose every junction is tangent has no corner to trim, so sweep it with the plain
        // Transformed mode -- no TrimShellCorner, no crash, and the result is identical because
        // there is nothing to mitre. Runs with a genuine corner keep RightCorner.
        bool anyCorner = false;
        for (size_t i = 0; i + 1 < count && !anyCorner; ++i) {
            const gp_Dir a = primFwdEnd(*prims[i], wireRadius);
            const gp_Dir b = primFwdStart(*prims[i + 1], wireRadius);
            if (a.Angle(b) > 0.05) anyCorner = true;   // same tangency threshold as the mitre path
        }
        ps.SetTransitionMode(anyCorner ? BRepBuilderAPI_RightCorner
                                       : BRepBuilderAPI_Transformed);
        ps.Add(prof);
        ps.Build();
        if (!ps.IsDone()) return TopoDS_Shape();
        if (!ps.MakeSolid()) return TopoDS_Shape();
        TopoDS_Shape shape = ps.Shape();

        // FLARE GUARD. A multi-edge PipeShell can silently balloon a whole span into
        // radial spikes when a near-degenerate edge (e.g. an 88-degree cone from a
        // layer-transition wrap) sits in the chain — the RightCorner transition machinery
        // corrupts it without failing. Reject any solid whose TRUE extent (optimal
        // bounds, NOT the loose control-point hull) exceeds the spine's own sample box by
        // more than the wire radius plus a margin; the caller then sweeps this span
        // per-primitive (each edge swept alone does not flare) and fuses.
        Bnd_Box spineBox;
        for (size_t i = 0; i < count; ++i)
            for (const auto& q : samplePrim(*prims[i], wireRadius)) spineBox.Add(q);
        if (!spineBox.IsVoid()) {
            Bnd_Box solidBox;
            BRepBndLib::AddOptimal(shape, solidBox, Standard_False, Standard_False);
            double sx0, sy0, sz0, sx1, sy1, sz1;
            spineBox.Get(sx0, sy0, sz0, sx1, sy1, sz1);
            double bx0, by0, bz0, bx1, by1, bz1;
            solidBox.Get(bx0, by0, bz0, bx1, by1, bz1);
            double excess = std::max({sx0 - bx0, sy0 - by0, sz0 - bz0,
                                      bx1 - sx1, by1 - sy1, bz1 - sz1});
            if (excess > 2.0 * wireRadius) return TopoDS_Shape();
        }

        return shape;
    } catch (const Standard_Failure&) {
        return TopoDS_Shape();
    } catch (const std::exception&) {
        return TopoDS_Shape();
    }
}

// Build ONE rectangular-wire solid for a single primitive, cross-section width(radial) x
// height(axial). Corners are REVOLVED (an annular wedge about the corner axis -- turns a corner
// tighter than the wire's half-width cleanly, exactly where a swept flat section would collapse);
// straights are prisms; gentle -Z transitions are piped (corrected Frenet). Used for rect/oblong
// columns, whose racetrack corners defeat a single swept pipe.
TopoDS_Shape rectPrimSolid(const Primitive& pr, double w, double h, const gp_Dir& axialStart,
                           const gp_Dir& axialEnd, double extendA = 0.0, double extendB = 0.0,
                           bool round = false, double radius = 0.0,
                           int splitOverride = -1) {
    // Section profile at a point: an exact circle for round/litz wire, else the oriented rectangle.
    // NOTE: do NOT split this into arcs. The analytic Cylinder/Torus surfaces a closed circular
    // profile produces are exactly what gmsh meshes best here; splitting them was measured to make
    // things WORSE -- the copper-only mesh, which builds fine, started failing with "Invalid
    // boundary mesh (overlapping facets)", and the coated one merely swapped one periodic-surface
    // complaint for another. Only the swept BLEND needs a split profile (see rectPrimSolid's
    // BLEND branch), because a SWEPT closed profile makes a B-spline surface rather than an
    // analytic one.
    // MVB_SPLIT_PROFILE: build the round section from ARCS rather than one closed circle edge.
    // A closed circular profile extrudes/revolves into a PERIODIC surface (cylinder / torus with a
    // seam), which gmsh cannot mesh directly ("Impossible to mesh periodic surface"). The boolean
    // union normally hides this by trimming those surfaces -- which is precisely why the union
    // looked indispensable -- but it pays for it with seam artefacts (self-overlapping faces,
    // sliver sheets). Splitting the profile removes periodicity AT THE SOURCE, so the primitives
    // are meshable on their own and the fragile union becomes optional (see MVB_NO_FUSE).
    const int splitPieces = splitOverride >= 0 ? splitOverride
                          : std::getenv("MVB_SPLIT_PROFILE")
                            ? std::atoi(std::getenv("MVB_SPLIT_PROFILE")) : 0;
    auto profile = [&](const gp_Pnt& c, const gp_Dir& tan, const gp_Dir& ax) {
        if (!round) return rectProfileWire(c, tan, ax, w, h);
        if (splitPieces > 0) {
            // EQUAL-AREA POLYGON, not split arcs: split cylindrical panels are still patches
            // of a periodic basis surface, and after a fragment re-maps their u-range across
            // the 2 pi seam gmsh refuses them ('Impossible to mesh periodic surface',
            // measured on 04_forward's lead panel and 02_flyback's descent leg). Twelve
            // PLANAR facets cannot be periodic, and the radius correction keeps the copper
            // area exact (r' = r / sqrt(N sin(2pi/N) / (2pi))).
            const int N = std::max(12, 2 * splitPieces);
            const double rEq = radius / std::sqrt(N * std::sin(kTwoPi / N) / kTwoPi);
            gp_Ax2 frame(c, tan);
            BRepBuilderAPI_MakePolygon poly;
            for (int k = 0; k < N; ++k) {
                const double a2 = kTwoPi * k / N;
                gp_XYZ pt = c.XYZ() + frame.XDirection().XYZ() * (rEq * std::cos(a2)) +
                            frame.YDirection().XYZ() * (rEq * std::sin(a2));
                poly.Add(gp_Pnt(pt));
            }
            poly.Close();
            return poly.Wire();
        }
        return wireProfileWire(c, tan, radius, 0);
    };
    try {
        if (pr.kind == Primitive::SEG) {
            gp_XYZ dir = pr.seg.b.XYZ() - pr.seg.a.XYZ();
            if (dir.Modulus() < 1e-12) return TopoDS_Shape();
            gp_Dir t(dir);
            // extendA/extendB grow the prism past a SHARP, elbow-less junction (a terminal lead's
            // own 90-degree bend) so the two lead prisms overlap and the fuse closes the mitre.
            gp_XYZ u = dir / dir.Modulus();
            gp_Pnt sa(pr.seg.a.XYZ() - u * extendA), sb(pr.seg.b.XYZ() + u * extendB);
            gp_XYZ ext = sb.XYZ() - sa.XYZ();
            // A round section is rotation-invariant, and a rect straight whose ends want the SAME
            // orientation, are plain prisms (SEG -> cylinder for round). Only a rect straight whose
            // ends differ (a toroidal chord advancing azimuthally to the next turn) is TWISTED by
            // lofting between the two oriented rectangles, so it meets its elbows at their own angle.
            if (round || axialStart.Angle(axialEnd) < 0.01) {
                TopoDS_Face face = BRepBuilderAPI_MakeFace(profile(sa, t, axialStart)).Face();
                return BRepPrimAPI_MakePrism(face, gp_Vec(ext)).Shape();
            }
            BRepOffsetAPI_ThruSections loft(Standard_True);  // solid
            loft.AddWire(rectProfileWire(sa, t, axialStart, w, h));
            loft.AddWire(rectProfileWire(sb, t, axialEnd, w, h));
            loft.Build();
            if (!loft.IsDone()) return TopoDS_Shape();
            return loft.Shape();
        }
        const gp_Dir& axial = axialStart;
        if (pr.kind == Primitive::ARC3) {
            double r = pr.arc.v0.Modulus();
            if (r < 1e-12 || std::abs(pr.arc.sweep) < 1e-12) return TopoDS_Shape();
            gp_Pnt start(pr.arc.c.XYZ() + pr.arc.v0);
            gp_XYZ tan = pr.arc.axis.Crossed(pr.arc.v0);  // start tangent (d/dt of the arc at t=0)
            if (tan.Modulus() < 1e-12) return TopoDS_Shape();
            // Revolve the section about the corner's bend axis: for round that is a torus segment,
            // for rect a clean annular wedge in the (radial x axial) bend plane (its section axis IS
            // the arc axis -- matching make_toroidal_quarter_swept_rectangle).
            TopoDS_Face face =
                BRepBuilderAPI_MakeFace(profile(start, gp_Dir(tan), gp_Dir(pr.arc.axis))).Face();
            gp_Dir revDir = pr.arc.sweep > 0 ? gp_Dir(pr.arc.axis) : gp_Dir(pr.arc.axis).Reversed();
            return BRepPrimAPI_MakeRevol(face, gp_Ax1(pr.arc.c, revDir), std::abs(pr.arc.sweep))
                .Shape();
        }
        if (!round && pr.kind == Primitive::SPIRAL && !pr.spiral.blend) {
            // RISING RECT CORNER (Alf, 18_stacked): the pipe-shell swept the section from a
            // SAMPLED chord tangent, so the solid began its rotation before the straight's
            // real end and its lateral surfaces matched neither neighbour. Loft EXACT
            // analytic sections instead: profiles on the helical corner at the TRUE tangent,
            // the first/last EXACTLY on the neighbouring prisms' end planes.
            const Spiral& sp = pr.spiral;
            const double azSpan = sp.az1 - sp.az0;
            const int K = std::max(5, (int)std::ceil(std::abs(azSpan) / (kPi / 12.0)) + 1);
            BRepOffsetAPI_ThruSections loft(Standard_True, Standard_False);
            for (int k2 = 0; k2 < K; ++k2) {
                const double t = (double)k2 / (K - 1);
                const double az = sp.az0 + azSpan * t;
                const double rr = sp.r0 + (sp.r1 - sp.r0) * t;
                const double yy = sp.y0 + (sp.y1 - sp.y0) * t;
                const gp_Pnt c(sp.cx + rr * std::cos(az), yy, sp.cz - rr * std::sin(az));
                gp_XYZ tv(-rr * std::sin(az) * azSpan + (sp.r1 - sp.r0) * std::cos(az),
                          sp.y1 - sp.y0,
                          -rr * std::cos(az) * azSpan - (sp.r1 - sp.r0) * std::sin(az));
                if (tv.Modulus() < 1e-15) return TopoDS_Shape();
                loft.AddWire(rectProfileWire(c, gp_Dir(tv), axial, w, h));
            }
            loft.Build();
            if (loft.IsDone() && !loft.Shape().IsNull()) return loft.Shape();
            return TopoDS_Shape();
        }
        // BLEND / SPIRAL (the gentle axial transition that carries the wire from one turn to the
        // next): pipe the section over the analytic edge. This is the ONLY swept primitive in an
        // otherwise analytic rect-column turn (every straight is a prism, every corner a revolve),
        // and building it with ONE unchecked sweep mode is what made the FEM one-body export
        // unmeshable: on e138 the corrected-Frenet pipe closed with a SELF-INTERSECTING seam wire
        // on each of the 6 transitions, so the union of the 56 primitives came out at exactly the
        // right copper volume but BRepCheck-invalid, no ShapeFix/UnifySameDomain repair could clear
        // it, and gmsh rejected it ("the 1D mesh seems not to be forming a closed loop"). Try the
        // sweep modes in order and keep the first BRepCheck accepts -- the same validity-checked
        // ladder sweepWire() already applies to whole-run pipes. A ROUND section is rotation-
        // invariant, so the simple MakePipe (exact surfaces where it can build them) goes first for
        // it; a rect section has exactly one correct frame and must not use it. (ABT #332)
        double sampleR = round ? radius : std::min(w, h) / 2.0;
        TopoDS_Edge e = primEdge(pr, sampleR);
        if (e.IsNull()) return TopoDS_Shape();
        auto pts = samplePrim(pr, sampleR);
        if (pts.size() < 2) return TopoDS_Shape();
        gp_Dir t0 = (pts[1].XYZ() - pts[0].XYZ()).Modulus() > 1e-12
                        ? gp_Dir(pts[1].XYZ() - pts[0].XYZ())
                        : gp_Dir(1, 0, 0);
        const TopoDS_Wire spine = BRepBuilderAPI_MakeWire(e).Wire();
        // Round section: sweep a SPLIT (two half-arc) profile so the pipe has no periodic seam
        // face -- see wireProfileWireSplit. A rect profile is already four separate edges.
        // 6 arcs, not 2: with only two half-arcs the swept patches are still parameterised over a
        // half period and gmsh reports "Impossible to mesh periodic surface"; 6 puts every patch
        // well inside one period. (Measured on e138: 2 and 3 pieces -> periodic-surface failures,
        // 4 -> the union collapses, 6 -> meshes.)
        const int blendPieces = std::getenv("MVB_BLEND_PIECES")
                                ? std::atoi(std::getenv("MVB_BLEND_PIECES")) : 6;
        const TopoDS_Wire blendProf =
            round ? wireProfileWireSplit(pts.front(), t0, radius, blendPieces)
                  : profile(pts.front(), t0, axial);
        TopoDS_Shape fallback;   // first body that built at all, valid or not (drawing tolerates it)
        for (int mode = 0; mode < 3; ++mode) {
            try {
                BRepOffsetAPI_MakePipeShell ps(spine);
                // NB: do NOT apply MVB_SWEEP_TOL here. Tightening the BLEND sweep merely trades
                // one failure for another on e138: the self-overlapping -Z ramp face becomes a
                // PERIODIC surface gmsh cannot mesh at all ("Impossible to mesh periodic surface").
                // The split profile already prevents periodicity at the default tolerance.
                if (mode == 0) ps.SetMode(Standard_False);       // corrected Frenet (torsion-stable)
                else if (mode == 1) ps.SetMode(Standard_True);   // Frenet
                else ps.SetDiscreteMode();                       // discrete
                ps.Add(blendProf);
                ps.Build();
                if (!ps.IsDone() || !ps.MakeSolid()) continue;
                TopoDS_Shape s = ps.Shape();
                if (BRepCheck_Analyzer(s).IsValid()) return s;
                if (fallback.IsNull()) fallback = s;
            } catch (const Standard_Failure&) {
                continue;
            }
        }
        return fallback;
    } catch (const Standard_Failure&) {
        return TopoDS_Shape();
    }
}

// Rect wire on a rect/oblong column: build every primitive as its own rect solid. When the turns
// are SEPARATED (the conducting/copper footprint), consecutive primitives touch only along the
// path, so fusing yields one clean spiral solid -- ideal for FEM. When the turns TOUCH (the outer/
// insulation footprint), fusing would merge the stacked flat faces into an electrically-wrong solid
// brick (current would short straight through instead of spiralling), so we detect the collapse (a
// merged brick loses almost all its faces) and keep the per-primitive COMPOUND, which shows every
// turn as its own solid and honours every MKF position exactly.
bool hasDegenerateSheetFace(const TopoDS_Shape& shape, double wireRadius);   // defined below

TopoDS_Shape emitRectColumn(const ConductorPath& path) {
    if (std::getenv("MVB_DIAG"))
        std::cerr << "[emitRectColumn] ENTER '" << path.name << "' prims=" << path.prims.size()
                  << " femReady=" << path.femReady << " roundProfile=" << path.roundProfile << "\n";
    // The straight/blend section's "axial" axis: the column axis Y for concentric columns, or the
    // toroid's AZIMUTHAL tangent (perpendicular to the poloidal plane) at the primitive's position.
    // (Corners take their axis from the arc itself, inside rectPrimSolid.)
    const bool round = path.roundProfile;
    const double radius = path.wireRadius;
    auto axialFor = [&](const gp_Pnt& c) -> gp_Dir {
        if (path.toroidal) {
            gp_XYZ az(-c.Z(), 0.0, c.X());  // tangent to the big ring (hole axis Y, ring in XZ)
            if (az.Modulus() > 1e-9) return gp_Dir(az);
        }
        return gp_Dir(0, 1, 0);
    };
    // A corner (ARC3) is revolved about its OWN bend axis, so its section's axial is that arc axis
    // for the whole elbow. A straight must therefore meet each neighbouring elbow on THAT elbow's
    // axis (not the raw azimuthal, which differs by a degree or so and steps the faces). So a
    // straight's end axial = the adjacent elbow's arc axis, flipped to the azimuthal hemisphere so
    // the twist loft goes the short way; a straight next to another straight (a lead) uses azimuthal.
    auto alignHemisphere = [](const gp_XYZ& axis, const gp_Dir& ref) {
        gp_Dir a(axis);
        return a.Dot(ref) >= 0.0 ? a : gp_Dir(a.Reversed());
    };
    // Tangent at a primitive end (for detecting a lead's own sharp bend).
    double sampR = std::min(path.wireWidth, path.wireHeight) / 2.0;
    auto entryDir = [&](const Primitive& p) {
        auto pts = samplePrim(p, sampR);
        return pts.size() >= 2 && (pts[1].XYZ() - pts[0].XYZ()).Modulus() > 1e-12
                   ? gp_Dir(pts[1].XYZ() - pts[0].XYZ()) : gp_Dir(1, 0, 0);
    };
    auto exitDir = [&](const Primitive& p) {
        auto pts = samplePrim(p, sampR);
        size_t n = pts.size();
        return n >= 2 && (pts[n - 1].XYZ() - pts[n - 2].XYZ()).Modulus() > 1e-12
                   ? gp_Dir(pts[n - 1].XYZ() - pts[n - 2].XYZ()) : gp_Dir(1, 0, 0);
    };
    auto segLen = [](const Primitive& p) {
        return p.kind == Primitive::SEG ? p.seg.a.Distance(p.seg.b) : 0.0;
    };
    // A terminal lead's own bend is two isLead SEGs with no elbow between them. Build a REAL rounded
    // elbow there -- the same revolved-rectangle corner the winding turns use -- by trimming both
    // lead prisms back by R*tan(angle/2) and inserting an ARC3 of radius R. Only isLead<->isLead
    // corners qualify, so every winding corner is untouched. rotationAxis[i] carries the elbow's
    // bend axis so the trimmed lead prism aligns its section to it (flush, as the turns do).
    std::vector<double> trimStart(path.prims.size(), 0.0), trimEnd(path.prims.size(), 0.0);
    std::vector<gp_Dir> axisAfter(path.prims.size(), gp_Dir(0, 1, 0));
    std::vector<gp_Dir> axisBefore(path.prims.size(), gp_Dir(0, 1, 0));
    std::vector<bool> hasAxisAfter(path.prims.size(), false), hasAxisBefore(path.prims.size(), false);
    std::vector<TopoDS_Shape> leadElbows;
    for (size_t i = 0; i + 1 < path.prims.size(); ++i) {
        const Primitive& A = path.prims[i];
        const Primitive& B = path.prims[i + 1];
        if (!(A.isLead && B.isLead)) continue;
        gp_Dir dA = exitDir(A), dB = entryDir(B);
        double ang = dA.Angle(dB);
        if (ang < 0.05) continue;
        gp_XYZ bendAxis = dA.XYZ().Crossed(dB.XYZ());
        if (bendAxis.Modulus() < 1e-12) continue;
        bendAxis.Normalize();
        double R = std::min({std::min(path.wireWidth, path.wireHeight), 0.4 * segLen(A),
                             0.4 * segLen(B)});
        if (R < 1e-9) continue;
        double trim = R * std::tan(ang / 2.0);
        gp_Pnt P = primEndpoints(A).second;
        gp_XYZ Ap = P.XYZ() - dA.XYZ() * trim;
        gp_XYZ nA = dB.XYZ() - dA.XYZ() * dB.XYZ().Dot(dA.XYZ());  // dB perp to dA -> toward centre
        if (nA.Modulus() < 1e-12) continue;
        nA.Normalize();
        Primitive elbow;
        elbow.kind = Primitive::ARC3;
        elbow.arc.c = gp_Pnt(Ap + nA * R);
        elbow.arc.axis = bendAxis;
        elbow.arc.v0 = Ap - (Ap + nA * R);  // centre -> Ap, magnitude R
        elbow.arc.sweep = ang;
        gp_XYZ Bp = P.XYZ() + dB.XYZ() * trim;
        if (primEndpoints(elbow).second.Distance(gp_Pnt(Bp)) > 1e-6)
            elbow.arc.axis = bendAxis * -1.0;
        TopoDS_Shape es = rectPrimSolid(elbow, path.wireWidth, path.wireHeight, gp_Dir(bendAxis),
                                        gp_Dir(bendAxis), 0.0, 0.0, round, radius);
        if (es.IsNull()) continue;
        leadElbows.push_back(es);
        trimEnd[i] = trim;
        trimStart[i + 1] = trim;
        axisAfter[i] = gp_Dir(elbow.arc.axis);
        hasAxisAfter[i] = true;
        axisBefore[i + 1] = gp_Dir(elbow.arc.axis);
        hasAxisBefore[i + 1] = true;
    }
    std::vector<TopoDS_Shape> solids;
    solids.reserve(path.prims.size());
    // Everything needed to REBUILD each prim with a split profile if the union fails and the
    // butt-joined chain must be meshed directly (closed-circle profiles are best for the FUSE
    // -- analytic quadrics -- but their full-circle seam edges are unmeshable in a fragment:
    // "distance 6.28319 [2 pi] between first and last node in 1D mesh of surface N").
    struct SplitArg { const Primitive* pr; gp_Dir axialA, axialB; double extA, extB; };
    std::vector<SplitArg> splitArgs;
    int compoundFaces = 0;
    for (size_t i = 0; i < path.prims.size(); ++i) {
        const Primitive& pr = path.prims[i];
        auto [pa, pb] = primEndpoints(pr);
        gp_Dir axialA = axialFor(pa), axialB = axialFor(pb);
        if (!round && pr.kind == Primitive::SEG) {  // a round section needs no orientation
            if (i > 0 && path.prims[i - 1].kind == Primitive::ARC3)
                axialA = alignHemisphere(path.prims[i - 1].arc.axis, axialA);
            if (i + 1 < path.prims.size() && path.prims[i + 1].kind == Primitive::ARC3)
                axialB = alignHemisphere(path.prims[i + 1].arc.axis, axialB);
            // Align the trimmed lead prism to the inserted elbow's axis (flush, as the turns).
            if (hasAxisBefore[i]) axialA = alignHemisphere(axisBefore[i].XYZ(), axialA);
            if (hasAxisAfter[i]) axialB = alignHemisphere(axisAfter[i].XYZ(), axialB);
        }
        // Negative extend == trim the prism back so the inserted lead elbow meets it flush.
        // NB do NOT grow the junctions to force a transversal boolean: measured on e138, any
        // overlap (0.05-0.3 x wireRadius) fragments the union instead of cleaning it, because the
        // grown straight pokes out of the neighbouring corner/blend.
        double extA = -trimStart[i], extB = -trimEnd[i];
        TopoDS_Shape s = rectPrimSolid(pr, path.wireWidth, path.wireHeight, axialA, axialB, extA,
                                       extB, round, radius, /*splitOverride=*/-1);
        if (s.IsNull()) {
            // NO silent drops (Alf): a failed primitive solid means the emitted copper is
            // missing from the conductor -- 17_cllc shipped its secondary WITHOUT the exit
            // lead through this very 'continue' and the watertight battery could not see it
            // (every present solid was valid; the absent one had no witness).
            throw std::runtime_error(
                "ConductorBuilder: rectPrimSolid failed for '" + pr.label + "' of '" +
                path.name + "' (trims " + std::to_string(extA) + "/" + std::to_string(extB) +
                " m). Refusing to emit a conductor with missing copper.");
        }
        if (std::getenv("MVB_DIAG") && !BRepCheck_Analyzer(s).IsValid())
            std::cerr << "[rect-column] INVALID prim solid [" << pr.label << "]\n";
        splitArgs.push_back({&pr, axialA, axialB, extA, extB});
        solids.push_back(s);
        for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next()) ++compoundFaces;
    }
    for (const auto& es : leadElbows) {
        solids.push_back(es);
        for (TopExp_Explorer e(es, TopAbs_FACE); e.More(); e.Next()) ++compoundFaces;
    }
    if (solids.empty()) {
        throw std::runtime_error(
            "ConductorBuilder: rectangular-wire rect/oblong column produced no solids for '" +
            path.name + "'");
    }
    const bool diag = std::getenv("MVB_DIAG") != nullptr;
    if (diag) {
        int nseg = 0, narc = 0, nspiral = 0, nblend = 0;
        for (const auto& pr : path.prims)
            switch (pr.kind) {
                case Primitive::SEG:    ++nseg; break;
                case Primitive::ARC3:   ++narc; break;
                case Primitive::SPIRAL: ++nspiral; break;
                case Primitive::BLEND:  ++nblend; break;
            }
        std::cerr << "[rect-column] " << path.prims.size() << " prims: SEG=" << nseg
                  << " ARC3=" << narc << " SPIRAL=" << nspiral << " BLEND=" << nblend
                  << " elbows=" << leadElbows.size() << "\n";
        // Which INPUT solids are already invalid: a union cannot be blamed for a bad argument.
        for (size_t i = 0; i < solids.size(); ++i) {
            if (BRepCheck_Analyzer(solids[i]).IsValid()) continue;
            const char* kind = i < path.prims.size()
                                   ? (path.prims[i].kind == Primitive::SEG    ? "SEG"
                                      : path.prims[i].kind == Primitive::ARC3 ? "ARC3"
                                      : path.prims[i].kind == Primitive::BLEND ? "BLEND" : "SPIRAL")
                                   : "elbow";
            std::cerr << "[rect-column]   input prim " << i << " (" << kind << ") INVALID\n";
        }
    }
    // EXPECTED COPPER: profile area x centreline length. The per-primitive solids overlap only at
    // their junctions, so a correct union lands within a few percent of this. This is the ONLY
    // check that catches a union which silently DROPPED bodies -- see acceptFused below.
    double spineLen = 0.0;
    for (const auto& pr : path.prims) {
        auto pts = samplePrim(pr, path.wireRadius);
        for (size_t k = 1; k < pts.size(); ++k) spineLen += pts[k].Distance(pts[k - 1]);
    }
    const double sectionArea =
        round ? kPi * path.wireRadius * path.wireRadius : path.wireWidth * path.wireHeight;
    const double expectedVol = sectionArea * spineLen;
    auto volumeOf = [](const TopoDS_Shape& s) {
        GProp_GProps g; BRepGProp::VolumeProperties(s, g); return g.Mass();
    };
    // What BRepCheck actually objects to, by sub-shape type and status code -- "invalid" alone
    // says nothing about whether the body is meshable.
    auto checkReport = [](const TopoDS_Shape& s) {
        BRepCheck_Analyzer an(s);
        std::map<std::string, int> counts;
        const char* kinds[] = {"VERTEX", "EDGE", "WIRE", "FACE", "SHELL", "SOLID"};
        const TopAbs_ShapeEnum types[] = {TopAbs_VERTEX, TopAbs_EDGE,  TopAbs_WIRE,
                                          TopAbs_FACE,   TopAbs_SHELL, TopAbs_SOLID};
        for (int k = 0; k < 6; ++k) {
            for (TopExp_Explorer e(s, types[k]); e.More(); e.Next()) {
                if (an.IsValid(e.Current())) continue;
                Handle(BRepCheck_Result) res = an.Result(e.Current());
                if (res.IsNull()) continue;
                for (const BRepCheck_Status& st : res->Status())
                    if (st != BRepCheck_NoError)
                        counts[std::string(kinds[k]) + ":" + std::to_string(static_cast<int>(st))]++;
            }
        }
        std::string out;
        for (const auto& kv : counts) out += " " + kv.first + "x" + std::to_string(kv.second);
        // Where the bad faces are, so a defect can be tied back to the primitive that made it.
        for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next()) {
            if (an.IsValid(e.Current())) continue;
            Bnd_Box bb; BRepBndLib::Add(e.Current(), bb);
            if (bb.IsVoid()) continue;
            double x0, y0, z0, x1, y1, z1; bb.Get(x0, y0, z0, x1, y1, z1);
            char buf[160];
            std::snprintf(buf, sizeof(buf), " @(%.2f,%.2f,%.2f)mm", 0.5 * (x0 + x1) * 1e3,
                          0.5 * (y0 + y1) * 1e3, 0.5 * (z0 + z1) * 1e3);
            out += buf;
        }
        return out.empty() ? std::string(" (none)") : out;
    };
    // Validate + finish a candidate union. Returns a null shape if the candidate is not a single,
    // valid, FULL-VOLUME conductor.
    //
    // THE VOLUME CHECK IS THE LOAD-BEARING ONE. BRepAlgoAPI_Fuse can report IsDone() and hand back
    // a perfectly VALID single solid that is just the first argument, having silently discarded
    // every tool: on the e138 spiral (56 primitives, round wire) the union returned 0.215 mm3 of an
    // expected 13.8 mm3 -- the lead stub alone. Nothing else catches that. nsol==1 is satisfied,
    // BRepCheck is satisfied, and the old face-count heuristic was short-circuited for round wire
    // ("a round section can't brick"), which is true but irrelevant -- it was the only guard, so a
    // 98%-of-the-copper loss sailed through into the mesher and the FEM solver, where it looked
    // like a solver bug (no azimuthal current, |I_net| ~ 3e-7 A). Never accept a union on topology
    // alone; measure the copper. (ABT #332)
    auto acceptFused = [&](TopoDS_Shape fused, const char* how, bool* collapsed = nullptr) -> TopoDS_Shape {
        if (fused.IsNull()) return TopoDS_Shape();
        int nsol = 0, fusedFaces = 0;
        for (TopExp_Explorer e(fused, TopAbs_SOLID); e.More(); e.Next()) ++nsol;
        for (TopExp_Explorer e(fused, TopAbs_FACE); e.More(); e.Next()) ++fusedFaces;
        const double v = volumeOf(fused);
        // A true rect spiral keeps most of its faces (each turn's flats survive); a brick-merge
        // (stacked flats fused away) keeps only a small fraction -> reject as a short. A ROUND
        // section can't brick -- turns touch only along a line -- so the face test does not apply
        // to it. Volume, below, is checked for BOTH.
        const bool notCollapsed = round || fusedFaces > compoundFaces / 2;
        if (collapsed) *collapsed = !notCollapsed;
        const bool okVol = expectedVol > 0.0 && v > 0.80 * expectedVol && v < 1.25 * expectedVol;
        // SLIVER GUARD. The union can weld the primitives into a solid of exactly the right copper
        // that nevertheless carries a near-degenerate SHEET face -- a wide, micron-thin patch left
        // where two co-located primitive faces did not quite merge. It is invisible to BRepCheck
        // and to the volume test, but a mesher cannot resolve it: gmsh reports "Invalid boundary
        // mesh (overlapping facets) on surface N surface N" because the sheet's two sides collide.
        // Measured on e138 after an OCCT rebuild: face 4.3 x 1.3 mm across but only 33 um thick,
        // which took the reference design from meshing in 5 s to not meshing at all. The same
        // detector already guards fuseAllSolids and the single-body sweep; the rect-column union
        // was the one path without it, so a slivered result was accepted and shipped downstream.
        // Rejecting here lets the PAIRWISE strategy (or the compound) be tried instead.
        const bool okSliver = !hasDegenerateSheetFace(fused, path.wireRadius);
        // SELF-INTERSECTION CHECK -- the missing acceptance criterion. BRepCheck_Analyzer tests
        // topological validity (closed shells, orientation, tolerances) and says nothing about
        // whether the body's own faces PASS THROUGH each other. A union of ~56 nearly-tangent
        // primitives can be BRepCheck-valid, carry exactly the right copper volume, and still
        // self-intersect; the first thing that notices is the mesher, reporting "Invalid boundary
        // mesh (overlapping facets) on surface N surface N" -- the same surface against itself --
        // which no element size can repair. Because the defect lives in OCCT's boolean
        // micro-behaviour it can appear or vanish across an OCCT rebuild: exactly the silent,
        // environment-dependent breakage that must never reach a solver. BOPAlgo's self-interference
        // test answers the question directly, so a self-intersecting union is treated as a FAILED
        // union, letting the pairwise strategy (then the compound) be tried instead.
        // REPORTED, NOT GATED. Making this an acceptance criterion was tried and rejected: it
        // turns away unions that OCCT calls self-intersecting but that build correct, meshable
        // geometry today (it broke the rectangular-column E-core zigzag regression). It is kept as
        // a diagnostic because it is the only direct read on the defect class that produces
        // "overlapping facets on surface N surface N" downstream -- and note the e138 failure does
        // NOT trip it, which is itself the useful finding: that body is geometrically sound and the
        // fault is in DISCRETISING a seam remnant far thinner than the element size.
        bool okSelfInt = true;
        try {
            BRepAlgoAPI_Check chk(fused, /*bRunParallel=*/Standard_False,
                                  /*bTestSelfInt=*/Standard_True);
            okSelfInt = chk.IsValid();
        } catch (const Standard_Failure&) {
            okSelfInt = false;   // the checker itself failing is not a pass
        }
        bool valid = BRepCheck_Analyzer(fused).IsValid();
        // The fuse of the analytic per-primitive solids routinely lands ONE solid that BRepCheck
        // rejects over small tolerance/seam defects at the junctions. That is a REPAIRABLE shape,
        // not a failed union -- and silently dropping to the compound here is what left the
        // "continuous" winding as N disconnected per-turn bodies, which gmsh then cannot mesh
        // ("1D mesh seems not to be forming a closed loop", ABT #332).
        if (nsol == 1 && notCollapsed && okVol && okSliver && !valid) {
            // REPAIR LADDER. The union of the analytic per-primitive solids lands ONE solid of the
            // right copper volume that BRepCheck rejects over junction seams -- on e138, six
            // SelfIntersectingWire / UnorientableShape faces, one per turn. That is NOT cosmetic:
            // gmsh's 1D stage chokes on exactly those wires ("the 1D mesh seems not to be forming a
            // closed loop"), so the body must be repaired, not waved through. Plain ShapeFix_Shape
            // does not touch self-intersections -- its wire tool needs FixSelfIntersectionMode
            // switched on explicitly -- and UnifySameDomain can dissolve the seam outright by
            // merging the co-domain cylinder/torus faces it sits between. Try both, and combined.
            struct Repair { const char* name; int mode; };
            static const Repair kLadder[] = {
                {"unify", 0},
                {"shapefix", 1},
                {"shapefix+selfint", 2},
                {"unify+shapefix+selfint", 3},
            };
            for (const Repair& r : kLadder) {
                TopoDS_Shape cand = fused;
                try {
                    if (r.mode == 0 || r.mode == 3) {
                        ShapeUpgrade_UnifySameDomain u(cand, Standard_True, Standard_True,
                                                       Standard_True);
                        u.Build();
                        cand = u.Shape();
                    }
                    if (r.mode >= 1) {
                        ShapeFix_Shape fixer(cand);
                        fixer.SetPrecision(1e-7);
                        fixer.SetMaxTolerance(1e-5);
                        if (r.mode >= 2) {
                            // Reach the wire tool and enable self-intersection repair (off by
                            // default) plus its dependencies.
                            Handle(ShapeFix_Wire) fw =
                                fixer.FixSolidTool()->FixShellTool()->FixFaceTool()->FixWireTool();
                            if (!fw.IsNull()) {
                                fw->FixSelfIntersectionMode() = 1;
                                fw->ModifyTopologyMode() = Standard_True;
                                fw->FixIntersectingEdgesMode() = 1;
                                fw->FixSelfIntersectingEdgeMode() = 1;
                            }
                        }
                        fixer.Perform();
                        cand = fixer.Shape();
                    }
                } catch (const Standard_Failure&) {
                    continue;
                }
                if (cand.IsNull()) continue;
                int rsol = 0;
                for (TopExp_Explorer e(cand, TopAbs_SOLID); e.More(); e.Next()) ++rsol;
                const double rv = volumeOf(cand);
                const bool rok = rsol == 1 && BRepCheck_Analyzer(cand).IsValid() &&
                                 rv > 0.80 * expectedVol && rv < 1.25 * expectedVol;
                if (diag) std::cerr << "[rect-column] repair '" << r.name << "' (" << how
                                    << "): nsol=" << rsol << " v=" << rv * 1e9 << "mm3 valid="
                                    << (rok ? 1 : 0)
                                    << (rok ? std::string() : "  BRepCheck:" + checkReport(cand))
                                    << "\n";
                if (rok) { fused = cand; valid = true; break; }
            }
        }
        if (nsol == 1 && notCollapsed && okVol && okSliver && valid) {
            // Merge co-domain faces, with a LINEAR TOLERANCE wide enough to absorb the seam
            // remnants the union leaves behind. Those remnants are the real meshing hazard: the
            // union can leave a strip only tens of microns wide where two primitives' surfaces
            // nearly coincide (measured on e138: a 4.3 x 1.3 mm face just 33 um across). It is
            // sound geometry -- OCCT's own self-intersection test passes -- but when the element
            // size is an order of magnitude coarser than the strip, the 2D mesher lays triangles
            // that fold over each other and gmsh rejects the surface against ITSELF ("overlapping
            // facets on surface N surface N"). At the default (Precision::Confusion, ~1e-7) such
            // faces are far too far apart to merge; a tolerance tied to the WIRE (a small fraction
            // of its radius) merges the remnant while staying far below any real feature.
            // MVB_UNIFY_TOL overrides (metres).
            ShapeUpgrade_UnifySameDomain unify(fused, Standard_True, Standard_True, Standard_True);
            // NB: do NOT widen UnifySameDomain's linear tolerance to try to absorb seam remnants.
            // Measured: at 0.25*wireRadius it merged nothing on e138 (face count unchanged at
            // 92/198, still unmeshable) while breaking the rectangular-column E-core zigzag
            // regression -- a tolerance that large lets legitimately distinct faces merge.
            unify.Build();
            TopoDS_Shape merged = unify.Shape();
            // UnifySameDomain merges co-domain faces; it must not change the copper.
            if (BRepCheck_Analyzer(merged).IsValid() && volumeOf(merged) > 0.80 * expectedVol) {
                if (diag) std::cerr << "[rect-column] FUSED spiral solid via " << how << " ("
                                    << fusedFaces << "/" << compoundFaces << " faces, v="
                                    << v * 1e9 << "mm3 exp=" << expectedVol * 1e9 << "mm3) from "
                                    << solids.size() << " prims\n";
                return merged;
            }
            if (diag) std::cerr << "[rect-column] UnifySameDomain invalid/lossy; returning the "
                                << "un-unified fused solid (still ONE connected body, " << how
                                << ")\n";
            return fused;
        }
        if (diag) std::cerr << "[rect-column] fuse REJECTED (" << how << "): nsol=" << nsol
                            << " notCollapsed=" << notCollapsed << " okVol=" << okVol
                            << " (v=" << v * 1e9 << "mm3 exp=" << expectedVol * 1e9 << "mm3)"
                            << " okSliver=" << okSliver << " okSelfInt=" << okSelfInt
                            << " valid=" << valid << " round=" << round << " faces=" << fusedFaces
                            << "/" << compoundFaces
                            << (valid ? std::string() : "  BRepCheck:" + checkReport(fused)) << "\n";
        return TopoDS_Shape();
    };
    // The per-primitive solids OVERLAP slightly at each junction (consecutive sections are
    // co-located but not bit-coincident), so a real boolean intersection is needed to weld them
    // -- GlueShift, which skips coincident-face splitting, leaves them unwelded. A small fuzzy
    // value absorbs the sub-micron section misalignment so the union stays clean. It must stay
    // SMALL: at 3e-6 the e138 union fragments into 17 solids, and no fuzzy value at all makes the
    // one-shot union weld the 56-primitive round spiral (measured 0 / 1e-9 / 1e-8 / 1e-7 -> one
    // INVALID solid; 3e-7 / 1e-6 -> the collapse described above) -- the fix for that is the
    // pairwise strategy below, not a bigger tolerance.
    const double fuzzy = std::getenv("MVB_FUSE_FUZZY")
                         ? std::atof(std::getenv("MVB_FUSE_FUZZY")) : 1e-7;
    // MVB_NO_FUSE: skip the boolean union entirely and return the butt-joined CHAIN of analytic
    // primitives. The union is the single most fragile step in the whole 3D path -- welding 56
    // nearly-tangent bodies routinely yields a solid that is topologically valid and volumetrically
    // correct yet carries seam artefacts no mesher can resolve (self-overlapping faces, sliver
    // sheets, periodic surfaces). The chain avoids it: consecutive primitives already share
    // co-located end faces, and gmsh's FRAGMENT makes them conforming, so the meshed conductor is
    // ONE connected region -- which is what FEM actually requires. A single BREP solid was never
    // the requirement; it was an assumption.
    if (std::getenv("MVB_NO_FUSE")) {
        if (diag) std::cerr << "[rect-column] MVB_NO_FUSE: skipping the boolean union, returning "
                            << solids.size() << " butt-joined primitives for fragment welding\n";
    } else {
    // MVB_FUSE_BISECT=1: binary-search the smallest prefix of the chain whose bulk union
    // FAILS acceptance, and name the primitive at the boundary -- one run instead of a
    // rebuild-per-guess. Diagnostic only; falls through to the normal strategies after.
    if (std::getenv("MVB_FUSE_BISECT") && solids.size() > 2) {
        auto tryPrefix = [&](size_t n) -> bool {
            try {
                TopTools_ListOfShape args, tools;
                args.Append(solids.front());
                for (size_t i = 1; i < n; ++i) tools.Append(solids[i]);
                BRepAlgoAPI_Fuse f;
                f.SetArguments(args);
                f.SetTools(tools);
                f.SetFuzzyValue(fuzzy);
                f.Build();
                if (!f.IsDone()) return false;
                GProp_GProps gp; BRepGProp::VolumeProperties(f.Shape(), gp);
                double expect = 0.0;
                for (size_t i = 0; i < n; ++i) {
                    GProp_GProps g2; BRepGProp::VolumeProperties(solids[i], g2);
                    expect += g2.Mass();
                }
                int nsol = 0;
                for (TopExp_Explorer e(f.Shape(), TopAbs_SOLID); e.More(); e.Next()) ++nsol;
                return nsol == 1 && gp.Mass() > 0.90 * expect &&
                       BRepCheck_Analyzer(f.Shape()).IsValid();
            } catch (const Standard_Failure&) { return false; }
        };
        size_t lo = 2, hi = solids.size();   // lo: works (assume); hi: fails (observed)
        if (!tryPrefix(hi)) {
            while (hi - lo > 1) {
                size_t mid = (lo + hi) / 2;
                if (tryPrefix(mid)) lo = mid; else hi = mid;
                std::cerr << "[fuse-bisect] prefix " << mid << ": "
                          << (lo == mid ? "OK" : "FAIL") << "\n";
            }
            std::cerr << "[fuse-bisect] first failing prefix = " << hi << "; boundary prim ["
                      << path.prims[std::min(hi - 1, path.prims.size() - 1)].label << "]\n";
        } else {
            std::cerr << "[fuse-bisect] full bulk fuse ACCEPTS under bisect criteria\n";
        }
    }
    // STRATEGY 1: one bulk union (arguments = first solid, tools = the rest). Fast, and it is what
    // welds a short rect chain, so it stays first.
    bool bulkCollapsed = false;
    try {
        TopTools_ListOfShape args, tools;
        args.Append(solids.front());
        for (size_t i = 1; i < solids.size(); ++i) tools.Append(solids[i]);
        BRepAlgoAPI_Fuse fuse;
        fuse.SetArguments(args);
        fuse.SetTools(tools);
        fuse.SetFuzzyValue(fuzzy);
        fuse.Build();
        if (fuse.IsDone()) {
            TopoDS_Shape ok = acceptFused(fuse.Shape(), "bulk", &bulkCollapsed);
            if (!ok.IsNull()) return ok;
        } else if (diag) {
            std::cerr << "[rect-column] bulk fuse not done\n";
        }
    } catch (const Standard_Failure&) {
        if (diag) std::cerr << "[rect-column] bulk fuse threw -> pairwise\n";
    }
    // STRATEGY 2: PAIRWISE union up a balanced tree. Every boolean sees exactly two bodies that
    // genuinely overlap, so BOPAlgo never has to reason about many arguments at once -- which is
    // where the one-shot union can silently drop tools on a long chain of near-tangent sections.
    // O(N) booleans instead of one, but each is trivial.
    //
    // NOT attempted when the bulk union COLLAPSED. A collapse means the turns genuinely TOUCH and
    // the union fused their stacked flats into an electrically-wrong brick (current would short
    // straight through instead of spiralling) -- there the per-turn compound is the CORRECT answer,
    // as the header comment says, and re-attempting the same union pairwise is both pointless and
    // unsafe: on the E70 5x1 mm stacked rect winding (272 primitives, union 1090 of 31474 mm3) it
    // segfaults inside BOPAlgo. Escalate only when the bulk union failed on VALIDITY or lost
    // material, which is the case the pairwise tree was added for.
    if (bulkCollapsed && diag)
        std::cerr << "[rect-column] bulk union collapsed the touching turns -> keeping the "
                  << "per-turn compound (no pairwise retry)\n";
    try {
        if (bulkCollapsed) throw Standard_Failure();
        std::vector<TopoDS_Shape> level = solids;
        auto lvol = [](const TopoDS_Shape& sh) {
            GProp_GProps gp; BRepGProp::VolumeProperties(sh, gp); return gp.Mass();
        };
        // A failed PAIR no longer abandons the whole tree: both operands carry to the next
        // level unfused and the reduction continues until no pair makes progress. Each merge
        // is volume-guarded (the OCC silent-drop class: on 02_flyback the 384-prim bulk
        // union returned 6.1 of 230.5 mm3).
        while (level.size() > 1) {
            std::vector<TopoDS_Shape> next;
            next.reserve(level.size() / 2 + 1);
            bool progressed = false;
            for (size_t i = 0; i + 1 < level.size(); i += 2) {
                const double va = lvol(level[i]), vb = lvol(level[i + 1]);
                TopoDS_Shape fusedPair;
                try {
                    TopTools_ListOfShape a, t;
                    a.Append(level[i]);
                    t.Append(level[i + 1]);
                    BRepAlgoAPI_Fuse f;
                    f.SetArguments(a);
                    f.SetTools(t);
                    f.SetFuzzyValue(fuzzy);
                    f.Build();
                    if (f.IsDone() && !f.Shape().IsNull() &&
                        lvol(f.Shape()) >= 0.98 * (va + vb))
                        fusedPair = f.Shape();
                } catch (const Standard_Failure&) {}
                if (!fusedPair.IsNull()) { next.push_back(fusedPair); progressed = true; }
                else { next.push_back(level[i]); next.push_back(level[i + 1]); }
            }
            if (level.size() % 2) next.push_back(level.back());
            if (!progressed) break;
            level.swap(next);
        }
        if (level.size() == 1) {
            TopoDS_Shape merged = acceptFused(level.front(), "pairwise");
            if (!merged.IsNull()) return merged;
        } else if (level.size() * 4 < solids.size()) {
            // PARTIAL weld that made real progress (>= 8x fewer pieces): hand the compound of
            // welded groups to the consumer's fragment instead of 100s of butt-joined prims
            // (measured on 02_flyback: 384 raw prims degenerate the downstream fragment --
            // "the 1D mesh seems not to be forming a closed loop" -- while a dozen groups
            // mesh). Guard: total copper conserved and every group BRepCheck-valid; else the
            // exact chain below stays the answer.
            // Single-prim leftovers (pieces that never welded -- identity-preserved through
            // the carry-over) are REBUILT with a split profile: an unwelded quadric prim's
            // closed side face is periodic and unmeshable (04_forward: 'Impossible to mesh
            // periodic surface 655' from a 38-group partial whose singles kept full cylinders).
            for (auto& piece : level) {
                for (size_t si = 0; si < solids.size(); ++si) {
                    if (!piece.IsSame(solids[si])) continue;
                    TopoDS_Shape rs = rectPrimSolid(*splitArgs[si].pr, path.wireWidth,
                                                    path.wireHeight, splitArgs[si].axialA,
                                                    splitArgs[si].axialB, splitArgs[si].extA,
                                                    splitArgs[si].extB, round, radius,
                                                    /*splitOverride=*/6);
                    if (!rs.IsNull()) piece = rs;
                    break;
                }
            }
            double total = 0.0; bool allValid = true;
            for (auto& piece : level) {
                if (!BRepCheck_Analyzer(piece).IsValid()) {
                    // A weld seam can leave a locally-invalid face; give OCC's healer one
                    // shot before abandoning the whole partial weld -- but only accept a
                    // repair that CONSERVES the copper (the healer must not become another
                    // silent geometry rewriter).
                    const double before = lvol(piece);
                    try {
                        ShapeFix_Shape fixer(piece);
                        fixer.Perform();
                        TopoDS_Shape healed = fixer.Shape();
                        if (BRepCheck_Analyzer(healed).IsValid() &&
                            std::abs(lvol(healed) - before) <= 0.005 * before)
                            piece = healed;
                    } catch (const Standard_Failure&) {}
                }
                total += lvol(piece);
                if (!BRepCheck_Analyzer(piece).IsValid()) { allValid = false; break; }
            }
            if (allValid && std::abs(total - expectedVol) <= 0.02 * expectedVol) {
                BRep_Builder bb; TopoDS_Compound comp; bb.MakeCompound(comp);
                for (const auto& piece : level) bb.Add(comp, piece);
                if (diag) std::cerr << "[rect-column] PARTIAL pairwise weld: " << level.size()
                                    << " groups from " << solids.size() << " prims (v="
                                    << total * 1e9 << "mm3)\n";
                return comp;
            }
            if (diag) std::cerr << "[rect-column] pairwise fuse: " << level.size()
                                << " left, allValid=" << allValid << " -> chain fallback\n";
        } else if (diag) {
            std::cerr << "[rect-column] pairwise fuse failed (" << level.size() << " left)\n";
        }
    } catch (const Standard_Failure&) {
        if (diag) std::cerr << "[rect-column] pairwise fuse threw -> compound fallback\n";
    }
    }
    // HARD GUARD (ABT #332): a FEM export must never receive a "continuous" winding that is
    // actually a set of DISCONNECTED per-turn solids. That is not the requested geometry (it is the
    // idealised-rings model wearing the real-winding label), it silently changes the physics, and
    // downstream it only surfaces as an unmeshable body. Drawing (femReady=false) may keep the
    // compound unchecked — it renders fine and never reaches a solver.
    //
    // A CHAIN of butt-joined per-primitive solids is NOT that failure. Consecutive primitives share
    // a co-located end face by construction, so the compound is one electrically continuous
    // conductor; gmsh's fragment welds coincident faces into a conforming assembly and the meshed
    // conductor comes out as a single connected region -- which is all the FEM needs. (It is also
    // the only geometry that meshes: the boolean union of these same primitives is BRepCheck-
    // invalid at every axial transition, and the one-piece B-spline sweep, though valid, takes
    // gmsh >20 min against 2.5 s for the analytic pieces.) So the guard checks CONNECTEDNESS of the
    // chain, not the solid count.
    if (path.femReady && solids.size() > 1) {
        // Link solids whose bounding boxes touch within a fraction of the wire radius; require one
        // connected component. Butt-joined neighbours overlap-to-touching, disjoint rings do not.
        const double touch = 0.25 * path.wireRadius;
        const size_t n = solids.size();
        std::vector<std::array<double, 6>> bx(n);
        for (size_t i = 0; i < n; ++i) {
            Bnd_Box bb; BRepBndLib::Add(solids[i], bb);
            if (bb.IsVoid()) { bx[i] = {0, 0, 0, 0, 0, 0}; continue; }
            bb.Get(bx[i][0], bx[i][1], bx[i][2], bx[i][3], bx[i][4], bx[i][5]);
        }
        std::vector<size_t> uf(n);
        for (size_t i = 0; i < n; ++i) uf[i] = i;
        std::function<size_t(size_t)> find = [&](size_t a) {
            while (uf[a] != a) { uf[a] = uf[uf[a]]; a = uf[a]; }
            return a;
        };
        for (size_t i = 0; i < n; ++i)
            for (size_t j = i + 1; j < n; ++j) {
                const bool sep = bx[i][0] > bx[j][3] + touch || bx[j][0] > bx[i][3] + touch ||
                                 bx[i][1] > bx[j][4] + touch || bx[j][1] > bx[i][4] + touch ||
                                 bx[i][2] > bx[j][5] + touch || bx[j][2] > bx[i][5] + touch;
                if (sep) continue;
                size_t ra = find(i), rb = find(j);
                if (ra != rb) uf[ra] = rb;
            }
        std::set<size_t> roots;
        for (size_t i = 0; i < n; ++i) roots.insert(find(i));
        if (roots.size() > 1) {
            throw std::runtime_error(
                "ConductorBuilder: real-winding FEM export for '" + path.name + "' produced " +
                std::to_string(solids.size()) + " primitives forming " +
                std::to_string(roots.size()) + " DISCONNECTED bodies (single-body sweep failed, "
                "the boolean union was rejected, and the pieces do not touch). Refusing to return "
                "a disconnected conductor for FEM. Re-run with MVB_DIAG=1 to see the rejection "
                "reason; MVB_FUSE_FUZZY tunes the union tolerance (default 1e-7).");
        }
        // STRATEGY 2.5: LINEAR chain-order group welding. The chain is welded one neighbour
        // at a time into maximal groups; a junction whose weld fails acceptance (single solid,
        // copper conserved, BRepCheck valid) simply STARTS A NEW GROUP. Neighbours genuinely
        // overlap (junction lenses), so most welds succeed; the hard boundaries -- lead
        // junctions, layer transitions whose tubes lie near-tangent over whole faces (the
        // 2-layer 02_flyback chain that defeats every bulk/tree union) -- become group seams.
        // Groups butt on clean single-circle quadric discs, which the consumer's fragment
        // glues without the phase-mismatch slivers of split-profile discs; a group left as a
        // SINGLE prim is rebuilt with a split profile (its untrimmed quadric side face would
        // be periodic and unmeshable).
        {
            auto lvol2 = [](const TopoDS_Shape& sh) {
                GProp_GProps gp; BRepGProp::VolumeProperties(sh, gp); return gp.Mass();
            };
            std::vector<std::pair<TopoDS_Shape, int>> groups;   // shape, prim count
            TopoDS_Shape acc = solids.front();
            double accMass = lvol2(acc);
            int accCount = 1;
            size_t accStart = 0;
            std::vector<size_t> singleIdx;   // splitArgs index of single-prim groups
            auto flush = [&](size_t nextStart) {
                if (accCount == 1) singleIdx.push_back(accStart);
                groups.push_back({acc, accCount});
                accStart = nextStart;
            };
            for (size_t i2 = 1; i2 < solids.size(); ++i2) {
                bool welded = false;
                const double mt = lvol2(solids[i2]);
                try {
                    BRepAlgoAPI_Fuse f(acc, solids[i2]);
                    f.SetFuzzyValue(fuzzy);
                    f.Build();
                    if (f.IsDone() && !f.Shape().IsNull()) {
                        int nsol = 0;
                        for (TopExp_Explorer e(f.Shape(), TopAbs_SOLID); e.More(); e.Next()) ++nsol;
                        const double m = lvol2(f.Shape());
                        if (nsol == 1 && m >= 0.97 * accMass + 0.5 * mt &&
                            BRepCheck_Analyzer(f.Shape()).IsValid()) {
                            acc = f.Shape(); accMass = m; ++accCount; welded = true;
                        }
                    }
                } catch (const Standard_Failure&) {}
                if (!welded) {
                    flush(i2);
                    acc = solids[i2]; accMass = mt; accCount = 1;
                }
            }
            flush(solids.size());
            if (groups.size() > 1 && groups.size() * 6 < solids.size()) {
                // rebuild the single-prim groups with split profiles
                bool okSingles = true;
                size_t gi = 0;
                for (auto& g : groups) {
                    if (g.second == 1) {
                        // find which splitArgs index: singles were recorded in order
                        (void)gi;
                    }
                }
                std::vector<TopoDS_Shape> rebuilt;
                for (size_t si : singleIdx) {
                    TopoDS_Shape rs = rectPrimSolid(*splitArgs[si].pr, path.wireWidth,
                                                    path.wireHeight, splitArgs[si].axialA,
                                                    splitArgs[si].axialB, splitArgs[si].extA,
                                                    splitArgs[si].extB, round, radius,
                                                    /*splitOverride=*/6);
                    if (rs.IsNull()) { okSingles = false; break; }
                    rebuilt.push_back(rs);
                }
                double total = 0.0;
                for (auto& g : groups) total += lvol2(g.first);
                if (okSingles && std::abs(total - expectedVol) <= 0.03 * expectedVol) {
                    BRep_Builder bb; TopoDS_Compound comp; bb.MakeCompound(comp);
                    size_t ri = 0;
                    size_t siPos = 0;
                    for (auto& g : groups) {
                        if (g.second == 1 && siPos < singleIdx.size() && ri < rebuilt.size()) {
                            bb.Add(comp, rebuilt[ri]); ++ri; ++siPos;
                        } else {
                            bb.Add(comp, g.first);
                        }
                    }
                    if (diag) std::cerr << "[rect-column] LINEAR group weld: " << groups.size()
                                        << " groups (" << singleIdx.size()
                                        << " split singles) from " << solids.size()
                                        << " prims, v=" << total * 1e9 << "mm3\n";
                    return comp;
                }
                if (diag) std::cerr << "[rect-column] linear weld rejected (groups="
                                    << groups.size() << " v=" << total * 1e9 << " exp="
                                    << expectedVol * 1e9 << "mm3 singlesOk=" << okSingles
                                    << ")\n";
            } else if (diag) {
                std::cerr << "[rect-column] linear weld insufficient (" << groups.size()
                          << " groups from " << solids.size() << " prims)\n";
            }
        }
        // STRATEGY 3: fuse the WRAP CHAIN alone and keep the LEADS as separate pieces. The
        // fuse-bisect diagnostic showed the union poisoned from the very first LEAD<->wrap
        // junction (prefix 3 of 384 already fails acceptance on 02_flyback), while wrap-only
        // chains fuse routinely. Leads butt the fused body with real perpendicular overlap
        // lenses, so the consumer's weld/fragment joins them -- the exact configuration that
        // meshed ETD34. Lead pieces are rebuilt with SPLIT profiles (an unfused quadric
        // cylinder's periodic side face is unmeshable).
        {
            std::vector<TopoDS_Shape> wrapSolids, leadSolids;
            bool rebuildOk = true;
            for (size_t i2 = 0; i2 < splitArgs.size(); ++i2) {
                const Primitive& pr2 = *splitArgs[i2].pr;
                if (pr2.isLead || pr2.isConnection) {
                    TopoDS_Shape rs = rectPrimSolid(pr2, path.wireWidth, path.wireHeight,
                                                    splitArgs[i2].axialA, splitArgs[i2].axialB,
                                                    splitArgs[i2].extA, splitArgs[i2].extB,
                                                    round, radius, /*splitOverride=*/6);
                    if (rs.IsNull()) { rebuildOk = false; break; }
                    leadSolids.push_back(rs);
                } else {
                    wrapSolids.push_back(solids[i2]);
                }
            }
            if (rebuildOk && !wrapSolids.empty() && !leadSolids.empty()) {
                try {
                    TopTools_ListOfShape args, tools;
                    args.Append(wrapSolids.front());
                    for (size_t i2 = 1; i2 < wrapSolids.size(); ++i2)
                        tools.Append(wrapSolids[i2]);
                    BRepAlgoAPI_Fuse f;
                    f.SetArguments(args);
                    f.SetTools(tools);
                    f.SetFuzzyValue(fuzzy);
                    f.Build();
                    if (f.IsDone()) {
                        GProp_GProps gp; BRepGProp::VolumeProperties(f.Shape(), gp);
                        double expect = 0.0;
                        for (const auto& ws : wrapSolids) {
                            GProp_GProps g2; BRepGProp::VolumeProperties(ws, g2);
                            expect += g2.Mass();
                        }
                        int nsol = 0;
                        for (TopExp_Explorer e(f.Shape(), TopAbs_SOLID); e.More(); e.Next()) ++nsol;
                        if (nsol == 1 && gp.Mass() > 0.90 * expect &&
                            !hasDegenerateSheetFace(f.Shape(), path.wireRadius) &&
                            BRepCheck_Analyzer(f.Shape()).IsValid()) {
                            BRep_Builder bb; TopoDS_Compound comp; bb.MakeCompound(comp);
                            bb.Add(comp, f.Shape());
                            for (const auto& ls : leadSolids) bb.Add(comp, ls);
                            if (diag) std::cerr << "[rect-column] WRAP-ONLY fuse accepted (v="
                                                << gp.Mass() * 1e9 << "mm3) + "
                                                << leadSolids.size()
                                                << " split lead piece(s) for consumer welding\n";
                            return comp;
                        }
                        if (diag) std::cerr << "[rect-column] wrap-only fuse rejected (nsol="
                                            << nsol << " v=" << gp.Mass() * 1e9 << " exp="
                                            << expect * 1e9 << "mm3)\n";
                    }
                } catch (const Standard_Failure& e2) {
                    if (diag) std::cerr << "[rect-column] wrap-only fuse threw ("
                                        << e2.GetMessageString() << ")\n";
                }
            }
        }
        // Rebuild every prim with a 6-arc SPLIT profile for the chain: the closed-circle
        // originals carry periodic side faces and full-circle seam edges that a fragment
        // cannot mesh (measured on 02_flyback: gmsh's 1D loop ends 2 pi from its start).
        {
            std::vector<TopoDS_Shape> resplit;
            resplit.reserve(splitArgs.size());
            bool okSplit = true;
            for (const auto& saPack : splitArgs) {
                TopoDS_Shape rs = rectPrimSolid(*saPack.pr, path.wireWidth, path.wireHeight,
                                                saPack.axialA, saPack.axialB, saPack.extA,
                                                saPack.extB, round, radius, /*splitOverride=*/6);
                if (rs.IsNull()) { okSplit = false; break; }
                resplit.push_back(rs);
            }
            if (okSplit && resplit.size() == solids.size()) solids.swap(resplit);
            else if (diag) std::cerr << "[rect-column] split-profile rebuild failed -- "
                                        "keeping quadric prims\n";
        }
        if (diag) std::cerr << "[rect-column] returning the CONNECTED butt-joined chain of "
                            << n << " primitives for FEM (fragment welds it downstream)\n";
    }
    BRep_Builder b;
    TopoDS_Compound comp;
    b.MakeCompound(comp);
    for (const auto& s : solids) b.Add(comp, s);
    return comp;
}

// Split the path into maximal continuous runs: at revisited points (the wire crossing
// itself), at closed rings, and at lead<->wrap boundaries (the only non-tangent,
// 90-degree junctions of the path — a pipe swept across them grows an unbounded mitre
// spike, while wrap chains are tangent-continuous by construction and sweep as one pipe).
std::vector<std::pair<size_t, size_t>> continuousRuns(const ConductorPath& path) {
    auto quantize = [](const gp_Pnt& p) {
        return std::make_tuple(static_cast<long long>(std::llround(p.X() * 1e9)),
                               static_cast<long long>(std::llround(p.Y() * 1e9)),
                               static_cast<long long>(std::llround(p.Z() * 1e9)));
    };
    std::map<std::tuple<long long, long long, long long>, int> visits;
    for (const auto& pr : path.prims) {
        auto [a, b] = primEndpoints(pr);
        ++visits[quantize(a)];
        ++visits[quantize(b)];
    }
    // A junction between consecutive primitives counts twice; more visits = crossover.
    std::vector<std::pair<size_t, size_t>> runs;
    auto isClosedRing = [&](const Primitive& pr) {
        auto [a, b] = primEndpoints(pr);
        return quantize(a) == quantize(b);
    };
    size_t start = 0;
    for (size_t i = 0; i < path.prims.size(); ++i) {
        auto [a, b] = primEndpoints(path.prims[i]);
        if (isClosedRing(path.prims[i])) {
            // A closed spine segfaults OCCT's pipe-shell; isolate it (swept by the exact
            // revolve in the piecewise fallback).
            if (i > start) runs.push_back({start, i});
            runs.push_back({i, i + 1});
            start = i + 1;
            continue;
        }
        // Split at lead<->wrap and inter-layer-link<->wrap boundaries alike: both meet the
        // wrap chain at 90-degree corners that would flare the pipe-shell, and both sweep
        // robustly as cylinders + sphere elbows on their own.
        bool leadBoundary = (i + 1 < path.prims.size()) &&
                            ((path.prims[i].isLead != path.prims[i + 1].isLead) ||
                             (path.prims[i].isConnection != path.prims[i + 1].isConnection));
        if (i + 1 < path.prims.size() && (visits[quantize(b)] > 2 || leadBoundary)) {
            if (std::getenv("MVB_DEBUG_RUNS")) {
                std::cerr << "SPLIT after prim " << i << " [" << path.prims[i].label
                          << "] visits=" << visits[quantize(b)] << " leadBoundary="
                          << leadBoundary << " at (" << b.X() << "," << b.Y() << ","
                          << b.Z() << ")\n";
            }
            runs.push_back({start, i + 1});
            start = i + 1;
        }
    }
    if (start < path.prims.size()) runs.push_back({start, path.prims.size()});
    return runs;
}

// Per-piece fallback: sweep each primitive on its own (exact revolve for arcs, exact
// cylinder + sphere elbows for segments, pipe for spirals), then fuse.
TopoDS_Shape sweepPiecewise(const Primitive* const* prims, size_t count, double wireRadius,
                            int wirePolygonSegments, const std::vector<gp_Pnt>& flatCaps) {
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    // A terminal end (a free end of the whole conductor) gets NO sphere cap: the swept
    // cylinder's flat end face is left exposed, so downstream FEM (Ansys, MFEM, ...) has a
    // planar surface to assign the current boundary condition on. Interior joints keep
    // their sphere so the wire envelope stays continuous across the elbows.
    auto isFlatCap = [&](const gp_Pnt& p) {
        for (const auto& f : flatCaps)
            if (p.Distance(f) < 1e-9) return true;
        return false;
    };

    for (size_t pi = 0; pi < count; ++pi) {
        const auto& pr = *prims[pi];
        if (pr.kind == Primitive::ARC3) {
            double radius = pr.arc.v0.Modulus();
            if (radius < 1e-12 || std::abs(pr.arc.sweep) < 1e-12) continue;
            gp_Pnt start(pr.arc.c.XYZ() + pr.arc.v0);
            gp_XYZ tangent = pr.arc.axis.Crossed(pr.arc.v0);
            TopoDS_Face prof = wireProfile(start, gp_Dir(tangent), wireRadius, 0);
            BRepPrimAPI_MakeRevol rev(prof, gp_Ax1(pr.arc.c, gp_Dir(pr.arc.axis)),
                                      pr.arc.sweep);
            if (!rev.IsDone() || rev.Shape().IsNull()) {
                throw std::runtime_error("ConductorBuilder: arc revolve failed for [" +
                                         pr.label + "]");
            }
            builder.Add(compound, rev.Shape());
            continue;
        }
        if (pr.kind == Primitive::SEG) {
            double len = pr.seg.a.Distance(pr.seg.b);
            if (len < 1e-12) continue;
            gp_Dir dir(gp_XYZ(pr.seg.b.XYZ() - pr.seg.a.XYZ()));
            // Prism over a SPLIT (6-arc) circular profile, NOT BRepPrimAPI_MakeCylinder: the
            // exact quadric cylinder's side face is a closed PERIODIC surface, and once this
            // piece is WELDED into the conductor mid-body (station-overlap leads), that face
            // survives with internal trims and gmsh refuses it ("Impossible to mesh periodic
            // surface", measured on ETD34's entrance lead). Six cylindrical panels carry the
            // identical geometry without periodicity.
            if (std::getenv("MVB_SEG_QUADRIC")) {   // diagnostic: the old exact cylinder
                builder.Add(compound,
                            BRepPrimAPI_MakeCylinder(gp_Ax2(pr.seg.a, dir), wireRadius, len)
                                .Shape());
            } else {
                TopoDS_Face segProf = BRepBuilderAPI_MakeFace(
                    wireProfileWireSplit(pr.seg.a, dir, wireRadius, 6)).Face();
                TopoDS_Shape prism = BRepPrimAPI_MakePrism(segProf, gp_Vec(dir) * len).Shape();
                if (std::getenv("MVB_DIAG")) {
                    GProp_GProps gp_; BRepGProp::VolumeProperties(prism, gp_);
                    const double expect = kPi * wireRadius * wireRadius * len;
                    if (std::abs(gp_.Mass() - expect) > 0.01 * expect)
                        std::cerr << "[seg-prism] VOLUME MISMATCH [" << pr.label << "] got "
                                  << gp_.Mass()*1e9 << " expect " << expect*1e9 << " mm3\n";
                }
                builder.Add(compound, prism);
            }
            // Sphere joints keep the wire envelope continuous across the plane elbows —
            // except at a terminal end, which stays a flat cylinder cap (FEM surface).
            if (!isFlatCap(pr.seg.a))
                builder.Add(compound, BRepPrimAPI_MakeSphere(pr.seg.a, wireRadius).Shape());
            if (!isFlatCap(pr.seg.b))
                builder.Add(compound, BRepPrimAPI_MakeSphere(pr.seg.b, wireRadius).Shape());
            continue;
        }
        TopoDS_Edge e = primEdge(pr, wireRadius);
        TopoDS_Shape pipe;
        if (!e.IsNull()) {
            try {
                TopoDS_Wire spine = BRepBuilderAPI_MakeWire(e).Wire();
                auto pts = samplePrim(pr, wireRadius);
                gp_Vec t0(pts[1].XYZ() - pts[0].XYZ());
                if (t0.Magnitude() > 1e-12) {
                    // Follows the requested faceting -- never a hardcoded round profile.
                    TopoDS_Wire prof = wireProfileWire(pts.front(), gp_Dir(t0),
                                                       wireRadius, wirePolygonSegments);
                    BRepOffsetAPI_MakePipeShell ps(spine);
                    ps.Add(prof);
                    ps.Build();
                    if (ps.IsDone() && ps.MakeSolid()) pipe = ps.Shape();
                }
            } catch (const Standard_Failure&) {
            }
        }
        if (!pipe.IsNull()) {
            builder.Add(compound, pipe);
        } else {
            auto pts = samplePrim(pr, wireRadius);
            for (size_t j = 0; j + 1 < pts.size(); ++j) {
                double len = pts[j].Distance(pts[j + 1]);
                if (len < 1e-12) continue;
                gp_Dir dir(gp_XYZ(pts[j + 1].XYZ() - pts[j].XYZ()));
                builder.Add(compound,
                            BRepPrimAPI_MakeCylinder(gp_Ax2(pts[j], dir), wireRadius, len)
                                .Shape());
                if (!isFlatCap(pts[j]))
                    builder.Add(compound, BRepPrimAPI_MakeSphere(pts[j], wireRadius).Shape());
            }
            if (!isFlatCap(pts.back()))
                builder.Add(compound, BRepPrimAPI_MakeSphere(pts.back(), wireRadius).Shape());
        }
    }
    return compound;
}

// Fuse all solids of a conductor's compound into as FEW bodies as possible, in ONE general
// boolean (BRepAlgoAPI_Fuse over the whole set via SetArguments/SetTools). A single BOP with
// OBB culling is both far faster than N sequential pairwise fuses (which is O(N) booleans on
// swept BSpline solids — minutes for a many-piece winding) AND more robust: OCCT resolves all
// the mutual intersections together, so overlapping chunk pipes + elbow/junction spheres merge
// into one solid while any piece that genuinely touches nothing stays its own solid. No glue:
// the pieces OVERLAP volumetrically (junction spheres straddle the pipe caps), which is a real
// intersection, not the coincident-face case glue is for. A gross-corruption volume guard
// rejects the known OCCT "ate a whole body" defect; on any failure the exact (unfused) compound
// is returned rather than losing copper — same guarantee as before, just reached far less often.
// A conductor face whose bounding extent is many wire-radii long but whose area is a sliver:
// the degenerate SEAM faces OCCT's boolean union leaves where two swept pipes (or a pipe and
// its junction sphere) merge -- the thin "line/cone" sheets seen on a multilayer winding. A
// genuine tube-wall face of the same extent has area ~ 2*pi*r*length, orders of magnitude
// larger; end caps and short faces (diag <= a few radii) are exempt.
bool hasDegenerateSheetFace(const TopoDS_Shape& shape, double wireRadius) {
    for (TopExp_Explorer fe(shape, TopAbs_FACE); fe.More(); fe.Next()) {
        Bnd_Box fb;
        BRepBndLib::Add(fe.Current(), fb);
        if (fb.IsVoid()) continue;
        double x0, y0, z0, x1, y1, z1;
        fb.Get(x0, y0, z0, x1, y1, z1);
        double diag = std::sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) +
                                (z1 - z0) * (z1 - z0));
        if (diag <= 4.0 * wireRadius) continue;
        GProp_GProps fp;
        BRepGProp::SurfaceProperties(fe.Current(), fp);
        if (fp.Mass() < 0.2 * wireRadius * diag) return true;
    }
    return false;
}

// Pairwise, coverage-guarded weld of a compound's solids into ONE body. A single N-way
// BRepAlgo fuse can silently DROP an input piece (measured: ETD34's exit lead, 2.75 mm3 --
// only 0.54% of the winding, invisible to total-mass guards); welding one piece at a time,
// biggest first, lets each step be verified (one solid out, mass conserved to the overlap
// lens) and REJECTED piece-wise. Pieces that will not weld stay separate in the returned
// compound -- the caller decides whether that is acceptable.
TopoDS_Shape weldSolidsPairwise(const TopoDS_Shape& compound, double fuzzy) {
    std::vector<TopoDS_Shape> pieces;
    for (TopExp_Explorer e(compound, TopAbs_SOLID); e.More(); e.Next())
        pieces.push_back(e.Current());
    if (pieces.size() <= 1) return compound;
    // Weld in the MILLIMETRE frame: OCC booleans, like MakePipeShell (see rawGrownSolid),
    // are unreliable on sub-millimetre features at metre scale -- measured on 03_buck's
    // tangential lead-corner lens, where the metre-frame fuse united the bodies but ATE
    // 37 mm3 of copper (caught by the volume guard below). Pure scaling is an exact
    // affine map both ways; fuzzy scales with the frame so the physical tolerance is
    // unchanged.
    gp_Trsf up, down;
    up.SetScale(gp_Pnt(0, 0, 0), 1000.0);
    down.SetScale(gp_Pnt(0, 0, 0), 1.0 / 1000.0);
    for (auto& p : pieces)
        p = BRepBuilderAPI_Transform(p, up, Standard_True).Shape();
    fuzzy *= 1000.0;
    auto backDown = [&](const TopoDS_Shape& s) {
        return BRepBuilderAPI_Transform(s, down, Standard_True).Shape();
    };
    auto massOf = [](const TopoDS_Shape& sh) {
        GProp_GProps gp; BRepGProp::VolumeProperties(sh, gp); return gp.Mass();
    };
    std::sort(pieces.begin(), pieces.end(),
              [&](const TopoDS_Shape& x, const TopoDS_Shape& y){ return massOf(x) > massOf(y); });
    TopoDS_Shape acc = pieces[0];
    double accMass = massOf(acc);
    std::vector<TopoDS_Shape> loose;
    const bool weldDiag = std::getenv("MVB_DIAG") != nullptr;
    auto weldOne = [&](const TopoDS_Shape& piece) -> bool {
        const double mt = massOf(piece);
        try {
            BRepAlgoAPI_Fuse f(acc, piece);
            if (fuzzy > 0.0) f.SetFuzzyValue(fuzzy);
            f.Build();
            if (!f.IsDone()) {
                if (weldDiag) std::cerr << "[weld] fuse not done (piece " << mt*1e9 << " mm3)\n";
                return false;
            }
            TopoDS_Shape out = f.Shape();
            int nsol = 0;
            for (TopExp_Explorer e(out, TopAbs_SOLID); e.More(); e.Next()) ++nsol;
            const double m = massOf(out);
            if (nsol == 1 && m >= 0.98 * accMass + 0.5 * mt && m >= 0.90 * (accMass + mt)
                && BRepCheck_Analyzer(out).IsValid()) {   // an invalid weld is a failed weld
                acc = out; accMass = m; return true;
            }
            if (weldDiag) std::cerr << "[weld] step rejected: nsol=" << nsol
                                    << " out=" << m*1e9 << " acc=" << accMass*1e9
                                    << " piece=" << mt*1e9 << " mm3\n";
        } catch (const Standard_Failure& e) {
            if (weldDiag) std::cerr << "[weld] threw (" << e.GetMessageString() << ")\n";
        }
        return false;
    };
    for (size_t i = 1; i < pieces.size(); ++i)
        if (!weldOne(pieces[i])) loose.push_back(pieces[i]);
    // second pass: a piece can weld once its coincident neighbours are in
    for (auto it = loose.begin(); it != loose.end();)
        if (weldOne(*it)) it = loose.erase(it); else ++it;
    if (loose.empty()) return backDown(acc);
    BRep_Builder bb; TopoDS_Compound out; bb.MakeCompound(out);
    bb.Add(out, acc);
    for (const auto& l : loose) bb.Add(out, l);
    return backDown(out);
}

TopoDS_Shape fuseAllSolids(const TopoDS_Shape& compound, double wireRadius) {
    std::vector<TopoDS_Shape> solids;
    for (TopExp_Explorer exp(compound, TopAbs_SOLID); exp.More(); exp.Next()) {
        solids.push_back(exp.Current());
    }
    if (solids.size() <= 1) return compound;
    auto volumeOf = [](const TopoDS_Shape& s) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(s, props);
        return props.Mass();
    };
    double vSum = 0.0;      // OVER-counts the true copper: junction spheres overlap the pipes
    double vMax = 0.0;      // the largest single piece — a valid fuse can't drop below this
    for (const auto& s : solids) {
        double v = volumeOf(s);
        vSum += v;
        vMax = std::max(vMax, v);
    }

    try {
        TopTools_ListOfShape args, tools;
        args.Append(solids.front());
        for (size_t i = 1; i < solids.size(); ++i) tools.Append(solids[i]);
        BRepAlgoAPI_Fuse fuse;
        fuse.SetArguments(args);
        fuse.SetTools(tools);
        fuse.SetUseOBB(true);         // bounding-box cull the mostly-disjoint pieces -> fast
        fuse.SetRunParallel(true);
        fuse.Build();
        if (fuse.IsDone() && !fuse.Shape().IsNull()) {
            double v = volumeOf(fuse.Shape());
            int n = 0;
            for (TopExp_Explorer e(fuse.Shape(), TopAbs_SOLID); e.More(); e.Next()) ++n;
            // A correct union sits in [vMax, vSum] (overlaps subtract, nothing is added). A
            // result well below vMax means OCCT swallowed a body; above vSum means it doubled
            // geometry (or a corrupt input reported a nonsense volume). Either way the boolean
            // is untrustworthy -> keep the exact compound.
            // Accept the union only when it welded pieces (fewer solids), the volume is sane,
            // AND it is clean: the boolean can leave degenerate seam-sliver faces where pipes
            // merge (the thin "line/cone" sheets on multilayer windings). Those break meshing
            // and render as stray sheets, so if any appear keep the exact unfused compound --
            // geometrically correct, one PRODUCT per piece -- rather than a corrupt single body.
            if (v >= 0.9 * vMax && v <= 1.02 * vSum &&
                n < static_cast<int>(solids.size()) &&
                !hasDegenerateSheetFace(fuse.Shape(), wireRadius)) {
                return fuse.Shape();
            }
        }
    } catch (const Standard_Failure& e) {
        std::cerr << "WARN ConductorBuilder: conductor fuse threw (" << e.GetMessageString()
                  << "); returning the exact unfused compound\n";
        return compound;
    }
    std::cerr << "WARN ConductorBuilder: conductor fuse did not weld the pieces; returning the "
                 "exact unfused compound (geometry correct, one PRODUCT per solid)\n";
    return compound;
}

// Drop degenerate (near-zero-volume) solids from a conductor result. OCCT's fuse can leave
// thin sheet/shell solids behind when it welds a dense multilayer winding: they carry no
// copper, render as stray sheets, and break tetrahedral meshing. Any real conductor piece --
// even a lead-junction sphere -- is orders of magnitude above the 1e-12 m^3 (1e-3 mm^3) floor.
TopoDS_Shape pruneDegenerateSolids(const TopoDS_Shape& shape) {
    constexpr double kMinSolidVolume = 1e-12;   // m^3
    std::vector<TopoDS_Shape> kept;
    int total = 0;
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        ++total;
        GProp_GProps props;
        BRepGProp::VolumeProperties(exp.Current(), props);
        if (std::abs(props.Mass()) >= kMinSolidVolume) kept.push_back(exp.Current());
    }
    if (kept.empty() || static_cast<int>(kept.size()) == total) return shape;
    if (kept.size() == 1) return kept.front();
    TopoDS_Compound out;
    BRep_Builder b;
    b.MakeCompound(out);
    for (const auto& s : kept) b.Add(out, s);
    return out;
}

// Fallback when the whole-run pipe sweep fails: sweep the run in per-ELECTRICAL-TURN chunks.
// A one-turn spine (~a handful of edges) sweeps reliably through MakePipeShell where the full
// multi-turn spine (160+ edges of a many-turn winding) trips it and forces the per-primitive
// path — which shatters the winding into hundreds of tiny cylinder/sphere solids that then
// overwhelm the boolean fuse. Chunking keeps a bad turn LOCAL: only the turn whose pipe fails
// drops to per-primitive, the rest stay clean one-turn pipes. Consecutive chunks are made to
// SHARE their boundary primitive (each chunk also sweeps the first primitive of the next turn),
// Adjacent chunks meet end-to-end at their shared boundary point; a wire-radius sphere there
// overlaps both pipe caps so the downstream fuse can weld the chunks into one body.
TopoDS_Shape sweepRunChunked(const Primitive* const* prims, size_t count, double wireRadius,
                             int wirePolygonSegments, const std::vector<gp_Pnt>& flatCaps) {
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    auto isFlatCap = [&](const gp_Pnt& p) {
        for (const auto& f : flatCaps)
            if (p.Distance(f) < 1e-9) return true;
        return false;
    };
    size_t b = 0;
    while (b < count) {
        size_t e = b + 1;
        while (e < count && prims[e]->turnOrdinal == prims[b]->turnOrdinal) ++e;
        TopoDS_Shape chunk = sweepRun(prims + b, e - b, wireRadius, wirePolygonSegments);
        if (chunk.IsNull()) {
            chunk = sweepPiecewise(prims + b, e - b, wireRadius, wirePolygonSegments, flatCaps);
        }
        if (!chunk.IsNull()) builder.Add(compound, chunk);
        if (e < count) {
            gp_Pnt j = primEndpoints(*prims[e - 1]).second;
            if (!isFlatCap(j))
                builder.Add(compound, BRepPrimAPI_MakeSphere(j, wireRadius).Shape());
        }
        b = e;
    }
    return compound;
}

// ---------------------------------------------------------------------------------------
// Conformal per-primitive toroid (mitre joints). A dense toroid's winding cannot become one
// solid (MakePipe self-intersects on the packed hole spine; a boolean fuse of the per-run pipes
// comes out invalid), and the plain per-primitive compound OVERLAPS at every joint (double
// material -> not a valid FEM assembly). Here each primitive is built as its OWN round solid,
// grown a hair past both endpoints, then sliced by the ANGLE-BISECTOR plane through each shared
// point. Two equal round tubes mitred on their bisector leave IDENTICAL elliptical faces (same
// centre, tilt and radius), so neighbouring solids ABUT on a coincident face instead of
// overlapping: a conformal, mesh-shareable assembly. Centrelines/crossings never move.
// ANALYTIC end tangents, not sampled chords. Chord directions from samplePrim run up to half
// a sample step's turn away from the true tangent (a 16-sample wrap chord is ~5-6 deg off),
// which made emitToroidConformal misread every TANGENT wrap junction as a corner and slice it
// on a chord-tilted bisector: the mitre faces then no longer coincide and the 'conformal'
// assembly falls apart into disconnected pieces (measured: PQ33 round column, 18 prims,
// 34 spurious cuts, 14 components).
static gp_Dir spiralTangent(const Spiral& sp, bool atStart) {
    const double daz = sp.az1 - sp.az0;
    const double az = atStart ? sp.az0 : sp.az1;
    double rp = 0.0, yp = 0.0;
    const double r = atStart ? sp.r0 : sp.r1;
    if (std::fabs(daz) > 1e-12 && !sp.blend) {
        rp = (sp.r1 - sp.r0) / daz;
        yp = (sp.y1 - sp.y0) / daz;
    }   // cosine-blend spirals have purely azimuthal end tangents (r' = y' = 0)
    // P(az) = (cx + r cos az, y, cz - r sin az)  [azPointC convention]
    gp_XYZ t(rp * std::cos(az) - r * std::sin(az), yp,
             -rp * std::sin(az) - r * std::cos(az));
    if (daz < 0) t *= -1.0;
    if (t.Modulus() < 1e-15) return gp_Dir(1, 0, 0);
    return gp_Dir(t);
}
static gp_Dir primFwdStart(const Primitive& p, double r) {
    switch (p.kind) {
        case Primitive::SEG: {
            gp_Vec v(p.seg.a, p.seg.b);
            if (v.Magnitude() > 1e-12) return gp_Dir(v);
            break;
        }
        case Primitive::ARC3: {
            gp_XYZ t = p.arc.axis.Crossed(p.arc.v0);
            if (p.arc.sweep < 0) t *= -1.0;
            if (t.Modulus() > 1e-15) return gp_Dir(t);
            break;
        }
        case Primitive::SPIRAL:
            return spiralTangent(p.spiral, /*atStart=*/true);
        case Primitive::BLEND:
            if (p.blendc.u.Modulus() > 1e-15) return gp_Dir(p.blendc.u);
            break;
    }
    auto pts = samplePrim(p, r);
    for (size_t i = 1; i < pts.size(); ++i) {
        gp_Vec v(pts.front(), pts[i]);
        if (v.Magnitude() > 1e-12) return gp_Dir(v);
    }
    return gp_Dir(1, 0, 0);
}
static gp_Dir primFwdEnd(const Primitive& p, double r) {
    switch (p.kind) {
        case Primitive::SEG: {
            gp_Vec v(p.seg.a, p.seg.b);
            if (v.Magnitude() > 1e-12) return gp_Dir(v);
            break;
        }
        case Primitive::ARC3: {
            gp_XYZ ve = rotateXYZ(p.arc.v0, p.arc.axis, p.arc.sweep);
            gp_XYZ t = p.arc.axis.Crossed(ve);
            if (p.arc.sweep < 0) t *= -1.0;
            if (t.Modulus() > 1e-15) return gp_Dir(t);
            break;
        }
        case Primitive::SPIRAL:
            return spiralTangent(p.spiral, /*atStart=*/false);
        case Primitive::BLEND:
            if (p.blendc.u.Modulus() > 1e-15) return gp_Dir(p.blendc.u);
            break;
    }
    auto pts = samplePrim(p, r);
    for (size_t i = pts.size(); i-- > 1;) {
        gp_Vec v(pts[i - 1], pts.back());
        if (v.Magnitude() > 1e-12) return gp_Dir(v);
    }
    return gp_Dir(1, 0, 0);
}

// A round solid for one primitive, grown by overA/overB past its start/end so a tilted bisector
// plane fully crosses the tube there (SEG -> longer cylinder, ARC3 -> wider revolve, else pipe).
// A tangent junction passes 0 -> the tube ends flush on a perpendicular cap, no growth, no cut.
// segments <= 0 keeps the EXACT round section (analytic cylinders/tori); segments > 0 emits the
// n-gon prism/revolve. This is not cosmetic: an exact round sweep produces PERIODIC surfaces
// (seam-carrying cylinders and tori), and gmsh's "Impossible to mesh periodic surface" is a
// known failure cluster on them -- the faceted section has no seam at all, which is why the
// pipeline has always meshed with segments=12. The conformal assembler must honour the same
// knob the retired sweep path did.
static TopoDS_Shape rawGrownSolid(const Primitive& pr, double r, double overA, double overB,
                                  int segments) {
    if (pr.kind == Primitive::SEG) {
        gp_Vec d(pr.seg.a, pr.seg.b);
        double len = d.Magnitude();
        if (len < 1e-12) return {};
        gp_Dir dir(d);
        gp_Pnt a = pr.seg.a.Translated(gp_Vec(dir) * (-overA));
        const double total = len + overA + overB;
        if (segments <= 0)
            return BRepPrimAPI_MakeCylinder(gp_Ax2(a, dir), r, total).Shape();
        TopoDS_Face prof = wireProfile(a, dir, r, segments);
        return BRepPrimAPI_MakePrism(prof, gp_Vec(dir) * total).Shape();
    }
    if (pr.kind == Primitive::ARC3) {
        double radius = pr.arc.v0.Modulus();
        if (radius < 1e-12 || std::abs(pr.arc.sweep) < 1e-12) return {};
        double sgn = pr.arc.sweep < 0 ? -1.0 : 1.0;
        double ddA = (overA / radius) * sgn, ddB = (overB / radius) * sgn;
        double total = pr.arc.sweep + ddA + ddB;
        if (std::abs(total) > 1.9 * kPi) return {};  // a near-full revolve would self-close: caller retries
        gp_XYZ v0e = rotateXYZ(pr.arc.v0, pr.arc.axis, -ddA);
        gp_Pnt start(pr.arc.c.XYZ() + v0e);
        gp_XYZ tangent = pr.arc.axis.Crossed(v0e);
        TopoDS_Face prof = wireProfile(start, gp_Dir(tangent), r, segments);
        BRepPrimAPI_MakeRevol rev(prof, gp_Ax1(pr.arc.c, gp_Dir(pr.arc.axis)), total);
        return rev.IsDone() ? rev.Shape() : TopoDS_Shape();
    }
    // BLEND / SPIRAL: swept as exact pipes. A pipe has no revolve/extrude axis to overhang
    // along, so mitre growth EXTENDS THE SPINE instead: a short straight run along the analytic
    // end tangent (G1-continuous, so MakePipeShell sweeps it as one smooth body) which the
    // bisector cut then trims back -- identical overhang semantics to the SEG/ARC3 branches.
    // Without this, a corner cut at a spiral/blend end carved a wedge with nothing to fill it
    // (measured: every lead<->wrap corner on PQ33 left a gap -> 5 disconnected components).
    auto pts = samplePrim(pr, r);
    if (pts.size() < 2) return {};
    TopoDS_Edge e = primEdge(pr, r);
    if (e.IsNull()) return {};
    try {
        const auto ends = primEndpoints(pr);
        const gp_Dir tA = primFwdStart(pr, r), tB = primFwdEnd(pr, r);
        BRepBuilderAPI_MakeWire mw;
        gp_Pnt spineStart = ends.first;
        if (overA > 0.0) {
            spineStart = gp_Pnt(ends.first.XYZ() - tA.XYZ() * overA);
            mw.Add(BRepBuilderAPI_MakeEdge(spineStart, ends.first).Edge());
        }
        mw.Add(e);
        if (overB > 0.0)
            mw.Add(BRepBuilderAPI_MakeEdge(ends.second,
                                           gp_Pnt(ends.second.XYZ() + tB.XYZ() * overB)).Edge());
        if (!mw.IsDone()) return {};
        TopoDS_Wire spine = mw.Wire();
        // ONE attempt, at the REQUESTED faceting -- no round-profile retry (Alf: no
        // fallbacks). A failed faceted sweep returns null and the caller throws naming the
        // primitive, so the geometry that produced it gets fixed instead of silently mixing
        // a round (periodic-surface, gmsh-hostile) piece into a faceted conductor.
        //
        // The sweep runs in a MILLIMETRE frame: OCC's MakePipeShell mis-frames
        // sub-millimetre faceted profiles in METRES (measured standalone: the identical
        // helix/profile that fails at x1 -- PipeNotDone on a 0.1 mm 12-gon, 10_emi/
        // 13_current_sense/20_iso/24_margin -- sweeps cleanly at x1000). Pure scaling is an
        // EXACT affine map both ways, so this is a numerical frame choice, not an
        // approximation. TODO: migrate the other MakePipeShell sites to the same frame as
        // they are touched.
        TopoDS_Wire prof = wireProfileWire(spineStart, tA, r, segments);
        gp_Trsf up, down;
        up.SetScale(gp_Pnt(0, 0, 0), 1000.0);
        down.SetScale(gp_Pnt(0, 0, 0), 1.0 / 1000.0);
        TopoDS_Wire spineMm =
            TopoDS::Wire(BRepBuilderAPI_Transform(spine, up, Standard_True).Shape());
        TopoDS_Wire profMm =
            TopoDS::Wire(BRepBuilderAPI_Transform(prof, up, Standard_True).Shape());
        BRepOffsetAPI_MakePipeShell ps(spineMm);
        ps.Add(profMm);
        ps.Build();
        if (ps.IsDone() && ps.MakeSolid())
            return BRepBuilderAPI_Transform(ps.Shape(), down, Standard_True).Shape();
        if (std::getenv("MVB_DIAG"))
            std::cerr << "[sweep-fail] MakePipeShell status=" << (int)ps.GetStatus()
                      << " isDone=" << ps.IsDone() << " segments=" << segments
                      << " r=" << r << "\n";
    } catch (const Standard_Failure& f) {
        if (std::getenv("MVB_DIAG"))
            std::cerr << "[sweep-fail] OCC exception: "
                      << (f.GetMessageString() ? f.GetMessageString() : "(null)") << "\n";
    }
    return {};
}

// Trim the grown overhang of `solid` that pokes past the mitre plane (through P, material side =
// keepDir): SUBTRACT a small knife box sitting on the discard side, localized to the junction.
// Neither of the two 'obvious' cuts is correct here:
//  - a global half-space keep is WRONG for near-closed wraps -- the wrap legitimately curves back
//    across the junction plane far from the joint, so the far side would be sliced off;
//  - the previous 'tight beam' knife (8r-wide box, depth spanning the solid, kept via Common) was
//    NOT a half-space: it discarded everything outside the beam laterally. Measured on PQ33: it
//    kept 17% of the first wrap and split it into 3 solids -> disconnected slivers.
// Locality bounds: the grown overhang lies within `grow`+r of P along -keepDir and within
// ~2.5r+`grow` of P laterally (bisector tilt <= ~50 deg), so a (3r+grow)-lateral x (grow+2r)-deep
// knife covers it with margin. A volume guard rejects any trim that removed more than a stub's
// worth of material (i.e. the knife reached a neighbouring pass of the same solid).
static TopoDS_Shape localMitreTrim(const TopoDS_Shape& solid, const gp_Pnt& P, const gp_Dir& keepDir,
                                   const gp_Dir& stubDir, double r, double grow) {
    // Oriented knife: the overhang stub extends along `stubDir` beyond the plane, so the box
    // follows its lateral DRIFT direction and stays tube-radius-thin on the perpendicular
    // (axial/pitch) axis -- a symmetric (3r+grow)-wide box nicked the neighbouring wrap pass
    // one pitch (~3r) away and detached a 0.37 mm fragment (measured, PQ33 junction 9->10).
    const gp_XYZ w = keepDir.XYZ() * -1.0;                              // into the discard side
    const double cosT = std::max(0.2, std::abs(stubDir.XYZ().Dot(w)));  // clamp: tilt <= ~78 deg
    const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    // Knife extent. It must cover the overhang and NOTHING else: a CURVED piece bends back
    // into the far part of the box and gets a rectangular notch gouged out of its keep side
    // (visible as square voids beside every dragback). Depth is therefore just the overhang plus
    // the tube radius -- what the tilted section actually needs -- and the perpendicular half is
    // barely wider than the tube.
    const double ell = r / cosT + 0.25 * r;   // mitre ellipse semi-extent + margin
    const double drift = (grow + r) * sinT;   // stub tip's lateral reach
    const double depth = grow + r;
    gp_XYZ u = stubDir.XYZ() - w * stubDir.XYZ().Dot(w);
    if (u.Modulus() > 1e-9) u.Normalize();
    else u = gp_Ax2(P, gp_Dir(w)).XDirection().XYZ();   // stub ~normal to plane: no drift
    const gp_XYZ v = w.Crossed(u);
    const double vHalf = 1.05 * r;
    gp_Pnt corner(P.XYZ() - u * ell - v * vHalf);
    TopoDS_Shape box = BRepPrimAPI_MakeBox(gp_Ax2(corner, gp_Dir(w), gp_Dir(u)),
                                           drift + 2.0 * ell, 2.0 * vHalf, depth).Shape();
    // BOUND THE KNIFE TO THE JUNCTION. A box has straight sides, so on a CURVED piece the tube
    // bends back into its far region and the box bites a rectangular notch out of material that
    // should have stayed -- the square voids beside every dragback. Intersecting the box with a
    // sphere centred on the joint keeps the cut strictly local: nothing further from the corner
    // than the overhang itself can be reached, whatever the piece does downstream.
    try {
        TopoDS_Shape ball =
            BRepPrimAPI_MakeSphere(P, grow + 1.5 * r).Shape();
        BRepAlgoAPI_Common local(box, ball);
        if (local.IsDone() && !local.Shape().IsNull()) box = local.Shape();
    } catch (const Standard_Failure&) {
    }
    try {
        GProp_GProps gpBefore, gpAfter;
        BRepGProp::VolumeProperties(solid, gpBefore);
        BRepAlgoAPI_Cut cut(solid, box);
        if (cut.IsDone() && !cut.Shape().IsNull()) {
            BRepGProp::VolumeProperties(cut.Shape(), gpAfter);
            const double removed = gpBefore.Mass() - gpAfter.Mass();
            // stub bound: the overhang is at most a full-radius tube of length grow+2r, doubled
            // for the tilted-ellipse wedge. More than that = the knife ate distant material.
            const double stubMax = 2.0 * kPi * r * r * (grow + 2.0 * r);
            if (removed <= stubMax && gpAfter.Mass() > 0.0) return cut.Shape();
        }
    } catch (const Standard_Failure&) {
    }
    return solid;
}

TopoDS_Shape emitToroidConformal(const std::vector<const Primitive*>& ptrs, double wireRadius,
                                 int segments) {
    const double over = 1.3 * wireRadius;   // covers a bisector tilt up to ~50 deg half-angle
    const double tanThresh = 0.05;          // rad (~3 deg): below this a junction is tangent -> no cut
    const bool diag = std::getenv("MVB_MITRE_DIAG") != nullptr;
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    size_t n = ptrs.size();
    if (n == 0) return compound;
    std::vector<gp_Dir> fs(n), fe(n);
    for (size_t i = 0; i < n; ++i) {
        fs[i] = primFwdStart(*ptrs[i], wireRadius);
        fe[i] = primFwdEnd(*ptrs[i], wireRadius);
    }
    int nCut = 0, nRepaired = 0, nInvalid = 0;
    std::vector<TopoDS_Shape> built;
    for (size_t i = 0; i < n; ++i) {
        // Bend at each end: a tangent joint (wrap SEG<->ARC) needs no growth and no boolean; only a
        // real corner (the entrance/exit leads) is grown and sliced on its angle-bisector plane.
        // Junctions may carry an ENDPOINT MISMATCH (layouts hand out lead ends offset from the wrap
        // start by up to a wire radius -- measured 0.36 mm on PQ33): the mitre plane goes through
        // the MIDPOINT of the two endpoints (same plane for both sides) and the growth is widened
        // by the mismatch so the overhang always reaches past the shared plane.
        const double dpS = (i > 0)
            ? primEndpoints(*ptrs[i - 1]).second.Distance(primEndpoints(*ptrs[i]).first) : 0.0;
        const double dpE = (i + 1 < n)
            ? primEndpoints(*ptrs[i]).second.Distance(primEndpoints(*ptrs[i + 1]).first) : 0.0;
        bool bentS = false;
        gp_Dir nS = fs[i];
        if (diag && i > 0) {
            static const char* kn[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
            std::cerr << "[mitre]   junction " << i - 1 << "->" << i << " " << kn[ptrs[i - 1]->kind]
                      << "->" << kn[ptrs[i]->kind] << " angle=" << fe[i - 1].Angle(fs[i]) * 180.0 / kPi
                      << " deg" << (fe[i - 1].Angle(fs[i]) > tanThresh ? " CORNER" : "")
                      << " ['" << ptrs[i - 1]->label << "' -> '" << ptrs[i]->label << "']\n";
        }
        if (i > 0 && fe[i - 1].Angle(fs[i]) > tanThresh) {
            gp_Vec s(fe[i - 1].XYZ());
            s += gp_Vec(fs[i].XYZ());
            if (s.Magnitude() > 1e-9) { nS = gp_Dir(s); bentS = true; }
        }
        bool bentE = false;
        gp_Dir nE = fe[i];
        if (i + 1 < n && fe[i].Angle(fs[i + 1]) > tanThresh) {
            gp_Vec s(fe[i].XYZ());
            s += gp_Vec(fs[i + 1].XYZ());
            if (s.Magnitude() > 1e-9) { nE = gp_Dir(s); bentE = true; }
        }
        // Growth per end: corners grow by over+mismatch (then trimmed on the mitre plane);
        // a TANGENT junction with an endpoint gap is bridged by growing the EARLIER prim's end
        // flush forward (no cut) -- one side only, so the bridge is never doubled.
        const double overS = bentS ? over + dpS : 0.0;
        const double overE = bentE ? over + dpE : (dpE > 1e-9 ? dpE : 0.0);
        TopoDS_Shape solid = rawGrownSolid(*ptrs[i], wireRadius, overS, overE, segments);
        if (solid.IsNull()) {  // ARC clamp (near-full revolve): fall back to a flush, uncut tube
            solid = rawGrownSolid(*ptrs[i], wireRadius, 0.0, 0.0, segments);
            bentS = bentE = false;
        }
        if (solid.IsNull()) {
            // Never skip silently: a missing primitive is missing copper, and its neighbours
            // may still touch each other so even the connectivity contract can stay satisfied
            // (measured: faceted blend sweeps returned null and every rect-column turn quietly
            // detached from the next).
            static const char* kindName[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
            const auto ends = primEndpoints(*ptrs[i]);
            throw std::runtime_error(
                "ConductorBuilder: the conformal assembler could not sweep primitive " +
                std::to_string(i) + " (" + kindName[ptrs[i]->kind] + ") '" + ptrs[i]->label +
                "' (" + std::to_string(ends.first.X()) + "," + std::to_string(ends.first.Y()) +
                "," + std::to_string(ends.first.Z()) + ")->(" + std::to_string(ends.second.X()) +
                "," + std::to_string(ends.second.Y()) + "," + std::to_string(ends.second.Z()) +
                "). Refusing to emit a conductor with missing copper.");
        }
        if (bentS) {
            const gp_Pnt J(0.5 * (primEndpoints(*ptrs[i - 1]).second.XYZ() +
                                  primEndpoints(*ptrs[i]).first.XYZ()));
            solid = localMitreTrim(solid, J, nS, gp_Dir(fs[i].XYZ() * -1.0), wireRadius, overS);
            ++nCut;
        }
        if (bentE) {
            const gp_Pnt J(0.5 * (primEndpoints(*ptrs[i]).second.XYZ() +
                                  primEndpoints(*ptrs[i + 1]).first.XYZ()));
            solid = localMitreTrim(solid, J, gp_Dir(nE.XYZ() * -1.0), fe[i], wireRadius, overE);
            ++nCut;
        }
        // STEP round-trip robustness: the boolean's cut curve on a torus can carry an edge tolerance
        // that a STEP export degrades into an invalid solid (valid in memory, invalid on reload).
        // ShapeFix tightens edges/tolerances in place so the written solid survives the round-trip.
        if (!solid.IsNull() && (bentS || bentE)) {
            try {
                ShapeFix_Shape fix(solid);
                fix.Perform();
                if (!fix.Shape().IsNull()) solid = fix.Shape();
            } catch (const Standard_Failure&) {
            }
        }
        // A slice that self-intersected or emptied the tube (very short/sharp primitive) is caught
        // here and rebuilt as a plain flush-capped tube -- always valid, at worst a hair of overlap
        // at that single joint. Never emit an invalid solid.
        if (solid.IsNull() || !BRepCheck_Analyzer(solid).IsValid()) {
            TopoDS_Shape flush = rawGrownSolid(*ptrs[i], wireRadius, 0.0, 0.0, segments);
            if (!flush.IsNull() && BRepCheck_Analyzer(flush).IsValid()) { solid = flush; ++nRepaired; }
            else {
                // NEVER drop a primitive silently: the copper it carries simply disappears from
                // the conductor, and downstream checks can miss it (a short corner's neighbours
                // may still touch, so even the connectivity contract stays satisfied -- that is
                // exactly how a missing dragback corner reached a STEP). Fail loudly instead,
                // naming the primitive so the geometry that produced it can be fixed.
                static const char* kindName[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
                const auto ends = primEndpoints(*ptrs[i]);
                std::ostringstream detail;
                detail << "ConductorBuilder: the conformal assembler cannot build a valid solid for "
                       << "primitive " << i << " (" << kindName[ptrs[i]->kind] << ") '"
                       << ptrs[i]->label << "' (" << ends.first.X() << "," << ends.first.Y() << ","
                       << ends.first.Z() << ")->(" << ends.second.X() << "," << ends.second.Y()
                       << "," << ends.second.Z() << ")";
                if (ptrs[i]->kind == Primitive::ARC3) {
                    detail << " sweep=" << ptrs[i]->arc.sweep
                           << " centrelineRadius=" << ptrs[i]->arc.v0.Modulus()
                           << " wireRadius=" << wireRadius;
                    if (ptrs[i]->arc.v0.Modulus() <= wireRadius * (1.0 + 1e-6)) {
                        detail << " -- the centreline radius does not exceed the wire radius, so the "
                                  "swept tube is a HORN TORUS touching its own axis of revolution "
                                  "(no valid solid exists; give the corner a larger bend radius)";
                    }
                }
                detail << ". Refusing to emit a conductor with missing copper.";
                throw std::runtime_error(detail.str());
            }
        }
        builder.Add(compound, solid);
        if (diag) {
            built.push_back(solid);
            GProp_GProps gpr;
            BRepGProp::VolumeProperties(solid, gpr);
            Bnd_Box bb;
            BRepBndLib::Add(solid, bb);
            double x0, y0, z0, x1, y1, z1;
            bb.Get(x0, y0, z0, x1, y1, z1);
            const double len = primEndpoints(*ptrs[i]).first.Distance(primEndpoints(*ptrs[i]).second);
            const double nominal = kPi * wireRadius * wireRadius * len;
            std::cerr << "[prim] " << i << " '" << ptrs[i]->label << "' vol=" << gpr.Mass()
                      << " nominal=" << nominal
                      << " ratio=" << (nominal > 0 ? gpr.Mass() / nominal : 0.0)
                      << " y=[" << y0 << "," << y1 << "]\n";
        }
    }
    // Per-solid analysis (union-find connectivity, volumes, interpenetration) is O(n^2) EXACT
    // BRepExtrema/classification -- ~15k queries on a 170-prim toroid, tens of minutes. Keep it
    // behind its own flag so plain MVB_MITRE_DIAG (junction angles, DROPPED prims, the counters)
    // stays usable during normal iteration.
    if (diag && std::getenv("MVB_MITRE_INTERPEN")) {
        // component structure over ALL pairs (mirrors the test's union-find), with volumes
        std::vector<int> uf(built.size());
        for (size_t i = 0; i < uf.size(); ++i) uf[i] = (int)i;
        std::function<int(int)> find = [&](int a) {
            while (uf[a] != a) a = uf[a] = uf[uf[a]];
            return a;
        };
        for (size_t i = 0; i < built.size(); ++i)
            for (size_t j = i + 1; j < built.size(); ++j) {
                try {
                    BRepExtrema_DistShapeShape d(built[i], built[j]);
                    if (d.IsDone() && d.Value() < 1e-6) uf[find((int)i)] = find((int)j);
                } catch (const Standard_Failure&) {
                }
            }
        for (size_t i = 0; i < built.size(); ++i) {
            GProp_GProps gp;
            BRepGProp::VolumeProperties(built[i], gp);
            int nsol = 0;
            std::string sub;
            for (TopExp_Explorer ex(built[i], TopAbs_SOLID); ex.More(); ex.Next()) {
                ++nsol;
                GProp_GProps gs;
                BRepGProp::VolumeProperties(ex.Current(), gs);
                char vb[32]; std::snprintf(vb, sizeof(vb), " %.4g", gs.Mass()); sub += vb;
            }
            std::cerr << "[mitre]   solid " << i << " comp=" << find((int)i) << " vol=" << gp.Mass()
                      << " nsolids=" << nsol << (nsol > 1 ? " subvols:" + sub : "") << "\n";
        }
        for (size_t i = 1; i < built.size(); ++i) {
            const double dp = primEndpoints(*ptrs[i - 1]).second.Distance(primEndpoints(*ptrs[i]).first);
            double ds = -1.0;
            try {
                BRepExtrema_DistShapeShape ext(built[i - 1], built[i]);
                if (ext.IsDone()) ds = ext.Value();
            } catch (const Standard_Failure&) {
            }
            if (dp > 1e-9 || ds > 1e-9)
                std::cerr << "[mitre]   junction " << i - 1 << "->" << i << " endpointGap=" << dp
                          << " solidGap=" << ds << "\n";
        }
        // Interpenetration sweep over ALL bbox-touching pairs (not just chain neighbours), the
        // same junction-grid classifier as the test battery, with LABELS -- post-prune solid
        // indices in test output cannot be mapped back to primitives. O(n^2) solid
        // classification dominates the run (tens of minutes on a 170-prim toroid), so it needs
        // its OWN flag: plain MVB_MITRE_DIAG stays fast enough for routine use.
        {
            std::vector<Bnd_Box> bxs(built.size());
            for (size_t i = 0; i < built.size(); ++i) BRepBndLib::Add(built[i], bxs[i]);
            for (size_t i = 0; i < built.size(); ++i)
                for (size_t j = i + 1; j < built.size(); ++j) {
                    if (bxs[i].Distance(bxs[j]) > 1e-9) continue;
                    double x0, y0, z0, x1, y1, z1, u0, v0, w0, u1, v1, w1;
                    bxs[i].Get(x0, y0, z0, x1, y1, z1);
                    bxs[j].Get(u0, v0, w0, u1, v1, w1);
                    const double ax = std::max(x0, u0), bx = std::min(x1, u1);
                    const double ay = std::max(y0, v0), by = std::min(y1, v1);
                    const double az = std::max(z0, w0), bz = std::min(z1, w1);
                    if (ax >= bx || ay >= by || az >= bz) continue;
                    // classify per SOLID: a built[] entry can be a multi-solid compound (a trim
                    // that fragmented), and BRepClass3d_SolidClassifier on a compound silently
                    // misclassifies -- exploding is what makes this sweep agree with the test
                    // battery's flat-solid enumeration.
                    auto insideAny = [](const TopoDS_Shape& sh, const gp_Pnt& pp) {
                        for (TopExp_Explorer ex(sh, TopAbs_SOLID); ex.More(); ex.Next()) {
                            BRepClass3d_SolidClassifier c(TopoDS::Solid(ex.Current()), pp, 1e-9);
                            if (c.State() == TopAbs_IN) return true;
                        }
                        return false;
                    };
                    int hits = 0;
                    gp_Pnt hit;
                    constexpr int N = 6;
                    for (int gx = 0; gx < N && hits == 0; ++gx)
                        for (int gy = 0; gy < N && hits == 0; ++gy)
                            for (int gz = 0; gz < N && hits == 0; ++gz) {
                                gp_Pnt pp(ax + (bx - ax) * (gx + 0.5) / N,
                                          ay + (by - ay) * (gy + 0.5) / N,
                                          az + (bz - az) * (gz + 0.5) / N);
                                if (insideAny(built[i], pp) && insideAny(built[j], pp)) {
                                    ++hits;
                                    hit = pp;
                                }
                            }
                    if (hits > 0)
                        std::cerr << "[mitre]   INTERPENETRATION prims " << i << "<->" << j
                                  << " ['" << ptrs[i]->label << "' vs '" << ptrs[j]->label
                                  << "'] at (" << hit.X() << "," << hit.Y() << "," << hit.Z()
                                  << ")\n";
                }
        }
    }
    if (diag) {
        std::cerr << "[mitre] prims=" << n << " boolean-cuts=" << nCut << " repaired=" << nRepaired
                  << " dropped-invalid=" << nInvalid << "\n";
    }
    return compound;
}

TopoDS_Shape emitConductor(const ConductorPath& path, int wirePolygonSegments) {
    if (std::getenv("MVB_DIAG"))
        std::cerr << "[emitConductor] '" << path.name << "' prims=" << path.prims.size()
                  << " useRectSolids=" << path.useRectSolids << " roundProfile="
                  << path.roundProfile << " femReady=" << path.femReady << " toroidal="
                  << path.toroidal << " singleBodyCapable=" << path.singleBodyCapable
                  << " isRect=" << path.isRectangular << "\n";
    // FEM ROUND-WIRE windings (EVERY column type: round/oblong/rect columns AND toroids):
    // the CONFORMAL MITRE ASSEMBLY, checked FIRST so no legacy strategy runs. See the block
    // comment further down (kept with the legacy paths) and docs/: research-converged
    // construction, no booleans on the winding.
    if (path.femReady && !path.isRectangular) {
        std::vector<const Primitive*> cptrs;
        cptrs.reserve(path.prims.size());
        for (const auto& pr : path.prims) cptrs.push_back(&pr);
        return pruneDegenerateSolids(
            emitToroidConformal(cptrs, path.wireRadius, wirePolygonSegments));
    }
    // Rect/oblong-column rectangular wire: the flat section can't sweep the racetrack corners, so
    // build every primitive as its own rect solid (prisms + revolved corners) and fuse.
    // MVB_RECT_SINGLE_BODY: for ROUND wire on a RECTANGULAR column, attempt the whole-path
    // single sweep FIRST -- one watertight MakePipeShell solid over the full G1 centerline, i.e.
    // no boolean union at all -- and only fall back to the per-primitive analytic path + fuse
    // (emitRectColumn) if the sweep is rejected. Rect columns were excluded from the single-body
    // path because pipe framing on the long racetrack straights was seen to displace the section
    // off the crossings ("a radius or two off"); that observation predates the centring check, so
    // rather than excluding the class we now let the ACCEPTANCE BATTERY decide per design: volume
    // match, BRepCheck validity, degenerate-sheet scan, and centrelineStaysCentred sampling every
    // wrap crossing (which is precisely the guard that would catch the historical displacement).
    // Rejected sweeps lose nothing -- the exact analytic compound is still built.
    const bool rectWholeSweep = path.useRectSolids && path.roundProfile && !path.toroidal &&
                                path.femReady && std::getenv("MVB_RECT_SINGLE_BODY") != nullptr;
    if (path.useRectSolids && !rectWholeSweep) return emitRectColumn(path);

    // Sweep each maximal continuous run as ONE pipe (the whole wrap chain when the wire
    // never crosses itself); lead runs emit as exact cylinders + sphere elbows. Fuse the
    // handful of pieces into a single solid.
    std::vector<const Primitive*> ptrs;
    ptrs.reserve(path.prims.size());
    for (const auto& pr : path.prims) ptrs.push_back(&pr);

    // MVB_SOLID_DUMP: the PRIMITIVE map. A CAD tree flattens the compound's piecewise runs, so
    // the index a viewer shows against "<conductor> N" is the primitive index -- print each one's
    // label and endpoints in mm so a solid picked by eye can be named in code.
    if (std::getenv("MVB_SOLID_DUMP")) {
        for (size_t i = 0; i < ptrs.size(); ++i) {
            auto [pa, pb] = primEndpoints(*ptrs[i]);
            static const char* kn[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
            std::cerr << "[prim] " << path.name << " " << i << " " << kn[ptrs[i]->kind]
                      << " '" << ptrs[i]->label << "' a=(" << pa.X() * 1e3 << ","
                      << pa.Y() * 1e3 << "," << pa.Z() * 1e3 << ") b=(" << pb.X() * 1e3 << ","
                      << pb.Y() * 1e3 << "," << pb.Z() * 1e3 << ") len_mm="
                      << pa.Distance(pb) * 1e3 << "\n";
        }
    }

    // The conductor's two free ends (the terminal faces): the start of the first
    // primitive (the entrance lead's outboard end) and the end of the last (the exit
    // lead's) — swept with a flat cap for FEM current-injection surfaces.
    std::vector<gp_Pnt> flatCaps;
    if (!path.prims.empty()) {
        flatCaps.push_back(primEndpoints(path.prims.front()).first);
        flatCaps.push_back(primEndpoints(path.prims.back()).second);
    }

    // SINGLE BODY (for meshing/FEM): sweep the ENTIRE conductor centerline as ONE pipe over a
    // G1 (tangent-continuous) spine -- buildFilletedWire keeps every lead/link/helix/racetrack
    // exact and rounds only the sharp junctions into physical wire bends, filleting only the
    // lead/link side of a lead<->wrap corner so the wrap endpoints (MKF crossings) stay exactly on
    // the centerline. With no C0 corner to bridge, MakePipeShell closes into a single watertight
    // solid -- no boolean union (so no seam slivers, no multi-solid compound), flat terminal caps
    // from MakeSolid's end sections. The rounded bends move the copper volume well under a percent.
    // Emitted for ROUND and OBLONG concentric columns (singleBodyCapable); rectangular columns and
    // toroids keep the exact per-run compound because MakePipeShell mis-frames their centerlines
    // (see ConductorPath::singleBodyCapable). Accept only a clean, valid, watertight single solid
    // whose volume matches the swept copper; otherwise fall back to the exact per-run compound.
    // Round/litz wire only sweeps a single body when FEM geometry is asked for; rectangular wire
    // has no valid round-profile compound fallback, so it MUST take the single-body path regardless.
    // MVB_NO_SINGLE_BODY: diagnostic switch -- skip the whole-spine single-body sweep and take
    // the per-run compound (+ fuse) fallback directly. Used to isolate sweep-crease artefacts
    // (the single pipe can fold onto itself at the wrap->lead fillet; measured 3.7 um between
    // adjacent patches on ETD34 at (0.09, 10.30, -7.96) mm) from the fused alternative.
    const bool skipSingleBody = std::getenv("MVB_NO_SINGLE_BODY") != nullptr;
    // FEM (femReady) ROUND-wire windings now default to the EXACT per-run compound, welded into
    // one region by the consumer's boolean fragment -- NOT the whole-spine single-body sweep.
    // Measured (STEP face-face distances, ETD34 + e138): the single pipe folds onto itself at
    // junctions -- 3.7 um between adjacent patches at the wrap->lead fillet (ETD34, and its
    // 'overlapping facets surface 41 surface 41' pre-dates this change), 4.7 um at the same
    // racetrack corner of every e138 turn -- artefacts of MakePipeShell's framing that no
    // element size can discretise, while the analytic compound keeps the layout's true
    // clearances. The acceptance battery cannot see these creases (volume, BRepCheck, centring
    // all pass), so the choice cannot be made per-design after the fact. Rect WIRE still MUST
    // take the single-body sweep (it has no round-profile fallback); MVB_SINGLE_BODY=1 forces
    // the old behaviour for experiments.
    // Scope: CONCENTRIC columns only. Toroids keep their existing flow (simple-MakePipe single
    // body, else the conformal mitre compound) -- their junction geometry is poloidal, not the
    // shallow lead/wrap wedge measured here. Single-turn conductors also keep the single body
    // (their fuse is clean; no junction wedge exists).
    const bool multiTurnPath = !path.prims.empty() &&
                               path.prims.back().turnOrdinal != path.prims.front().turnOrdinal;
    const bool roundFemCompound = path.femReady && !path.isRectangular && !path.toroidal &&
                                  multiTurnPath && std::getenv("MVB_SINGLE_BODY") == nullptr;
    // FEM ROUND-WIRE windings (EVERY column type, toroids included): the CONFORMAL
    // MITRE ASSEMBLY -- each primitive its own exact analytic solid, tangent junctions
    // abutting on flush perpendicular discs, corner junctions grown and sliced on the
    // angle-bisector plane so neighbours share IDENTICAL elliptical faces. No fuse, no
    // weld, no overlap anywhere: geometric coincidence is exact by construction, which is
    // the input class OCCT's boolean spec and gmsh's fragment actually support (and what
    // gmsh >= 4.13's OCCBooleanGlue accelerates). This retires the whole strategy zoo
    // (whole-spine sweeps, pairwise welds, wrap overshoots) for round wire: the 2026-08
    // research review (FiQuS Pancake3D, OCCT boolean spec, Cubit imprint/merge) converged
    // on exactly this construction, and the toroid path had already proven it in-tree.
    // FEM RECT-WIRE windings on round/oblong columns: per-run compound + weld, like the round
    // path -- NOT the whole-spine single body. The single body's wrap->lead junction fillet
    // grazes the top wrap (measured on 03_buck: a 2.24-rad fragment-imprint band on the thin
    // side face whose boundary tetgen cannot recover -- edge 1/61 on curve 513). Per-run
    // sweeps end AT the crossings, the wrap overshoot supplies the weld lens, and the leads
    // are exact mitred prisms.
    const bool rectWireCompound = path.femReady && path.isRectangular && !path.toroidal &&
                                  path.singleBodyCapable && multiTurnPath &&
                                  std::getenv("MVB_SINGLE_BODY") == nullptr;
    if (rectWireCompound) {
        const bool diagR = std::getenv("MVB_DIAG") != nullptr;
        BRep_Builder bb; TopoDS_Compound comp; bb.MakeCompound(comp);
        bool okAll = true;
        double total = 0.0;
        auto vol = [](const TopoDS_Shape& sh) {
            GProp_GProps gp; BRepGProp::VolumeProperties(sh, gp); return gp.Mass();
        };
        // A lead CORNER (the quarter helix joining the radial run tangentially to the wrap --
        // Alf's 03 directive) is emitted inside the lead run but chains G1 to the wrap, so it
        // is swept WITH the wrap in the same pipe shell: the corner<->wrap seam then never
        // exists as a boolean junction. (Sweeping the corner separately and overshooting it
        // along the wrap was measured to fail: the extension hugs the wrap with coplanar
        // top/bottom faces, and OCC's fuse cannot weld that tangential lens.) The only weld
        // left at the corner is the prism's 0.35 r overshoot poking transversally into it.
        const auto runsAll = continuousRuns(path);
        std::vector<char> claimed(ptrs.size(), 0);
        auto leadRunAt = [&](size_t k) {
            return ptrs[runsAll[k].first]->isLead || ptrs[runsAll[k].first]->isConnection;
        };
        for (size_t k = 0; k < runsAll.size(); ++k) {
            if (leadRunAt(k)) continue;
            if (k > 0 && leadRunAt(k - 1))
                for (size_t i = runsAll[k - 1].second; i-- > runsAll[k - 1].first;) {
                    if (ptrs[i]->kind == Primitive::SPIRAL) claimed[i] = 1; else break;
                }
            if (k + 1 < runsAll.size() && leadRunAt(k + 1))
                for (size_t i = runsAll[k + 1].first; i < runsAll[k + 1].second; ++i) {
                    if (ptrs[i]->kind == Primitive::SPIRAL) claimed[i] = 1; else break;
                }
        }
        // Simple run-corner-wraps-corner-run path: with the lead corners the WHOLE path is
        // G1 (every junction tangent by construction), so sweep it as ONE pipe shell -- no
        // boolean junction anywhere. This is Alf's 03 directive taken to its conclusion:
        // the prism overshoot + proud height previously poked flaps out of the corner
        // ("weird joints"), and the buried end faces sat inside the winding OD. One sweep
        // has no seams: the corner flows into the connection, whose flat start face is the
        // corner's end section at the turns' outer diameter (+ the minimum-bend margin).
        const bool cornersBothEnds = runsAll.size() == 3 &&
                                     leadRunAt(0) && !leadRunAt(1) && leadRunAt(2) &&
                                     claimed[runsAll[0].second - 1] &&
                                     claimed[runsAll[2].first];
        if (cornersBothEnds) {
            try {
                BRepBuilderAPI_MakeWire wm;
                bool okW = true;
                for (size_t i = 0; i < ptrs.size(); ++i) {
                    TopoDS_Edge pe = primEdge(*ptrs[i], path.wireRadius);
                    if (pe.IsNull()) { okW = false; break; }
                    wm.Add(pe);
                    if (!wm.IsDone()) { okW = false; break; }
                }
                if (okW) {
                    auto pts = samplePrim(*ptrs[0], path.wireRadius);
                    gp_Dir t0(pts[1].XYZ() - pts[0].XYZ());
                    TopoDS_Shape whole = sweepWire(wm.Wire(), pts.front(), t0,
                                                   path.wireRadius, 0, /*rectangular=*/true,
                                                   path.wireWidth, path.wireHeight,
                                                   gp_Dir(0, 1, 0));
                    if (!whole.IsNull()) {
                        double spineLen = 0.0;
                        for (size_t i = 0; i < ptrs.size(); ++i) {
                            auto sp = samplePrim(*ptrs[i], path.wireRadius);
                            for (size_t j = 1; j < sp.size(); ++j)
                                spineLen += sp[j].Distance(sp[j - 1]);
                        }
                        const double expected = path.wireWidth * path.wireHeight * spineLen;
                        const double v = vol(whole);
                        int nsol = 0;
                        for (TopExp_Explorer e2(whole, TopAbs_SOLID); e2.More(); e2.Next())
                            ++nsol;
                        if (nsol == 1 && expected > 0.0 && v > 0.9 * expected &&
                            v < 1.1 * expected && BRepCheck_Analyzer(whole).IsValid()) {
                            if (diagR) std::cerr << "[rect-wire] G1 whole-path sweep: 1 solid, "
                                                 << "v=" << v * 1e9 << "mm3 (expected "
                                                 << expected * 1e9 << ")\n";
                            return whole;
                        }
                        if (diagR) std::cerr << "[rect-wire] G1 whole-path sweep rejected: "
                                             << "nsol=" << nsol << " v=" << v * 1e9
                                             << " expected=" << expected * 1e9 << "\n";
                    } else if (diagR) {
                        std::cerr << "[rect-wire] G1 whole-path sweep returned null\n";
                    }
                }
            } catch (const Standard_Failure& f) {
                if (diagR) std::cerr << "[rect-wire] G1 whole-path sweep threw: "
                                     << (f.GetMessageString() ? f.GetMessageString() : "(null)")
                                     << "\n";
            }
        }
        for (size_t k = 0; k < runsAll.size(); ++k) {
            const auto [b, e] = runsAll[k];
            const bool leadRun = leadRunAt(k);
            if (leadRun) {
                // straight L-route: one oriented prism per SEG, extended 0.35 r past interior
                // junction ends so consecutive prisms (and the wrap behind) overlap for the weld
                for (size_t i = b; i < e; ++i) {
                    if (claimed[i]) continue;   // corner helix -> swept with the wrap run
                    if (ptrs[i]->kind != Primitive::SEG) { okAll = false; break; }
                    gp_XYZ d = ptrs[i]->seg.b.XYZ() - ptrs[i]->seg.a.XYZ();
                    const double len = d.Modulus();
                    if (len < 1e-12) continue;
                    d /= len;
                    const double ovA = (i > b || b > 0) ? 0.35 * path.wireRadius : 0.0;
                    const double ovB = (i + 1 < e || e < path.prims.size()) ? 0.35 * path.wireRadius : 0.0;
                    gp_Pnt a2(ptrs[i]->seg.a.XYZ() - d * ovA);
                    try {
                        // 10% PROUD in height: the lead crosses the wrap with both rect
                        // sections axially aligned, so exact-height prisms meet the wrap in
                        // COPLANAR top/bottom faces -- the boolean's tangent T-junction line
                        // is exactly the edge tetgen cannot recover (measured on 03_buck:
                        // edge 3/3 on curve 34, at every ladder size down to 49 um). A lead
                        // 5% proud each side turns every intersection transversal; the extra
                        // copper is confined to the short lead runs.
                        TopoDS_Face prof = BRepBuilderAPI_MakeFace(
                            rectProfileWire(a2, gp_Dir(d), gp_Dir(0, 1, 0),
                                            path.wireWidth, 1.10 * path.wireHeight)).Face();
                        TopoDS_Shape prism =
                            BRepPrimAPI_MakePrism(prof, gp_Vec(d) * (len + ovA + ovB)).Shape();
                        if (prism.IsNull()) { okAll = false; break; }
                        bb.Add(comp, prism); total += vol(prism);
                    } catch (const Standard_Failure&) { okAll = false; break; }
                }
            } else {
                try {
                    // wrap wire = [claimed corners of the previous lead run] + wraps +
                    // [claimed corners of the next lead run], all G1-chained by construction
                    std::vector<size_t> idxs;
                    if (k > 0) {   // TRAILING claimed corners of the previous lead run
                        size_t lo = runsAll[k - 1].second;
                        while (lo > runsAll[k - 1].first && claimed[lo - 1]) --lo;
                        for (size_t i = lo; i < runsAll[k - 1].second; ++i) idxs.push_back(i);
                    }
                    for (size_t i = b; i < e; ++i) idxs.push_back(i);
                    if (k + 1 < runsAll.size())   // LEADING claimed corners of the next lead run
                        for (size_t i = runsAll[k + 1].first;
                             i < runsAll[k + 1].second && claimed[i]; ++i)
                            idxs.push_back(i);
                    BRepBuilderAPI_MakeWire wm;
                    for (size_t i : idxs) {
                        TopoDS_Edge pe = primEdge(*ptrs[i], path.wireRadius);
                        if (!pe.IsNull()) wm.Add(pe);
                        if (!wm.IsDone()) { okAll = false; break; }
                    }
                    if (!okAll || !wm.IsDone()) { okAll = false; break; }
                    auto pts = samplePrim(*ptrs[idxs.front()], path.wireRadius);
                    if (pts.size() < 2) { okAll = false; break; }
                    gp_Dir t0(pts[1].XYZ() - pts[0].XYZ());
                    TopoDS_Shape run = sweepWire(wm.Wire(), pts.front(), t0, path.wireRadius, 0,
                                                 /*rectangular=*/true, path.wireWidth,
                                                 path.wireHeight, gp_Dir(0, 1, 0));
                    if (run.IsNull()) { okAll = false; break; }
                    bb.Add(comp, run); total += vol(run);
                } catch (const Standard_Failure&) { okAll = false; break; }
            }
            if (!okAll) break;
        }
        if (okAll) {
            TopoDS_Shape welded = weldSolidsPairwise(comp, 1e-7);
            int nsol = 0;
            for (TopExp_Explorer e2(welded, TopAbs_SOLID); e2.More(); e2.Next()) ++nsol;
            if (diagR) std::cerr << "[rect-wire] per-run compound: " << nsol
                                 << " solid(s) after weld, v=" << vol(welded) * 1e9
                                 << "mm3 (runs total " << total * 1e9 << ")\n";
            if (diagR && nsol > 1) {
                int si = 0;
                for (TopExp_Explorer e2(welded, TopAbs_SOLID); e2.More(); e2.Next(), ++si) {
                    Bnd_Box bx; BRepBndLib::Add(e2.Current(), bx);
                    double x0,y0,z0,x1,y1,z1; bx.Get(x0,y0,z0,x1,y1,z1);
                    std::fprintf(stderr, "[rect-wire]   solid %d v=%.3fmm3 "
                                 "x[%.2f,%.2f] y[%.2f,%.2f] z[%.2f,%.2f]\n", si,
                                 vol(e2.Current())*1e9, x0*1e3,x1*1e3, y0*1e3,y1*1e3, z0*1e3,z1*1e3);
                }
            }
            return pruneDegenerateSolids(welded);
        }
        if (diagR) std::cerr << "[rect-wire] per-run compound failed -- falling back to the "
                                "whole-spine single body\n";
    }
    if (!skipSingleBody && !roundFemCompound &&
        (path.singleBodyCapable || rectWholeSweep) && (path.isRectangular || path.femReady)) {
        const bool diag = std::getenv("MVB_DIAG") != nullptr;
        TopoDS_Shape whole;
        TopoDS_Wire spine = buildFilletedWire(ptrs.data(), ptrs.size(), path.wireRadius);
        if (diag && spine.IsNull())
            std::cerr << "[single-body] spine NULL (untrimmable corner / edge build failed)\n";
        if (diag) std::cerr << "[single-body] wireRadius=" << path.wireRadius
                            << " condRadius=" << path.condRadius
                            << " wireWidth=" << path.wireWidth << " wireHeight=" << path.wireHeight
                            << " isRect=" << path.isRectangular << " prims=" << ptrs.size() << "\n";
        if (diag) {
            auto kn = [](int k) {
                return k == Primitive::SEG ? "SEG" : k == Primitive::ARC3 ? "ARC3"
                     : k == Primitive::SPIRAL ? "SPIRAL" : "BLEND";
            };
            for (size_t i = 0; i < ptrs.size(); ++i) {
                auto pts = samplePrim(*ptrs[i], path.wireRadius);
                gp_Pnt a = pts.front(), b = pts.back();
                double gap = -1.0;
                if (i + 1 < ptrs.size())
                    gap = b.Distance(samplePrim(*ptrs[i + 1], path.wireRadius).front());
                std::fprintf(stderr,
                             "[prim %2zu] %-6s a=(%9.4f,%9.4f,%9.4f) b=(%9.4f,%9.4f,%9.4f) "
                             "gap->next=%.6f mm  %s\n",
                             i, kn(ptrs[i]->kind), a.X() * 1e3, a.Y() * 1e3, a.Z() * 1e3,
                             b.X() * 1e3, b.Y() * 1e3, b.Z() * 1e3, gap * 1e3,
                             ptrs[i]->label.c_str());
            }
        }
        if (!spine.IsNull()) {
            auto p0pts = samplePrim(*ptrs[0], path.wireRadius);
            gp_Dir t0 = (p0pts.size() >= 2 &&
                         (p0pts[1].XYZ() - p0pts[0].XYZ()).Modulus() > 1e-12)
                            ? gp_Dir(p0pts[1].XYZ() - p0pts[0].XYZ())
                            : gp_Dir(1, 0, 0);
            whole = sweepWire(spine, p0pts.front(), t0, path.wireRadius, /*exact profile*/ 0,
                              path.isRectangular, path.wireWidth, path.wireHeight,
                              /*axialAxis=*/gp_Dir(0, 1, 0),
                              // Simple MakePipe for toroids AND the rect-column whole sweep: on
                              // both, MakePipeShell's corrected-Frenet framing drifts the section
                              // off the spine (measured 18 um on a 56-prim racetrack spine -- enough
                              // to close a 19.5 um inter-wrap clearance to 1.5 um and make the body
                              // unmeshable); MakePipe's framing keeps the section centred.
                              /*preferSimplePipe=*/path.toroidal,
                              /*tryFixedBinormal=*/rectWholeSweep);
        }
        if (!whole.IsNull()) {
            int nsol = 0;
            for (TopExp_Explorer e(whole, TopAbs_SOLID); e.More(); e.Next()) ++nsol;
            double spineLen = 0.0;
            for (size_t i = 0; i < ptrs.size(); ++i) {
                auto pts = samplePrim(*ptrs[i], path.wireRadius);
                for (size_t k = 1; k < pts.size(); ++k) spineLen += pts[k].Distance(pts[k - 1]);
            }
            // Swept copper cross-section area: rectangle width*height, else circle.
            double area = path.isRectangular ? path.wireWidth * path.wireHeight
                                             : kPi * path.wireRadius * path.wireRadius;
            double expected = area * spineLen;
            GProp_GProps gp;
            BRepGProp::VolumeProperties(whole, gp);
            double v = gp.Mass();
            bool okVol = expected > 0.0 && v > 0.9 * expected && v < 1.1 * expected;
            // The degenerate-sheet-face heuristic keys on the round wireRadius to catch boolean-fuse
            // slivers; a rectangular wire has legitimately thin flat faces and never fuses, so it is
            // guarded by BRepCheck validity + the exact volume match instead.
            bool okDeg = path.isRectangular || !hasDegenerateSheetFace(whole, path.wireRadius);
            bool okChk = BRepCheck_Analyzer(whole).IsValid();
            // Toroid sweeps additionally must keep the section centred on the spine: a valid, right-
            // volume body whose section drifted off an interior crossing is geometrically WRONG.
            // Sample the wrap (non-lead) centrelines -- lead tips sit at the end caps (~0 clearance)
            // and would false-reject -- and require every one to clear 0.85*wireRadius.
            bool okCentre = true;
            if (path.toroidal || rectWholeSweep) {
                std::vector<gp_Pnt> interior;
                for (size_t i = 0; i < ptrs.size(); ++i) {
                    if (ptrs[i]->isLead) continue;
                    auto pts = samplePrim(*ptrs[i], path.wireRadius);
                    interior.insert(interior.end(), pts.begin(), pts.end());
                }
                okCentre = centrelineStaysCentred(whole, interior, path.wireRadius);
            }
            if (nsol == 1 && okVol && okDeg && okChk && okCentre) {
                if (diag) std::cerr << "[single-body] ACCEPTED v=" << v << " exp=" << expected
                                    << "\n";
                return whole;
            }
            if (diag) std::cerr << "[single-body] rejected nsol=" << nsol << " okVol=" << okVol
                                << " (v=" << v << " exp=" << expected << ") okDeg=" << okDeg
                                << " okChk=" << okChk << " okCentre=" << okCentre << "\n";
        } else if (diag && !spine.IsNull()) {
            std::cerr << "[single-body] sweep NULL (MakePipeShell failed on the G1 spine)\n";
        }
        // The rect-column whole-sweep attempt did not produce an acceptable body: take the exact
        // per-primitive analytic path it would otherwise have taken. Nothing is lost by trying.
        if (rectWholeSweep) return emitRectColumn(path);
    }

    // The per-run compound below sweeps ROUND profiles; a rectangular winding must come out of the
    // single-body path or not at all -- never a silently round-profiled body. Surface the failure.
    if (path.isRectangular) {
        throw std::runtime_error(
            "ConductorBuilder: rectangular-wire single-body sweep failed for '" + path.name +
            "' (fixed-binormal MakePipeShell did not close a valid solid); refusing to fall back to "
            "a round-profile compound");
    }

    // FEM: a conformal (mitre-jointed) toroid compound -- a valid, non-overlapping, mesh-shareable
    // assembly for the dense toroids where the single-solid sweep can't close. Drawing (femReady
    // false) skips it and keeps the fast overlapping per-run compound below.
    if (path.toroidal && path.femReady) {
        return pruneDegenerateSolids(emitToroidConformal(ptrs, path.wireRadius, wirePolygonSegments));
    }

    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    size_t solidIndex = 0;
    for (auto [b, e] : continuousRuns(path)) {
        bool closedRing = (e - b == 1) && ptrs[b]->kind != Primitive::SEG &&
                          [&] {
                              auto [pa, pb] = primEndpoints(*ptrs[b]);
                              return pa.Distance(pb) < 1e-9;
                          }();
        // runs never mix lead/link and wrap primitives (continuousRuns splits at both)
        bool leadRun = ptrs[b]->isLead || ptrs[b]->isConnection;
        TopoDS_Shape run;
        if (std::getenv("MVB_DEBUG_RUNS"))
            std::cerr << "[sweep] " << path.name << " run[" << b << "," << e << ") first='"
                      << ptrs[b]->label << "' last='" << ptrs[e - 1]->label << "'"
                      << (leadRun ? " LEAD" : "") << (closedRing ? " RING" : "") << std::endl;
        if (!closedRing && !leadRun) {
            run = sweepRun(ptrs.data() + b, e - b, path.wireRadius, wirePolygonSegments);
        }
        if (std::getenv("MVB_DEBUG_RUNS") && run.IsNull() && !leadRun) {
            std::cerr << "PIPE FAILED for run of " << (e - b) << " prims starting ["
                      << ptrs[b]->label << "]\n";
        }
        if (run.IsNull() && !leadRun && !closedRing) {
            // The whole-run pipe failed: retry per electrical turn (short spines pipe
            // reliably) before falling all the way to per-primitive. Keeps the winding a
            // handful of weldable pipe chunks instead of hundreds of loose primitives.
            run = sweepRunChunked(ptrs.data() + b, e - b, path.wireRadius,
                                  wirePolygonSegments, flatCaps);
        }
        if (run.IsNull()) {
            // Per-primitive sweeps for this span (arcs revolve exactly).
            run = sweepPiecewise(ptrs.data() + b, e - b, path.wireRadius,
                                 wirePolygonSegments, flatCaps);
        }
        if (std::getenv("MVB_DEBUG_RUNS") && !run.IsNull()) {
            Bnd_Box bb;
            // Optimal bounds: the default Bnd_Box of a swept BSpline surface is its
            // control-point hull, up to ~2.5x larger than the actual surface.
            BRepBndLib::AddOptimal(run, bb, Standard_False, Standard_False);
            double x0, y0, z0, x1, y1, z1;
            bb.Get(x0, y0, z0, x1, y1, z1);
            Bnd_Box sb;
            for (size_t k = b; k < e; ++k)
                for (const auto& q : samplePrim(*ptrs[k], path.wireRadius)) sb.Add(q);
            double sx0, sy0, sz0, sx1, sy1, sz1;
            sb.Get(sx0, sy0, sz0, sx1, sy1, sz1);
            double excess = std::max({sx0 - x0, sy0 - y0, sz0 - z0,
                                      x1 - sx1, y1 - sy1, z1 - sz1}) - path.wireRadius;
            if (excess > 1e-4) {
                std::cerr << "RUN FLARE [" << b << "," << e << ") first=["
                          << ptrs[b]->label << "] last=[" << ptrs[e - 1]->label
                          << "] excess " << excess << " m\n";
            }
        }
        // MVB_SOLID_DUMP: map the VIEWER's solid index back to the primitives that made it.
        // The drawing compound emits one solid per continuous RUN, not per primitive, so
        // "<conductor> 55" in a CAD tree is run 55 -- unmappable without this.
        if (!run.IsNull() && std::getenv("MVB_SOLID_DUMP")) {
            // Explode the run: a piecewise sweep returns a COMPOUND of one solid per primitive,
            // and a CAD tree lists those individually. Numbering them here in emission order is
            // what makes "<conductor> 58" in a viewer resolvable to a label.
            size_t k = 0;
            for (TopExp_Explorer ex(run, TopAbs_SOLID); ex.More(); ex.Next(), ++k) {
                Bnd_Box bb;
                BRepBndLib::AddOptimal(ex.Current(), bb, Standard_False, Standard_False);
                double x0, y0, z0, x1, y1, z1;
                bb.Get(x0, y0, z0, x1, y1, z1);
                // Attribute by GEOMETRY, never by position: a piecewise run emits a cylinder
                // AND sphere elbows per primitive, so solid k is not primitive b+k. Match the
                // solid's bbox centre to the nearest primitive centreline.
                const gp_Pnt mid(0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1));
                const Primitive* src = ptrs[b];
                double bestD = 1e30;
                for (size_t q = b; q < e; ++q) {
                    for (const auto& sp : samplePrim(*ptrs[q], path.wireRadius)) {
                        const double d = sp.Distance(mid);
                        if (d < bestD) { bestD = d; src = ptrs[q]; }
                    }
                }
                GProp_GProps gp;
                BRepGProp::VolumeProperties(ex.Current(), gp);
                std::cerr << "[solid] " << path.name << " " << solidIndex++
                          << " vol_mm3=" << gp.Mass() * 1e9
                          << " run[" << b << "," << e << ") '" << src->label
                          << "' bbox_mm x=[" << x0 * 1e3 << "," << x1 * 1e3
                          << "] y=[" << y0 * 1e3 << "," << y1 * 1e3
                          << "] z=[" << z0 * 1e3 << "," << z1 * 1e3 << "]\n";
            }
        }
        if (!run.IsNull()) builder.Add(compound, run);
    }
    // A real multi-turn winding's swept pipes always leave degenerate seam faces when the
    // boolean unions them (the fuse-guard proves this and rejects the result), so running the
    // costly BOP just to fall back to the compound is wasted work -- minutes on a 32-turn coil.
    // Skip straight to the exact clean compound when the conductor spans more than one turn;
    // fuse only simple single-turn conductors, where the union is both cheap and clean.
    bool multiTurn = !path.prims.empty() &&
                     path.prims.back().turnOrdinal != path.prims.front().turnOrdinal;
    // FEM geometry becomes ONE solid: single-turn conductors via the (cheap, clean) N-way fuse;
    // multi-turn conductors via the PAIRWISE coverage-guarded weld -- the N-way fuse of a whole
    // winding either leaves seam slivers or silently drops lead pieces, but welding the runs one
    // piece at a time with each step verified builds the same single body robustly (the station-
    // overlap leads guarantee real volumetric overlap at every junction). Drawing keeps the fast
    // compound.
    TopoDS_Shape result = compound;
    if (path.femReady) {
        result = multiTurn ? weldSolidsPairwise(compound, 1e-7)
                           : fuseAllSolids(compound, path.wireRadius);
    }
    result = pruneDegenerateSolids(result);

    // MVB_COVERAGE_CHECK: is the conductor FULLY CONDUCTING? Walk every primitive's centreline
    // and classify each sample against the emitted solids. A sample that lies in no solid is a
    // gap in the copper -- a run whose sweep failed and whose fallback also dropped it, a piece
    // pruned as degenerate, or a chunk the pipe silently lost. Nothing else in the pipeline
    // notices: the compound is still "valid", it is just missing wire.
    if (std::getenv("MVB_COVERAGE_CHECK")) {
        // BRepClass3d_SolidClassifier is non-copyable (its explorer is), so hold them by pointer.
        std::vector<std::unique_ptr<BRepClass3d_SolidClassifier>> cls;
        for (TopExp_Explorer ex(result, TopAbs_SOLID); ex.More(); ex.Next())
            cls.push_back(std::make_unique<BRepClass3d_SolidClassifier>(ex.Current()));
        size_t miss = 0, total = 0;
        std::string worst;
        for (const auto& pr : path.prims) {
            for (const auto& q : samplePrim(pr, path.wireRadius)) {
                ++total;
                bool in = false;
                for (auto& c : cls) {
                    c->Perform(q, 1e-9);
                    if (c->State() == TopAbs_IN || c->State() == TopAbs_ON) { in = true; break; }
                }
                if (!in) {
                    if (!miss)
                        worst = pr.label + " at (" + std::to_string(q.X() * 1e3) + "," +
                                std::to_string(q.Y() * 1e3) + "," + std::to_string(q.Z() * 1e3) + ") mm";
                    ++miss;
                }
            }
        }
        std::cerr << "[coverage] " << path.name << ": " << (miss ? "GAP" : "ok") << " -- " << miss
                  << "/" << total << " centreline samples outside the copper"
                  << (miss ? ("; first: " + worst) : "") << "\n";
    }
    return result;
}

// ---------------------------------------------------------------------------------------
// Connection replay: MKF's drawn ConnectionReservedSpace rectangles (layer == "") are the
// authoritative routes — the pink (terminal) and blue (link) boxes of the Painter SVG.
using RSpace = OpenMagnetics::ConnectionReservedSpace;

struct PlanePt {
    double x = 0, y = 0;   // MKF 2D window coordinates (x = radial, y = axial)
};

bool rectIsVertical(const RSpace& s) { return s.dimensions.at(1) > s.dimensions.at(0); }

// Terminal waypoints from MKF's terminal rect group (1 rect = leave at own level;
// stub + run = L-route along the window edge). `station` = the connecting turn.
std::vector<PlanePt> terminalWaypoints(const std::vector<const RSpace*>& group,
                                       const PlanePt& station, const std::string& who) {
    if (group.empty()) {
        throw std::runtime_error("ConductorBuilder: no drawn terminal lead for " + who +
                                 " in MKF's connection reserved spaces");
    }
    const RSpace* run = nullptr;
    for (const RSpace* s : group) {
        if (!rectIsVertical(*s)) run = s;
    }
    if (!run) {
        throw std::runtime_error("ConductorBuilder: terminal lead group for " + who +
                                 " has no horizontal run rect");
    }
    // The run rect spans [turnX - w/2, borderX + w/2], symmetric about its centre, so
    // borderX = 2*centre.x - turnX; the run's level is the rect's own y (with no stub
    // MKF still routes at the blocked edge slot, within half a wire of the station) —
    // all MKF data, nothing invented.
    // Border = the drawn run box's own FAR EDGE minus half a wire (box height = one wire
    // OD). Station-independent -- the old mirror (2*centre - station) assumed the station
    // at the box's near edge and INVERTED the lead when a dragback-displaced attach sat
    // past the centre (17_cllc secondary exit: attach ridden +5.3 mm -> "border" landed
    // inside the winding and the exit lead was mangled into the ring band). For an
    // undisplaced station the two derivations are algebraically identical. A lead never
    // runs inward: the border is clamped outward of the attach.
    double borderX = run->coordinates.at(0) + 0.5 * run->dimensions.at(0) -
                     0.5 * run->dimensions.at(1);
    borderX = std::max(borderX, station.x);
    if (std::getenv("MVB_DIAG"))
        std::cerr << "[terminalWaypoints] " << who << ": station=(" << station.x * 1e3 << ","
                  << station.y * 1e3 << ") run rect c=(" << run->coordinates.at(0) * 1e3 << ","
                  << run->coordinates.at(1) * 1e3 << ") dims=(" << run->dimensions.at(0) * 1e3
                  << "x" << run->dimensions.at(1) * 1e3 << ") -> borderX=" << borderX * 1e3
                  << "\n";
    double edgeY = run->coordinates.at(1);
    // SVG semantics (Alf, 2026-08-07): a pink box that COVERS the wire surface at the
    // turn is NOT a segment -- it only marks the connection of the turn to another
    // segment (03: the horizontal run). So a vertical connection exists only when the
    // attach row sits genuinely OUTSIDE the run box's own height band; inside it, the
    // route is straight at the turn's row and any edgeY-station.y offset is layout
    // noise (03_buck: a sub-micron offset built a degenerate 3-point route; 06_llc:
    // a 0.13 mm one collapsed into a diagonal run). Threshold = MKF's own drawn box.
    if (std::abs(edgeY - station.y) <= 0.5 * run->dimensions.at(1)) {
        return {{station.x, station.y}, {borderX, station.y}};
    }
    return {{station.x, station.y}, {station.x, edgeY}, {borderX, edgeY}};
}

// Split one conductor's drawn terminal rects (MKF emission order: the ENTRANCE lead's
// rects first, then the EXIT lead's, each as [vertical stub?, horizontal run] -- the run
// closes its group) into the entrance and exit lead groups. With exactly two runs, every
// stub pairs to the run whose centre is nearest: pure emission order misgroups whenever
// MKF draws a route run-first instead of stub-first (09_planar's CONTIGUOUS secondary:
// [H V H V] split into three sequential groups). ONE drawn lead is still usable: it
// becomes the entrance and the exit group returns EMPTY -- the caller synthesizes the
// minimal straight-out exit, LOUDLY, because that is builder-invented routing and the
// missing lead is an MKF gap to fix upstream (seen on 13_current_sense_er95).
std::pair<std::vector<const RSpace*>, std::vector<const RSpace*>>
splitTerminalGroups(const std::vector<const RSpace*>& terminalRects, const std::string& who) {
    std::vector<std::vector<const RSpace*>> groups;
    std::vector<const RSpace*> hs;
    for (const RSpace* s : terminalRects)
        if (!rectIsVertical(*s)) hs.push_back(s);
    if (hs.size() == 2) {
        groups.assign(2, {});
        groups[0].push_back(hs[0]);
        groups[1].push_back(hs[1]);
        for (const RSpace* s : terminalRects) {
            if (!rectIsVertical(*s)) continue;
            double d[2];
            for (int g = 0; g < 2; ++g)
                d[g] = std::hypot(s->coordinates.at(0) - hs[g]->coordinates.at(0),
                                  s->coordinates.at(1) - hs[g]->coordinates.at(1));
            groups[d[0] <= d[1] ? 0 : 1].push_back(s);
        }
    } else {
        groups.emplace_back();
        for (const RSpace* s : terminalRects) {
            groups.back().push_back(s);
            if (!rectIsVertical(*s)) groups.emplace_back();   // run closes the group
        }
        if (!groups.empty() && groups.back().empty()) groups.pop_back();
    }
    if (groups.size() == 1) {
        std::cerr << "[ConductorBuilder] " << who << ": MKF drew only ONE terminal lead; "
                     "using it as the entrance and synthesizing a straight-out exit at the "
                     "last turn's level\n";
        return {groups[0], {}};
    }
    if (groups.size() != 2) {
        std::string dump;
        for (const RSpace* s : terminalRects) {
            char buf[128];
            std::snprintf(buf, sizeof buf, " [%.3g x %.3g at (%.3g, %.3g) %s]",
                          s->dimensions.at(0) * 1e3, s->dimensions.at(1) * 1e3,
                          s->coordinates.at(0) * 1e3, s->coordinates.at(1) * 1e3,
                          rectIsVertical(*s) ? "V" : "H");
            dump += buf;
        }
        throw std::runtime_error(
            "ConductorBuilder: expected 2 terminal lead groups (entrance, exit) for " + who +
            ", got " + std::to_string(groups.size()) +
            ". Terminal rects (w x h at (x,y), V=vertical):" + dump);
    }
    return {groups[0], groups[1]};
}

// ---------------------------------------------------------------------------------------
// Wrap planners. Each appends the connecting geometry of ONE wrap — crossing k-1 to
// crossing k — for its column shape. Under MKF's real-winding model the turnsDescription
// holds the N+1 window crossings of an N-turn winding (the first entry is the beginning
// of the first turn), so no closed loops exist anywhere.

// One placed turn ring of ANY conductor: centre (r, y) in the window half-plane and its
// CONDUCTING (copper) radius. The Z-transition end-run planner routes around these.
struct RingInv {
    double r = 0, y = 0, rw = 0;
};

// Is the move from crossing s to crossing n a Z-RETURN (needs the end-run route) rather than
// an ordinary helix wrap or a serpentine U layer-link? A U-link is a predominantly-radial
// step (handled separately). A LAYER CHANGE that also travels more than one wrap pitch
// axially is a return (the serpentine link's axial move stays under a pitch); a SAME-radius
// jump (section-interleaved order) is a return only when it far exceeds the conductor's own
// median pitch — sparse windings legitimately advance several wire-ODs per ordinary wrap.
// zOrderAdvance: the conductor's SIGNED within-layer advance when every same-radius step moves
// the same way (a z-order winding); 0.0 when the direction alternates (classic serpentine).
// Within-layer pitch statistics PER RADIUS BAND. A conductor's layers do not share a pitch:
// a dense inner layer packs at one wire OD while a sparse outer layer legitimately advances
// several ODs per turn (fence-post SPREAD over a half-full layer). Judging every transition
// against the CONDUCTOR's median pitch therefore misreads the sparse layer's ordinary wraps as
// Z-returns and emits vertical descents where the wire should wrap (measured on the litz
// fixture: layer 1's 4 turns, pitch 2.55 mm against layer 0's 0.69 mm median, came out as three
// straight drops instead of three rings). Each band is judged against its OWN pitch.
struct PitchBand {
    double radius = 0.0;
    double medianPitch = 0.0;
    double advance = 0.0;   // signed, 0 when the layer's direction alternates
};
std::vector<PitchBand> computePitchBands(const std::vector<PlanePt>& stations, double wireRadius) {
    std::vector<std::pair<double, std::vector<double>>> raw;   // (band radius, |dy| samples)
    std::vector<std::pair<int, int>> signs;                    // (pos, neg) per band
    for (size_t i = 0; i + 1 < stations.size(); ++i) {
        const PlanePt& a = stations[i];
        const PlanePt& b = stations[i + 1];
        if (std::abs(b.x - a.x) > wireRadius) continue;   // layer change, not a within-layer step
        size_t band = raw.size();
        for (size_t k = 0; k < raw.size(); ++k) {
            if (std::abs(raw[k].first - a.x) <= wireRadius) { band = k; break; }
        }
        if (band == raw.size()) { raw.push_back({a.x, {}}); signs.push_back({0, 0}); }
        const double dy = b.y - a.y;
        raw[band].second.push_back(std::abs(dy));
        (dy >= 0.0 ? signs[band].first : signs[band].second)++;
    }
    std::vector<PitchBand> bands;
    for (size_t k = 0; k < raw.size(); ++k) {
        auto samples = raw[k].second;
        // A band needs at least TWO within-layer steps before its median means anything: with
        // a single sample the band adopts that very step as its own "normal pitch", so any
        // jump -- including a genuine dragback -- looks ordinary (measured on 14_dab, whose
        // outer secondary layer holds two turns 15 mm apart: the return came out as a 15 mm
        // rise helix). One-sample bands fall back to the conductor-wide median.
        if (samples.size() < 2) continue;
        std::sort(samples.begin(), samples.end());
        PitchBand band;
        band.radius = raw[k].first;
        band.medianPitch = samples[samples.size() / 2];
        if (signs[k].first == 0 || signs[k].second == 0)
            band.advance = (signs[k].first > 0 ? band.medianPitch : -band.medianPitch);
        bands.push_back(band);
    }
    return bands;
}
// The band a station belongs to; falls back to the whole-conductor values when a radius has no
// within-layer step of its own (a single-turn layer).
PitchBand bandAt(const std::vector<PitchBand>& bands, double radius, double wireRadius,
                 double fallbackPitch, double fallbackAdvance) {
    for (const auto& band : bands)
        if (std::abs(band.radius - radius) <= wireRadius) return band;
    return {radius, fallbackPitch, fallbackAdvance};
}

bool isZReturn(const PlanePt& s, const PlanePt& n, double wireRadius, double medianPitch,
               double zOrderAdvance = 0.0) {
    if (std::abs(n.x - s.x) > wireRadius && std::abs(n.y - s.y) <= std::abs(n.x - s.x))
        return false;   // serpentine U-link
    const bool layerChange = std::abs(n.x - s.x) > wireRadius;
    // A layer change moving AGAINST a consistently-advancing winding's own direction is a
    // return no matter how short: the diagonal cone it would otherwise get sweeps backwards
    // through the band the sibling parallels occupy (23_interleaved: s0's 1.2 mm back-move
    // crossed s1's helix 0.38 mm off centre).
    if (layerChange && zOrderAdvance != 0.0 &&
        (n.y - s.y) * (zOrderAdvance > 0.0 ? 1.0 : -1.0) < -2.0 * wireRadius)
        return true;
    const double thr = layerChange ? std::max(medianPitch, 2.0 * wireRadius)
                                   : std::max(6.0 * wireRadius, 1.6 * medianPitch);
    return std::abs(n.y - s.y) > thr;
}

// ROUND column: within a layer, one full 360-degree cylindrical spiral about the column
// axis. A SERPENTINE (U) layer transition connects the same end of two adjacent layers, so
// the two crossings share the connection-plane azimuth and axial height but sit at different
// radii — the wire does not wrap around again, it steps straight out radially to the next
// layer. That radial step is exactly MKF's blue inter-layer link box: a short straight wire
// leaving along -Z until it reaches the layer it connects to, rather than a full conical
// wrap around the whole top of the winding. Detect it as a predominantly-radial move (radius
// change over one wire radius, and larger than the axial change).
//
// A Z-style transition instead flies back from the end of one layer to the FAR end of the
// next. Modelling it as a diagonal cone (linear r/y over one revolution) INTERPENETRATES the
// layer it climbs over: with layers exactly one wire-OD apart, a diagonal at intermediate
// radius passes closer to the inner layer's rings than the copper sum wherever its axial
// offset from a ring is under the pitch (provably infeasible in-between — the clearance band
// is narrower than the pitch for any real winding). Since turn positions are never moved,
// the only exact-clearance realization is the END-RUN: leave the winding past its axial end,
// travel in the free annulus OUTSIDE every placed ring, and re-enter the target layer from
// beyond the far end. That is also how a real crossover is routed when the turns cannot
// locally bulge. corridorIdx staggers the end-run lanes of a winding's parallel strands
// (which transition in lockstep and would otherwise coincide).
// A dragback the wrap must ride OVER: its azimuth, and the height the wire is lifted there.
struct WrapBump {
    double azimuth = 0.0;   // the dragback's azimuth; the raised half-turn is centred on it
    // How far the wire is RAISED, as a distance. Not a count of dragbacks: the wire that laid
    // the dragback sets the clearance needed, and different windings use different wire.
    double distance = 0.0;
};

// The tallest azimuth COLUMN in a bump list. Dragbacks SHARING an azimuth stack radially --
// their clearances add (a return lies in its destination layer's space, the step over it
// displaces the next layer, and so on outward) -- while returns in DIFFERENT fan slots sit
// side by side, so a wire riding over the whole fan only needs to clear the tallest stack.
// Every junction raise (wrap ends, dragback chain ends, the exit lead's lift) MUST come from
// this one helper on the same bump list, so abutting pieces land on exactly the same offset.
std::pair<double, double> tallestBumpColumn(const std::vector<WrapBump>& bumps) {
    std::vector<std::pair<double, double>> cols;   // (azimuth, summed distance)
    for (const auto& bmp : bumps) {
        bool found = false;
        for (auto& c : cols) {
            if (std::abs(std::remainder(bmp.azimuth - c.first, kTwoPi)) < 1e-9) {
                c.second += bmp.distance;
                found = true;
                break;
            }
        }
        if (!found) cols.push_back({bmp.azimuth, bmp.distance});
    }
    double raise = 0.0, az = 0.0;
    for (const auto& c : cols)
        if (c.second > raise) { raise = c.second; az = c.first; }
    return {raise, az};
}

// Emit a helical sweep (r,y linear in azimuth) that RIDES OVER the given dragbacks: plain helix,
// a cosine-blended ramp out, a blended ramp back in, plain helix. The blend's end tangents are
// purely azimuthal, so the pieces meet tangentially and the conformal assembler sees no corner.
void appendBumpedSweep(ConductorPath& path, double r0, double y0, double azStart, double r1,
                       double y1, double azEnd, const std::vector<WrapBump>& bumps,
                       double wireRadius, const std::string& label, size_t ordinal,
                       bool isConnection) {
    auto radiusAt = [&](double az) {
        const double t = (az - azStart) / (azEnd - azStart);
        return r0 + (r1 - r0) * t;
    };
    auto heightAt = [&](double az) {
        const double t = (az - azStart) / (azEnd - azStart);
        return y0 + (y1 - y0) * t;
    };
    auto pushPiece = [&](double a0, double a1, double ra, double rb, double ya, double yb,
                         bool blend, const char* suffix, double cz = 0.0) {
        if (a1 - a0 < 1e-12) return;
        Primitive pr;
        // A piece at constant radius AND constant height is a plain circular arc, not a helix:
        // the spiral machinery degenerates on it (zero pitch leaves a zero-length line in the
        // cylinder's parameter space, and the pipe sweep then fails outright). The dragback's
        // run-out and run-in are exactly this case whenever no bump falls inside them.
        if (std::abs(rb - ra) < 1e-12 && std::abs(yb - ya) < 1e-12) {
            pr.kind = Primitive::ARC3;
            pr.arc.c = gp_Pnt(0, ya, cz);
            pr.arc.axis = gp_XYZ(0, 1, 0);
            pr.arc.v0 = azPointC(0, cz, ra, ya, a0).XYZ() - pr.arc.c.XYZ();
            pr.arc.sweep = a1 - a0;
        }
        else {
            pr.kind = Primitive::SPIRAL;
            pr.spiral = {0, cz, ra, ya, a0, rb, yb, a1};
            pr.spiral.blend = blend;
        }
        pr.label = label + suffix;
        pr.turnOrdinal = ordinal;
        pr.isConnection = isConnection;
        path.prims.push_back(std::move(pr));
    };
    // The verticals live in a TIGHT FAN of a few degrees around the station plane (each vertical
    // takes the feasible azimuth closest to the plane, see the vertical-fan pre-scan), so ONE
    // raised treatment centred on the plane clears every one of them: same-azimuth columns
    // stack (clearances add), columns side by side do not, and the ride-over only needs the
    // tallest stack.
    const double raise = tallestBumpColumn(bumps).first;
    if (raise <= 0.0) {
        // A full revolution at CONSTANT radius and height (the first turn of a U layer -- see
        // appendRoundWrap) is a CLOSED circle: swept whole it makes a periodic surface with a
        // seam, which the conformal assembler cannot build ("could not sweep primitive"). Two
        // half-arcs sweep as ordinary patches and meet tangentially. Bumped rings are already
        // split into head/level/tail below, so only the unbumped case needs this.
        if (azEnd - azStart > kPi + 1e-9 && std::abs(r1 - r0) < 1e-12 &&
            std::abs(y1 - y0) < 1e-12) {
            const double azMid = 0.5 * (azStart + azEnd);
            pushPiece(azStart, azMid, r0, r0, y0, y0, false, "");
            pushPiece(azMid, azEnd, r0, r0, y0, y0, false, "");
            return;
        }
        pushPiece(azStart, azEnd, r0, r1, y0, y1, false, "");
        return;
    }
    // ALF'S RULE: the bump affects the -Z HALF ONLY. The +Z half is always the plain
    // 180-degree arc of the wrap's own circle, from XY plane to XY plane. The raised -Z half
    // is THE SAME CIRCLE, at the wrap's own radius, translated bodily by -raise in z (centre
    // cz = -raise) -- and the two halves are bridged at the XY-plane crossings by a STRAIGHT
    // SEGMENT whose length IS the bump dimension at that turn (Alf, 2026-08-07).
    //
    // Why the straight riser and not a radius-grown circle (the dropped-centre rule this
    // replaces, radius hypot(r, raise), which met the plane at (+-r, 0) with no bridge):
    //   * the wire keeps its EXACT radius all the way round, so it stays concentric with the
    //     column instead of bulging outward by hypot(r,raise) - r over the whole raised half;
    //   * the dip at the station side is EXACTLY the raise, not raise + (hypot(r,raise) - r);
    //   * the junctions are exactly TANGENT: at x = -+r the arc's tangent is purely axial
    //     (d/daz of (r cos az, y, -r sin az) is (0, 0, -+r) there), which is the riser's own
    //     direction -- so this model is G1 where the grown-radius one left a raise/r kink for
    //     the assembler to mitre.
    const double c = raise;
    // A sweep short enough to live entirely inside the -Z half is emitted raised end to end:
    // both its ends are station-side, already inside the raised region, so no riser is due.
    if (azEnd - azStart <= kPi + 1e-9) {
        pushPiece(azStart, azEnd, radiusAt(azStart), radiusAt(azEnd), y0, y1, false,
                  " (over dragback)", -c);
        return;
    }
    // Full wrap: raised head [station .. XY plane at x=-r], level +Z half [exactly pi..2pi in
    // the wrap's own frame], raised tail [XY plane at x=+r .. station]. The boundaries are the
    // GLOBAL XY-plane crossings, never offsets from this wrap's own span, so every layer's
    // raised half covers the same azimuths.
    double aHead = kPlaneAz + kPi / 2.0;          // z = 0, x = -r
    double aTail = kPlaneAz + kTwoPi - kPi / 2.0; // z = 0, x = +r
    while (aHead <= azStart + 1e-9) aHead += kTwoPi;
    while (aTail >= azEnd - 1e-9) aTail -= kTwoPi;
    if (!(aHead < aTail)) {   // fan wider than the half-turn between the plane crossings
        throw std::runtime_error("ConductorBuilder: bump region boundaries inverted for '" +
                                 label + "' (vertical fan too wide for the ride-over model)");
    }
    // The riser at a plane crossing: straight along z, from the raised half's endpoint to the
    // level half's endpoint (or back), length exactly `raise`.
    auto pushRiser = [&](double az, double rAt, double yAt, bool up) {
        const gp_Pnt raised = azPointC(0, -c, rAt, yAt, az);
        const gp_Pnt level = azPointC(0, 0.0, rAt, yAt, az);
        Primitive pr;
        pr.kind = Primitive::SEG;
        pr.seg = up ? Seg{raised, level} : Seg{level, raised};
        pr.label = label + " (bump riser)";
        pr.turnOrdinal = ordinal;
        pr.isConnection = isConnection;
        path.prims.push_back(std::move(pr));
    };
    pushPiece(azStart, aHead, radiusAt(azStart), radiusAt(aHead), heightAt(azStart),
              heightAt(aHead), false, " (over dragback)", -c);
    pushRiser(aHead, radiusAt(aHead), heightAt(aHead), /*up=*/true);
    pushPiece(aHead, aTail, radiusAt(aHead), radiusAt(aTail), heightAt(aHead), heightAt(aTail),
              false, "");
    pushRiser(aTail, radiusAt(aTail), heightAt(aTail), /*up=*/false);
    pushPiece(aTail, azEnd, radiusAt(aTail), radiusAt(azEnd), heightAt(aTail), y1, false,
              " (over dragback)", -c);
}

// azS / azE are the CROSSING azimuths at the wrap's two stations: kPlaneAz normally, a fan
// slot when the crossing belongs to a dragback or terminal lead spread off the plane.
void appendRoundWrap(ConductorPath& path, const PlanePt& s, const PlanePt& n,
                     double wireRadius, const std::string& label, size_t ordinal,
                     const std::vector<WrapBump>& bumps = {}, double azS = kPlaneAz,
                     double azE = kPlaneAz, const std::vector<WrapBump>& bumpsEnd = {}) {
    // U (SERPENTINE) LAYER LINK -- Alf, 2026-08-07, 14_dab. Layers wound in U (this one bottom
    // to top, the next top to bottom) connect DIFFERENTLY from a dragback, and differently from
    // a cone: the wire leaves the last turn's crossing on a straight HORIZONTAL segment, radial
    // and at CONSTANT HEIGHT, until it reaches the next layer's radius; from there it continues
    // as a NORMAL round turn -- that first turn of the new layer also at constant height. The
    // height only begins to change at the layer's SECOND turn.
    //
    // So the transition carries the segment AND one full revolution: it is a turn, exactly as
    // MKF counts it. Emitting only the segment dropped a whole turn of copper (and of
    // inductance) per layer change; emitting a cone over the revolution instead ("what the fuck
    // is this monstrosity") gave the right turn count with the wrong shape.
    if (std::abs(n.x - s.x) > wireRadius && std::abs(n.y - s.y) <= std::abs(n.x - s.x)) {
        const double rsRaise = tallestBumpColumn(bumps).first;
        const double reRaise = tallestBumpColumn(bumpsEnd).first;
        // The segment is TANGENTIAL to the arc it leaves -- at the crossing the turn's own
        // tangent is the -X direction -- and runs just far enough for the tangent to reach the
        // next layer's radius. A tangent from radius r0 meets radius r1 after L = sqrt(r1^2 -
        // r0^2), having advanced dAz = atan2(L, r0) in azimuth; the wire leaves its arc
        // smoothly (the junction is exactly tangent) instead of stepping sideways out of it,
        // which a radial segment does.
        const double L = std::sqrt(std::max(0.0, n.x * n.x - s.x * s.x));
        const double dAz = std::atan2(L, s.x);
        Primitive step;
        step.kind = Primitive::SEG;
        step.seg = {azPointC(0, -rsRaise, s.x, s.y, azS),
                    azPointC(0, -reRaise, n.x, s.y, azS + dAz)};
        step.label = label + " (layer link)";
        step.turnOrdinal = ordinal;
        step.isConnection = true;
        path.prims.push_back(std::move(step));
        // The new layer's FIRST turn: a normal round turn, from where the tangent landed, at
        // the height the segment arrived at (the height only changes from the SECOND turn).
        appendBumpedSweep(path, n.x, s.y, azS + dAz, n.x, n.y, azE + kTwoPi, bumpsEnd,
                          wireRadius, label, ordinal, false);
        return;
    }
    // The wrap spans EXACTLY one turn, crossing to crossing. No overshoot: stretching the first
    // and last wrap past the connection plane existed only so a boolean fuse would find
    // overlapping material at the lead, and the conformal assembler fuses nothing -- it was
    // leaving copper sticking out into the lead corridor. Geometry is never extended to satisfy
    // a downstream algorithm; if two pieces must join, they are constructed to meet.
    appendBumpedSweep(path, s.x, s.y, azS, n.x, n.y, azE + kTwoPi, bumps, wireRadius,
                      label, ordinal, /*isConnection=*/false);
}

// ---- Z-return END-RUN planning (round columns) ------------------------------------------
// A Z layer/section return cannot be drawn as a diagonal cone (it interpenetrates the layer
// it climbs over — with layers one wire-OD apart the clearance band is provably narrower
// than the pitch), and no fixed lane assignment survives every winding topology (sibling
// helices, leads and other returns each forbid different azimuth windows per fixture). So
// returns are planned LAST, against the full built obstacle field: for each pending return,
// candidate end-run geometries (azimuth lanes x arc directions) are generated and tested
// with the same sampling/clearance rule as the collision gate, and the first conflict-free
// candidate is inserted. Turns are never moved; the return rides its own source/target ring
// to the lane (exempt adjacent-ordinal contact), runs radially/axially in the free space
// outside every ring, and re-enters past the winding's end.
// ROUNDED CORNERS for a round-wire polyline: emit the run as straight legs joined by exact
// tangent fillet arcs, the same treatment every toroidal corner gets (physical wire bends; a
// sharp butt leaves the conformal assembler a 90-degree junction, and two tubes merely
// overlapping there is what a mitred corner is meant to replace).
//
// Fillet radius must EXCEED the wire radius: a corner whose CENTRELINE radius equals the wire
// radius sweeps a HORN TORUS, whose tube touches its own axis of revolution -- OCC cannot build
// a valid solid from it. Where a leg is too short to host a big enough fillet, that corner is
// left SHARP rather than rounded with an invalid radius.
void appendFilletedPolyline(std::vector<Primitive>& out, const std::vector<gp_Pnt>& raw,
                            double wireRadius, const std::string& labelPrefix, size_t ordinal,
                            bool isLead, bool isConnection) {
    std::vector<gp_Pnt> p;
    for (const auto& q : raw)
        if (p.empty() || p.back().Distance(q) > 1e-12) p.push_back(q);
    if (p.size() < 2) return;
    const size_t n = p.size();
    const double bTarget = 1.5 * wireRadius;      // > wireRadius: never a horn torus
    const double bFloor = 1.05 * wireRadius;      // below this a fillet is not buildable

    std::vector<gp_XYZ> dir(n - 1);
    std::vector<double> len(n - 1);
    for (size_t i = 0; i + 1 < n; ++i) {
        gp_XYZ d = p[i + 1].XYZ() - p[i].XYZ();
        len[i] = d.Modulus();
        dir[i] = d / len[i];
    }
    // Per-corner bend radius and its tangent length, then a pass that shrinks any pair sharing
    // a leg until both fit; a corner whose radius falls under the floor is dropped (kept sharp).
    std::vector<double> radius(n, 0.0), tanLen(n, 0.0), halfTan(n, 0.0);
    for (size_t i = 1; i + 1 < n; ++i) {
        const double cosT = std::max(-1.0, std::min(1.0, dir[i - 1].Dot(dir[i])));
        if (cosT > 1.0 - 1e-9) continue;                       // collinear: no corner
        halfTan[i] = std::sqrt((1.0 - cosT) / (1.0 + cosT));   // tan(turn/2)
        radius[i] = bTarget;
        tanLen[i] = bTarget * halfTan[i];
    }
    // A fillet may never consume most of a leg: the straight stub that survives is what carries
    // the run into its endpoint, and a MKF turn crossing sits exactly at a lead's end -- eat the
    // leg and the arc curves away from the crossing, leaving it outside the copper (measured on
    // the E-core zigzag and EP oblong leads). Each corner takes at most 40% of its shorter
    // neighbouring leg, so every leg keeps at least 20% of its length straight.
    for (size_t i = 1; i + 1 < n; ++i) {
        if (halfTan[i] <= 0.0) continue;
        const double cap = 0.4 * std::min(len[i - 1], len[i]);
        if (tanLen[i] > cap) {
            tanLen[i] = cap;
            radius[i] = cap / halfTan[i];
        }
    }
    for (size_t i = 1; i + 1 < n; ++i)
        if (radius[i] < bFloor) { radius[i] = 0.0; tanLen[i] = 0.0; }

    auto pushSeg = [&](const gp_Pnt& a, const gp_Pnt& b, size_t idx) {
        if (a.Distance(b) < 1e-12) return;
        Primitive pr;
        pr.kind = Primitive::SEG;
        pr.seg = {a, b};
        pr.label = labelPrefix + " seg " + std::to_string(idx);
        pr.turnOrdinal = ordinal;
        pr.isLead = isLead;
        pr.isConnection = isConnection;
        out.push_back(std::move(pr));
    };
    for (size_t i = 0; i + 1 < n; ++i) {
        const gp_Pnt a(p[i].XYZ() + dir[i] * tanLen[i]);
        const gp_Pnt b(p[i + 1].XYZ() - dir[i] * tanLen[i + 1]);
        pushSeg(a, b, i);
        if (i + 2 >= n || tanLen[i + 1] <= 0.0) continue;
        const size_t c = i + 1;
        const gp_Pnt tA(p[c].XYZ() - dir[i] * tanLen[c]);
        const double cosT = std::max(-1.0, std::min(1.0, dir[i].Dot(dir[c])));
        gp_XYZ normal = dir[c] - dir[i] * cosT;
        const double nm = normal.Modulus();
        if (nm < 1e-12) continue;
        normal /= nm;
        gp_XYZ axis = dir[i].Crossed(dir[c]);
        const double am = axis.Modulus();
        if (am < 1e-12) continue;
        Primitive pr;
        pr.kind = Primitive::ARC3;
        pr.arc.c = gp_Pnt(tA.XYZ() + normal * radius[c]);
        pr.arc.axis = axis / am;
        pr.arc.v0 = tA.XYZ() - pr.arc.c.XYZ();
        pr.arc.sweep = std::acos(cosT);
        pr.label = labelPrefix + " corner " + std::to_string(c);
        pr.turnOrdinal = ordinal;
        pr.isLead = isLead;
        pr.isConnection = isConnection;
        out.push_back(std::move(pr));
    }
}

struct PendingZ {
    size_t pathIdx = 0, insertAt = 0, ordinal = 0;
    PlanePt s, n;
    std::string label;
};

// PHYSICAL Z DRAGBACK (round columns). When a layer finishes, the wire must get back to the
// start of the next layer. A real winder does not detour around the whole coil: the wire keeps
// wrapping to the far side, steps out by one wire OD onto the layer it just finished, runs
// AXIALLY back along it -- that axial run is the dragback -- and wraps in to the next layer's
// first station. The layer wound afterwards then rides over that dragback, which is what
// appendRoundWrap's bumps model.
//
// The dragback sits at its FAN SLOT azimuth (kPlaneAz when it is the only vertical there): the
// vertical-fan pre-scan spreads the window's verticals a few degrees apart around the station
// plane so sibling parallels' returns -- which used to share the plane and run through each
// other -- each get their own angle, exactly as a real winder lays them side by side.
void appendZDragback(ConductorPath& path, const PlanePt& s, const PlanePt& n, double wireRadius,
                     double azD, const std::string& label, size_t ordinal,
                     const std::vector<WrapBump>& bumpsOut,
                     const std::vector<WrapBump>& bumpsIn) {
    // The return is CENTRED ON ITS SLOT's axial plane: the wrap arrives at its station on that
    // plane, steps out parallel to the radial direction onto the next layer's radius, drops
    // parallel to Y, and steps out once more onto the first station of the layer that continues.
    // No lateral offset and no stub, so every leg is exactly axis-parallel and the two corners
    // between them are plain 90-degree junctions the assembler mitres at 45 degrees.
    // The raise comes from tallestBumpColumn on the SAME list the adjoining wrap used, so the
    // chain and the wrap meet at one identical point.
    // The raised positions follow appendBumpedSweep's dropped-centre rule (cz = -raise, radius
    // same radius, z dropped by raise) on the SAME bump lists, so chain and wraps meet at identical
    // points.
    const double chainRaise = tallestBumpColumn(bumpsOut).first;
    const gp_Pnt pS = azPointC(0, -chainRaise, s.x, s.y, azD);
    // The turn this return FEEDS rides over it, exactly as every turn at or outside the return's
    // radius does: its first station therefore sits one wire OD further out than the bare layer
    // radius. bumpsIn comes from the WINDOW-WIDE pre-scan and already CONTAINS this return (it is
    // taken at the destination radius, which is the return's own radius), so its diameter must not
    // be added again.
    const double destRaise = tallestBumpColumn(bumpsIn).first;
    // Run-out at the top, descent at the next layer's radius, run-in at the bottom -- the bottom
    // step is the MIRROR of the top one and closes the chain onto the turn it feeds. It must land
    // on that turn's ACTUAL start: the conformal assembler bridges an endpoint mismatch by GROWING
    // the piece, so a descent stopping at the bare station (one OD short of the lifted turn) grew
    // 1.58 mm further down, straight through the entrance lead.
    std::vector<gp_Pnt> direct{
        pS,
        azPointC(0, -chainRaise, n.x, s.y, azD),
        azPointC(0, -chainRaise, n.x, n.y, azD),
        azPointC(0, -destRaise,  n.x, n.y, azD)};
    appendFilletedPolyline(path.prims, direct, wireRadius, label + " (dragback)", ordinal,
                           /*isLead=*/false, /*isConnection=*/true);
}

std::vector<Primitive> buildZEndRun(const PlanePt& s, const PlanePt& n, double azA, double azOff,
                                    bool outFwd, bool homeFwd, double rOut, double yEnd,
                                    const std::string& label, size_t ordinal) {
    std::vector<Primitive> prims;
    const double fwdSweep = std::fmod(azOff - azA + 2.0 * kTwoPi, kTwoPi);   // A -> azOff going +
    auto addArc = [&](double r, double y, double azFrom, bool forward, double sweep,
                      const char* what) {
        if (sweep < 1e-9) return;
        Primitive arc;
        arc.kind = Primitive::ARC3;
        arc.arc.c = gp_Pnt(0, y, 0);
        arc.arc.axis = forward ? gp_XYZ(0, 1, 0) : gp_XYZ(0, -1, 0);  // +Y advances azimuth
        gp_Pnt start = azPointC(0, 0, r, y, azFrom);
        arc.arc.v0 = start.XYZ() - arc.arc.c.XYZ();
        arc.arc.sweep = sweep;
        arc.label = label + " (Z end-run " + std::string(what) + ")";
        arc.turnOrdinal = ordinal;
        arc.isConnection = true;
        prims.push_back(std::move(arc));
    };
    auto addSeg = [&](const PlanePt& a, const PlanePt& b, const char* what) {
        gp_Pnt pa = azPointC(0, 0, a.x, a.y, azOff);
        gp_Pnt pb = azPointC(0, 0, b.x, b.y, azOff);
        if (pa.Distance(pb) < 1e-12) return;
        Primitive step;
        step.kind = Primitive::SEG;
        step.seg = {pa, pb};
        step.label = label + " (Z end-run " + std::string(what) + ")";
        step.turnOrdinal = ordinal;
        step.isConnection = true;
        prims.push_back(std::move(step));
    };
    // out-arc: ride the source ring from the plane to the lane, in either direction.
    addArc(s.x, s.y, azA, outFwd, outFwd ? fwdSweep : (kTwoPi - fwdSweep), "out-arc");
    addSeg({s.x, s.y}, {rOut, s.y}, "ramp");
    addSeg({rOut, s.y}, {rOut, yEnd}, "descent");
    addSeg({rOut, yEnd}, {n.x, yEnd}, "duck");
    addSeg({n.x, yEnd}, {n.x, n.y}, "rise");
    // home-arc: ride the target ring from the lane back to the plane.
    addArc(n.x, n.y, azOff, homeFwd, homeFwd ? (kTwoPi - fwdSweep) : fwdSweep, "home-arc");
    return prims;
}

void planZEndRuns(std::vector<ConductorPath>& paths, std::vector<PendingZ>& pending,
                  const std::vector<RingInv>& rings) {
    if (pending.empty() || rings.empty()) return;
    double maxEdge = -1e30, yLoEdge = 1e30, yHiEdge = -1e30;
    for (const RingInv& g : rings) {
        maxEdge = std::max(maxEdge, g.r + g.rw);
        yLoEdge = std::min(yLoEdge, g.y - g.rw);
        yHiEdge = std::max(yHiEdge, g.y + g.rw);
    }
    // Sampled polyline cache of every existing primitive (obstacle field), kept current as
    // end-runs are inserted.
    std::vector<std::vector<std::vector<gp_Pnt>>> polys(paths.size());
    for (size_t p = 0; p < paths.size(); ++p) {
        polys[p].reserve(paths[p].prims.size());
        for (const auto& pr : paths[p].prims)
            polys[p].push_back(samplePrim(pr, paths[p].wireRadius));
    }
    std::vector<size_t> inserted(paths.size(), 0);   // per-path insertion offset
    for (auto& pz : pending) {
        ConductorPath& P = paths[pz.pathIdx];
        const double crw = P.condRadius > 0 ? P.condRadius : P.wireRadius;
        const double rOut = maxEdge + 1.2 * crw;
        const bool viaBottom = std::abs(pz.n.y - yLoEdge) <= std::abs(pz.n.y - yHiEdge);
        const double yEnd = viaBottom ? (yLoEdge - 1.2 * crw) : (yHiEdge + 1.2 * crw);
        // Candidate lanes x arc directions, scored by worst clearance slack vs everything built.
        double bestSlack = -1e30;
        std::vector<Primitive> best;
        const double azA = kPlaneAz + P.seamRot;   // this path's rotated connection plane
        for (int ai = 0; ai < 16 && bestSlack < 0.0; ++ai) {
            const double azOff = azA + 0.35 + ai * (kTwoPi - 0.7) / 15.0;
            for (int dir = 0; dir < 4 && bestSlack < 0.0; ++dir) {
                auto cand = buildZEndRun(pz.s, pz.n, azA, azOff, dir & 1, dir & 2, rOut, yEnd,
                                         pz.label, pz.ordinal);
                double slack = 1e30;
                std::string worst;
                for (size_t k = 0; k < cand.size(); ++k) {
                    auto cp = samplePrim(cand[k], P.wireRadius);
                    for (size_t q = 0; q < paths.size(); ++q) {
                        const auto& Q = paths[q];
                        const double qrw = Q.condRadius > 0 ? Q.condRadius : Q.wireRadius;
                        const double minGap = gateMinSeparation(crw, qrw);
                        for (size_t j = 0; j < Q.prims.size(); ++j) {
                            if (q == pz.pathIdx) {
                                const auto& pj = Q.prims[j];
                                const size_t d = pj.turnOrdinal > pz.ordinal
                                                     ? pj.turnOrdinal - pz.ordinal
                                                     : pz.ordinal - pj.turnOrdinal;
                                if (d <= 1) continue;   // own adjacent wraps/leads share the crossings
                            }
                            const double sl = polyPolyDistance(cp, polys[q][j]) - minGap;
                            if (sl < slack) {
                                slack = sl;
                                if (std::getenv("MVB_DIAG"))
                                    worst = cand[k].label + " vs " + Q.name + " [" +
                                            Q.prims[j].label + "]";
                            }
                        }
                    }
                }
                if (slack > bestSlack) {
                    bestSlack = slack;
                    best = std::move(cand);
                    if (std::getenv("MVB_DIAG") && !worst.empty())
                        std::cerr << "[z-endrun]   cand az=" << (azOff - azA) << " dirs=" << dir
                                  << " slack=" << slack * 1e6 << "um worst: " << worst << "\n";
                }
            }
        }
        if (std::getenv("MVB_DIAG"))
            std::cerr << "[z-endrun] " << pz.label << ": bestSlack=" << bestSlack * 1e6
                      << "um " << (bestSlack >= 0 ? "CLEAR" : "CONFLICT") << "\n";
        // Insert the best candidate (a conflict-free one when found; otherwise the least-bad —
        // checkCollisions stays the final arbiter and reports it honestly).
        const size_t at = pz.insertAt + inserted[pz.pathIdx];
        P.prims.insert(P.prims.begin() + at, best.begin(), best.end());
        std::vector<std::vector<gp_Pnt>> bp;
        for (const auto& pr : best) bp.push_back(samplePrim(pr, P.wireRadius));
        polys[pz.pathIdx].insert(polys[pz.pathIdx].begin() + at, bp.begin(), bp.end());
        inserted[pz.pathIdx] += best.size();
    }
}

// RECTANGULAR column: the racetrack of build_concentric_rect_column_turn — straights at
// exactly the MAS radial clearance on all four faces, corner arcs of radius = clearance
// (grown to the wire's minimum bend radius with the arc axis pulled inward when the wire
// sits closer than that; identical rule and throw as TurnBuilder). The whole racetrack
// runs at the SOURCE crossing's station; the transition to the next crossing is the
// final straight of the -Z face (the user-approved rectangular-face zigzag: three sides
// and the corners stay on-station, one face carries the move).
struct RectStation {
    double y = 0, zPos = 0, xPos = 0, cornerR = 0, segX = 0, segZ = 0;
};
RectStation rectStation(const PlanePt& p, double halfW, double halfD, double minBend,
                        const std::string& who) {
    double clearance = p.x - halfW;
    if (clearance <= 0.0) {
        throw std::runtime_error(
            "ConductorBuilder: crossing radial position " + std::to_string(p.x) +
            " m lies inside the column (half-width " + std::to_string(halfW) +
            " m) for " + who + " — inconsistent MAS turn/bobbin data");
    }
    double cornerR = clearance;
    double inset = 0.0;
    if (cornerR < minBend) {
        inset = minBend - clearance;
        cornerR = minBend;
        if (inset > std::min(halfW, halfD)) {
            throw std::runtime_error(
                "ConductorBuilder: minimum bend radius " + std::to_string(minBend) +
                " m cannot be accommodated (clearance " + std::to_string(clearance) +
                " m, column half-dims " + std::to_string(halfW) + "/" +
                std::to_string(halfD) + " m) for " + who +
                " — adjacent corner arcs would cross");
    }
    }
    return {p.y, halfD + clearance, halfW + clearance, cornerR, halfW - inset,
            halfD - inset};
}

// ride0 / ride1: the -Z displacement of this and the next crossing (one laid OD per
// dragback level at or inside the depth -- see rectRideFor). For a CROSS-LAYER return
// (isReturn): chainRide = displacement of the DESCENT's face (returns of deeper levels
// only), destRide = displacement of the destination crossing (own level included), xSlot =
// this conductor's descent lane on the face.
void appendRectWrap(ConductorPath& path, const RectStation& s0, const RectStation& s1,
                    const std::string& label, size_t ordinal, double wireRadius,
                    double ride0, double rideBack0, bool isReturn, double chainRide,
                    double destRide, double xSlot,
                    // When the NEXT transition is a dragback, its descent x-slot: the rising
                    // turn ENDS there and corners straight into the step-out (Alf, 17_cllc:
                    // no run past the crossing and back -- the U-fold pieces must not exist).
                    double stopAtX = std::numeric_limits<double>::quiet_NaN(),
                    // First transition feeding an entrance lead corner: the first straight
                    // BEGINS this far along its own direction (the corner occupies [0, startAtX]).
                    double startAtX = std::numeric_limits<double>::quiet_NaN()) {
    auto pushSeg = [&](const gp_Pnt& a, const gp_Pnt& b, const char* what) {
        if (a.Distance(b) < 1e-12) return;
        Primitive pr;
        pr.kind = Primitive::SEG;
        pr.seg = {a, b};
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    auto pushCorner = [&](double cxx, double czz, double azStart, double yNow,
                          const char* what) {
        Primitive pr;
        pr.kind = Primitive::ARC3;
        pr.arc.c = gp_Pnt(cxx, yNow, czz);
        pr.arc.axis = gp_XYZ(0, 1, 0);
        pr.arc.v0 = gp_XYZ(s0.cornerR * std::cos(azStart), 0,
                           -s0.cornerR * std::sin(azStart));
        pr.arc.sweep = kPi / 2.0;
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    const double y = s0.y;
    // Displaced -Z geometry (the dragback reservation): the -Z face and its corners sit
    // `ride` further out; the lateral +-X faces EXTEND to reach them. The +Z side never
    // moves.
    const double zN0 = s0.zPos + ride0;         // this crossing's -Z face depth
    const double cZ0 = s0.segZ + ride0;         // -Z corner centres' z
    const double zP0 = s0.zPos + rideBack0;     // the BACK (+Z) face, displaced by the
    const double cP0 = s0.segZ + rideBack0;     // OPPOSITE side's dragback levels
    if (!isReturn) {
        if (std::abs(s1.zPos - s0.zPos) > 1e-12) {
            throw std::runtime_error(
                "ConductorBuilder: cross-layer transition of " + label +
                " reached the rising-turn branch (detection drift vs the rectReturns "
                "pre-scan)");
        }
        // ALF'S RISING RECT TURN (2026-08-06): a same-depth transition distributes its pitch
        // over the turn's WHOLE path length -- straights AND corner arcs -- each piece
        // climbing proportionally to its length, so the climb rate is uniform and every
        // junction is tangent-continuous. The pitch is the turnsDescription delta (next
        // crossing's y minus this one's), so parallels wind as PARALLEL rising ribbons a
        // full row apart and can never cross. The EXTENDED lateral length from any dragback
        // reservation joins this same distribution (Alf).
        const double q = 0.5 * kPi * s0.cornerR;
        // A rising turn that feeds a dragback ends AT the descent slot, so the pitch is
        // distributed over the path it actually travels (crossing-to-slot shortening).
        const double endX = std::isnan(stopAtX) ? 0.0 : stopAtX;
        const double begX = std::isnan(startAtX) ? 0.0 : startAtX;
        if (!std::isnan(stopAtX) && std::abs(stopAtX) > s0.segX)
            throw std::runtime_error(
                "ConductorBuilder: dragback x-slot " + std::to_string(stopAtX) +
                " m lies outside the -Z face straight of " + label);
        if (begX < 0.0 || begX > s0.segX)
            throw std::runtime_error(
                "ConductorBuilder: entrance corner offset " + std::to_string(begX) +
                " m lies outside the -Z face straight of " + label);
        const double L = 4.0 * s0.segX + 4.0 * s0.segZ + 2.0 * ride0 + 2.0 * rideBack0 +
                         4.0 * q - endX - begX;
        if (L < 1e-12) {
            throw std::runtime_error(
                "ConductorBuilder: rect turn of " + label +
                " has no path length to distribute its pitch over");
        }
        const double dy = s1.y - s0.y;
        double arc = 0.0;
        auto yAt = [&](double at) { return s0.y + dy * at / L; };
        auto riseSeg = [&](double ax, double az2, double bx, double bz2, double len,
                           const char* what) {
            pushSeg(gp_Pnt(ax, yAt(arc), az2), gp_Pnt(bx, yAt(arc + len), bz2), what);
            arc += len;
        };
        auto riseCorner = [&](double cxx, double czz, double azStart, const char* what) {
            Primitive pr;
            pr.kind = Primitive::SPIRAL;
            pr.spiral = {cxx, czz, s0.cornerR, yAt(arc), azStart,
                         s0.cornerR, yAt(arc + q), azStart + kPi / 2.0};
            pr.label = label + std::string(" ") + what;
            pr.turnOrdinal = ordinal;
            path.prims.push_back(std::move(pr));
            arc += q;
        };
        riseSeg(-begX, -zN0, -s0.segX, -zN0, s0.segX - begX, "face -Z out");
        riseCorner(-s0.segX, -cZ0, kPi / 2.0, "corner -X-Z");
        riseSeg(-s0.xPos, -cZ0, -s0.xPos, +cP0, 2.0 * s0.segZ + ride0 + rideBack0,
                "face -X");
        riseCorner(-s0.segX, +cP0, kPi, "corner -X+Z");
        riseSeg(-s0.segX, +zP0, +s0.segX, +zP0, 2.0 * s0.segX, "face +Z");
        riseCorner(+s0.segX, +cP0, 3.0 * kPi / 2.0, "corner +X+Z");
        riseSeg(+s0.xPos, +cP0, +s0.xPos, -cZ0, 2.0 * s0.segZ + ride0 + rideBack0,
                "face +X");
        riseCorner(+s0.segX, -cZ0, 0.0, "corner +X-Z");
        riseSeg(+s0.segX, -zN0, endX, -zN0, s0.segX - endX, "face -Z in");
        return;
    }
    // CROSS-LAYER DRAGBACK (Alf, 2026-08-06): the CHAIN ONLY, no revolution -- exactly like
    // the round column's appendZDragback. The rising wrap already delivered the wire to this
    // crossing at the full pitch; an extra flat revolution here retraced the row the wrap
    // had just climbed to, 0.21 mm off its own final approach (measured on 02_flyback,
    // gate-blind because the pieces are ordinal-adjacent). The return runs ON the -Z face:
    // along the face to this conductor's x slot, step OUT to the descent depth (the
    // destination layer's face, displaced by any DEEPER dragback levels), descend axially
    // there, half a step out onto the destination's own displaced face, and along the row
    // to the next crossing. Everything wound at or outside the destination depth rides one
    // OD further out (the reservation), so the descent always has its own free face.
    {
        const double zDesc = s1.zPos + chainRide;   // the descent's face depth
        const double zDest = s1.zPos + destRide;    // the destination crossing's depth
        // The chain STARTS at the slot: the preceding rising turn already ends there
        // (stopAtX) and corners into the step-out -- the old run from the face crossing
        // to the slot doubled back over the just-arrived straight (Alf, 17_cllc:
        // "Secondary parallel 110 should not exist"). Only a chain with NO preceding
        // wrap (a return as the conductor's first transition) keeps the crossing run.
        std::vector<gp_Pnt> pts;
        if (ordinal == 0) pts.push_back(gp_Pnt(0, y, -zN0));
        pts.insert(pts.end(), {gp_Pnt(xSlot, y, -zN0),
                                gp_Pnt(xSlot, y, -zDesc),
                                gp_Pnt(xSlot, s1.y, -zDesc),
                                gp_Pnt(xSlot, s1.y, -zDest),
                                gp_Pnt(0, s1.y, -zDest)});
        appendFilletedPolyline(path.prims, pts, wireRadius, label + " (dragback)", ordinal,
                               /*isLead=*/false, /*isConnection=*/true);
    }
}

// OBLONG column (stadium): straights on the two flat +-X faces at x = +-radial, half caps
// of radius = radial about the end centres (0, +-straightHalf) — the same decomposition
// as build_concentric_oblong_turn. The -Z cap's second quarter is the transition: a
// spiral about the cap centre from az 0 to the crossing apex at az pi/2 (the flat faces
// stay on-station; the round end carries the move).
void appendOblongWrap(ConductorPath& path, const PlanePt& s, const PlanePt& n,
                      double straightHalf, const std::string& label, size_t ordinal) {
    auto pushSeg = [&](const gp_Pnt& a, const gp_Pnt& b, const char* what) {
        if (a.Distance(b) < 1e-12) return;
        Primitive pr;
        pr.kind = Primitive::SEG;
        pr.seg = {a, b};
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    auto pushCap = [&](double cz, double azStart, double sweep, const char* what) {
        Primitive pr;
        pr.kind = Primitive::ARC3;
        pr.arc.c = gp_Pnt(0, s.y, cz);
        pr.arc.axis = gp_XYZ(0, 1, 0);
        pr.arc.v0 = gp_XYZ(s.x * std::cos(azStart), 0, -s.x * std::sin(azStart));
        pr.arc.sweep = sweep;
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    const double sh = straightHalf;
    pushCap(-sh, kPi / 2.0, kPi / 2.0, "cap -Z out");
    pushSeg(gp_Pnt(-s.x, s.y, -sh), gp_Pnt(-s.x, s.y, +sh), "face -X");
    pushCap(+sh, kPi, kPi, "cap +Z");
    pushSeg(gp_Pnt(+s.x, s.y, +sh), gp_Pnt(+s.x, s.y, -sh), "face +X");
    auto pushRamp = [&](double r0, double y0, double az0, double r1, double y1,
                        double az1, const char* what) {
        Primitive ramp;
        ramp.kind = Primitive::SPIRAL;
        ramp.spiral = {0, -sh, r0, y0, az0, r1, y1, az1, /*blend=*/true};
        ramp.label = label + std::string(" ") + what;
        ramp.turnOrdinal = ordinal;
        path.prims.push_back(std::move(ramp));
    };
    if (std::abs(n.x - s.x) > 1e-12 && std::abs(n.y - s.y) > 1e-12) {
        // Same two-phase rule as the rectangular ramp: radial move first (clearing the
        // crossed layer), axial move at the new clearance second.
        pushRamp(s.x, s.y, 0.0, n.x, s.y, kPi / 4.0, "cap -Z ramp (radial)");
        pushRamp(n.x, s.y, kPi / 4.0, n.x, n.y, kPi / 2.0, "cap -Z ramp (axial)");
    } else {
        pushRamp(s.x, s.y, 0.0, n.x, n.y, kPi / 2.0, "cap -Z ramp");
    }
}

// ---------------------------------------------------------------------------------------
// TOROIDAL wraps. Build frame: toroid hole axis along +Y, hole plane at y = 0 (the
// assembly is later counter-rotated by -pi/2 about X, mapping build (x,y,z) to final
// (x,z,-y), so a MAS hole-plane point (cx, cy, 0) is build (cx, 0, cy)). One wrap is a
// poloidal loop around the core cross-section: up the inner tube at the inner crossing,
// over the top face along the straight chord between corner exits (the exact
// construction of build_toroidal_turn), down the outer tube through the OUTER crossing
// (additionalCoordinates — MAS data, hit exactly), back under the bottom face to the
// NEXT inner crossing. Tube lengths mirror TurnBuilder's derived clearances:
// tube = halfDepth + layerOffset, layerOffset = max(0, wwRadialHeight - innerRadial -
// wireRadius) — the transition run uses the deeper of the two adjacent crossings'
// clearances so it clears the intervening rings.
struct ToroCross {
    gp_XY pin;    // inner crossing, build-frame horizontal (x, z) = MAS (x, y)
    gp_XY pout;   // outer crossing
    double tube;  // vertical tube length above/below the hole plane
};

// Bottom-half BASE running depth of the wrap between two crossings — the clearance under
// the core face shared by all rings. Ring-to-ring separation is added on top as an
// explicit stagger (see appendToroWrap RING STAGGER); the per-ring tube length already
// differs via layerOffset, but that clamps to 0 once the window is full, so the explicit
// stagger is what actually nests the rings.
double toroWrapDepth(const ToroCross& c0, const ToroCross& c1, double /*bend*/) {
    return std::max(c0.tube, c1.tube);
}

void appendToroWrap(ConductorPath& path, const ToroCross& c0, const ToroCross& c1,
                    double bend, double bottomExtraDepth, const std::string& label,
                    size_t ordinal) {
    auto P = [](const gp_XY& h, double y) { return gp_Pnt(h.X(), y, h.Y()); };
    auto pushSeg = [&](const gp_Pnt& a, const gp_Pnt& b, const char* what) {
        if (a.Distance(b) < 1e-12) return;
        Primitive pr;
        pr.kind = Primitive::SEG;
        pr.seg = {a, b};
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    auto pushArc = [&](const gp_Pnt& center, const gp_XYZ& axis, const gp_XYZ& v0,
                       const char* what) {
        Primitive pr;
        pr.kind = Primitive::ARC3;
        pr.arc.c = center;
        pr.arc.axis = axis;
        pr.arc.v0 = v0;
        pr.arc.sweep = kPi / 2.0;
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    const gp_XYZ yHat(0, 1, 0);
    const double b = bend;

    // Top half: turn k-1's inner leg to its outer leg.
    gp_XY dH = c0.pout - c0.pin;
    double lTop = dH.Modulus();
    if (lTop <= 2.0 * b) {
        throw std::runtime_error("ConductorBuilder: toroidal crossings of " + label +
                                 " are closer than two bend radii (" + std::to_string(lTop) +
                                 " m) — no room for the face run");
    }
    dH.Divide(lTop);
    gp_XYZ d3(dH.X(), 0, dH.Y());
    double t0 = c0.tube, rh0 = t0 + b;

    pushSeg(P(c0.pin, 0), P(c0.pin, t0), "inner tube up");
    pushArc(P(c0.pin + dH * b, t0), yHat.Crossed(d3), d3 * (-b), "top inner corner");
    pushSeg(P(c0.pin + dH * b, rh0), P(c0.pout - dH * b, rh0), "top chord");
    pushArc(P(c0.pout - dH * b, t0), yHat.Crossed(d3), yHat * b, "top outer corner");

    // Bottom half: turn k-1's outer leg to turn k's inner leg. Both corners swing PURELY
    // RADIALLY (parallel among packed neighbours — a corner swung toward the per-turn
    // chord target crosses the neighbouring tubes' corners), and the direction change
    // happens down at the running depth as an arc-line-arc of the same bend radius, so
    // the whole half stays tangent-continuous.
    //
    // RING STAGGER (bottomExtraDepth): an OUTER ring's return sweeps radially inward
    // across the radii the INNER rings occupy; at equal depth its inward corner-swing
    // clips the inner ring's return (MKF packs different rings' rim crossings under one
    // OD apart — no corner clearance, MKF ABT #231). Each ring's whole bottom return
    // therefore tucks one wire OD deeper below the core than the ring inside it, so the
    // outer sweeps pass UNDER the inner returns and nest instead of colliding. The top
    // half is already auto-staggered (its over-core height = the tube length, which
    // grows for inner rings), so only the bottom needs this. Single-ring windings pass
    // bottomExtraDepth = 0 and are unchanged.
    double tb = toroWrapDepth(c0, c1, b) + bottomExtraDepth, rhb = tb + b;
    gp_XY eH = c1.pin - c0.pout;
    double lBot = eH.Modulus();
    if (lBot <= 2.0 * b) {
        throw std::runtime_error("ConductorBuilder: toroidal return crossings of " + label +
                                 " are closer than two bend radii (" + std::to_string(lBot) +
                                 " m) — no room for the face run");
    }
    eH.Divide(lBot);
    gp_XYZ e3(eH.X(), 0, eH.Y());

    // Mirror of the TOP half at the (ring-staggered) running depth: the corners swing
    // TOWARD the straight chord (not radially), exactly like the top's, so the chord runs
    // parallel among a ring's packed neighbours. Different rings never meet because the
    // ring stagger puts each ring's whole return one OD deeper than the ring inside it
    // (the radial-swing + Dubins routing this replaces was for same-depth rings and
    // clipped tightly-packed neighbours).
    pushSeg(P(c0.pout, t0), P(c0.pout, -tb), "outer tube down");
    pushArc(P(c0.pout + eH * b, -tb), e3.Crossed(yHat), e3 * (-b), "bottom outer corner");
    pushSeg(P(c0.pout + eH * b, -rhb), P(c1.pin - eH * b, -rhb), "bottom chord");
    pushArc(P(c1.pin - eH * b, -tb), e3.Crossed(yHat), yHat * (-b), "bottom inner corner");
    pushSeg(P(c1.pin, -tb), P(c1.pin, 0), "inner tube up to crossing");
}

// RING-TRANSITION wrap (the Z connection), dragback in the FREE BAND the ring stagger already
// leaves between the two rings' returns (level 2k). The run itself is NOT the old straight
// diagonal -- that passed within a wire OD of the destination ring's inner tubes as it cut
// across the face. Instead it follows the CORE'S CENTRAL CIRCLE: a radial leg out of the source
// turn's outer descent, a fillet corner onto the centre-radius arc, the azimuthal arc along the
// core's mid-annulus (max distance from BOTH rings' inner and outer tubes), a fillet corner off
// it, and a radial leg arriving straight in front of the destination turn. All five pieces are
// exact tangent constructions (line-circle fillets), so the conformal mitre assembly sees only
// tangent junctions.
void appendToroTransitionBand(ConductorPath& path, const ToroCross& c0, const ToroCross& c1,
                              double bend, double extraDepth, double rMid,
                              const std::string& label, size_t ordinal) {
    auto P = [](const gp_XY& h, double y) { return gp_Pnt(h.X(), y, h.Y()); };
    auto pushSeg = [&](const gp_Pnt& a, const gp_Pnt& b, const char* what) {
        if (a.Distance(b) < 1e-12) return;
        Primitive pr;
        pr.kind = Primitive::SEG;
        pr.seg = {a, b};
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    auto pushArc90 = [&](const gp_Pnt& center, const gp_XYZ& axis, const gp_XYZ& v0,
                         const char* what) {
        Primitive pr;
        pr.kind = Primitive::ARC3;
        pr.arc.c = center;
        pr.arc.axis = axis;
        pr.arc.v0 = v0;
        pr.arc.sweep = kPi / 2.0;
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    // Horizontal arc from startXY to endXY about centerXY at height y; the signed sweep about
    // +Y is derived from the endpoints (always the < pi solution).
    auto pushArcH = [&](const gp_XY& centerXY, double y, const gp_XY& startXY,
                        const gp_XY& endXY, const char* what) {
        gp_XYZ v0(startXY.X() - centerXY.X(), 0, startXY.Y() - centerXY.Y());
        gp_XYZ v1(endXY.X() - centerXY.X(), 0, endXY.Y() - centerXY.Y());
        double sweep = std::atan2(v0.Crossed(v1).Dot(gp_XYZ(0, 1, 0)), v0.Dot(v1));
        if (std::abs(sweep) < 1e-12) return;
        Primitive pr;
        pr.kind = Primitive::ARC3;
        pr.arc.c = gp_Pnt(centerXY.X(), y, centerXY.Y());
        // ALWAYS a positive sweep about the matching axis: OCC's revolve/edge builders refuse
        // negative angles, so a clockwise arc emitted with sweep<0 produced an invalid solid
        // that the conformal assembler DROPPED -- the visibly missing dragback corner.
        pr.arc.axis = sweep < 0 ? gp_XYZ(0, -1, 0) : gp_XYZ(0, 1, 0);
        pr.arc.v0 = v0;
        pr.arc.sweep = std::abs(sweep);
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    auto pol = [](double r, double a) { return gp_XY(r * std::cos(a), r * std::sin(a)); };
    const gp_XYZ yHat(0, 1, 0);
    const double b = bend;
    // Fillet bend radius. It must be STRICTLY GREATER than the wire radius: a corner whose
    // centreline radius equals the wire radius sweeps a HORN TORUS -- the tube touches its own
    // axis of revolution, so OCC returns a self-touching (invalid) solid and the conformal
    // assembler drops it, which is exactly the "missing dragback corner" seen in the STEP.
    // 1.5x is also the more physical figure: real wire cannot bend to its own radius without
    // crushing the inner fibre.
    const double bf = 1.5 * bend;

    // ---- source turn's top half (identical to a normal wrap's) ----
    gp_XY dH = c0.pout - c0.pin;
    double lTop = dH.Modulus();
    if (lTop <= 2.0 * b) {
        throw std::runtime_error("ConductorBuilder: toroidal crossings of " + label +
                                 " are closer than two bend radii (" + std::to_string(lTop) +
                                 " m) — no room for the face run");
    }
    dH.Divide(lTop);
    gp_XYZ d3(dH.X(), 0, dH.Y());
    const double t0 = c0.tube, rh0 = t0 + b;
    pushSeg(P(c0.pin, 0), P(c0.pin, t0), "inner tube up");
    pushArc90(P(c0.pin + dH * b, t0), yHat.Crossed(d3), d3 * (-b), "top inner corner");
    pushSeg(P(c0.pin + dH * b, rh0), P(c0.pout - dH * b, rh0), "top chord");
    pushArc90(P(c0.pout - dH * b, t0), yHat.Crossed(d3), yHat * b, "top outer corner");

    // ---- band depth and sweep sense ----
    const double tb = toroWrapDepth(c0, c1, b) + extraDepth, rhb = tb + b;
    const double rOut0 = c0.pout.Modulus(), rIn1 = c1.pin.Modulus();
    const double az0 = std::atan2(c0.pout.Y(), c0.pout.X());
    const double az1 = std::atan2(c1.pin.Y(), c1.pin.X());
    const double dAz = std::remainder(az1 - az0, 2.0 * kPi);
    const double sg = dAz < 0 ? -1.0 : 1.0;
    // Fillet geometry: the out-leg meets the centre circle from OUTSIDE (fillet centre at
    // rMid+b), the in-leg leaves it towards the hole (fillet centre at rMid-b). Tangency radii
    // on the radial lines and tangency angles on the circle are exact.
    const double d1 = std::asin(bf / (rMid + bf)), rt1 = std::sqrt((rMid + bf) * (rMid + bf) - bf * bf);
    const double d2 = std::asin(bf / (rMid - bf)), rt2 = std::sqrt((rMid - bf) * (rMid - bf) - bf * bf);
    if (rOut0 - b <= rt1 || rIn1 + b >= rt2 || std::abs(dAz) <= d1 + d2)
        throw std::runtime_error(
            "ConductorBuilder: ring-transition dragback of " + label +
            " cannot follow the core's central radius (r=" + std::to_string(rMid) +
            " m): the source outer / destination inner crossings leave no room for the "
            "radial legs and fillets");

    // ---- descend at the source azimuth ----
    gp_XY e0xy = c0.pout;
    e0xy.Divide(-rOut0);   // radial inward unit at az0
    gp_XYZ e0(e0xy.X(), 0, e0xy.Y());
    pushSeg(P(c0.pout, t0), P(c0.pout, -tb), "outer tube down");
    pushArc90(P(c0.pout + e0xy * b, -tb), e0.Crossed(yHat), e0 * (-b), "bottom outer corner");

    // ---- the 5-piece dragback along the core's central circle ----
    pushSeg(P(c0.pout + e0xy * b, -rhb), P(pol(rt1, az0), -rhb), "dragback out leg");
    pushArcH(pol(rMid + bf, az0 + sg * d1), -rhb, pol(rt1, az0), pol(rMid, az0 + sg * d1),
             "dragback outer fillet");
    pushArcH(gp_XY(0, 0), -rhb, pol(rMid, az0 + sg * d1), pol(rMid, az1 - sg * d2),
             "dragback central arc");
    pushArcH(pol(rMid - bf, az1 - sg * d2), -rhb, pol(rMid, az1 - sg * d2), pol(rt2, az1),
             "dragback inner fillet");
    gp_XY eIxy = c1.pin;
    eIxy.Divide(-rIn1);   // radial inward unit at az1 (direction of travel)
    gp_XYZ eI(eIxy.X(), 0, eIxy.Y());
    pushSeg(P(pol(rt2, az1), -rhb), P(c1.pin - eIxy * b, -rhb), "dragback in leg");

    // ---- arrive at the destination turn ----
    pushArc90(P(c1.pin - eIxy * b, -tb), eI.Crossed(yHat), yHat * (-b), "bottom inner corner");
    pushSeg(P(c1.pin, -tb), P(c1.pin, 0), "inner tube up to crossing");
}


// Toroidal terminal lead rect: MKF draws ONE radial rect per lead in the hole plane,
// spanning [crossing radius, radialBorder] along the crossing's angle — its near edge IS
// the crossing (up to MKF's 1e-9 serialisation rounding). Decode its outward unit
// direction and reserved length; the caller picks the realizable route (MKF ABT #230:
// the reservation ignores both the window bore and other rings' turns).
struct ToroLeadRect {
    gp_XY dir;        // unit, from the crossing outward
    double length;    // to the rect's far EDGE (the reserved envelope end)
};
ToroLeadRect toroLeadRect(const RSpace& s, const gp_XY& crossing, const std::string& who) {
    double th = s.rotation * kPi / 180.0;
    gp_XY u(std::cos(th), std::sin(th));
    gp_XY cen(s.coordinates.at(0), s.coordinates.at(1));
    gp_XY e1 = cen - u * (s.dimensions.at(0) / 2.0);
    gp_XY e2 = cen + u * (s.dimensions.at(0) / 2.0);
    double d1 = (e1 - crossing).Modulus();
    double d2 = (e2 - crossing).Modulus();
    // MKF rounds rect coordinates to 1e-9 and the rotation to 1e-6 deg; anything beyond
    // 1e-7 m means the rect does not belong to this crossing.
    if (std::min(d1, d2) > 1e-7) {
        throw std::runtime_error(
            "ConductorBuilder: toroidal terminal rect for " + who +
            " does not touch its crossing (nearest end " +
            std::to_string(std::min(d1, d2)) + " m away) — MKF reserved-space mismatch");
    }
    gp_XY farEdge = (d1 < d2) ? e2 : e1;
    gp_XY dir = farEdge - crossing;
    double rectLen = dir.Modulus();
    if (rectLen > 1e-12) dir.Divide(rectLen);
    double crossR = crossing.Modulus();
    if (rectLen > 1e-12 &&
        (crossR < 1e-12 || std::abs(dir.Dot(crossing) / crossR - 1.0) > 1e-6)) {
        throw std::runtime_error(
            "ConductorBuilder: toroidal terminal rect for " + who +
            " is not radial from its crossing — unexpected MKF reserved-space layout");
    }
    return {dir, rectLen};
}

double pointSegDistance2d(const gp_XY& p, const gp_XY& a, const gp_XY& b) {
    gp_XY d = b - a;
    double len2 = d.Dot(d);
    if (len2 < 1e-30) return (p - a).Modulus();
    double t = std::clamp((p - a).Dot(d) / len2, 0.0, 1.0);
    return (p - (a + d * t)).Modulus();
}

} // anonymous namespace

// ---------------------------------------------------------------------------------------
template <typename CoilT, typename WireT>
std::vector<NamedShape> buildAllImpl(const CoilT& coil,
                                     const MAS::CoreBobbinProcessedDescription& bobbinPd,
                                     bool isToroidal,
                                     std::vector<RSpace> allSpaces,
                                     const ConductorBuilder::Options& opts,
                                     std::vector<ConductorBuilder::PathPolyline>* polyOut = nullptr) {
    auto turnsOpt = coil.get_turns_description();
    if (!turnsOpt || turnsOpt->empty()) {
        throw std::runtime_error(
            "ConductorBuilder: coil has no turnsDescription — the magnetic must be wound "
            "by MKF before building real-winding conductors");
    }

    // Drawn connection routes only (layer == ""): the pink/blue boxes of the SVG.
    std::vector<RSpace> drawn;
    for (auto& s : allSpaces) {
        if (s.layer.empty()) drawn.push_back(std::move(s));
    }

    // Wire per winding.
    std::map<std::string, MAS::Wire> wireMap;
    for (const auto& winding : coil.get_functional_description()) {
        const auto& wireVar = winding.get_wire();
        if (std::holds_alternative<std::string>(wireVar)) {
            wireMap[winding.get_name()] =
                OpenMagnetics::find_wire_by_name(std::get<std::string>(wireVar));
        } else {
            wireMap[winding.get_name()] = std::get<WireT>(wireVar);
        }
    }

    // Group turns per (winding, parallel) preserving turnsDescription order (= electrical
    // order, MKF's contract).
    struct ConductorTurns {
        std::string winding;
        int64_t parallel;
        std::vector<const MAS::Turn*> turns;
    };
    std::vector<ConductorTurns> conductors;
    for (const auto& winding : coil.get_functional_description()) {
        for (int64_t k = 0; k < winding.get_number_parallels(); ++k) {
            conductors.push_back({winding.get_name(), k, {}});
        }
    }
    for (const auto& turn : *turnsOpt) {
        auto it = std::find_if(conductors.begin(), conductors.end(), [&](const ConductorTurns& c) {
            return c.winding == turn.get_winding() && c.parallel == turn.get_parallel();
        });
        if (it == conductors.end()) {
            throw std::runtime_error("ConductorBuilder: turn '" + turn.get_name() +
                                     "' references winding '" + turn.get_winding() +
                                     "' parallel " + std::to_string(turn.get_parallel()) +
                                     " which is not declared in coil.functionalDescription");
        }
        it->turns.push_back(&turn);
    }
    for (const auto& c : conductors) {
        if (c.turns.empty()) {
            throw std::runtime_error("ConductorBuilder: winding '" + c.winding + "' parallel " +
                                     std::to_string(c.parallel) +
                                     " has no turns in turnsDescription");
        }
    }

    // Column geometry (concentric) / toroidal window data.
    MAS::ColumnShape columnShape = bobbinPd.get_column_shape();
    // IRREGULAR columns (EFD / EPX-style poles: flat-sided profiles that are none of the three
    // analytic shapes): approximate the WRAP PATH with the column's bounding RECTANGLE
    // (columnWidth x columnDepth). This is a stated geometric approximation, not a silent one --
    // the turn positions themselves come from MKF's turnsDescription and are never moved; only the
    // racetrack route BETWEEN crossings changes, and the bounding rect ENCLOSES the real pole, so
    // the approximated wrap can only be farther from the core, never inside it. Without this every
    // EFD-family design fails real-winding export outright.
    if (columnShape == MAS::ColumnShape::IRREGULAR && !isToroidal &&
        bobbinPd.get_column_width()) {
        std::cerr << "[ConductorBuilder] IRREGULAR column: approximating the wrap path with the "
                     "bounding rectangle (" << *bobbinPd.get_column_width() * 1e3 << " x "
                  << bobbinPd.get_column_depth() * 1e3 << " mm half-dims); turn positions are "
                     "MKF's and unchanged\n";
        columnShape = MAS::ColumnShape::RECTANGULAR;
    }
    const double halfD = bobbinPd.get_column_depth();
    double halfW = 0.0;
    if (!isToroidal &&
        (columnShape == MAS::ColumnShape::RECTANGULAR ||
         columnShape == MAS::ColumnShape::OBLONG)) {
        auto wOpt = bobbinPd.get_column_width();
        if (!wOpt) {
            throw std::runtime_error(
                "ConductorBuilder: bobbin processedDescription has no columnWidth "
                "(required for rectangular/oblong column racetracks)");
        }
        halfW = *wOpt;
    }
    if (!isToroidal && columnShape != MAS::ColumnShape::ROUND &&
        columnShape != MAS::ColumnShape::RECTANGULAR &&
        columnShape != MAS::ColumnShape::OBLONG) {
        throw std::runtime_error(
            "ConductorBuilder: unsupported column shape for real-winding conductors "
            "(only round/rectangular/oblong concentric columns and toroids are modelled)");
    }
    // Oblong with no straight section degenerates to a round column (same rule as
    // build_concentric_oblong_turn).
    const double oblongHalf = halfD - halfW;
    const bool effectivelyRound =
        columnShape == MAS::ColumnShape::ROUND ||
        (columnShape == MAS::ColumnShape::OBLONG && oblongHalf <= 0.0);
    // MKF's 2D window x maps to the connection plane at z = -(x + zoff): the crossing
    // sits one clearance past the column on the -Z side, which for rectangular/oblong
    // columns is columnDepth - columnWidth deeper than the radial position itself.
    const double zoff =
        (isToroidal || effectivelyRound) ? 0.0 : (halfD - halfW);

    double wwRadialHeight = 0.0;
    if (isToroidal) {
        const auto& wws = bobbinPd.get_winding_windows();
        if (wws.empty() || !wws[0].get_radial_height()) {
            throw std::runtime_error(
                "ConductorBuilder: toroidal bobbin has no winding-window radialHeight "
                "(required for the conductor's face-run clearance)");
        }
        wwRadialHeight = wws[0].get_radial_height().value();
    }

    auto station = [](const MAS::Turn* t) -> PlanePt {
        const auto& c = t->get_coordinates();
        if (c.size() < 2) {
            throw std::runtime_error("ConductorBuilder: turn '" + t->get_name() +
                                     "' has fewer than 2 coordinates");
        }
        return {c[0], c[1]};
    };

    // All inner crossings with their wire radii — the corridor test for toroidal
    // in-plane leads (an MKF lead rect drawn through another ring's tubes is
    // unrealizable in the hole plane; the wire then leaves axially instead).
    std::vector<std::pair<gp_XY, double>> allToroCrossings;
    if (isToroidal) {
        for (const auto& ct : conductors) {
            const MAS::Wire& w = wireMap.at(ct.winding);
            for (const MAS::Turn* t : ct.turns) {
                const auto& c = t->get_coordinates();
                if (c.size() < 2) continue;
                auto [ww, wh] = TurnBuilder::wireDimensions(w, *t, opts.paintCoating);
                allToroCrossings.push_back({gp_XY(c[0], c[1]), std::min(ww, wh) / 2.0});
            }
        }
    }

    // Every placed ring of EVERY conductor with its conducting radius (window half-plane) —
    // the obstacle inventory the round-column Z-transition end-run routes around.
    std::vector<RingInv> allRings;
    if (!isToroidal) {
        for (const auto& ct : conductors) {
            const MAS::Wire& w = wireMap.at(ct.winding);
            for (const MAS::Turn* t : ct.turns) {
                const auto& c = t->get_coordinates();
                if (c.size() < 2) continue;
                auto [cw2, ch2] = TurnBuilder::wireDimensions(w, *t, /*paintCoating=*/false);
                allRings.push_back({c[0], c[1], std::min(cw2, ch2) / 2.0});
            }
        }
    }

    std::vector<ConductorPath> paths;
    paths.reserve(conductors.size());
    std::vector<PendingZ> pendingZ;   // Z-returns, planned after all conductors are built

    // Per-winding index (stable, MAS order) for the seam-azimuth stagger. MKF draws every winding's
    // seam/leads on the SAME reference plane, so different windings' terminal leads collide there.
    // Rotating each winding rigidly about the column axis by a distinct angle sends its leads out at a
    // different azimuth (rings unchanged, magnetics unchanged) so they no longer overlap. The first
    // winding keeps angle 0, so single-winding inductors are byte-for-byte unchanged.
    std::map<std::string,int> windingIdx;
    for (const auto& ct : conductors) windingIdx.emplace(ct.winding, (int)windingIdx.size());
    const int nWindings = (int)windingIdx.size();
    // Connection FACE per winding (Alf, 2026-08-07: "use isolation side"): windings on the
    // PRIMARY isolation side connect on the reference face (-Z, face 0); every other
    // isolation side (secondary, tertiary, ...) connects on the opposite face (+Z, face 1).
    // The split is the SAFETY-ISOLATION barrier -- mains-side terminals leave one side of
    // the component, output-side terminals the other -- not an arbitrary alternation by
    // winding order. A rect/oblong column maps onto itself only under 180 deg, so exactly
    // these two faces exist; round columns keep the single fan plane (aimed at the core
    // opening) and ignore this.
    std::map<std::string, int> windingFace;
    for (const auto& w : coil.get_functional_description())
        windingFace[w.get_name()] =
            (w.get_isolation_side() == MAS::IsolationSide::PRIMARY) ? 0 : 1;
    // Some cores open on ONE side only (Alf, 2026-08-07: EP/EPX-style plates wrap one Z
    // face). Before honouring the isolation-side faces, probe the core solids over both
    // +-Z connection corridors (the slot span x the winding's height band x the lead's
    // z reach). A blocked face takes no connections: every winding is forced to the open
    // face -- the collision gate arbitrates the crowding loudly, never the core. Both
    // corridors blocked -> refuse. This must run BEFORE emission: the dragback space
    // reservation and the lead attachment key on the face during prim construction.
    if (!isToroidal && !opts.coreObstacles.empty() &&
        (columnShape == MAS::ColumnShape::RECTANGULAR ||
         columnShape == MAS::ColumnShape::OBLONG)) {
        double maxR = 0.0, maxY = 0.0, maxOD = 0.0;
        for (const auto& ct : conductors) {
            if (ct.turns.empty()) continue;
            const MAS::Wire& w = wireMap.at(ct.winding);
            auto [ww, wh] = TurnBuilder::wireDimensions(w, *ct.turns.front(), true);
            maxOD = std::max(maxOD, std::max(ww, wh));
            for (const MAS::Turn* t : ct.turns) {
                maxR = std::max(maxR, station(t).x);
                maxY = std::max(maxY, std::abs(station(t).y));
            }
        }
        const double zStart = halfD + maxOD;      // just outside the winding's outer face
        // Probe out to the CORE's own z extent (its C dimension), never just the lead-tip
        // reach (Alf, 2026-08-07): a one-sided plate sits BEYOND the copper tip, and the
        // connection must exit the core envelope to be usable.
        double zCore = 0.0;
        for (const auto& obst : opts.coreObstacles) {
            Bnd_Box bx; BRepBndLib::Add(obst, bx);
            double x0, y0, z0, x1, y1, z1; bx.Get(x0, y0, z0, x1, y1, z1);
            zCore = std::max({zCore, std::abs(z0), std::abs(z1)});
        }
        const double zEnd = std::max(maxR + 4.0 * maxOD, zCore);
        const double xs = static_cast<double>(conductors.size()) * maxOD;
        auto corridorBlocked = [&](double sgn) {
            for (const auto& obst : opts.coreObstacles)
                for (TopExp_Explorer se(obst, TopAbs_SOLID); se.More(); se.Next()) {
                    BRepClass3d_SolidClassifier cls(se.Current());
                    for (int zi = 0; zi <= 8; ++zi)
                        for (int xi = -2; xi <= 2; ++xi)
                            for (int yi = -2; yi <= 2; ++yi) {
                                cls.Perform(gp_Pnt(xs * xi / 2.0, maxY * yi / 2.0,
                                                   sgn * (zStart + (zEnd - zStart) * zi / 8.0)),
                                            1e-7);
                                if (cls.State() == TopAbs_IN || cls.State() == TopAbs_ON)
                                    return true;
                            }
                }
            return false;
        };
        const bool blockNeg = corridorBlocked(-1.0);   // face 0: -Z, the reference face
        const bool blockPos = corridorBlocked(+1.0);   // face 1: +Z (seam rotated 180 deg)
        if (blockNeg && blockPos)
            throw std::runtime_error(
                "ConductorBuilder: the core blocks BOTH +-Z connection corridors "
                "(probed z=[" + std::to_string(zStart * 1e3) + "," +
                std::to_string(zEnd * 1e3) + "] mm) -- no rect-column terminal exit "
                "is possible");
        if (blockNeg || blockPos) {
            const int open = blockNeg ? 1 : 0;
            for (auto& [wname, f] : windingFace)
                if (f != open) {
                    std::cerr << "[ConductorBuilder] core blocks the "
                              << (blockNeg ? "-Z" : "+Z") << " connection corridor: winding '"
                              << wname << "' forced from its isolation-side face to the open "
                              << (open == 1 ? "+Z" : "-Z") << " face\n";
                    f = open;
                }
        }
    }

    // RECT-column DRAGBACKS (Alf, 2026-08-06), the round-column model transliterated:
    // a cross-layer return runs on the CONNECTION (-Z) face -- step out to the destination
    // layer's depth, descend axially there, land at the next crossing -- and RESERVES the
    // space of a layer: every turn at or outside the return's destination depth is DISPLACED
    // one laid wire OD further out on its -Z side (the -Z face moves out, the lateral +-X
    // faces extend to reach it; the extra lateral length joins the rising turn's pitch
    // distribution). Same-depth transitions need no routing at all (the rising turn).
    // Parallels' descents take staggered x slots on the face, one slot per conductor
    // (a conductor's own returns stack in depth, exactly like the round columns' layer
    // reuse). This replaces the out-of-window end-run entirely.
    struct RectReturn {
        size_t ci, trans;
        double srcZ, dstZ, diam, xSlot;
        int side;   // which z-side of SPACE its winding's connection face is on
                    // (0 = -Z, 1 = +Z after the per-winding 180-degree seam rotation)
    };
    std::vector<RectReturn> rectReturns;
    // (dstZ, diam) per DISTINCT level, PER SIDE OF SPACE: a return only reserves space on
    // the face it runs on. Every conductor's geometry then displaces on BOTH its z-sides --
    // its connection side by its own side's rides, its BACK side by the opposite side's
    // (the 180-degree-rotated windings interleave, and an undisplaced back face sat 77 um
    // from a displaced connection face on 02_flyback).
    std::vector<std::pair<double, double>> rectRideLevels[2];
    if (!isToroidal && columnShape == MAS::ColumnShape::RECTANGULAR) {
        std::map<size_t, double> slotOf;   // conductor -> x slot
        double maxDiam = 0.0;
        for (size_t cv = 0; cv < conductors.size(); ++cv) {
            const auto& ct = conductors[cv];
            const MAS::Wire& w = wireMap.at(ct.winding);
            auto [cw2, ch2] =
                TurnBuilder::wireDimensions(w, *ct.turns.front(), /*paintCoating=*/true);
            const double diam = std::min(cw2, ch2);
            auto [ew, eh] = TurnBuilder::wireDimensions(w, *ct.turns.front(), opts.paintCoating);
            const bool rectW = w.get_type() == MAS::WireType::RECTANGULAR ||
                               w.get_type() == MAS::WireType::PLANAR;
            const double rw = rectW ? 0.5 * std::hypot(ew, eh) : 0.5 * std::min(ew, eh);
            const double mb = rw * 1.02;
            if (std::getenv("MVB_DIAG")) {
                for (size_t i = 0; i < ct.turns.size(); ++i) {
                    const RectStation st =
                        rectStation(station(ct.turns[i]), halfW, halfD, mb, ct.winding);
                    std::cerr << "[rect] ci=" << cv << " turn " << i << " zPos=" << st.zPos
                              << " y=" << st.y << " segX=" << st.segX << " segZ=" << st.segZ
                              << "\n";
                }
            }
            for (size_t i = 0; i + 1 < ct.turns.size(); ++i) {
                const RectStation a =
                    rectStation(station(ct.turns[i]), halfW, halfD, mb, ct.winding);
                const RectStation b =
                    rectStation(station(ct.turns[i + 1]), halfW, halfD, mb, ct.winding);
                // ANY depth change is a return -- including a pure layer climb with no axial
                // move (its descent just has zero length). Requiring a y-move too let a
                // dy=0 layer change fall into the rising-turn branch, which ignores depth
                // (02_flyback's secondary: mismatched corners a wire OD apart).
                if (std::abs(b.zPos - a.zPos) > 1e-12) {
                    if (!slotOf.count(cv)) slotOf[cv] = 0.0;   // slot assigned below
                    rectReturns.push_back({cv, i, a.zPos, b.zPos, diam, 0.0,
                                           windingFace.at(ct.winding)});
                    maxDiam = std::max(maxDiam, diam);
                }
            }
        }
        // One x slot per conductor, centre-out around the crossing (x = 0), a coated OD
        // apart (layout criterion: insulation touching).
        int k = 0;
        for (auto& [cv, slot] : slotOf) {
            const int step = (k + 1) / 2;
            slot = (k % 2 == 1 ? +1.0 : -1.0) * step * maxDiam;
            if (k == 0) slot = 0.0;
            ++k;
        }
        for (auto& r : rectReturns) r.xSlot = slotOf.at(r.ci);
        // Distinct destination depths, ascending, PER SIDE. Each level displaces everything
        // at or beyond it (on that side) by the level's largest laid wire OD.
        for (const auto& r : rectReturns) {
            auto& levels = rectRideLevels[r.side];
            bool merged = false;
            for (auto& lv : levels) {
                if (std::abs(lv.first - r.dstZ) <= 0.5 * r.diam) {
                    lv.second = std::max(lv.second, r.diam);
                    merged = true;
                    break;
                }
            }
            if (!merged) levels.push_back({r.dstZ, r.diam});
        }
        std::sort(rectRideLevels[0].begin(), rectRideLevels[0].end());
        std::sort(rectRideLevels[1].begin(), rectRideLevels[1].end());
    }
    // Displacement of a turn's face on the given SIDE of space: one OD per ride level at or
    // inside its depth on that side.
    auto rectRideFor = [&](double zPos, int side) {
        double ride = 0.0;
        for (const auto& lv : rectRideLevels[side])
            if (lv.first <= zPos + 0.5 * lv.second) ride += lv.second;
        return ride;
    };

    // Parallel strands per winding: the tangential spread has to divide the core's opening
    // between every strand of every winding, so it needs the counts up front.
    std::map<std::string,int> nParallels;
    for (const auto& ct : conductors)
        nParallels[ct.winding] = std::max(nParallels[ct.winding], (int)ct.parallel + 1);
    // The spread itself is budgeted per conductor from the core's measured opening (see the
    // seam-angle block below); ~2/3 of the circle is only its unconstrained upper bound.

    // EVERY dragback in the window, from EVERY conductor, collected BEFORE anything is built.
    // A return is a physical obstruction in the winding build: whatever is wound at or outside
    // its radius rides over it, whether or not it belongs to the same winding or the same
    // parallel. Collecting them per-conductor (the old behaviour) meant a return laid by
    // Primary parallel 0 was invisible to parallel 2 and to the Secondary, which wound straight
    // through it (06_llc_xfmr_eq4128: 487 um between parallel 0's descent and parallel 2's
    // wrap, against 500 um of copper). The radius rule below still keeps ordering physical --
    // only turns at or outside a return's radius ride over it, and those are exactly the ones
    // wound after it.
    // THE VERTICAL FAN (Alf, 2026-08-05). Every vertical run in the window -- terminal
    // connections and dragback returns -- used to sit exactly on the station plane, where they
    // collided: both parallels of a winding dragging back through each other (06/14/23), sibling
    // entrance leads overlapping (11/24), two windings' leads too close (13). They are now
    // spread in a TIGHT fan around the plane: each vertical takes the smallest azimuth offset
    // that clears every other vertical it could actually touch. Verticals in the SAME radial
    // band need distinct angles; verticals of different bands (returns of different layers,
    // leads at well-separated heights) reuse the same angle -- except a lead RUNS THROUGH every
    // layer radius on its way out, so leads are also cleared against any return whose axial span
    // covers the lead's height.
    struct LaidDragback { double az; double radius; double diam; };
    std::vector<LaidDragback> laidDragbacks;
    std::map<std::pair<size_t, size_t>, double> dragAzOf;   // (conductor, transition) -> azimuth
    std::map<size_t, double> leadAzIn, leadAzOut;           // conductor -> lead azimuth
    double fanWidth = 0.0;
    if (effectivelyRound) {
        struct Vert {          // one azimuth-slot request
            int kind;          // 0 = terminal lead GROUP, 1 = dragback
            size_t ci;         // conductor index (dragbacks)
            size_t trans;      // transition index (dragbacks)
            bool entrance;     // leads
            double r;          // innermost radius of the vertical run
            double y0, y1;     // axial extent of the run (y0 <= y1)
            double rw;         // COATED wire radius -- slot clearance is physical
            double rwBare;     // bare-copper radius -- what the collision gate enforces
            // Leads: one vert per (conductor, side). `segs` is the MKF-drawn 2D lead
            // route in the (r, y) half-plane (vertical stub + horizontal run -- Alf,
            // 2026-08-07: the connection positions are MKF's, never invented), used for
            // exact capsule feasibility. ys = the RUN row(s); attachY = the crossing row
            // the lead attaches at; wname groups parallels of one winding so the packer
            // keeps them in ADJACENT slots.
            std::vector<double> ys;
            std::vector<size_t> cis;
            std::vector<std::array<double, 4>> segs;   // (r0, y0, r1, y1)
            double destAdvance = 0.0;  // dragbacks: the climb of the layer they land on
            double attachY = 0.0;
            std::string wname;
        };
        // Min distance between two 2D segment sets -- the exact clearance criterion for
        // two lead routes (or a route and a dragback descent) sharing an azimuth plane.
        auto seg2Dist = [](const std::array<double, 4>& a, const std::array<double, 4>& b) {
            auto clamp01 = [](double t) { return t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t); };
            const double ax = a[2] - a[0], ay = a[3] - a[1];
            const double bx = b[2] - b[0], by = b[3] - b[1];
            double best = 1e30;
            // sample-free segment/segment distance via endpoint-projection (exact for the
            // axis-aligned L-routes involved here; conservative-tight otherwise)
            auto ptSeg = [&](double px, double py, const std::array<double, 4>& s2,
                             double sx, double sy) {
                const double L2 = sx * sx + sy * sy;
                const double t = L2 > 0.0
                    ? clamp01(((px - s2[0]) * sx + (py - s2[1]) * sy) / L2) : 0.0;
                return std::hypot(px - (s2[0] + sx * t), py - (s2[1] + sy * t));
            };
            best = std::min(best, ptSeg(a[0], a[1], b, bx, by));
            best = std::min(best, ptSeg(a[2], a[3], b, bx, by));
            best = std::min(best, ptSeg(b[0], b[1], a, ax, ay));
            best = std::min(best, ptSeg(b[2], b[3], a, ax, ay));
            // proper crossing check
            auto cross = [](double ux, double uy, double vx, double vy) {
                return ux * vy - uy * vx;
            };
            const double d1 = cross(ax, ay, b[0] - a[0], b[1] - a[1]);
            const double d2 = cross(ax, ay, b[2] - a[0], b[3] - a[1]);
            const double d3 = cross(bx, by, a[0] - b[0], a[1] - b[1]);
            const double d4 = cross(bx, by, a[2] - b[0], a[3] - b[1]);
            if (d1 * d2 < 0.0 && d3 * d4 < 0.0) best = 0.0;
            return best;
        };
        auto routeDist = [&](const std::vector<std::array<double, 4>>& A,
                             const std::vector<std::array<double, 4>>& B) {
            double best = 1e30;
            for (const auto& sa : A)
                for (const auto& sb : B) best = std::min(best, seg2Dist(sa, sb));
            return best;
        };
        std::vector<Vert> verts;
        struct WireRow {       // one station row, with its conductor's signed advance
            double r, y, adv, rw;   // rw = BARE radius (the gate's criterion)
            size_t ci;
            // Azimuth validity of this row's copper in the fan frame (crossings sit at
            // c = 0): the FIRST station's wrap exists only for c >= 0, the LAST only for
            // c <= 0 -- the winding ends there. Without the bounds, the last row's drift
            // band extrapolated PHANTOM copper past the final crossing, covered the whole
            // reach, and the fan dropped a constraint it could in fact satisfy at c > 0
            // (11_pushpull: Secondary 1 par 0's stub is clear just past par 1's last
            // crossing, where par 1's ended wrap no longer exists).
            double cLo = -1e30, cHi = 1e30;
        };
        std::vector<WireRow> rows;   // NB: WireRow.rw carries the BARE radius: rows are MKF's
                                     // geometry, and lead corridors past them are planned at
                                     // the full bare separation (see the drift model below)
        for (size_t cv = 0; cv < conductors.size(); ++cv) {
            const auto& ct = conductors[cv];
            const MAS::Wire& w = wireMap.at(ct.winding);
            auto [cw, chh] = TurnBuilder::wireDimensions(w, *ct.turns.front(),
                                                         /*paintCoating=*/true);
            const double rwCoat = 0.5 * std::min(cw, chh);
            auto [bw, bh] = TurnBuilder::wireDimensions(w, *ct.turns.front(),
                                                        /*paintCoating=*/false);
            const double rwBare = 0.5 * std::min(bw, bh);
            // Detection mirrors the emission loop EXACTLY (same wire footprint, same median
            // pitch), so every dragback the loop emits has a fan azimuth waiting for it.
            auto [ew, eh] = TurnBuilder::wireDimensions(w, *ct.turns.front(), opts.paintCoating);
            const bool rectW = w.get_type() == MAS::WireType::RECTANGULAR ||
                               w.get_type() == MAS::WireType::PLANAR;
            const double rwEmit = rectW ? 0.5 * std::hypot(ew, eh) : 0.5 * std::min(ew, eh);
            double mp = 0.0, adv = 0.0;
            {
                std::vector<double> pitches;
                int pos = 0, neg = 0;
                for (size_t i = 0; i + 1 < ct.turns.size(); ++i) {
                    PlanePt a = station(ct.turns[i]), b = station(ct.turns[i + 1]);
                    if (std::abs(b.x - a.x) <= rwEmit) {
                        const double dy = b.y - a.y;
                        pitches.push_back(std::abs(dy));
                        (dy >= 0.0 ? pos : neg)++;
                    }
                }
                if (!pitches.empty()) {
                    std::sort(pitches.begin(), pitches.end());
                    mp = pitches[pitches.size() / 2];
                    if (pos == 0 || neg == 0) adv = (pos > 0 ? mp : -mp);
                }
            }
            // Mirrors the emission loop: every station is emitted (no trailing prune), so
            // every transition claims its fan slot and every station is an obstacle row.
            const size_t nEmitP = ct.turns.size();
            std::vector<PlanePt> stationsP;
            for (const MAS::Turn* t : ct.turns) stationsP.push_back(station(t));
            const auto bandsP = computePitchBands(stationsP, rwEmit);
            for (size_t i = 0; i + 1 < nEmitP; ++i) {
                PlanePt a = station(ct.turns[i]), b = station(ct.turns[i + 1]);
                const PitchBand bandP = bandAt(bandsP, a.x, rwEmit, mp, adv);
                if (isZReturn(a, b, rwEmit, bandP.medianPitch, bandP.advance)) {
                    Vert dv;
                    dv.kind = 1; dv.ci = cv; dv.trans = i; dv.entrance = false;
                    dv.r = b.x; dv.y0 = std::min(a.y, b.y); dv.y1 = std::max(a.y, b.y);
                    dv.rw = rwCoat; dv.rwBare = rwBare;
                    // the climb of the layer this return LANDS on (its own band, not the
                    // conductor's net direction -- a serpentine conductor has none)
                    dv.destAdvance = bandAt(bandsP, b.x, rwEmit, 0.0, 0.0).advance;
                    verts.push_back(std::move(dv));
                }
            }
            const PlanePt pf = station(ct.turns.front()), pl = station(ct.turns[nEmitP - 1]);
            // MKF-DRAWN LEAD ROUTES (Alf, 2026-08-07): the fan models EXACTLY the
            // waypoints the emission will sweep -- MKF's drawn vertical connection and
            // run, never an invented straight-out. One vert per (conductor, side);
            // parallels of a winding keep ADJACENT slots via the packer's anchoring.
            std::vector<const RSpace*> tRects;
            for (const auto& sp2 : drawn)
                if (sp2.winding == ct.winding && sp2.parallel == ct.parallel && sp2.isTerminal)
                    tRects.push_back(&sp2);
            const std::string whoV =
                ct.winding + " parallel " + std::to_string(ct.parallel);
            auto [egrp, xgrp] = splitTerminalGroups(tRects, whoV);
            auto routeVert = [&](const std::vector<PlanePt>& wpv, bool entrance,
                                 const PlanePt& attach) {
                Vert v;
                v.kind = 0; v.ci = cv; v.trans = 0; v.entrance = entrance;
                v.rw = rwCoat; v.rwBare = rwBare;
                v.r = 1e30; v.y0 = 1e30; v.y1 = -1e30;
                for (const auto& pw : wpv) {
                    v.r = std::min(v.r, pw.x);
                    v.y0 = std::min(v.y0, pw.y);
                    v.y1 = std::max(v.y1, pw.y);
                }
                for (size_t i2 = 0; i2 + 1 < wpv.size(); ++i2)
                    v.segs.push_back({wpv[i2].x, wpv[i2].y, wpv[i2 + 1].x, wpv[i2 + 1].y});
                v.ys = {wpv.back().y};
                v.cis = {cv};
                v.attachY = attach.y;
                v.wname = ct.winding;
                verts.push_back(std::move(v));
            };
            if (std::getenv("MVB_LEAD_DIAG")) {
                auto dump = [&](const char* side, const std::vector<const RSpace*>& g,
                                const std::vector<PlanePt>& wpv) {
                    std::fprintf(stderr, "[lead-route] %s %s:", whoV.c_str(), side);
                    for (const RSpace* r2 : g)
                        std::fprintf(stderr, " rect(%.4f x %.4f @ %.4f,%.4f)",
                                     r2->dimensions.at(0) * 1e3, r2->dimensions.at(1) * 1e3,
                                     r2->coordinates.at(0) * 1e3, r2->coordinates.at(1) * 1e3);
                    std::fprintf(stderr, " -> wp:");
                    for (const auto& pw : wpv)
                        std::fprintf(stderr, " (%.4f,%.4f)", pw.x * 1e3, pw.y * 1e3);
                    std::fprintf(stderr, "\n");
                };
                dump("IN", egrp, terminalWaypoints(egrp, pf, whoV + " entrance"));
                if (!xgrp.empty())
                    dump("OUT", xgrp, terminalWaypoints(xgrp, pl, whoV + " exit"));
            }
            routeVert(terminalWaypoints(egrp, pf, whoV + " entrance"), true, pf);
            if (!xgrp.empty())
                routeVert(terminalWaypoints(xgrp, pl, whoV + " exit"), false, pl);
            else   // one drawn lead: synthesized straight-out exit (see splitTerminalGroups)
                routeVert({{pl.x, pl.y}, {pl.x + 1.0, pl.y}}, false, pl);
            for (size_t i = 0; i < nEmitP; ++i) {
                const PlanePt st = station(ct.turns[i]);
                rows.push_back({st.x, st.y, adv, rwBare, cv,
                                i == 0 ? 0.0 : -1e30,
                                i + 1 == nEmitP ? 0.0 : 1e30});
            }
        }
        // Angular clearance two verticals need between them, taken at the innermost radius where
        // both exist; 0 when they can never touch. The LAYOUT criterion is the COATED radii
        // sum -- wires may touch their insulation, never overlap it, exactly how MKF packs
        // adjacent turns. It sits strictly above the collision gate's fault criterion
        // (gateMinSeparation on bare copper), so fan-placed geometry can never sit AT the
        // gate's threshold where sampling noise flips the verdict.
        auto need = [&](const Vert& A, const Vert& B) -> double {
            const double dSep = A.rw + B.rw;
            if (A.kind == 0 && B.kind == 0) {
                // Two MKF-drawn lead routes: exact 2D capsule distance between the drawn
                // polylines. Routes MKF reserved disjoint (different windings' corridors)
                // coexist at ONE azimuth; only genuinely overlapping routes -- e.g. two
                // parallels' stubs sharing a layer band (11_pushpull's Secondary 1) --
                // take distinct, adjacent slots.
                if (routeDist(A.segs, B.segs) > dSep) return 0.0;
            } else if (A.kind != B.kind) {
                const Vert& L = A.kind == 0 ? A : B;
                const Vert& D = A.kind == 0 ? B : A;
                // The lead's drawn route against the dragback's vertical descent.
                std::vector<std::array<double, 4>> dseg{{D.r, D.y0, D.r, D.y1}};
                if (routeDist(L.segs, dseg) > dSep) return 0.0;
            } else {
                // Two dragbacks: a conductor's OWN returns deliberately stack in one column
                // (layer reuse -- the snap below puts them there), so they never need an
                // angle between them. Returns of DIFFERENT conductors always take distinct
                // angles: sharing a column stacks their raises into the other winding's
                // ride-over (23_interleaved: the secondary's exit rode 1.97 mm -- its own
                // return piled onto the primary's column -- while its entrance rode 1.07).
                if (A.ci == B.ci) return 0.0;
            }
            // Perpendicular clearance, not the chord: a vertical's RUN-OUT is a radial line
            // through the neighbour's band, so the closest approach between two slots dAz
            // apart is r*sin(dAz) at the innermost radius both occupy -- the chord formula
            // left a same-band dragback pair exactly 0.3 um under the gate.
            const double rr = std::max(std::min(A.r, B.r), dSep);
            return std::asin(std::min(1.0, dSep / rr));
        };
        // CENTRE-OUT symmetric packing (Alf, 2026-08-05): every vertical takes the FEASIBLE
        // azimuth CLOSEST TO THE STATION PLANE (0 deg), so the fan hugs the plane and only
        // spreads as far as real conflicts force it -- the bump treatment (one raised region,
        // steps at the XY-plane crossings) is centred there, and a vertical that drifts far
        // off-plane would slide out from under it. Dragbacks are placed FIRST (they live
        // inside the winding, under every bump); the leads pack around them, ordered by height
        // so close-height pairs meet their real clearance and far pairs share angles freely.
        std::stable_sort(verts.begin(), verts.end(), [](const Vert& a, const Vert& b) {
            if (a.kind != b.kind) return a.kind > b.kind;
            if (a.kind == 0) {
                // SHORTEST routes place first: a straight/short lead is feasible anywhere
                // near the plane, while a long stub must dodge the wrap tails that the
                // already-placed leads PIN at their crossings (a lead owns its crossing).
                // 11_pushpull S1: par 1's 0.65 mm down-stub takes the plane; par 0's
                // 2.9 mm stub then packs just past par 1's tail -- placed the other way
                // round, par 1 anchored next to par 0 and dragged its final wrap tail
                // straight onto par 0's stub.
                const double ea = a.y1 - a.y0, eb = b.y1 - b.y0;
                if (std::abs(ea - eb) > 1e-12) return ea < eb;
            }
            return a.y0 < b.y0;
        });
        // Forbidden azimuth intervals per LEAD: a lead runs radially at its own row, crossing
        // every ring band outside its attachment, and the wire rising from a nearby crossed
        // row DRIFTS in y with azimuth (its wrap's pitch). Where that drift brings the wire
        // within clearance of the lead's row, the azimuth is forbidden -- this is what put
        // 06_llc's P2 entrance 0.402 mm under P0's outer rows at the plane, and squeezed
        // 23's secondary entrances between their siblings' rising first wraps.
        struct AzInterval { double lo, hi; };
        // Validity bound of the fan's linearized models: half-way to the XY-plane crossings,
        // where the ride-over halves change regime. The fan-vs-core-opening guard and the
        // bump-boundary guard both throw far earlier for any real design.
        const double kFanReach = kPi / 4.0;
        std::vector<std::vector<AzInterval>> forbid(verts.size());
        for (size_t k = 0; k < verts.size(); ++k) {
            const Vert& L = verts[k];
            if (L.kind != 0) continue;
            for (const auto& R : rows) {
                // Rows radially inside the attachment are never crossed; the raises preserve
                // radial order (a taller stack lies under everything further out), so the
                // bare-radius test stays valid for the lifted lead too.
                if (R.r + R.rw < L.r - L.rw - 1e-12) continue;
                if (std::abs(R.adv) < 1e-12) continue;   // no consistent advance: gate arbitrates
                // Rows are MKF's own geometry: the corridor past them is planned at the FULL
                // BARE separation -- the gate's sag allowance (kMaxSagFraction) then remains
                // as the buffer between plan and verdict, with no invented factor. (The
                // coated criterion here forbade corridors the gate accepts and drove groups
                // to the reach edge.)
                const double clr = L.rwBare + R.rw;
                // ANY member's attachment row rides WITH the group's slot (the wrap it feeds
                // starts at that azimuth), so it constrains nothing in absolute azimuth --
                // the drift model only applies to rows whose crossings sit elsewhere. Without
                // this the group forbade itself and flew to the reach edge.
                const bool attachRow = R.ci == L.cis[0] &&
                                       std::abs(R.y - L.attachY) < 0.5 * (L.rwBare + R.rw);
                if (attachRow) continue;
                // Forbidden windows come from BOTH parts of the drawn route: the RUN row
                // (the radial run crosses the ring band at that height) and each VERTICAL
                // STUB (the stub occupies its whole y-interval at its own radius -- a
                // sibling's rising wrap sweeps through it; 11_pushpull: Secondary 1 par 0's
                // exit stub against par 1's 4.85 mm-pitch final wrap).
                std::vector<std::pair<double, double>> bands;
                for (size_t m = 0; m < L.ys.size(); ++m)
                    bands.push_back({L.ys[m] - clr, L.ys[m] + clr});
                for (const auto& sg : L.segs) {
                    if (std::abs(sg[0] - sg[2]) > 1e-12) continue;      // vertical stubs only
                    if (std::abs(R.r - sg[0]) > clr) continue;          // row off the stub radius
                    bands.push_back({std::min(sg[1], sg[3]) - clr, std::max(sg[1], sg[3]) + clr});
                }
                for (const auto& [bLo, bHi] : bands) {
                    const double a1 = (bLo - R.y) * kTwoPi / R.adv;
                    const double a2 = (bHi - R.y) * kTwoPi / R.adv;
                    AzInterval iv{std::min(a1, a2), std::max(a1, a2)};
                    // Clip to where this row's copper actually exists. A CLAMPED endpoint
                    // is the winding's END CROSSING, not a clearance boundary -- the end
                    // section still sits there, so when its y lies inside the forbidden
                    // band the interval extends an angular clearance PAST the crossing
                    // (point-obstacle model, the same asin criterion as need()).
                    const bool endIn = R.y >= bLo && R.y <= bHi;
                    const double dAng =
                        std::asin(std::min(1.0, clr / std::max(R.r, clr)));
                    if (iv.lo < R.cLo) iv.lo = endIn ? R.cLo - dAng : R.cLo;
                    if (iv.hi > R.cHi) iv.hi = endIn ? R.cHi + dAng : R.cHi;
                    if (iv.lo >= iv.hi) continue;
                    if (iv.hi < -kFanReach || iv.lo > kFanReach) continue;
                    // The linear drift model only means anything within the fan's reach: a
                    // pair whose forbidden band covers the whole reach cannot be fixed by
                    // azimuth at all -- leave it to the gate rather than launching the lead
                    // to a nonsense angle. Intervals are kept UNCLAMPED: truncating one at
                    // the reach edge fabricated a "feasible" boundary slot there while the
                    // real forbidden band continued past it.
                    if (iv.lo <= -kFanReach && iv.hi >= kFanReach) continue;
                    forbid[k].push_back(iv);
                }
            }
        }
        std::vector<double> az(verts.size(), 0.0);
        for (size_t k = 0; k < verts.size(); ++k) {
            // HARD constraints: the vertical-vs-vertical needs -- geometry the fan itself
            // places, non-negotiable. SOFT: the row-drift windows -- a linearized model of
            // MKF's wraps on the layout criterion, deliberately stricter than the gate; when
            // NO azimuth satisfies both, the vertical takes the tightest hard-feasible slot
            // and the collision gate keeps the final word (binding the model launched groups
            // to the reach edge for corridors the gate accepts).
            auto hardOk = [&](double c) {
                for (size_t j = 0; j < k; ++j)
                    if (std::abs(c - az[j]) + 1e-12 < need(verts[j], verts[k])) return false;
                return true;
            };
            auto softOk = [&](double c) {
                for (const auto& iv : forbid[k])
                    if (c > iv.lo + 1e-12 && c < iv.hi - 1e-12) return false;
                return true;
            };
            bool have = false;
            double best = 0.0;
            // Same-winding parallels' leads go TOGETHER (Alf, 11_pushpull: Secondary 1
            // parallel 0 next to parallel 1, not among the other windings): once one
            // member of a (winding, side) group is placed, its siblings pack NEXT TO
            // that anchor -- the proximity metric becomes distance to the anchor
            // instead of distance to the plane.
            double anchor = 0.0;
            bool anchored = false;
            if (verts[k].kind == 0)
                for (size_t j = 0; j < k && !anchored; ++j)
                    if (verts[j].kind == 0 && verts[j].wname == verts[k].wname &&
                        verts[j].entrance == verts[k].entrance) {
                        anchor = az[j];
                        anchored = true;
                    }
            // Alf, 2026-08-07 (14_dab): "they must have the same x coord" -- every parallel of
            // a winding leaves on the SAME azimuth, so the connections are aligned and every
            // strand keeps the same length. The anchor is therefore BINDING, not a preference:
            // once one member of a (winding, side) group is placed the rest take that exact
            // slot. If that is infeasible the collision gate says so -- the leads are never
            // silently fanned apart (Primary parallel 2 sat 2.6 deg off its siblings' plane).
            if (anchored) {
                az[k] = anchor;
                continue;
            }
            auto consider = [&](double c, bool withSoft) {
                if (!hardOk(c) || (withSoft && !softOk(c))) return;
                // Closest to the anchor (plane when unanchored) wins; ties go positive.
                const double ref = anchored ? std::abs(c - anchor) : std::abs(c);
                const double refBest = anchored ? std::abs(best - anchor) : std::abs(best);
                if (!have || ref < refBest - 1e-12 || (ref <= refBest + 1e-12 && c > best)) {
                    best = c;
                    have = true;
                }
            };
            // A conductor's returns SNAP onto the azimuth its earlier returns already took
            // whenever that is feasible: different-layer returns of one conductor stack in ONE
            // column (the layer-reuse rule, and what tallestBumpColumn's raise accounting
            // needs). Without the snap they landed fractions of a degree apart -- neither
            // stacked nor separated, and the ride-over raise under-counted the pile.
            for (int tier = 0; tier < 2 && !have; ++tier) {
                if (verts[k].kind != 1) break;
                for (size_t j = 0; j < k; ++j)
                    if (verts[j].kind == 1 && verts[j].ci == verts[k].ci)
                        consider(az[j], tier == 0);
            }
            std::vector<double> cand{0.0};
            if (anchored) cand.push_back(anchor);   // disjoint-route sibling: same azimuth
            for (size_t j = 0; j < k; ++j) {
                const double nd = need(verts[j], verts[k]);
                if (nd > 0.0) {
                    cand.push_back(az[j] + nd);
                    cand.push_back(az[j] - nd);
                }
            }
            for (const auto& iv : forbid[k]) {
                // Only TRUE interval endpoints inside the reach are meaningful slots.
                if (std::abs(iv.lo) <= kFanReach) cand.push_back(iv.lo);
                if (std::abs(iv.hi) <= kFanReach) cand.push_back(iv.hi);
            }
            for (int tier = 0; tier < 2 && !have; ++tier)
                for (double c : cand) consider(c, tier == 0);
            // max_j(az_j + need_jk) always satisfies the hard tier, so a slot always exists.
            az[k] = best;
        }
        // CLIMB-AWARE SLOT ORDER for sibling dragbacks (Alf's 14_dab). The slots are correct in
        // SEPARATION but their ORDER matters once the layer they land on climbs steeply. Two
        // parallels whose crossings sit dAz apart have their helices phase-shifted, so their
        // constant separation becomes  -(dY) + pitch*dAz/2pi : the climb either adds to the
        // station spacing or eats it. Assign the SAME set of slots so azimuth runs OPPOSITE to
        // the destination-station height -- then the shift always ADDS. (14_dab's secondary
        // returns onto a 15 mm-pitch layer: 2.29 deg of slot spacing ate 0.095 mm of an 0.855 mm
        // budget, leaving 0.755 against a 0.784 requirement.)
        // Permuting an already-validated set of slots cannot change any pairwise separation --
        // need() between two dragbacks depends only on their radii and wire, not on which
        // conductor holds which slot. Restricted to bands where each conductor has exactly ONE
        // return, so the consecutive-return snap below still holds.
        {
            std::vector<size_t> dragIdx;
            for (size_t k = 0; k < verts.size(); ++k)
                if (verts[k].kind == 1) dragIdx.push_back(k);
            std::vector<char> done(dragIdx.size(), 0);
            for (size_t a0 = 0; a0 < dragIdx.size(); ++a0) {
                if (done[a0]) continue;
                std::vector<size_t> group{dragIdx[a0]};
                done[a0] = 1;
                for (size_t b0 = a0 + 1; b0 < dragIdx.size(); ++b0) {
                    if (done[b0] ||
                        std::abs(verts[dragIdx[b0]].r - verts[dragIdx[a0]].r) >
                            verts[dragIdx[a0]].rw)
                        continue;
                    group.push_back(dragIdx[b0]);
                    done[b0] = 1;
                }
                if (group.size() < 2) continue;
                std::set<size_t> cis;
                for (size_t g : group) cis.insert(verts[g].ci);
                if (cis.size() != group.size()) continue;
                // the climb of the layer these returns land on, from their own band
                double advance = 0.0;
                for (size_t g : group)
                    if (std::abs(verts[g].destAdvance) > std::abs(advance))
                        advance = verts[g].destAdvance;
                if (std::abs(advance) < 1e-12) continue;
                std::vector<double> slots;
                for (size_t g : group) slots.push_back(az[g]);
                std::sort(slots.begin(), slots.end());
                // destination station height (a return lands at its low end when it descends)
                std::sort(group.begin(), group.end(), [&](size_t x, size_t y) {
                    return advance > 0.0 ? verts[x].y0 > verts[y].y0
                                         : verts[x].y0 < verts[y].y0;
                });
                for (size_t g = 0; g < group.size(); ++g) az[group[g]] = slots[g];
            }
        }
        double lo = 0.0, hi = 0.0;
        for (double a : az) { lo = std::min(lo, a); hi = std::max(hi, a); }
        fanWidth = hi - lo;
        for (size_t k = 0; k < verts.size(); ++k) {
            const double a = kPlaneAz + az[k];
            if (verts[k].kind == 1) {
                dragAzOf[{verts[k].ci, verts[k].trans}] = a;
            } else {
                // The group's slot binds EVERY parallel of the winding on this side.
                for (size_t m : verts[k].cis) {
                    if (verts[k].entrance) leadAzIn[m] = a;
                    else                   leadAzOut[m] = a;
                }
            }
        }
        // CONSECUTIVE dragback transitions of one conductor (a single-turn layer) chain end to
        // start: force them onto one azimuth so the chain stays connected.
        for (auto& [key, a] : dragAzOf) {
            if (key.second == 0) continue;
            auto prev = dragAzOf.find({key.first, key.second - 1});
            if (prev != dragAzOf.end()) a = prev->second;
        }
        for (const auto& v : verts)
            if (v.kind == 1)
                laidDragbacks.push_back({dragAzOf.at({v.ci, v.trans}), v.r, 2.0 * v.rw});
        if (std::getenv("MVB_DRAG_DIAG")) {
            std::cerr << "[fan] " << verts.size() << " verticals, width "
                      << fanWidth * 180.0 / kPi << " deg\n";
            for (size_t k = 0; k < verts.size(); ++k)
                for (const auto& iv : forbid[k])
                    std::cerr << "[fan-forbid] k=" << k << " kind=" << verts[k].kind
                              << " [" << iv.lo * 180.0 / kPi << "," << iv.hi * 180.0 / kPi
                              << "] deg\n";
            for (size_t k = 0; k < verts.size(); ++k)
                std::cerr << "[fan]   " << (verts[k].kind ? "drag" : (verts[k].entrance
                          ? "lead-in" : "lead-out")) << " ci=" << verts[k].ci
                          << " trans=" << verts[k].trans << " r=" << verts[k].r * 1e3
                          << "mm y=[" << verts[k].y0 * 1e3 << "," << verts[k].y1 * 1e3
                          << "]mm az=" << az[k] * 180.0 / kPi
                          << " deg off-plane\n";
        }
    }
    // Tallest vertical-fan column overall: every lead must clear it on its way out, so the
    // common tip radius sits this much further beyond the winding.
    double fanMaxRaise = 0.0;
    {
        std::vector<std::pair<double, double>> cols;
        for (const auto& ld : laidDragbacks) {
            bool found = false;
            for (auto& c : cols) {
                if (std::abs(std::remainder(ld.az - c.first, kTwoPi)) < 1e-9) {
                    c.second += ld.diam;
                    found = true;
                    break;
                }
            }
            if (!found) cols.push_back({ld.az, ld.diam});
        }
        for (const auto& c : cols) fanMaxRaise = std::max(fanMaxRaise, c.second);
    }

    for (size_t ci = 0; ci < conductors.size(); ++ci) {
        const auto& ct = conductors[ci];
        const MAS::Wire& wire = wireMap.at(ct.winding);
        const MAS::WireType wireType = wire.get_type();
        const bool rectWire = wireType == MAS::WireType::RECTANGULAR ||
                              wireType == MAS::WireType::PLANAR;
        // Rectangular/planar wire: ROUND column -> one body via fixed binormal; RECT/OBLONG column
        // or TOROID -> per-primitive rect solids fused (emitRectColumn), with the section oriented
        // by the local bend/spacing axis (column Y, or the toroid's azimuthal tangent). Only FOIL
        // (a single wide sheet, a different construction) still throws.
        if (wireType == MAS::WireType::FOIL) {
            throw std::runtime_error(
                "ConductorBuilder: real-winding does not support FOIL wire (a single wide sheet "
                "turn is a different construction) — winding '" + ct.winding + "'");
        }
        // Resolve through the first turn: turn.dimensions carries the OUTER footprint
        // (the context-less overload would fall back to the round wire's CONDUCTING
        // diameter for paintCoating=true). For rectangular wire dimensions[0]=width=RADIAL,
        // dimensions[1]=height=AXIAL (matches TurnBuilder::build_rect_profile).
        auto [wireW, wireH] =
            TurnBuilder::wireDimensions(wire, *ct.turns.front(), opts.paintCoating);
        // Envelope/bend radius: the round wire's is its radius; a rectangle's is its half-DIAGONAL
        // (the largest centre-to-corner reach, so the collision corridor and minimum bend stay
        // conservative for the rotated section).
        double wireRadius = rectWire ? 0.5 * std::hypot(wireW, wireH)
                                     : std::min(wireW, wireH) / 2.0;
        // Bare-copper footprint (paintCoating=false) for the collision gate — see ConductorPath::cond*.
        auto [copW, copH] = TurnBuilder::wireDimensions(wire, *ct.turns.front(), /*paintCoating=*/false);
        // Same minimum-bend rule as build_concentric_rect_column_turn: a swept corner
        // self-intersects when the arc radius is below the profile's radial half-extent.
        double minBend = wireRadius * 1.02;

        ConductorPath path;
        path.name = ct.winding + " parallel " + std::to_string(ct.parallel);
        path.wireRadius = wireRadius;
        path.condRadius = rectWire ? 0.5 * std::hypot(copW, copH) : std::min(copW, copH) / 2.0;
        path.condWidth = copW; path.condHeight = copH;
        path.isRectangular = rectWire;
        path.wireWidth = wireW;
        path.wireHeight = wireH;
        path.femReady = opts.femReady;
        if (rectWire) {
            // Round (or straight-less oblong) column: the section sweeps cleanly with a fixed
            // binormal -> one body. A rect/oblong column with real corners, OR a toroid: per-
            // primitive rect solids + fuse (emitRectColumn) -- the wide flat section can't sweep the
            // corners, and the toroid's poloidal elbows carry the section on the azimuthal axis.
            path.singleBodyCapable = !isToroidal && effectivelyRound;
            path.useRectSolids = isToroidal || !effectivelyRound;
            path.toroidal = isToroidal;
        } else {
            // Round / litz wire. Round & oblong columns sweep as one clean single body. A
            // RECTANGULAR column goes through the per-primitive ANALYTIC path (SEG -> cylinder,
            // ARC3 -> torus segment) and fuses to ONE FEM-ready solid -- the per-run swept-pipe
            // compound leaves un-fusable seam slivers there. A round TOROID also becomes ONE solid,
            // swept with the simple MakePipe (its framing keeps the section centred on the high-
            // torsion hole-threading spine, where corrected Frenet drifts a crossing off-centre); if
            // its spread spine is too complex for MakePipe to close, emitConductor drops to the exact
            // per-run compound (still crossing-exact, just multi-solid).
            path.singleBodyCapable = isToroidal ||
                                     (columnShape == MAS::ColumnShape::ROUND ||
                                      columnShape == MAS::ColumnShape::OBLONG);
            path.toroidal = isToroidal;
            bool rectColumnAnalytic = !isToroidal && columnShape == MAS::ColumnShape::RECTANGULAR;
            path.useRectSolids = rectColumnAnalytic;
            path.roundProfile = rectColumnAnalytic;
        }

        const auto& turns = ct.turns;

        // This conductor's drawn terminal rects, in MKF emission order (verified in
        // get_connection_reserved_spaces): the ENTRANCE lead's rects first, then the
        // EXIT lead's. Non-terminal (blue) link rects stay 2D documentation: the wrap
        // between the crossings IS the 3D transition.
        std::vector<const RSpace*> terminalRects;
        for (const auto& s : drawn) {
            if (s.winding != ct.winding || s.parallel != ct.parallel) continue;
            if (s.isTerminal) terminalRects.push_back(&s);
        }

        if (isToroidal) {
            if (terminalRects.size() != 2) {
                throw std::runtime_error(
                    "ConductorBuilder: expected 2 toroidal terminal lead rects "
                    "(entrance, exit) for " + path.name + ", got " +
                    std::to_string(terminalRects.size()));
            }
            double od = 2.0 * wireRadius;

            // Raw crossing data straight from MAS (TurnBuilder's clearance rules, hole
            // plane at y=0). The INNER crossing is the turn's exact position and is never
            // altered.
            auto toroCrossRaw = [&](const MAS::Turn* t) -> ToroCross {
                const auto& c = t->get_coordinates();
                if (c.size() < 2) {
                    throw std::runtime_error("ConductorBuilder: turn '" + t->get_name() +
                                             "' has fewer than 2 coordinates");
                }
                auto add = t->get_additional_coordinates();
                if (!(add && !add->empty() && (*add)[0].size() >= 2)) {
                    throw std::runtime_error(
                        "ConductorBuilder: toroidal turn '" + t->get_name() +
                        "' has no additionalCoordinates (outer XY-plane crossing); "
                        "refusing to invent the outer crossing");
                }
                gp_XY pin(c[0], c[1]);
                gp_XY pout((*add)[0][0], (*add)[0][1]);
                double layerOffset =
                    std::max(0.0, (wwRadialHeight - pin.Modulus()) - wireRadius);
                return {pin, pout, halfD + layerOffset};
            };

            // Ring index (0 = the ring against the core's inner wall): the rings stack at
            // the INNER crossings (each further ring one OD toward the hole centre); ring
            // 0 has the LARGEST inner radius. The OUTER crossings can't tell rings apart —
            // MKF packs them ~0.5 OD apart at the outer wall (non-physical: on the outer
            // face successive layers must stack outward by a full OD, MKF ABT #231).
            double maxInnerRadial = 0.0, minRawOuter = std::numeric_limits<double>::max();
            for (const MAS::Turn* t : turns) {
                ToroCross rc = toroCrossRaw(t);
                maxInnerRadial = std::max(maxInnerRadial, rc.pin.Modulus());
                minRawOuter = std::min(minRawOuter, rc.pout.Modulus());
            }
            auto ringIndexOfInner = [&](const gp_XY& pin) -> int {
                return std::max(0, static_cast<int>(
                                       std::llround((maxInnerRadial - pin.Modulus()) / od)));
            };

            // Outer crossings are MKF's raw values, TRUSTED: since the ABT #231 fix MKF
            // places them at the inner crossing's azimuth (no angular search) and now
            // COMPACTS the radius continuously -- each outer-ring crossing rests tangent
            // in the V-groove between the inner ring's crossings (or at the base radius
            // where the face is free), exactly like the physical wire under tension. The
            // historical builder-side re-placement at ring0-outer + k*OD is retired: it
            // forced a full-OD stack that erased MKF's interleaving. Bad layouts are still
            // caught loudly by the collision gate and the verification battery.
            auto toroCross = [&](const MAS::Turn* t) -> ToroCross { return toroCrossRaw(t); };
            (void)minRawOuter;

            ToroCross first = toroCross(turns.front());
            ToroCross last = toroCross(turns.back());
            int maxRingIndex = 0;
            for (const MAS::Turn* t : turns)
                maxRingIndex = std::max(maxRingIndex, ringIndexOfInner(toroCrossRaw(t).pin));

            // MVB_TORO_DIAG: per-turn toroidal layout audit (ring classification, raw-vs-
            // corrected outer crossings, tube heights) -- the data needed to judge ring
            // compacting and transition routing by numbers instead of eyeballing a STEP.
            if (std::getenv("MVB_TORO_DIAG")) {
                std::cerr << "[toro] " << path.name << ": " << turns.size() << " turns, od=" << od
                          << " maxInnerRadial=" << maxInnerRadial << " minRawOuter=" << minRawOuter
                          << " maxRingIndex=" << maxRingIndex << "\n";
                for (size_t i = 0; i < turns.size(); ++i) {
                    ToroCross raw = toroCrossRaw(turns[i]);
                    ToroCross cor = toroCross(turns[i]);
                    const double azIn = std::atan2(raw.pin.Y(), raw.pin.X()) * 180.0 / kPi;
                    const double azOutRaw = std::atan2(raw.pout.Y(), raw.pout.X()) * 180.0 / kPi;
                    const double azOutCor = std::atan2(cor.pout.Y(), cor.pout.X()) * 180.0 / kPi;
                    std::cerr << "[toro]   turn " << i << " '" << turns[i]->get_name() << "'"
                              << " ring=" << ringIndexOfInner(raw.pin)
                              << " rIn=" << raw.pin.Modulus() << " azIn=" << azIn
                              << " rOutRaw=" << raw.pout.Modulus() << " azOutRaw=" << azOutRaw
                              << " rOutCor=" << cor.pout.Modulus() << " azOutCor=" << azOutCor
                              << " tube=" << cor.tube << "\n";
                }
            }

            // Depth stagger (in ODs) for wrap i's bottom return. An intra-ring wrap tucks
            // under the rings inside it (ringIndex ODs deep). A RING-TRANSITION wrap whose
            // return jumps a large azimuth across the hole (MKF winds a spread section
            // U-order: the next ring starts back at the far end of the arc, so the outer
            // crossing and the next inner crossing are up to ~a section-angle apart) is a
            // straight chord across the middle of the winding — it must run BELOW every
            // ring's returns (maxRingIndex + 1 ODs deep) so it passes under them all
            // instead of clipping them. Short-jump transitions stay at their inner ring.
            auto wrapDepthOds = [&](size_t i) -> int {
                const ToroCross c0 = toroCross(turns[i]);
                const ToroCross c1 = toroCross(turns[i + 1]);
                int r0 = ringIndexOfInner(c0.pin), r1 = ringIndexOfInner(c1.pin);
                double aOut = std::atan2(c0.pout.Y(), c0.pout.X());
                double aIn = std::atan2(c1.pin.Y(), c1.pin.X());
                double jump = std::abs(std::remainder(aIn - aOut, kTwoPi));
                // Ring-change wraps take appendToroTransitionBand (dragback in the free band
                // at level 2k along the core's central radius); this same-ring depth rule is
                // only consulted for intra-ring wraps.
                if (r0 != r1 && jump > kPi / 4.0) return std::max(0, std::max(r0, r1) - 1);
                return std::max(r0, r1);
            };

            // Deepest chord level of this conductor (including the ring stagger): the
            // axial lead route runs one wire OD beyond it (contact with its own chords,
            // cleared by the gate's contact allowance).
            double deepestChord = 0.0;
            for (size_t i = 0; i + 1 < turns.size(); ++i) {
                deepestChord = std::max(
                    deepestChord,
                    toroWrapDepth(toroCross(turns[i]), toroCross(turns[i + 1]), wireRadius) +
                        wrapDepthOds(i) * od + wireRadius);
            }
            const double leadLevel = deepestChord + 2.0 * wireRadius;

            // Mean azimuth of this conductor's inner crossings (the middle of its arc).
            double sumSin = 0, sumCos = 0;
            for (const MAS::Turn* t : turns) {
                gp_XY p = toroCross(t).pin;
                sumSin += p.Y();
                sumCos += p.X();
            }
            double meanAng = std::atan2(sumSin, sumCos);

            // The winding's over-core top and under-core bottom envelope (the top chord
            // sits at rh0 = tube + bend; the deepest return at rhb), and the outermost
            // crossing — the lead's radial run clears both and reaches past the outer
            // diameter.
            double windingTop = 0.0, windingBot = 0.0, maxOuterR = 0.0;
            for (size_t i = 0; i + 1 < turns.size(); ++i) {
                ToroCross c0 = toroCross(turns[i]), c1 = toroCross(turns[i + 1]);
                windingTop = std::max(windingTop, c0.tube + wireRadius);
                windingBot = std::max(windingBot,
                                      toroWrapDepth(c0, c1, wireRadius) + wrapDepthOds(i) * od +
                                          wireRadius);
            }
            for (const MAS::Turn* t : turns)
                maxOuterR = std::max(maxOuterR, toroCross(t).pout.Modulus());
            (void)meanAng; (void)allToroCrossings; (void)toroLeadRect;

            // Toroidal terminal lead: the standard 90-degree terminal. From the crossing
            // (where the turn starts) the wire continues AXIALLY out of the hole in the
            // poloidal direction it is wound (this builder winds inner-tube-UP at every
            // inner crossing, so the exit — continuing the last turn — leaves over the TOP
            // (+Z), and the entrance — feeding the first turn, which then goes up — comes
            // from BELOW (-Z)), then turns 90 degrees and runs RADIALLY outward, over the
            // core face, past the outer diameter. The axial leg is continuous with the
            // adjoining turn's inner tube; the radial leg clears the winding envelope.
            // ROUTED terminal leads: candidate routes are VALIDATED against the winding's own
            // wrap obstacle field before emission -- the chord lanes (turn k outer -> turn k+1
            // inner) sweep every azimuth the winding spans, including the terminal's own, so no
            // single fixed route is safe for every density:
            //   A) the classic 90-degree drop at the rim -- kept wherever it is clear (sparse
            //      windings; preserves the historical geometry);
            //   B) an inward slant to the free hole-centre corridor -- dense/full-circle
            //      windings whose chord lane crosses the terminal azimuth (measured 0.33-0.42 mm
            //      centreline clearance vs a 1.4-2.0 mm envelope with route A on the CMC and
            //      dense-toroid fixtures); physically this matches a real part, where over-wound
            //      layers press the terminal lead inward.
            // If NEITHER clears, the layout genuinely reserves no lead corridor (2-layer packed
            // holes leave no od-wide gap anywhere) -> THROW; turn positions are never moved and
            // lead space is MKF's to reserve (MKF ABT #187).
            double minCrossR = std::numeric_limits<double>::max();
            for (const MAS::Turn* t : turns)
                minCrossR = std::min(minCrossR, toroCross(t).pin.Modulus());
            // Terminal-lead obstacle field = the EMITTED wrap primitives themselves, sampled.
            // Leads are routed AFTER the wraps are emitted so every candidate is validated
            // against the real geometry, whatever its topology -- the previous hand-maintained
            // obstacle model (analytic chords/tubes) went stale the moment the ring transition
            // changed shape. Exemption is geometric: primitives touching the lead's own
            // crossing point are its continuation. Other windings on the shared core are still
            // conservative vertical capsules at their crossings (their prims are not visible
            // from this conductor's builder; the assembly collision gate remains downstream).
            auto emitToroLead = [&](const ToroCross& cross, bool isExit, size_t ordinal,
                                    const std::string& who) {
                double crossR = cross.pin.Modulus();
                if (crossR < 1e-9) return;
                double level = isExit ? (windingTop + od) : -(windingBot + od);
                double beyondR = maxOuterR + 3.0 * od;
                gp_XY dir = cross.pin;
                dir.Divide(crossR);
                gp_Pnt pCross(cross.pin.X(), 0, cross.pin.Y());
                gp_Pnt pOut(dir.X() * beyondR, level, dir.Y() * beyondR);
                // route = pOut -> elbow -> pCross (entrance) or reversed (exit); validate the
                // two legs against every emitted primitive except those that TOUCH pCross
                // (the lead is their continuation).
                struct Obst {
                    std::vector<gp_Pnt> pts;
                    double r;
                    std::string what;
                    size_t turnIdx;
                };
                std::vector<Obst> obst;
                for (const auto& pr : path.prims) {
                    auto [pa, pb] = primEndpoints(pr);
                    if (pa.Distance(pCross) < 1e-9 || pb.Distance(pCross) < 1e-9) continue;
                    obst.push_back({samplePrim(pr, wireRadius), wireRadius, pr.label,
                                    pr.turnOrdinal});
                }
                for (const auto& [pxy, r] : allToroCrossings) {
                    bool own = false;
                    for (const MAS::Turn* t : turns)
                        if ((toroCrossRaw(t).pin - pxy).Modulus() < 1e-9) { own = true; break; }
                    if (own) continue;
                    obst.push_back({{gp_Pnt(pxy.X(), -(windingBot + od), pxy.Y()),
                                     gp_Pnt(pxy.X(), windingTop + od, pxy.Y())},
                                    r, "other-winding tube", std::numeric_limits<size_t>::max()});
                }
                std::string worstWhat;
                size_t worstTurn = 0;
                int worstLeg = -1;
                // Clearance basis = the collision gate's own criterion: the exact BARE-radii
                // sum, credited with the sampled obstacles' sag bound (their chords cut inside
                // true arcs by their sagitta). There is exactly ONE definition of a collision
                // in this builder (checkCollisions).
                const double bareOwn = path.condRadius > 0 ? path.condRadius : wireRadius;
                auto routeWorst2 = [&](const gp_Pnt& elbow, const gp_Pnt& out) {
                    double worst = std::numeric_limits<double>::max();
                    const gp_Pnt* poly[3] = {&out, &elbow, &pCross};
                    for (int k = 0; k < 2; ++k)
                        for (const auto& o : obst)
                            for (size_t q = 0; q + 1 < o.pts.size(); ++q) {
                                const double d = segSegDistance(*poly[k], *poly[k + 1],
                                                                o.pts[q], o.pts[q + 1]);
                                const double bareObst = std::min(o.r, bareOwn + (o.r - wireRadius));
                                const double required = gateMinSeparation(bareOwn, bareObst) -
                                                        samplingSag(wireRadius) -
                                                        samplingSag(o.r);
                                if (d - required < worst) {
                                    worst = d - required;
                                    worstWhat = o.what;
                                    worstTurn = o.turnIdx;
                                    worstLeg = k;
                                }
                            }
                    return worst;   // >= 0 means clear under the gate's criterion
                };
                auto routeWorst = [&](const gp_Pnt& elbow) { return routeWorst2(elbow, pOut); };
                const gp_Pnt elbowA(cross.pin.X(), level, cross.pin.Y());
                // MVB_LEAD_NO_VALIDATE=1: DIAGNOSTIC ONLY -- emit the classic 90-degree drop
                // without route validation, so a layout the router refuses can still be
                // exported to STEP and inspected by eye. Never a production path: the emitted
                // copper may interpenetrate the winding (that is exactly what one goes to look
                // at). The collision gate is downstream and will still throw unless the caller
                // also inspects with the gate relaxed.
                const bool leadNoValidate = std::getenv("MVB_LEAD_NO_VALIDATE") != nullptr;
                // RECT wire keeps the classic 90-degree drop unconditionally: the round-envelope
                // obstacle model (centreline distance vs radius sums) over-rejects flat sections
                // that legitimately nest closer than their circumscribed circles (measured 0.49 mm
                // false interference on the rect-toroid fixture); the downstream collision gate
                // checks rect pairs with the correct axial/in-plane split and still guards the
                // emitted geometry.
                if (path.isRectangular || leadNoValidate) {
                    auto pushLeadSegR = [&](const gp_Pnt& a, const gp_Pnt& b, const char* wh) {
                        if (a.Distance(b) < 1e-12) return;
                        Primitive pr;
                        pr.kind = Primitive::SEG;
                        pr.seg = {a, b};
                        pr.label = who + std::string(" ") + wh;
                        pr.turnOrdinal = ordinal;
                        pr.isLead = true;
                        path.prims.push_back(std::move(pr));
                    };
                    if (isExit) {
                        pushLeadSegR(pCross, elbowA, "lead axial");
                        pushLeadSegR(elbowA, pOut, "lead radial");
                    } else {
                        pushLeadSegR(pOut, elbowA, "lead radial");
                        pushLeadSegR(elbowA, pCross, "lead axial");
                    }
                    return;
                }
                const bool toroDiag = std::getenv("MVB_TORO_DIAG") != nullptr;
                const double worstA = routeWorst(elbowA);
                if (toroDiag)
                    std::cerr << "[toro]   lead '" << who << "' routeA worst=" << worstA
                              << " culprit=" << worstWhat << " turn=" << worstTurn
                              << " leg=" << (worstLeg == 0 ? "radial" : "axial") << "\n";
                gp_Pnt elbow;
                bool routed = false;
                double bestShort = worstA;   // least-bad shortfall, for the error message
                if (worstA >= 0.0) {
                    elbow = elbowA;
                    routed = true;
                } else {
                    // Slant LADDER, depth x azimuth: clearance past rim-hugging chords/tubes
                    // varies with slant depth (32 um short at one depth clears a hair deeper),
                    // and the free gap sits BETWEEN neighbouring tubes, so the route may also
                    // rotate by fractions of the wire pitch (a real part's lead is dressed into
                    // exactly that gap; measured 9 um short at the crossing azimuth on the
                    // full-single-layer T40 that a third of a pitch clears).
                    const double rMax = std::min(0.55 * crossR, minCrossR - 1.5 * od);
                    const double azPitch = od / crossR;   // one wire od at the crossing radius
                    for (double f : {1.0, 0.75, 0.5, 0.3, 0.15}) {
                        const double rDrop = rMax * f;
                        if (rDrop <= wireRadius) break;
                        for (double g : {0.0, -1.0 / 3, 1.0 / 3, -2.0 / 3, 2.0 / 3}) {
                            const double ca = std::cos(g * azPitch), sa = std::sin(g * azPitch);
                            const gp_XY dirR(dir.X() * ca - dir.Y() * sa,
                                             dir.X() * sa + dir.Y() * ca);
                            const gp_Pnt elbowB(dirR.X() * rDrop, level, dirR.Y() * rDrop);
                            const gp_Pnt outB(dirR.X() * beyondR, level, dirR.Y() * beyondR);
                            const double worstB = routeWorst2(elbowB, outB);
                            if (toroDiag)
                                std::cerr << "[toro]     slant f=" << f << " g=" << g
                                          << " worst=" << worstB << " culprit=" << worstWhat
                                          << " turn=" << worstTurn << "\n";
                            bestShort = std::max(bestShort, worstB);
                            if (worstB >= 0.0) {
                                elbow = elbowB;
                                pOut = outB;
                                routed = true;
                                break;
                            }
                        }
                        if (routed) break;
                    }
                }
                if (!routed)
                    throw std::runtime_error(
                        "ConductorBuilder: no clear terminal-lead route for '" + who +
                        "': the straight 90-degree drop interferes with the winding by " +
                        std::to_string(-worstA) +
                        " m and every hole-centre slant depth by >= " + std::to_string(-bestShort) +
                        " m (last-checked worst obstacle: " + worstWhat + " of turn " +
                        std::to_string(worstTurn) + " vs " +
                        (worstLeg == 0 ? "radial" : "axial") +
                        " leg; bare-copper envelopes). The layout reserves no lead corridor "
                        "(dense/multi-layer toroid hole); turn positions are never moved -- "
                        "fix the winding data or the MKF blocking (MKF ABT #187).");
                auto pushLeadSeg = [&](const gp_Pnt& a, const gp_Pnt& b, const char* wh) {
                    if (a.Distance(b) < 1e-12) return;
                    Primitive pr;
                    pr.kind = Primitive::SEG;
                    pr.seg = {a, b};
                    pr.label = who + std::string(" ") + wh;
                    pr.turnOrdinal = ordinal;
                    pr.isLead = true;
                    path.prims.push_back(std::move(pr));
                };
                // ROUND corner at the elbow, like every wrap corner: shorten both legs by the
                // fillet's tangent length and bend through an exact arc (works for the slanted
                // drop candidates too -- the fillet is generic in the two leg directions). The
                // old bare two-seg corner left a 90-degree junction for the mitre machinery,
                // whose trim on the short lead stubs sometimes fell back to flush overlapping
                // tubes instead of a clean joint.
                auto pushLeadRun = [&](const gp_Pnt& start, const gp_Pnt& mid, const gp_Pnt& end,
                                       const char* whA, const char* whB) {
                    gp_Vec inVec(start, mid), outVec(mid, end);
                    const double lenIn = inVec.Magnitude(), lenOut = outVec.Magnitude();
                    if (lenIn < 1e-12 || lenOut < 1e-12) {
                        pushLeadSeg(start, mid, whA);
                        pushLeadSeg(mid, end, whB);
                        return;
                    }
                    const gp_XYZ u = inVec.XYZ() / lenIn, v = outVec.XYZ() / lenOut;
                    const double cosTurn = u.Dot(v);
                    // > wireRadius, never == : an equal-radius corner is a horn torus (the tube
                    // touches its own revolution axis) and OCC rejects the solid.
                    const double b = 1.5 * wireRadius;
                    const double tangentLen =
                        b * std::sqrt(std::max(0.0, (1.0 - cosTurn) / (1.0 + cosTurn)));
                    if (cosTurn > 1.0 - 1e-9 || tangentLen > 0.49 * std::min(lenIn, lenOut)) {
                        pushLeadSeg(start, mid, whA);   // straight or no room: plain segs
                        pushLeadSeg(mid, end, whB);
                        return;
                    }
                    const gp_Pnt tangentA(mid.XYZ() - u * tangentLen);
                    const gp_Pnt tangentB(mid.XYZ() + v * tangentLen);
                    gp_XYZ normal = v - u * cosTurn;
                    normal /= normal.Modulus();
                    const gp_Pnt center(tangentA.XYZ() + normal * b);
                    pushLeadSeg(start, tangentA, whA);
                    {
                        Primitive pr;
                        pr.kind = Primitive::ARC3;
                        pr.arc.c = center;
                        gp_XYZ axis = u.Crossed(v);
                        pr.arc.axis = axis / axis.Modulus();
                        pr.arc.v0 = tangentA.XYZ() - center.XYZ();
                        pr.arc.sweep = std::acos(std::clamp(cosTurn, -1.0, 1.0));
                        pr.label = who + std::string(" lead corner");
                        pr.turnOrdinal = ordinal;
                        pr.isLead = true;
                        path.prims.push_back(std::move(pr));
                    }
                    pushLeadSeg(tangentB, end, whB);
                };
                if (isExit) {   // crossing -> up out of the hole -> radial out
                    pushLeadRun(pCross, elbow, pOut, "lead axial", "lead radial");
                } else {        // radial in -> down into the hole -> crossing (feeds turn 0)
                    pushLeadRun(pOut, elbow, pCross, "lead radial", "lead axial");
                }
            };

            // WRAPS FIRST, leads after: the lead router validates candidates against the
            // emitted primitives, so the winding must exist before its terminals are dressed.
            // A ring change is the Z-connection transition (dragback over the top); everything
            // else is a normal poloidal wrap.
            for (size_t i = 0; i + 1 < turns.size(); ++i) {
                const ToroCross c0 = toroCross(turns[i]), c1 = toroCross(turns[i + 1]);
                const int r0 = ringIndexOfInner(toroCrossRaw(turns[i]).pin);
                const int r1 = ringIndexOfInner(toroCrossRaw(turns[i + 1]).pin);
                const std::string wlabel = "wrap '" + turns[i]->get_name() + "' -> '" +
                                           turns[i + 1]->get_name() + "'";
                if (std::getenv("MVB_TORO_DIAG"))
                    std::cerr << "[toro]   wrap " << i << "->" << i + 1
                              << (r0 != r1 ? " TRANSITION" : "") << " topH=" << c0.tube + wireRadius
                              << " azOut=" << std::atan2(c0.pout.Y(), c0.pout.X()) * 180.0 / kPi
                              << " azNextIn=" << std::atan2(c1.pin.Y(), c1.pin.X()) * 180.0 / kPi
                              << "\n";
                if (r0 != r1)
                    appendToroTransitionBand(path, c0, c1, wireRadius,
                                             std::max(0, std::max(r0, r1) - 1) * od,
                                             0.5 * (maxInnerRadial + minRawOuter), wlabel, i);
                else
                    appendToroWrap(path, c0, c1, wireRadius, wrapDepthOds(i) * od, wlabel, i);
            }

            const size_t wrapPrimCount = path.prims.size();
            emitToroLead(first, /*isExit=*/false, 0, path.name + " entrance");
            // The entrance is the chain's HEAD: rotate its prims (appended above) to the front
            // so primitive order stays electrically continuous end-to-end.
            std::rotate(path.prims.begin(), path.prims.begin() + wrapPrimCount, path.prims.end());

            emitToroLead(last, /*isExit=*/true, turns.size() - 1, path.name + " exit");

            paths.push_back(std::move(path));
            continue;
        }

        // ---- concentric columns -------------------------------------------------------
        PlanePt first = station(turns.front());
        PlanePt last = station(turns.back());

        // Terminal grouping follows MKF's emission order: per (winding, parallel) the
        // ENTRANCE lead's rects are pushed first, then the EXIT lead's, each as
        // [vertical stub?, horizontal run] — the run closes its group. Proximity-based
        // grouping is ambiguous when a crossing sits within a wire of the window edge.
        std::vector<const RSpace*> entranceGroup, exitGroup;
        std::tie(entranceGroup, exitGroup) = splitTerminalGroups(terminalRects, path.name);

        // czRaise: how far the turn this lead attaches to has been lifted by the dragbacks it
        // rides over. The lead must meet the wire WHERE IT ACTUALLY IS, not at the nominal
        // station, or it is laid straight through the bumped turn.
        // How far the terminal leads stick OUT past the winding, in the connection
        // plane's radial (-Z) direction — the winding-window opening side. Like the
        // toroidal lead ("radially outward, past the outer diameter and beyond"), the
        // concentric lead runs out past the outermost turn so the flat terminal face sits
        // clear of the coil for FEM. maxTurnRadius is the outermost crossing.
        double maxTurnRadius = 0.0;
        for (const MAS::Turn* t : turns)
            maxTurnRadius = std::max(maxTurnRadius, station(t).x);
        // The tip plane sits beyond the winding AND beyond the tallest vertical-fan column, so
        // a lifted lead still runs OUTWARD all the way. Entrance and exit share the plane, so
        // a conductor's two terminals finish on the same face.
        // RECT columns: the terminal attach can sit on a dragback-DISPLACED face (rectRideFor),
        // so the common tip plane must clear the displaced attach too, or the tip lands INSIDE
        // the attach (17_cllc: exit attach ridden to 24.9 mm vs a 24.2 mm tip plane).
        double maxRide = 0.0;
        if (columnShape == MAS::ColumnShape::RECTANGULAR)
            for (int side2 = 0; side2 < 2; ++side2)
                for (const auto& [lvlZ, diam] : rectRideLevels[side2])
                    maxRide += diam;   // conservative: all reservation levels stacked
        const double leadTipRadius =
            maxTurnRadius + 4.0 * (2.0 * wireRadius) + fanMaxRaise + maxRide;

        auto pushPlaneSegs = [&](std::vector<PlanePt> wp, const std::string& what,
                                 size_t ordinal, bool stationAtFront, double liftRaise = 0.0,
                                 double azLead = kPlaneAz) {
            // Absorb intermediate waypoints closer than the wire radius to their
            // neighbour: a jog shorter than the wire's own radius lies entirely inside
            // the pipe body of the adjacent edge (and inside MKF's drawn rectangle,
            // whose height is one full wire OD), while sub-radius spine edges crash
            // OCCT's pipe-shell corner rounding. Endpoints (the exact stations/border)
            // are always kept.
            std::vector<PlanePt> kept;
            kept.push_back(wp.front());
            // Absorption bound: the RADIUS inside which a jog truly disappears into the wire
            // body. For round wire that is wireRadius; for RECT/PLANAR wire wireRadius is the
            // half-DIAGONAL (1.09 mm on 09_planar's flat trace), and absorbing an L-corner
            // within it turned the route into a diagonal that crossed the other lead's edge
            // row (measured: entrance seg crossing the exit row, centreline distance 0). The
            // conservative correct bound is the THIN half-dimension.
            const double absorbTol = rectWire
                ? 0.5 * std::min(path.wireWidth, path.wireHeight) : wireRadius;
            if (std::getenv("MVB_DIAG")) {
                std::cerr << "[lead-wp] " << path.name << " " << what << " absorbTol="
                          << absorbTol << " wr=" << wireRadius << " wp:";
                for (const auto& q : wp) std::cerr << " (" << q.x << "," << q.y << ")";
                std::cerr << "\n";
            }
            for (size_t i = 1; i + 1 < wp.size(); ++i) {
                if (std::hypot(wp[i].x - kept.back().x, wp[i].y - kept.back().y) <
                        absorbTol ||
                    std::hypot(wp[i].x - wp.back().x, wp[i].y - wp.back().y) < absorbTol) {
                    continue;
                }
                kept.push_back(wp[i]);
            }
            kept.push_back(wp.back());
            // ROUND wire gets rounded corners (the toroidal treatment); RECT wire keeps the
            // mitred butt its own machinery expects.
            std::vector<gp_Pnt> leadPts;
            leadPts.reserve(kept.size() + 1);
            for (const auto& q : kept) leadPts.push_back(azPointC(0, 0, q.x + zoff, q.y, azLead));
            // The lead runs STRAIGHT OUT at the level the turn finishes (Alf): every point,
            // tip included, carries the dropped-centre lift the adjoining wrap ends with
            // (cz = -raise at the SAME radius, from tallestBumpColumn on the same bump
            // list), so the attachment matches the wrap exactly and the run never dips back
            // to the bare radius -- the un-dip read as a zigzag at the exits. The TIP's
            // radius is chosen so it lands exactly on the conductor's common tip PLANE
            // (z = -leadTipRadius), whatever the lift or fan slot: entrance and exit still
            // finish on the same outermost face.
            if (!rectWire && leadPts.size() > 1) {
                const size_t tip = stationAtFront ? leadPts.size() - 1 : 0;
                const size_t nb = stationAtFront ? tip - 1 : 1;
                for (size_t i = 0; i < leadPts.size(); ++i) {
                    if (i == tip) continue;
                    if (liftRaise > 0.0)
                        leadPts[i] = azPointC(0, -liftRaise, kept[i].x + zoff,
                                              kept[i].y, azLead);
                }
                // The terminal RUN is a totally straight wire PARALLEL TO THE Z AXIS
                // (Alf), whatever fan slot its stub occupies: it keeps the last
                // waypoint's world X and lands on the conductor-common tip plane. The
                // old radial-at-azimuth run skewed by the fan angle (visible ~5 deg on
                // 11_pushpull's Secondary 1 parallel 0).
                leadPts[tip] = gp_Pnt(leadPts[nb].X(), kept[tip].y,
                                      -(leadTipRadius + zoff));
            }
            // The lead lies in its fan slot's axial plane and runs straight out radially, like
            // the dragback.
            if (!rectWire) {
                appendFilletedPolyline(path.prims, leadPts, wireRadius, what, ordinal,
                                       /*isLead=*/true, /*isConnection=*/false);
            }
            else {
                for (size_t i = 0; i + 1 < leadPts.size(); ++i) {
                    if (leadPts[i].Distance(leadPts[i + 1]) < 1e-12) continue;
                    Primitive pr;
                    pr.kind = Primitive::SEG;
                    pr.seg = {leadPts[i], leadPts[i + 1]};
                    pr.label = what + " seg " + std::to_string(i);
                    pr.turnOrdinal = ordinal;
                    pr.isLead = true;
                    path.prims.push_back(std::move(pr));
                }
            }
        };

        auto extendBorder = [&](std::vector<PlanePt>& wp) {
            // The BORDER is always the route's LAST waypoint (terminalWaypoints order:
            // station first; called before the entrance reversal). Push IT to the common
            // tip plane. Picking "the farthest-out waypoint" instead extended the STATION
            // itself whenever a dragback-displaced attach sat beyond the drawn border
            // (17_cllc secondary exit: attach 19.17 vs border 19.05), detaching the lead
            // from its turn -- both ends landed at the border.
            wp.back().x = std::max(wp.back().x, leadTipRadius);
        };

        // The conductor's typical within-layer axial advance per wrap: median |dy| over its
        // same-radius consecutive crossings. Calibrates the Z-transition detection so sparse
        // windings (pitch of several wire-ODs) keep their ordinary helix wraps.
        double medianPitch = 0.0, zOrderAdvance = 0.0;
        {
            std::vector<double> pitches;
            int pos = 0, neg = 0;
            for (size_t i = 0; i + 1 < turns.size(); ++i) {
                PlanePt a = station(turns[i]), b = station(turns[i + 1]);
                if (std::abs(b.x - a.x) <= wireRadius) {
                    const double dy = b.y - a.y;
                    pitches.push_back(std::abs(dy));
                    (dy >= 0.0 ? pos : neg)++;
                }
            }
            if (!pitches.empty()) {
                std::sort(pitches.begin(), pitches.end());
                medianPitch = pitches[pitches.size() / 2];
                if (pos == 0 || neg == 0) zOrderAdvance = (pos > 0 ? medianPitch : -medianPitch);
            }
        }
        // NO TRAILING PRUNE. Every station MKF describes is visited by the copper. A prune of
        // trailing return/link transitions used to live here (my over-application of Alf's
        // 23_interleaved note about zigzagging LEAD pieces, which the lead work itself fixed:
        // whole-lead lift, straight parallel-Z runs, MKF-drawn routes). It silently dropped
        // REAL TURNS -- the litz fixture lost a whole 4-turn outer layer, 23's secondary lost
        // turns 2-3 -- so the emitted body was a DIFFERENT MAGNETIC from the one MKF drew
        // (Alf, 2026-08-07: "litz_fixture.svg and COLLIDING_litz_fixture.step are not the same
        // design"). Worse, the terminal then attached to the pruned-back turn while still
        // consuming the run row MKF computed for the true last turn, which drove the exit lead
        // straight down through a dozen of its own turns.
        const size_t nEmit = turns.size();
        // Wraps between consecutive crossings. Z-returns are only RECORDED here (the connecting
        // end-run is planned after every conductor is built, against the full obstacle field —
        // see planZEndRuns); ordinary wraps and serpentine links are appended directly.
        // Z DRAGBACKS: one azimuth per transition, from the vertical-fan pre-scan -- assigned
        // up front so the wraps, emitted in the same pass, already know which dragbacks they
        // must ride over. The pre-scan mirrors this detection exactly; a transition it did not
        // see is a code drift, not a layout condition, so it throws.
        std::map<size_t, double> zDragbackAzimuth;
        std::vector<PlanePt> stationsE;
        for (const MAS::Turn* t : turns) stationsE.push_back(station(t));
        const auto bandsE = computePitchBands(stationsE, wireRadius);
        auto zReturnAt = [&](size_t i) {
            const PlanePt a = station(turns[i]), b = station(turns[i + 1]);
            const PitchBand band = bandAt(bandsE, a.x, wireRadius, medianPitch, zOrderAdvance);
            return isZReturn(a, b, wireRadius, band.medianPitch, band.advance);
        };
        if (effectivelyRound) {
            for (size_t i = 0; i + 1 < nEmit; ++i) {
                if (!zReturnAt(i)) continue;
                auto it = dragAzOf.find({ci, i});
                if (it == dragAzOf.end()) {
                    throw std::runtime_error(
                        "ConductorBuilder: dragback transition " + std::to_string(i) + " of " +
                        path.name + " has no azimuth in the vertical-fan pre-scan (detection "
                        "drift between the pre-scan and the emission loop)");
                }
                zDragbackAzimuth[i] = it->second;
            }
        }
        // Dragbacks (azimuth + the radius their axial run occupies) come from the WINDOW-WIDE
        // pre-scan above, so a turn rides over every return in the window at or inside its own
        // radius -- not just the ones its own conductor laid.
        auto bumpsForTurn = [&](double turnRadius) {
            std::vector<WrapBump> out;
            for (const auto& ld : laidDragbacks) {
                // EVERY layer at or outside a return rides over it. Returns SHARING an azimuth
                // stack (their raises add -- see tallestBumpColumn); returns in different fan
                // slots sit side by side. A return lies in its destination layer's space, so
                // that layer steps over it; the step puts the wire where the layer outside would
                // be, so that one steps over too, and so on outward.
                if (turnRadius < ld.radius - wireRadius) continue;
                // Raise by the dragback wire's own COATED diameter: that is what has to be
                // cleared, and different windings use different wire.
                out.push_back({ld.az, ld.diam});
            }
            return out;
        };
        // Crossing azimuths: kPlaneAz normally; the fan slot where the crossing belongs to a
        // dragback (both its ends) or to a terminal lead (the first/last crossing).
        std::vector<double> crossAz(nEmit, kPlaneAz);
        double azEntrance = kPlaneAz, azExit = kPlaneAz;
        if (effectivelyRound && !turns.empty()) {
            if (auto it = leadAzIn.find(ci); it != leadAzIn.end()) azEntrance = it->second;
            if (auto it = leadAzOut.find(ci); it != leadAzOut.end()) azExit = it->second;
            crossAz.front() = azEntrance;
            crossAz.back() = azExit;
            for (const auto& [m, azv] : zDragbackAzimuth) {
                crossAz[m] = azv;
                crossAz[m + 1] = azv;
            }
            // A LAYER LINK never changes azimuth: it hands its own azimuth FORWARD to the
            // next crossing, so the radial step and whatever follows (usually the exit lead)
            // stay collinear -- one straight run out to the tip plane, instead of stepping
            // out, walking an arc along the ring to a different slot, and going out again
            // (the zigzag Alf flagged on 23's secondary exits). A link feeding straight into
            // a dragback keeps that dragback's own slot (the link's azimuthal run covers the
            // difference there).
            for (size_t i = 0; i + 1 < nEmit; ++i) {
                const PlanePt a = station(turns[i]), b = station(turns[i + 1]);
                const bool isLink = std::abs(b.x - a.x) > wireRadius &&
                                    std::abs(b.y - a.y) <= std::abs(b.x - a.x);
                if (isLink && !zDragbackAzimuth.count(i + 1)) crossAz[i + 1] = crossAz[i];
            }
            // A dragback or link owns its crossings outright -- the lead then joins AT that
            // azimuth rather than pulling the shared crossing off it.
            if (turns.size() > 1) {
                azEntrance = crossAz.front();
                azExit = crossAz.back();
            }
        }

        // Entrance: MKF's drawn route, walked from the border TO the first station, then
        // extended straight out so the terminal sticks clear of the coil. Emitted at the
        // entrance's fan azimuth (the first crossing's), so the first wrap starts exactly
        // where the lead delivers the wire -- INCLUDING the wrap's raise: an interleaved
        // winding's first turn sits outside another winding's returns, so its first wrap
        // starts ridden-up and the lead must attach at that lifted radius (23_interleaved:
        // the bare secondary entrance ran level with the primary's raised head, 0.31 mm off).
        // RECT wire on a round column: a real CORNER connects the turn to its lead (Alf:
        // 03's ribbon butted the lead at 90 degrees with a bare notch). Quarter arc about
        // the vertical axis at the crossing, tangent to the wrap's azimuthal direction AND
        // to the radial run; the run continues parallel, offset one bend radius
        // tangentially. The rect single-body sweep rounds the section through it natively.
        // slope: the adjoining wrap's dy per unit arc length at the crossing. The corner and
        // the straight run BOTH carry it, so every junction is EXACTLY tangent -- a flat
        // corner met the rising wrap at a ~2 deg kink, and the spine filleter's 1.25*r trim
        // then ate the whole quarter arc into a degenerate edge (BRepAdaptor: No geometry).
        auto rectLeadCorner = [&](const PlanePt& st, double azL, bool isExit, double slope,
                                  const std::string& what) {
            const gp_XYZ Rhat(std::cos(azL), 0.0, -std::sin(azL));
            const gp_XYZ That(-std::sin(azL), 0.0, -std::cos(azL));
            const double Rc = std::max(minBend, 1.02 * 0.5 * path.wireWidth);
            const double q = 0.5 * kPi * Rc;
            const gp_XYZ C2 = azPointC(0, 0, st.x, st.y, azL).XYZ() + Rhat * Rc;
            const double runLen = leadTipRadius - st.x;
            // y along the wire: the crossing keeps the station row; corner and run continue
            // at the wrap's slope (entrance: the wire RISES into the wrap; exit: it keeps
            // climbing away).
            const double yP = st.y;
            const double yQ = isExit ? yP + slope * q : yP - slope * q;
            const double yTip = isExit ? yQ + slope * runLen : yQ - slope * runLen;
            Primitive arc;
            arc.kind = Primitive::SPIRAL;
            if (isExit) {
                // P (az azL+pi about C) -> Q (az azL+pi/2), azimuth decreasing.
                arc.spiral = {C2.X(), C2.Z(), Rc, yP, azL + kPi, Rc, yQ, azL + kPi / 2.0};
            } else {
                // Q (az azL-pi/2) -> P (az azL-pi), azimuth decreasing.
                arc.spiral = {C2.X(), C2.Z(), Rc, yQ, azL - kPi / 2.0, Rc, yP, azL - kPi};
            }
            arc.label = path.name + " " + what + " corner";
            arc.turnOrdinal = isExit ? turns.size() - 1 : 0;
            arc.isLead = true;
            const gp_Pnt Q(isExit ? C2 + That * Rc + gp_XYZ(0, yQ - st.y, 0)
                                  : C2 - That * Rc + gp_XYZ(0, yQ - st.y, 0));
            const gp_Pnt tipPt(Q.XYZ() + Rhat * runLen + gp_XYZ(0, yTip - yQ, 0));
            Primitive run;
            run.kind = Primitive::SEG;
            run.seg = isExit ? Seg{Q, tipPt} : Seg{tipPt, Q};
            run.label = path.name + " " + what + " run";
            run.turnOrdinal = arc.turnOrdinal;
            run.isLead = true;
            if (isExit) {
                path.prims.push_back(std::move(arc));
                path.prims.push_back(std::move(run));
            } else {
                path.prims.push_back(std::move(run));
                path.prims.push_back(std::move(arc));
            }
        };
        {
            std::vector<PlanePt> wp;
            if (effectivelyRound && rectWire) {
                // The tangent corner replaces the straight radial attach; it is only
                // defined for MKF routes WITHOUT a vertical connection. If MKF draws an
                // L here, the corner-through-a-stub geometry is unspecified -- refuse.
                if (terminalWaypoints(entranceGroup, first, path.name + " entrance").size() != 2)
                    throw std::runtime_error(
                        "ConductorBuilder: MKF drew a vertical connection on " + path.name +
                        "'s rect-wire entrance lead -- the tangent lead corner through an "
                        "L-route is not implemented");
                const double sIn =
                    turns.size() > 1
                        ? (station(turns[1]).y - station(turns[0]).y) /
                              (kTwoPi * std::max(station(turns[0]).x, 1e-9))
                        : 0.0;
                rectLeadCorner(first, azEntrance, /*isExit=*/false, sIn, "entrance lead");
            } else {
                // The MKF-DRAWN entrance route, for EVERY column shape (Alf, 2026-08-07:
                // the connection positions are MKF/MAS data -- the drawn vertical
                // connection and run row exactly as the Painter SVG shows them; azimuth
                // is the only dimension the 3D fan adds).
                PlanePt fLead = first;
                // Rect columns: the crossing sits on the DISPLACED -Z face (dragback
                // reservation), so the lead attaches there too.
                if (columnShape == MAS::ColumnShape::RECTANGULAR) {
                    const RectStation rsF =
                        rectStation(first, halfW, halfD, minBend, path.name);
                    const double rideF = rectRideFor(rsF.zPos, windingFace.at(ct.winding));
                    fLead.x += rideF;
                // RECT-WIRE LEAD CORNER (Alf, 18_stacked: "the connections are not
                // properly connected to the end of turn"): the ribbon cannot butt the
                // terminal run at 90 degrees -- the face straight ends a corner radius
                // short (see the appendRectWrap call) and this in-plane quarter wedge
                // turns it into the lead, conformal to both end planes.
                auto rectLeadCornerPrim = [&](const RectStation& rsA, double rideA,
                                              double rowY, bool isExit, size_t ord) {
                    const double cR = rsA.cornerR;
                    const double zA = rsA.zPos + rideA;
                    Primitive cnr;
                    cnr.kind = Primitive::ARC3;
                    cnr.arc.axis = gp_XYZ(0, -1, 0);
                    cnr.arc.sweep = kPi / 2.0;
                    if (isExit) {
                        cnr.arc.c = gp_Pnt(cR, rowY, -(zA + cR));
                        cnr.arc.v0 = gp_XYZ(0, 0, cR);
                    } else {
                        cnr.arc.c = gp_Pnt(-cR, rowY, -(zA + cR));
                        cnr.arc.v0 = gp_XYZ(cR, 0, 0);
                    }
                    cnr.label = path.name + (isExit ? " exit" : " entrance") + " lead corner";
                    cnr.turnOrdinal = ord;
                    cnr.isLead = true;
                    path.prims.push_back(std::move(cnr));
                };
                    // Corner only when the first transition is a rising turn that was
                    // actually shortened for it (mirrors the appendRectWrap call).
                    bool ret0 = false;
                    for (const auto& r : rectReturns)
                        if (r.ci == ci && r.trans == 0) ret0 = true;
                    if (rectWire && turns.size() > 1 && !ret0) {
                        rectLeadCornerPrim(rsF, rideF, first.y, /*isExit=*/false, 0);
                        fLead.x += rsF.cornerR;
                    }
                }
                wp = terminalWaypoints(entranceGroup, fLead, path.name + " entrance");
            }
            if (!wp.empty()) {
                extendBorder(wp);
                std::reverse(wp.begin(), wp.end());
                const double entrRaise =
                    tallestBumpColumn(bumpsForTurn(station(turns.front()).x)).first;
                pushPlaneSegs(wp, "entrance lead", 0, /*stationAtFront=*/false, entrRaise,
                              azEntrance);
            }
        }
        const bool dragDiag = std::getenv("MVB_DRAG_DIAG") != nullptr;
        for (size_t i = 0; i + 1 < nEmit; ++i) {
            PlanePt s = station(turns[i]);
            PlanePt nxt = station(turns[i + 1]);
            std::string label = "wrap '" + turns[i]->get_name() + "' -> '" +
                                turns[i + 1]->get_name() + "'";
            if (effectivelyRound && dragDiag) {
                double raise = 0.0;
                for (const auto& bmp : bumpsForTurn(s.x)) raise += bmp.distance;
                std::cerr << "[drag] turn " << i << " r=" << s.x << " y=" << s.y
                          << (zDragbackAzimuth.count(i) ? " DRAGBACK" : "")
                          << " bumpLength=" << raise << " laid=" << laidDragbacks.size() << "\n";
            }
            if (effectivelyRound) {
                auto zit = zDragbackAzimuth.find(i);
                if (zit != zDragbackAzimuth.end()) {
                    // NB: no push_back here -- every return in the window, this one included, is
                    // already in laidDragbacks from the pre-scan.
                    appendZDragback(path, s, nxt, wireRadius, zit->second, label, i,
                                    bumpsForTurn(s.x), bumpsForTurn(nxt.x));
                    continue;
                }
                const bool weldable = path.femReady && turns.size() > 1 &&
                                      (!rectWire || std::getenv("MVB_RECT_OVERSHOOT"));
                // A wrap overshoots wherever the chain BREAKS at its end: terminal leads
                // (first/last transition), serpentine U-links, and Z end-runs all butt the
                // wrap in a flat cap at the station -- every such junction needs the lens.
                appendRoundWrap(path, s, nxt, wireRadius, label, i, bumpsForTurn(s.x),
                                crossAz[i], crossAz[i + 1], bumpsForTurn(nxt.x));
            } else if (columnShape == MAS::ColumnShape::RECTANGULAR) {
                {
                    const RectStation rs0 = rectStation(s, halfW, halfD, minBend, path.name);
                    const RectStation rs1 = rectStation(nxt, halfW, halfD, minBend, path.name);
                    const int side = windingFace.at(ct.winding);
                    const RectReturn* ret = nullptr;
                    const RectReturn* nextRet = nullptr;
                    for (const auto& r : rectReturns) {
                        if (r.ci == ci && r.trans == i) ret = &r;
                        if (r.ci == ci && r.trans == i + 1) nextRet = &r;
                    }
                    double chainRide = 0.0, destRide = 0.0, xSlot = 0.0;
                    if (ret) {
                        destRide = rectRideFor(rs1.zPos, side);
                        // The descent lies on the destination face displaced only by the
                        // levels INSIDE it -- its own level's reservation is the space the
                        // descent itself occupies.
                        chainRide = destRide - ret->diam;
                        xSlot = ret->xSlot;
                        if (chainRide < -1e-12) {
                            throw std::runtime_error(
                                "ConductorBuilder: dragback level accounting negative for " +
                                label);
                        }
                    }
                    // RECT-WIRE LEAD CORNERS (Alf, 18_stacked): the first/last face
                    // straight ends a corner radius short of the crossing; the in-plane
                    // quarter wedge (emitted with the lead) turns into the terminal run.
                    const double kNaN = std::numeric_limits<double>::quiet_NaN();
                    const bool rw = path.isRectangular;
                    const double stopX = nextRet ? nextRet->xSlot
                                       : (rw && i + 2 == nEmit && !ret ? rs1.cornerR : kNaN);
                    const double startX = (rw && i == 0 && !ret) ? rs0.cornerR : kNaN;
                    appendRectWrap(path, rs0, rs1, label, i, wireRadius,
                                   rectRideFor(rs0.zPos, side),
                                   rectRideFor(rs0.zPos, 1 - side), ret != nullptr,
                                   chainRide, destRide, xSlot, stopX, startX);
                }
            } else {   // OBLONG with a real straight section
                appendOblongWrap(path, s, nxt, oblongHalf, label, i);
            }
        }

        // Exit: MKF's drawn route from the last station out to the border, extended out.
        // With no drawn exit group (MKF emitted one lead; see above), the synthetic route is the
        // minimal one: from the last station straight out radially at its own level.
        {
            std::vector<PlanePt> wp;
            if (effectivelyRound && rectWire) {
                if (!exitGroup.empty() &&
                    terminalWaypoints(exitGroup, last, path.name + " exit").size() != 2)
                    throw std::runtime_error(
                        "ConductorBuilder: MKF drew a vertical connection on " + path.name +
                        "'s rect-wire exit lead -- the tangent lead corner through an "
                        "L-route is not implemented");
                const double sOut =
                    nEmit > 1
                        ? (station(turns[nEmit - 1]).y - station(turns[nEmit - 2]).y) /
                              (kTwoPi * std::max(station(turns[nEmit - 1]).x, 1e-9))
                        : 0.0;
                rectLeadCorner(last, azExit, /*isExit=*/true, sOut, "exit lead");
            } else if (!exitGroup.empty()) {
                // The MKF-DRAWN exit route, every column shape (Alf, 2026-08-07 -- see
                // the entrance).
                PlanePt lLead = last;
                // Displaced -Z crossing, exactly like the entrance.
                if (columnShape == MAS::ColumnShape::RECTANGULAR) {
                    const RectStation rsL =
                        rectStation(last, halfW, halfD, minBend, path.name);
                    const double rideL = rectRideFor(rsL.zPos, windingFace.at(ct.winding));
                    lLead.x += rideL;
                // RECT-WIRE LEAD CORNER (Alf, 18_stacked: "the connections are not
                // properly connected to the end of turn"): the ribbon cannot butt the
                // terminal run at 90 degrees -- the face straight ends a corner radius
                // short (see the appendRectWrap call) and this in-plane quarter wedge
                // turns it into the lead, conformal to both end planes.
                auto rectLeadCornerPrim = [&](const RectStation& rsA, double rideA,
                                              double rowY, bool isExit, size_t ord) {
                    const double cR = rsA.cornerR;
                    const double zA = rsA.zPos + rideA;
                    Primitive cnr;
                    cnr.kind = Primitive::ARC3;
                    cnr.arc.axis = gp_XYZ(0, -1, 0);
                    cnr.arc.sweep = kPi / 2.0;
                    if (isExit) {
                        cnr.arc.c = gp_Pnt(cR, rowY, -(zA + cR));
                        cnr.arc.v0 = gp_XYZ(0, 0, cR);
                    } else {
                        cnr.arc.c = gp_Pnt(-cR, rowY, -(zA + cR));
                        cnr.arc.v0 = gp_XYZ(cR, 0, 0);
                    }
                    cnr.label = path.name + (isExit ? " exit" : " entrance") + " lead corner";
                    cnr.turnOrdinal = ord;
                    cnr.isLead = true;
                    path.prims.push_back(std::move(cnr));
                };
                    bool retLast = false;
                    for (const auto& r : rectReturns)
                        if (r.ci == ci && r.trans + 2 == nEmit) retLast = true;
                    if (rectWire && nEmit > 1 && !retLast) {
                        rectLeadCornerPrim(rsL, rideL, last.y, /*isExit=*/true, nEmit - 1);
                        lLead.x += rsL.cornerR;
                    }
                }
                wp = terminalWaypoints(exitGroup, lLead, path.name + " exit");
            } else {
                // MKF drew only one lead (see splitTerminalGroups): synthesized minimal
                // straight-out exit at the last turn's own row.
                wp.push_back(last);
                wp.push_back({leadTipRadius, last.y});
            }
            if (!wp.empty()) {
                extendBorder(wp);
                // The exit lead leaves the OUTERMOST turn, which rides over every dragback
                // beneath it, so it is z-lifted by the fan's tallest column -- computed by
                // the SAME helper on the SAME bump list as the last wrap's end raise, so the
                // two meet at one identical point.
                const double exitRaise =
                    tallestBumpColumn(bumpsForTurn(station(turns[nEmit - 1]).x)).first;
                pushPlaneSegs(wp, "exit lead", nEmit - 1, /*stationAtFront=*/true, exitRaise,
                              azExit);
            }
        }

        // Stagger this winding's seam azimuth so its leads clear the other windings' leads (all drawn
        // on the shared reference plane, ABT #240). A rigid rotation about the column axis leaves every
        // ring where it was, so it is valid wherever the rotation maps the COLUMN onto itself:
        //   - ROUND column: any angle -> spread the windings continuously over ~2/3 of the circle;
        //   - RECT / OBLONG column: 180deg-symmetric -> alternate windings by 180deg (opposite faces),
        //     which separates a 2-winding transformer (the common case) and never misaligns the core.
        // On a ROUND column, PARALLEL strands additionally get a SMALL per-parallel rotation:
        // MKF's 2D data fixes only (r, y) — the 3D azimuth is a modelling choice, and the rings
        // are circles (rotation-invariant), so rotating a strand's whole path separates its
        // leads/links/end-runs TANGENTIALLY from its sibling strands (physically: the strands
        // exit side by side) without touching any turn position.
        // The winding/parallel stagger is applied AFTER the core probe: how far the leads may be
        // spread depends on how wide the core's opening actually is (see below).
        double seamAngle = effectivelyRound ? 0.0 : (windingFace.at(ct.winding) * kPi);
        // Aim the whole seam sector at the core's window OPENING (Options::leadExitAzimuth):
        // a rigid rotation about the column axis, exactly like the stagger. ROUND columns
        // rotate freely; RECT/OBLONG columns only map onto themselves under 180 degrees, so
        // snap to the nearer of {0, pi} there.
        // ROUND columns only: a rect/oblong racetrack maps onto itself solely under 180
        // degrees, and flipping the seam side broke the EP-stadium and E-zigzag fixtures
        // while no validated design needs it -- E-family cores already open at -Z.
        if (std::getenv("MVB_DIAG"))
            std::cerr << "[lead-aim] gate: effRound=" << effectivelyRound
                      << " obstacles=" << opts.coreObstacles.size() << "\n";
        double freeArc = kTwoPi;   // angular width of the core's opening; kTwoPi = unconstrained
        if (effectivelyRound && (!opts.coreObstacles.empty() || !std::isnan(opts.leadExitAzimuth))) {
            double exitAz = opts.leadExitAzimuth;
            if (!opts.coreObstacles.empty()) {
                // Classify the CORE around the lead-tip circle: the widest arc the core does
                // not occupy is the real window opening. Sample a small axial band (the leads
                // run near the winding's ends and its middle).
                double maxTurnR = 0.0, maxAbsY = 0.0;
                for (const MAS::Turn* t : turns) {
                    maxTurnR = std::max(maxTurnR, station(t).x);
                    maxAbsY = std::max(maxAbsY, std::abs(station(t).y));
                }
                // Probe the WHOLE radial run of the lead, not just its tip: on a PQ the tip
                // radius clears the core entirely (fully free) while the run itself passes
                // straight through a plate at intermediate radii.
                // Take the outer bound from the COPPER THAT WAS ACTUALLY EMITTED, never from
                // leadTipRadius's formula: extendBorder keeps MKF's own border waypoint whenever
                // it reaches further out, so the real tip can sit well beyond
                // maxTurnRadius + 8*wireRadius (measured on 23_interleaved_llc_pq3530: tip at
                // r = 16.0 mm against a formula value of 14.3 mm). Probing only to the formula
                // samples the inside of the winding window, reports every azimuth free, aims
                // nothing -- and the lead runs on into the core plate.
                double rLead = 0.0;
                for (const auto& pr : path.prims) {
                    if (!pr.isLead) continue;
                    auto [pa, pb] = primEndpoints(pr);
                    rLead = std::max({rLead, std::hypot(pa.X(), pa.Z()), std::hypot(pb.X(), pb.Z())});
                }
                // Probe out to the CORE's own outer envelope, never just the emitted
                // copper's tip (Alf, 2026-08-07, 10_emi/EP13): the connection must EXIT
                // the core, so an azimuth is free only if the whole ray out PAST the core
                // clears it. An EP-style wrapping plate sits BEYOND the lead-tip radius --
                // stopping the probe at the copper reported every azimuth free and left
                // the fan on EP13's CLOSED -Z face.
                double rCore = 0.0;
                for (const auto& obst : opts.coreObstacles) {
                    Bnd_Box bx; BRepBndLib::Add(obst, bx);
                    double x0, y0, z0, x1, y1, z1; bx.Get(x0, y0, z0, x1, y1, z1);
                    rCore = std::max({rCore, std::hypot(x0, z0), std::hypot(x0, z1),
                                      std::hypot(x1, z0), std::hypot(x1, z1)});
                }
                const double rTip =
                    std::max({maxTurnR + 4.0 * (2.0 * wireRadius), rLead, rCore});
                const double r0 = maxTurnR + wireRadius;
                const double yProbe[3] = {0.0, 0.6 * maxAbsY, -0.6 * maxAbsY};
                const int N = 360, NR = 16;   // fine enough radially to never step over a core plate
                std::vector<char> free_(N, 1);
                for (const auto& obst : opts.coreObstacles) {
                    for (TopExp_Explorer se(obst, TopAbs_SOLID); se.More(); se.Next()) {
                        BRepClass3d_SolidClassifier cls(se.Current());
                        for (int k = 0; k < N; ++k) {
                            if (!free_[k]) continue;
                            const double az = kTwoPi * k / N;
                            for (int ri = 0; ri <= NR && free_[k]; ++ri) {
                                const double rp = r0 + (rTip - r0) * ri / NR;
                                for (double yp : yProbe) {
                                    cls.Perform(azPointC(0, 0, rp, yp, az), 1e-7);
                                    if (cls.State() == TopAbs_IN || cls.State() == TopAbs_ON) {
                                        free_[k] = 0;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                // widest free arc (circular)
                int bestLen = 0, bestStart = 0, curLen = 0, curStart = 0;
                for (int k = 0; k < 2 * N; ++k) {
                    if (free_[k % N]) {
                        if (curLen == 0) curStart = k;
                        if (++curLen > bestLen) { bestLen = curLen; bestStart = curStart; }
                        if (curLen >= 2 * N - 1) break;   // fully free
                    } else curLen = 0;
                }
                if (bestLen > 0 && bestLen < 2 * N - 1) {
                    exitAz = kTwoPi * (bestStart + 0.5 * bestLen) / N;
                    freeArc = kTwoPi * bestLen / N;
                    if (std::getenv("MVB_DIAG"))
                        std::cerr << "[lead-aim] " << path.name << ": core-free arc "
                                  << (360.0 * bestLen / N) << " deg, exit az "
                                  << (exitAz * 180.0 / kPi) << " deg (rTip=" << rTip * 1e3
                                  << "mm)\n";
                } else if (bestLen == 0) {
                    // FULLY BLOCKED: the core occupies every azimuth over the lead's radial run,
                    // so there is no opening to aim at and any azimuth terminates the lead inside
                    // the core. Refuse loudly -- silently keeping the default azimuth is what let
                    // 23_interleaved_llc's secondaries exit at 145-150 deg straight into the PQ
                    // plate and still export.
                    throw std::runtime_error(
                        "ConductorBuilder: no core-free exit azimuth for " + path.name +
                        "'s terminal leads -- the core occupies every azimuth over the lead's "
                        "radial run r=[" + std::to_string(r0 * 1e3) + "," +
                        std::to_string(rTip * 1e3) + "] mm. Set Options::leadExitAzimuth "
                        "explicitly, or shorten the lead so it stays inside the window.");
                } else if (std::getenv("MVB_DIAG")) {
                    // Genuinely free at every azimuth: nothing to aim at, keep the default.
                    std::cerr << "[lead-aim] " << path.name << ": no aim needed -- core free at "
                              << "every azimuth over r=[" << r0 * 1e3 << "," << rTip * 1e3
                              << "]mm\n";
                }
            }
            // The vertical fan swings with the seam as one rigid piece: if it is wider than the
            // core's free opening, some lead necessarily runs into the core no matter where the
            // aim points -- refuse loudly rather than export copper through a core plate.
            if (freeArc < kTwoPi && fanWidth > freeArc) {
                throw std::runtime_error(
                    "ConductorBuilder: the vertical-connection fan (" +
                    std::to_string(fanWidth * 180.0 / kPi) + " deg) is wider than the core's "
                    "free opening (" + std::to_string(freeArc * 180.0 / kPi) + " deg) for " +
                    path.name + " -- the window's connections cannot all clear the core");
            }
            if (!std::isnan(exitAz)) seamAngle += exitAz - kPlaneAz;
        }
        // NO TANGENTIAL SPREAD (Alf, 2026-08-04): every winding's and every parallel's terminals
        // enter and leave on the SAME station plane, straight out along -Z. The spread existed
        // only to stop different conductors' leads from sharing that plane, but it moved leads out
        // of the core's opening (23_interleaved_llc_pq3530: secondaries at 145-150 deg, inside the
        // PQ plate) and it makes the geometry far harder to read. Lead-lead clearance is the
        // layout's job; where the leads genuinely collide the gate now says so instead of the
        // spread hiding it.
        rotatePathAboutY(path, seamAngle);
        path.seamRot = seamAngle;
        paths.push_back(std::move(path));
    }

    // MVB_LEAD_CORE_CHECK: verify what the aim only PREDICTED -- classify the emitted lead
    // copper against the core solids and report any sample inside one. The aim probes a coarse
    // azimuth/radius/height lattice before the leads exist; this checks the actual centreline
    // after the seam rotation, which is the claim that matters.
    if (std::getenv("MVB_LEAD_CORE_CHECK") && !opts.coreObstacles.empty()) {
        int hits = 0;
        for (const auto& p : paths) {
            for (const auto& pr : p.prims) {
                if (!pr.isLead) continue;
                for (const auto& q : samplePrim(pr, p.wireRadius)) {
                    for (const auto& obst : opts.coreObstacles) {
                        for (TopExp_Explorer se(obst, TopAbs_SOLID); se.More(); se.Next()) {
                            BRepClass3d_SolidClassifier cls(se.Current());
                            cls.Perform(q, 1e-7);
                            if (cls.State() == TopAbs_IN) {
                                std::cerr << "[lead-core] " << p.name << " '" << pr.label
                                          << "' sample (" << q.X() * 1e3 << "," << q.Y() * 1e3
                                          << "," << q.Z() * 1e3 << ") mm is INSIDE the core\n";
                                ++hits;
                            }
                        }
                    }
                }
            }
        }
        std::cerr << "[lead-core] VERDICT: " << (hits ? "FAIL" : "PASS") << " (" << hits
                  << " lead samples inside the core)\n";
    }

    // Deferred Z-return end-runs: planned against the FULL obstacle field (every conductor's
    // helices, links and leads, plus already-planned end-runs), choosing a conflict-free
    // azimuth lane per return.
    planZEndRuns(paths, pendingZ, allRings);

    // Clean zero-net-progress retraces before anything consumes the paths, so the collision gate,
    // the sweep and the compound all see the same simple path.
    for (auto& p : paths) {
        const std::size_t n = dropRedundantExcursions(p);
        if (n && std::getenv("MVB_DIAG"))
            std::cerr << "[path] " << p.name << ": dropped " << n
                      << " retrace primitive(s) (out-and-back spurs)\n";
    }

    checkCollisions(paths);

    // Polyline capture mode: everything above ran (assembly, seam aiming, end-run planning,
    // the collision gate) but NO solid is built -- return the sampled centrelines instead.
    if (polyOut) {
        for (const auto& p : paths) {
            ConductorBuilder::PathPolyline pl;
            pl.name = p.name;
            pl.wireRadius = p.wireRadius;
            pl.isRectangular = p.isRectangular;
            pl.wireWidth = p.wireWidth;
            pl.wireHeight = p.wireHeight;
            for (const auto& pr : p.prims) {
                auto pts = samplePrim(pr, p.wireRadius);
                if (pts.size() < 2) continue;
                std::vector<std::array<double, 3>> seg;
                seg.reserve(pts.size());
                for (const auto& q : pts) seg.push_back({q.X(), q.Y(), q.Z()});
                pl.prims.push_back(std::move(seg));
            }
            if (!p.prims.empty()) {
                auto [a0, b0] = primEndpoints(p.prims.front());
                auto [a1, b1] = primEndpoints(p.prims.back());
                pl.end0 = {a0.X(), a0.Y(), a0.Z()};
                pl.end1 = {b1.X(), b1.Y(), b1.Z()};
                auto sf = samplePrim(p.prims.front(), p.wireRadius);
                auto sl = samplePrim(p.prims.back(), p.wireRadius);
                if (sf.size() >= 2) {
                    gp_XYZ d = sf.front().XYZ() - sf[1].XYZ();     // outward at the entrance tip
                    if (d.Modulus() > 1e-12) { d /= d.Modulus(); pl.dir0 = {d.X(), d.Y(), d.Z()}; }
                }
                if (sl.size() >= 2) {
                    gp_XYZ d = sl.back().XYZ() - sl[sl.size() - 2].XYZ();  // outward at the exit tip
                    if (d.Modulus() > 1e-12) { d /= d.Modulus(); pl.dir1 = {d.X(), d.Y(), d.Z()}; }
                }
            }
            polyOut->push_back(std::move(pl));
        }
        return {};
    }

    std::vector<NamedShape> out;
    out.reserve(paths.size());
    for (const auto& p : paths) {
        TopoDS_Shape cond = emitConductor(p, opts.wirePolygonSegments);
        out.push_back({cond, p.name});
        // PORT SURFACES for full-wave / FEM: the two free ends of a continuous round conductor are
        // its only PLANAR faces (the swept lateral surface is a cylinder/torus/BSpline, the elbows
        // are spheres/tori) -- MVB++ caps them flat precisely so a solver can put a port BC there.
        // Emit each as its own named face shape "<name> terminal <k>" so the mesher tags a named
        // surface group. Works for BOTH the single-body conductor and the multi-solid compound:
        // the free ends are known from the PATH (not guessed from "first planar face"), and the
        // lead-tip caps are planar full-cross-section faces in either representation. (The old
        // nsol == 1 gate assumed a compound has "no single pair of terminals" -- it does: the
        // same two path free ends.)
        if (!opts.femReady || p.isRectangular) continue;
        int nsol = 0;
        for (TopExp_Explorer e(cond, TopAbs_SOLID); e.More(); e.Next()) ++nsol;
        if (nsol < 1 || p.prims.empty()) continue;
        // The two free ends are KNOWN from the path -- take them from it, never from "the first two
        // planar faces found". That guess only held for the boolean-free single-body sweep, whose
        // lateral surfaces are all curved. The analytic per-primitive conductor is welded by a
        // boolean that leaves planar faces at junctions, so the guess returned the entrance-lead tip
        // and a 9e-18 m3 half-disc sliver at the lead/first-turn joint. A solver driven from THAT
        // pair pushes its whole current through 1.7 mm of lead and leaves the six turns dead -- the
        // FEM then reports |I_net| ~ 3e-7 A and a meaningless R_ac/R_dc, which reads as a solver
        // bug. Match each free end to the planar face centred on it, and demand a full wire
        // cross-section so a junction sliver can never be picked. (ABT #332)
        const gp_Pnt freeEnds[2] = {primEndpoints(p.prims.front()).first,
                                    primEndpoints(p.prims.back()).second};
        const double wantArea = kPi * p.wireRadius * p.wireRadius;
        for (int k = 0; k < 2; ++k) {
            TopoDS_Face bestFace;
            double bestD = 0.75 * p.wireRadius;
            for (TopExp_Explorer fe(cond, TopAbs_FACE); fe.More(); fe.Next()) {
                const TopoDS_Face& f = TopoDS::Face(fe.Current());
                if (BRepAdaptor_Surface(f).GetType() != GeomAbs_Plane) continue;
                GProp_GProps fp;
                BRepGProp::SurfaceProperties(f, fp);
                if (fp.Mass() < 0.80 * wantArea || fp.Mass() > 1.25 * wantArea) continue;
                const double d = fp.CentreOfMass().Distance(freeEnds[k]);
                if (d < bestD) { bestD = d; bestFace = f; }
            }
            if (bestFace.IsNull()) {
                std::cerr << "WARN ConductorBuilder: no terminal cap face found at free end " << k
                          << " of '" << p.name << "' -- the FEM port surface will be missing\n";
                continue;
            }
            out.push_back({bestFace, p.name + " terminal " + std::to_string(k)});
        }
    }
    return out;
}

std::vector<NamedShape> ConductorBuilder::buildAll(
    const MAS::Coil& coil, const MAS::CoreBobbinProcessedDescription& bobbin, bool isToroidal,
    const Options& opts) {
    // The reserved-space computation lives on OpenMagnetics::Coil; rebuild one (without
    // winding — the descriptions are already present) to obtain the drawn routes.
    nlohmann::json cj;
    to_json(cj, coil);
    OpenMagnetics::Coil omCoil(cj, /*windInConstructor=*/false);
    auto spaces = omCoil.get_connection_reserved_spaces();
    return buildAllImpl<MAS::Coil, MAS::Wire>(coil, bobbin, isToroidal, std::move(spaces), opts);
}

std::vector<NamedShape> ConductorBuilder::buildAll(
    const OpenMagnetics::Coil& coil, const MAS::CoreBobbinProcessedDescription& bobbin,
    bool isToroidal, const Options& opts) {
    OpenMagnetics::Coil coilCopy = coil;
    auto spaces = coilCopy.get_connection_reserved_spaces();
    return buildAllImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(coil, bobbin, isToroidal,
                                                                  std::move(spaces), opts);
}

std::vector<ConductorBuilder::PathPolyline> ConductorBuilder::buildAllPaths(
    const OpenMagnetics::Coil& coil, const MAS::CoreBobbinProcessedDescription& bobbin,
    bool isToroidal, const Options& opts) {
    OpenMagnetics::Coil coilCopy = coil;
    auto spaces = coilCopy.get_connection_reserved_spaces();
    std::vector<PathPolyline> out;
    buildAllImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(coil, bobbin, isToroidal,
                                                           std::move(spaces), opts, &out);
    return out;
}

} // namespace mvb
