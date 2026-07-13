#include "mvb/ConductorBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
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
#include <BRepAlgoAPI_Fuse.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BOPAlgo_GlueEnum.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopExp_Explorer.hxx>
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

constexpr double kPi = std::numbers::pi;
constexpr double kTwoPi = 2.0 * std::numbers::pi;
// Adjacent MKF slots are exactly one wire OD apart, so contact is normal; only
// penetration beyond numeric fuzz is a collision.
constexpr double kContactTol = 1e-7;
// Curved primitives are collision-checked as sampled polylines; the sampling step is
// chosen so each chord sags inward by at most this fraction of the wire radius. The same
// bound is granted back as contact allowance in the gate.
constexpr double kMaxSagFraction = 0.02;
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
};
struct ConductorPath {
    std::string name;
    double wireRadius = 0.0;
    std::vector<Primitive> prims;
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
    double maxSag = std::max(1e-6, kMaxSagFraction * wireRadius);
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
        double maxSag = std::max(1e-6, kMaxSagFraction * wireRadius);
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
void checkCollisions(const std::vector<ConductorPath>& paths) {
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
            double sagAllowance = kMaxSagFraction * (A.wireRadius + B.wireRadius);
            double minGap = A.wireRadius + B.wireRadius - sagAllowance - kContactTol;
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
                    if (d < minGap) {
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
                          << (A.wireRadius + B.wireRadius)
                          << " m (wire envelopes overlap). Turn positions are never moved; "
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
    gp_Ax2 plane(center, normal);
    if (segments <= 0) {
        gp_Circ circ(plane, radius);
        return BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire();
    }
    BRepBuilderAPI_MakePolygon poly;
    gp_Dir dx = plane.XDirection();
    gp_Dir dy = plane.YDirection();
    double offset = kPi / segments;
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
                           std::max(1e-6, kMaxSagFraction * wireRadius)));
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
        // Exact circular profile: one swept surface per edge instead of one per polygon
        // facet (a 16-gon profile makes multi-wrap sweeps take minutes), and the true
        // wire cross-section is round anyway.
        TopoDS_Wire prof = wireProfileWire(firstPts.front(), gp_Dir(t0), wireRadius, 0);

        BRepOffsetAPI_MakePipeShell ps(spine);
        // Mitre corners: RoundCorner's trim machinery needs edges longer than its
        // rounding radius and segfaults on the short MKF link segments
        // (BRepFill_TrimShellCorner::Perform, OCCT 7.9).
        ps.SetTransitionMode(BRepBuilderAPI_RightCorner);
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
        bool leadBoundary = (i + 1 < path.prims.size()) &&
                            (path.prims[i].isLead != path.prims[i + 1].isLead);
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
                            int wirePolygonSegments) {
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

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
            builder.Add(compound,
                        BRepPrimAPI_MakeCylinder(gp_Ax2(pr.seg.a, dir), wireRadius, len)
                            .Shape());
            // Sphere joints keep the wire envelope continuous across the plane elbows.
            builder.Add(compound, BRepPrimAPI_MakeSphere(pr.seg.a, wireRadius).Shape());
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
                    TopoDS_Wire prof = wireProfileWire(pts.front(), gp_Dir(t0),
                                                       wireRadius, 0);
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
                builder.Add(compound, BRepPrimAPI_MakeSphere(pts[j], wireRadius).Shape());
            }
            builder.Add(compound, BRepPrimAPI_MakeSphere(pts.back(), wireRadius).Shape());
        }
    }
    return compound;
}

// Fuse all solids of a compound — the sequential accumulator with a VOLUME GUARD.
// Adjacent layers of a winding touch along helical curves at exactly one wire OD — the
// most degenerate configuration for OCCT booleans, and fuses have been observed to
// silently EAT whole bodies on such inputs. A fuse is only accepted when its volume
// matches the summed piece volume (junction overlaps are well under 1%); otherwise the
// per-run compound is returned — geometrically exact, one PRODUCT in STEP.
TopoDS_Shape fuseAllSolids(const TopoDS_Shape& compound) {
    TopTools_ListOfShape solids;
    for (TopExp_Explorer exp(compound, TopAbs_SOLID); exp.More(); exp.Next()) {
        solids.Append(exp.Current());
    }
    if (solids.Extent() <= 1) return compound;
    auto volumeOf = [](const TopoDS_Shape& s) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(s, props);
        return props.Mass();
    };
    double vPieces = volumeOf(compound);
    try {
        TopTools_ListIteratorOfListOfShape it(solids);
        TopoDS_Shape acc = it.Value();
        it.Next();
        for (; it.More(); it.Next()) {
            BRepAlgoAPI_Fuse fuse(acc, it.Value());
            if (!fuse.IsDone() || fuse.Shape().IsNull()) {
                std::cerr << "WARN ConductorBuilder: fuse step not done; returning compound\n";
                return compound;
            }
            acc = fuse.Shape();
        }
        double vFused = volumeOf(acc);
        // Junction overlaps make vFused slightly SMALLER than the piece sum; losing any
        // whole piece removes far more than the overlaps ever could.
        constexpr double kMaxFuseVolumeLossFraction = 0.01;
        if (vFused < (1.0 - kMaxFuseVolumeLossFraction) * vPieces) {
            std::cerr << "WARN ConductorBuilder: fuse lost " << (vPieces - vFused)
                      << " m^3 of copper (OCCT tangent-contact defect); returning the "
                         "per-run compound instead\n";
            return compound;
        }
        return acc;
    } catch (const Standard_Failure& e) {
        std::cerr << "WARN ConductorBuilder: fuse threw (" << e.GetMessageString()
                  << "); returning compound\n";
        return compound;
    }
}

TopoDS_Shape emitConductor(const ConductorPath& path, int wirePolygonSegments) {
    // Sweep each maximal continuous run as ONE pipe (the whole wrap chain when the wire
    // never crosses itself); lead runs emit as exact cylinders + sphere elbows. Fuse the
    // handful of pieces into a single solid.
    std::vector<const Primitive*> ptrs;
    ptrs.reserve(path.prims.size());
    for (const auto& pr : path.prims) ptrs.push_back(&pr);

    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    for (auto [b, e] : continuousRuns(path)) {
        bool closedRing = (e - b == 1) && ptrs[b]->kind != Primitive::SEG &&
                          [&] {
                              auto [pa, pb] = primEndpoints(*ptrs[b]);
                              return pa.Distance(pb) < 1e-9;
                          }();
        bool leadRun = ptrs[b]->isLead;   // runs never mix lead and wrap primitives
        TopoDS_Shape run;
        if (!closedRing && !leadRun) {
            run = sweepRun(ptrs.data() + b, e - b, path.wireRadius, wirePolygonSegments);
        }
        if (std::getenv("MVB_DEBUG_RUNS") && run.IsNull() && !leadRun) {
            std::cerr << "PIPE FAILED for run of " << (e - b) << " prims starting ["
                      << ptrs[b]->label << "]\n";
        }
        if (run.IsNull()) {
            // Per-primitive sweeps for this span (arcs revolve exactly).
            run = sweepPiecewise(ptrs.data() + b, e - b, path.wireRadius,
                                 wirePolygonSegments);
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
        if (!run.IsNull()) builder.Add(compound, run);
    }
    return fuseAllSolids(compound);
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
    double borderX = 2.0 * run->coordinates.at(0) - station.x;
    double edgeY = run->coordinates.at(1);
    if (std::abs(edgeY - station.y) < 1e-12) {
        return {{station.x, station.y}, {borderX, station.y}};
    }
    return {{station.x, station.y}, {station.x, edgeY}, {borderX, edgeY}};
}

// ---------------------------------------------------------------------------------------
// Wrap planners. Each appends the connecting geometry of ONE wrap — crossing k-1 to
// crossing k — for its column shape. Under MKF's real-winding model the turnsDescription
// holds the N+1 window crossings of an N-turn winding (the first entry is the beginning
// of the first turn), so no closed loops exist anywhere.

// ROUND column: one full 360-degree spiral about the column axis — cylindrical within a
// layer, conical across a layer transition (radius and height both linear in azimuth).
void appendRoundWrap(ConductorPath& path, const PlanePt& s, const PlanePt& n,
                     const std::string& label, size_t ordinal) {
    Primitive wrap;
    wrap.kind = Primitive::SPIRAL;
    wrap.spiral = {0, 0, s.x, s.y, kPlaneAz, n.x, n.y, kPlaneAz + kTwoPi};
    wrap.label = label;
    wrap.turnOrdinal = ordinal;
    path.prims.push_back(std::move(wrap));
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

void appendRectWrap(ConductorPath& path, const RectStation& s0, const RectStation& s1,
                    const std::string& label, size_t ordinal) {
    auto pushSeg = [&](const gp_Pnt& a, const gp_Pnt& b, const char* what) {
        if (a.Distance(b) < 1e-12) return;
        Primitive pr;
        pr.kind = Primitive::SEG;
        pr.seg = {a, b};
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    auto pushCorner = [&](double cxx, double czz, double azStart, const char* what) {
        Primitive pr;
        pr.kind = Primitive::ARC3;
        pr.arc.c = gp_Pnt(cxx, s0.y, czz);
        pr.arc.axis = gp_XYZ(0, 1, 0);
        pr.arc.v0 = gp_XYZ(s0.cornerR * std::cos(azStart), 0,
                           -s0.cornerR * std::sin(azStart));
        pr.arc.sweep = kPi / 2.0;
        pr.label = label + std::string(" ") + what;
        pr.turnOrdinal = ordinal;
        path.prims.push_back(std::move(pr));
    };
    const double y = s0.y;
    // Crossing k-1 at the -Z face centre, moving -X (same chirality as the round wrap:
    // increasing azimuth from kPlaneAz).
    pushSeg(gp_Pnt(0, y, -s0.zPos), gp_Pnt(-s0.segX, y, -s0.zPos), "face -Z out");
    pushCorner(-s0.segX, -s0.segZ, kPi / 2.0, "corner -X-Z");
    pushSeg(gp_Pnt(-s0.xPos, y, -s0.segZ), gp_Pnt(-s0.xPos, y, +s0.segZ), "face -X");
    pushCorner(-s0.segX, +s0.segZ, kPi, "corner -X+Z");
    pushSeg(gp_Pnt(-s0.segX, y, +s0.zPos), gp_Pnt(+s0.segX, y, +s0.zPos), "face +Z");
    pushCorner(+s0.segX, +s0.segZ, 3.0 * kPi / 2.0, "corner +X+Z");
    pushSeg(gp_Pnt(+s0.xPos, y, +s0.segZ), gp_Pnt(+s0.xPos, y, -s0.segZ), "face +X");
    pushCorner(+s0.segX, -s0.segZ, 0.0, "corner +X-Z");
    // The transition: the final stretch of the -Z face carries the move to the next
    // crossing as S-blends whose tangent is -X at BOTH ends — tangent-continuous into
    // the corner behind and into the next wrap ahead, and passing exactly through the
    // crossing (a straight ramp would kink at both junctions).
    if (s0.segX < 1e-9) {
        throw std::runtime_error(
            "ConductorBuilder: corner bends consumed the whole -Z transition face for " +
            label + " — no room for the crossing ramp");
    }
    auto pushBlend = [&](const gp_Pnt& a, const gp_Pnt& b, const char* what) {
        if (a.Distance(b) < 1e-12) return;
        Primitive ramp;
        ramp.kind = Primitive::BLEND;
        ramp.blendc = {a, b, gp_XYZ(-1, 0, 0)};
        ramp.label = label + std::string(" ") + what;
        ramp.turnOrdinal = ordinal;
        path.prims.push_back(std::move(ramp));
    };
    if (std::abs(s1.zPos - s0.zPos) > 1e-12 && std::abs(s1.y - s0.y) > 1e-12) {
        // Layer transition that also moves axially (a Z-order return descends the whole
        // window). One diagonal blend would cut BETWEEN the layers, clipping every ramp
        // of the layer being crossed; a real winder pops the wire out to the new layer's
        // clearance first and only then runs it down at that clearance. Two planar
        // blends: radial move in the outer half of the face, axial move in the inner
        // half — the junction tangent is -X on both sides.
        gp_Pnt mid(+s0.segX / 2.0, y, -s1.zPos);
        pushBlend(gp_Pnt(+s0.segX, y, -s0.zPos), mid, "ramp -Z (radial)");
        pushBlend(mid, gp_Pnt(0, s1.y, -s1.zPos), "ramp -Z (axial)");
    } else {
        pushBlend(gp_Pnt(+s0.segX, y, -s0.zPos), gp_Pnt(0, s1.y, -s1.zPos), "ramp -Z");
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
                                     const ConductorBuilder::Options& opts) {
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
    const MAS::ColumnShape columnShape = bobbinPd.get_column_shape();
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

    auto planePoint = [&](double x2d, double y2d) {
        return azPointC(0, 0, x2d + zoff, y2d, kPlaneAz);
    };

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

    std::vector<ConductorPath> paths;
    paths.reserve(conductors.size());

    for (size_t ci = 0; ci < conductors.size(); ++ci) {
        const auto& ct = conductors[ci];
        const MAS::Wire& wire = wireMap.at(ct.winding);
        if (wire.get_type() == MAS::WireType::RECTANGULAR ||
            wire.get_type() == MAS::WireType::PLANAR ||
            wire.get_type() == MAS::WireType::FOIL) {
            throw std::runtime_error(
                "ConductorBuilder: rectangular/planar/foil wire is not supported for "
                "real-winding conductors yet (the lead cross-section orientation through "
                "the exit bends is undefined) — winding '" + ct.winding + "'");
        }
        // Resolve through the first turn: turn.dimensions carries the OUTER footprint
        // (the context-less overload would fall back to the round wire's CONDUCTING
        // diameter for paintCoating=true).
        auto [wireW, wireH] =
            TurnBuilder::wireDimensions(wire, *ct.turns.front(), opts.paintCoating);
        double wireRadius = std::min(wireW, wireH) / 2.0;
        // Same minimum-bend rule as build_concentric_rect_column_turn: a swept corner
        // self-intersects when the arc radius is below the profile's radial half-extent.
        double minBend = wireRadius * 1.02;

        ConductorPath path;
        path.name = ct.winding + " parallel " + std::to_string(ct.parallel);
        path.wireRadius = wireRadius;

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

            // Corrected crossing: an OUTER ring's outer crossing is re-placed at the
            // physical radial stack (ring 0's outer radius + ringIndex ODs), keeping
            // MKF's outer ANGLE — the inner crossing (the turn's real position) and ring 0
            // are untouched, so single-ring and single-layer builds are identical to
            // before. This is the ONE deviation from MKF's toroidal geometry, and only
            // because MKF's outer radius is self-overlapping (ABT #231); without it a
            // multilayer toroid cannot be realized in 3D without the wires intersecting.
            auto toroCross = [&](const MAS::Turn* t) -> ToroCross {
                ToroCross rc = toroCrossRaw(t);
                int k = ringIndexOfInner(rc.pin);
                if (k > 0) {
                    double ang = std::atan2(rc.pout.Y(), rc.pout.X());
                    double R = minRawOuter + k * od;
                    rc.pout = gp_XY(R * std::cos(ang), R * std::sin(ang));
                }
                return rc;
            };

            ToroCross first = toroCross(turns.front());
            ToroCross last = toroCross(turns.back());
            int maxRingIndex = 0;
            for (const MAS::Turn* t : turns)
                maxRingIndex = std::max(maxRingIndex, ringIndexOfInner(toroCrossRaw(t).pin));

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
                if (r0 != r1 && jump > kPi / 4.0) return maxRingIndex + 1;
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

            // Toroidal terminal lead. The inner crossing sits against the hole wall with
            // NO radial room (its envelope already touches the bore), so a radial lead —
            // MKF's pink box points radially outward, straight into the core, ABT #230 —
            // is unrealizable, and an axial riser collides with the neighbouring section's
            // riser. The physical lead runs TANGENT to the inner diameter: an arc at the
            // crossing radius, in the hole plane, peeling off along the wall into the
            // empty gap between this section and the next (away from the section's mean
            // azimuth). It is capped short of the neighbouring winding's nearest crossing
            // so two facing leads never meet.
            auto emitToroLead = [&](const ToroCross& cross, size_t ordinal,
                                    const std::string& who) {
                double crossR = cross.pin.Modulus();
                if (crossR < 1e-9) return;
                double crossAng = std::atan2(cross.pin.Y(), cross.pin.X());
                double delta = std::remainder(crossAng - meanAng, kTwoPi);
                double dirSign = (std::abs(delta) < 1e-9) ? 1.0 : (delta > 0 ? 1.0 : -1.0);

                // Angular room to the nearest OTHER crossing in the sweep direction.
                double gapAng = kTwoPi / 3.0;
                for (const auto& [pt, wr2] : allToroCrossings) {
                    if ((pt - cross.pin).Modulus() < 1e-9) continue;
                    double dd = std::remainder(std::atan2(pt.Y(), pt.X()) - crossAng, kTwoPi);
                    if (dirSign > 0 && dd > 1e-6) gapAng = std::min(gapAng, dd);
                    else if (dirSign < 0 && dd < -1e-6) gapAng = std::min(gapAng, -dd);
                }
                // Desired ~2.5 OD of arc, but capped so that when the facing lead (the
                // neighbouring section's, sweeping in symmetrically) also stops here, the
                // two ends stay a full OD apart: each takes (gap - 1.3 OD)/2 at most.
                double odAng = od / crossR;
                double arcAng = std::min(2.5 * odAng, std::max(0.0, (gapAng - 1.3 * odAng) / 2.0));
                if (arcAng < 1e-3) return;   // no room for a lead here

                Primitive pr;
                pr.kind = Primitive::ARC3;
                pr.arc.c = gp_Pnt(0, 0, 0);
                pr.arc.axis = gp_XYZ(0, -dirSign, 0);   // rotate the crossing about the
                pr.arc.v0 = gp_XYZ(cross.pin.X(), 0, cross.pin.Y());  // toroid axis
                pr.arc.sweep = arcAng;                  // away from the section centre
                pr.label = who + " lead";
                pr.turnOrdinal = ordinal;
                pr.isLead = true;
                path.prims.push_back(std::move(pr));
            };
            // MKF must have drawn both terminal leads (validates the winding is complete);
            // their radial geometry is replaced by the tangential arc above.
            (void)toroLeadRect;

            emitToroLead(first, 0, path.name + " entrance");

            for (size_t i = 0; i + 1 < turns.size(); ++i) {
                appendToroWrap(path, toroCross(turns[i]), toroCross(turns[i + 1]),
                               wireRadius, wrapDepthOds(i) * od,
                               "wrap '" + turns[i]->get_name() + "' -> '" +
                                   turns[i + 1]->get_name() + "'",
                               i);
            }

            emitToroLead(last, turns.size() - 1, path.name + " exit");

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
        {
            std::vector<std::vector<const RSpace*>> groups(1);
            for (const RSpace* s : terminalRects) {
                groups.back().push_back(s);
                if (!rectIsVertical(*s)) groups.emplace_back();   // run closes the group
            }
            if (!groups.empty() && groups.back().empty()) groups.pop_back();
            if (groups.size() != 2) {
                throw std::runtime_error(
                    "ConductorBuilder: expected 2 terminal lead groups (entrance, exit) "
                    "for " + path.name + ", got " + std::to_string(groups.size()));
            }
            entranceGroup = groups[0];
            exitGroup = groups[1];
        }

        auto pushPlaneSegs = [&](std::vector<PlanePt> wp, const std::string& what,
                                 size_t ordinal) {
            // Absorb intermediate waypoints closer than the wire radius to their
            // neighbour: a jog shorter than the wire's own radius lies entirely inside
            // the pipe body of the adjacent edge (and inside MKF's drawn rectangle,
            // whose height is one full wire OD), while sub-radius spine edges crash
            // OCCT's pipe-shell corner rounding. Endpoints (the exact stations/border)
            // are always kept.
            std::vector<PlanePt> kept;
            kept.push_back(wp.front());
            for (size_t i = 1; i + 1 < wp.size(); ++i) {
                if (std::hypot(wp[i].x - kept.back().x, wp[i].y - kept.back().y) <
                        wireRadius ||
                    std::hypot(wp[i].x - wp.back().x, wp[i].y - wp.back().y) < wireRadius) {
                    continue;
                }
                kept.push_back(wp[i]);
            }
            kept.push_back(wp.back());
            for (size_t i = 0; i + 1 < kept.size(); ++i) {
                if (std::hypot(kept[i + 1].x - kept[i].x, kept[i + 1].y - kept[i].y) < 1e-12) {
                    continue;
                }
                Primitive pr;
                pr.kind = Primitive::SEG;
                pr.seg = {planePoint(kept[i].x, kept[i].y),
                          planePoint(kept[i + 1].x, kept[i + 1].y)};
                pr.label = what + " seg " + std::to_string(i);
                pr.turnOrdinal = ordinal;
                pr.isLead = true;
                path.prims.push_back(std::move(pr));
            }
        };

        // Entrance: MKF's drawn route, walked from the border TO the first station.
        {
            auto wp = terminalWaypoints(entranceGroup, first, path.name + " entrance");
            std::reverse(wp.begin(), wp.end());
            pushPlaneSegs(wp, "entrance lead", 0);
        }

        // Wraps between consecutive crossings.
        for (size_t i = 0; i + 1 < turns.size(); ++i) {
            PlanePt s = station(turns[i]);
            PlanePt nxt = station(turns[i + 1]);
            std::string label = "wrap '" + turns[i]->get_name() + "' -> '" +
                                turns[i + 1]->get_name() + "'";
            if (effectivelyRound) {
                appendRoundWrap(path, s, nxt, label, i);
            } else if (columnShape == MAS::ColumnShape::RECTANGULAR) {
                appendRectWrap(path,
                               rectStation(s, halfW, halfD, minBend, path.name),
                               rectStation(nxt, halfW, halfD, minBend, path.name),
                               label, i);
            } else {   // OBLONG with a real straight section
                appendOblongWrap(path, s, nxt, oblongHalf, label, i);
            }
        }

        // Exit: MKF's drawn route from the last station out to the border.
        {
            auto wp = terminalWaypoints(exitGroup, last, path.name + " exit");
            pushPlaneSegs(wp, "exit lead", turns.size() - 1);
        }

        paths.push_back(std::move(path));
    }

    checkCollisions(paths);

    std::vector<NamedShape> out;
    out.reserve(paths.size());
    for (const auto& p : paths) {
        out.push_back({emitConductor(p, opts.wirePolygonSegments), p.name});
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

} // namespace mvb
