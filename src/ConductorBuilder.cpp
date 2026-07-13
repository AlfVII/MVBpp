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
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>
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
#include <BRepAlgoAPI_Fuse.hxx>
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
#include <limits>
#include <map>
#include <numbers>
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
// plane"). MKF's 2D winding-window cross-section maps into 3D at this azimuth:
// 2D (x = radial, y = axial) -> 3D (0, y, -x). Every turn starts/ends here; every drawn
// ConnectionReservedSpace rectangle (the pink/blue boxes of the Painter SVG) is replayed
// here verbatim — nothing is invented.
constexpr double kPlaneAz = kPi / 2.0;

// Azimuth convention (OCCT right-handed rotation about +Y, the column axis):
//   pos(r, y, az) = (r cos az, y, -r sin az)
gp_Pnt azPoint(double r, double y, double az) {
    return gp_Pnt(r * std::cos(az), y, -r * std::sin(az));
}
gp_Pnt planePoint(double radial, double axial) {
    return azPoint(radial, axial, kPlaneAz);
}

// ---------------------------------------------------------------------------------------
// Path model. ARC = open ring at a station (full wrap, no height change); SPIRAL = helical
// wrap or bulge piece; SEG = straight connection segment in the YZ plane (from MKF's
// drawn reserved-space rectangles).
struct Arc {
    double r = 0, y = 0;
    double azStart = 0, azSweep = 0;
};
struct Seg {
    gp_Pnt a, b;
};
struct Spiral {
    double r0 = 0, y0 = 0, az0 = 0;
    double r1 = 0, y1 = 0, az1 = 0;
};
struct Primitive {
    enum Kind { ARC, SEG, SPIRAL } kind = ARC;
    Arc arc{};
    Seg seg{};
    Spiral spiral{};
    std::string label;
    // Electrical turn ordinal this primitive belongs to (entrance lead = first turn's,
    // exit lead = last turn's, links = the source turn's). A continuous wire legitimately
    // contacts itself only between CONSECUTIVE turns; the gate exempts same-conductor
    // pairs with |ordinal diff| <= 1 and checks everything farther apart.
    size_t turnOrdinal = 0;
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
    if (p.kind == Primitive::ARC) {
        const Arc& a = p.arc;
        int n = curveSampleCount(a.r, a.azSweep, wireRadius);
        std::vector<gp_Pnt> pts;
        pts.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            double az = a.azStart + a.azSweep * i / (n - 1);
            pts.push_back(azPoint(a.r, a.y, az));
        }
        return pts;
    }
    const Spiral& sp = p.spiral;
    int n = curveSampleCount(std::max(sp.r0, sp.r1), std::abs(sp.az1 - sp.az0), wireRadius);
    std::vector<gp_Pnt> pts;
    pts.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        double t = static_cast<double>(i) / (n - 1);
        pts.push_back(azPoint(sp.r0 + (sp.r1 - sp.r0) * t, sp.y0 + (sp.y1 - sp.y0) * t,
                              sp.az0 + (sp.az1 - sp.az0) * t));
    }
    return pts;
}

double primDistance(const Primitive& pa, double ra, const Primitive& pb, double rb) {
    auto polyA = samplePrim(pa, ra);
    auto polyB = samplePrim(pb, rb);
    double best = std::numeric_limits<double>::max();
    for (size_t i = 0; i + 1 < polyA.size(); ++i)
        for (size_t j = 0; j + 1 < polyB.size(); ++j)
            best = std::min(best, segSegDistance(polyA[i], polyA[i + 1], polyB[j], polyB[j + 1]));
    return best;
}

bool arcsClearlyApart(const Arc& a, const Arc& b, double minGap) {
    return std::hypot(a.y - b.y, a.r - b.r) >= minGap;
}

std::pair<gp_Pnt, gp_Pnt> primEndpoints(const Primitive& p) {
    if (p.kind == Primitive::ARC) {
        return {azPoint(p.arc.r, p.arc.y, p.arc.azStart),
                azPoint(p.arc.r, p.arc.y, p.arc.azStart + p.arc.azSweep)};
    }
    if (p.kind == Primitive::SPIRAL) {
        return {azPoint(p.spiral.r0, p.spiral.y0, p.spiral.az0),
                azPoint(p.spiral.r1, p.spiral.y1, p.spiral.az1)};
    }
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
                    if (pa.kind == Primitive::ARC && pb.kind == Primitive::ARC &&
                        arcsClearlyApart(pa.arc, pb.arc, minGap)) {
                        continue;
                    }
                    if (ci == cj && shareEndpoint(pa, pb)) continue;
                    double d = primDistance(pa, A.wireRadius, pb, B.wireRadius);
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
//  - SPIRAL at constant radius and ARC rings: TRUE HELIX — a 2D line in the
//    (U = azimuth, V = height) parametric space of a cylinder about the column axis
//    (the canonical OCCT construction: Geom2d line on Geom_CylindricalSurface +
//    BRepLib::BuildCurves3d).
//  - SPIRAL with varying radius (crossover bulges): BSpline through the sampled
//    centreline.
TopoDS_Edge primEdge(const Primitive& pr, double wireRadius) {
    if (pr.kind == Primitive::SEG) {
        if (pr.seg.a.Distance(pr.seg.b) < 1e-12) return TopoDS_Edge();
        return BRepBuilderAPI_MakeEdge(pr.seg.a, pr.seg.b).Edge();
    }
    bool constantRadius =
        pr.kind == Primitive::ARC ||
        (pr.kind == Primitive::SPIRAL && std::abs(pr.spiral.r1 - pr.spiral.r0) < 1e-12);
    if (constantRadius) {
        double r = (pr.kind == Primitive::ARC) ? pr.arc.r : pr.spiral.r0;
        gp_Pnt2d p0, p1;
        if (pr.kind == Primitive::ARC) {
            p0 = gp_Pnt2d(pr.arc.azStart, pr.arc.y);
            p1 = gp_Pnt2d(pr.arc.azStart + pr.arc.azSweep, pr.arc.y);
        } else {
            p0 = gp_Pnt2d(pr.spiral.az0, pr.spiral.y0);
            p1 = gp_Pnt2d(pr.spiral.az1, pr.spiral.y1);
        }
        try {
            // Cylinder about +Y through the origin with U measured from +X toward -Z —
            // exactly the azPoint() convention: P(U,V) = (R cos U, V, -R sin U).
            gp_Ax3 axis(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0), gp_Dir(1, 0, 0));
            Handle(Geom_CylindricalSurface) cyl = new Geom_CylindricalSurface(axis, r);
            Handle(Geom2d_TrimmedCurve) seg2d = GCE2d_MakeSegment(p0, p1).Value();
            TopoDS_Edge e = BRepBuilderAPI_MakeEdge(seg2d, cyl).Edge();
            BRepLib::BuildCurves3d(e);
            return e;
        } catch (const Standard_Failure&) {
            return TopoDS_Edge();
        }
    }
    // Radius-varying spiral: BSpline through the sampled centreline.
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
        return ps.Shape();
    } catch (const Standard_Failure&) {
        return TopoDS_Shape();
    } catch (const std::exception&) {
        return TopoDS_Shape();
    }
}

// Split the path into maximal continuous runs at revisited points (the wire crossing
// itself at ring junctions). Everything within a run sweeps as one pipe; runs meet only
// at the crossover points and are fused afterwards.
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
    // Closed rings (start == end) are isolated into their own single-primitive runs:
    // a spine containing a closed loop is self-touching by construction, and OCCT's
    // pipe-shell segfaults on it — the isolated ring is swept by the exact revolve in
    // the piecewise fallback instead.
    std::vector<std::pair<size_t, size_t>> runs;
    auto isClosedRing = [&](const Primitive& pr) {
        auto [a, b] = primEndpoints(pr);
        return quantize(a) == quantize(b);
    };
    size_t start = 0;
    for (size_t i = 0; i < path.prims.size(); ++i) {
        auto [a, b] = primEndpoints(path.prims[i]);
        if (isClosedRing(path.prims[i])) {
            if (i > start) runs.push_back({start, i});
            runs.push_back({i, i + 1});
            start = i + 1;
            continue;
        }
        // Split at crossover junctions AND at SEG<->SPIRAL kind changes: a pipe swept
        // across the 90-degree lead-to-wrap corner grows an unbounded mitre spike
        // (observed: centimetre-scale flare). Pure-spiral chains are tangent-continuous
        // (no corners), and pure-SEG chains emit as exact cylinders with sphere elbows.
        bool kindChange = (i + 1 < path.prims.size()) &&
                          (path.prims[i].kind == Primitive::SEG) !=
                              (path.prims[i + 1].kind == Primitive::SEG);
        if (i + 1 < path.prims.size() && (visits[quantize(b)] > 2 || kindChange)) {
            runs.push_back({start, i + 1});
            start = i + 1;
        }
    }
    if (start < path.prims.size()) runs.push_back({start, path.prims.size()});
    return runs;
}

// Per-piece fallback: sweep each primitive on its own (pipe along its single edge,
// revolve for rings), then fuse everything into one solid.
TopoDS_Shape sweepPiecewise(const Primitive* const* prims, size_t count, double wireRadius,
                            int wirePolygonSegments) {
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    for (size_t pi = 0; pi < count; ++pi) {
        const auto& pr = *prims[pi];
        if (pr.kind == Primitive::ARC) {
            gp_Pnt c = azPoint(pr.arc.r, pr.arc.y, pr.arc.azStart);
            gp_Dir tangent(-std::sin(pr.arc.azStart), 0.0, -std::cos(pr.arc.azStart));
            TopoDS_Face prof = wireProfile(c, tangent, wireRadius, 0);
            BRepPrimAPI_MakeRevol rev(prof, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)),
                                      pr.arc.azSweep);
            if (!rev.IsDone() || rev.Shape().IsNull()) {
                throw std::runtime_error("ConductorBuilder: ring revolve failed for [" +
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

// Fuse all solids of a compound in ONE boolean operation (arguments + tools, parallel
// mode) — the sequential pairwise accumulator is O(n^2) in faces and grinds for a
// 34-wrap conductor. Returns the compound unchanged (with a warning) on failure.
TopoDS_Shape fuseAllSolids(const TopoDS_Shape& compound) {
    TopTools_ListOfShape solids;
    for (TopExp_Explorer exp(compound, TopAbs_SOLID); exp.More(); exp.Next()) {
        solids.Append(exp.Current());
    }
    if (solids.Extent() <= 1) return compound;
    // VOLUME-GUARDED fuse. Adjacent layers of a winding touch along helical curves at
    // exactly one wire OD — the most degenerate configuration for OCCT booleans, and
    // both the multi-tool and the sequential fuse have been observed to silently EAT
    // whole bodies (rings, layers) on such inputs. A fuse is only accepted when its
    // volume matches the summed piece volume (junction overlaps are well under 1%);
    // otherwise the per-run compound is returned — geometrically exact, one PRODUCT in
    // STEP, a dozen meaningful solids instead of one.
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
    // Sweep each maximal continuous run as ONE pipe (the whole conductor when the wire
    // never crosses itself); runs meet only at the crossover points. Fuse the handful of
    // pieces into a single solid in one boolean.
    std::vector<const Primitive*> ptrs;
    ptrs.reserve(path.prims.size());
    for (const auto& pr : path.prims) ptrs.push_back(&pr);

    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    for (auto [b, e] : continuousRuns(path)) {
        // Single-primitive closed rings go straight to the exact revolve (a closed spine
        // segfaults OCCT's pipe-shell); everything else sweeps as one pipe.
        bool closedRing = (e - b == 1) && ptrs[b]->kind == Primitive::ARC;
        bool segChain = ptrs[b]->kind == Primitive::SEG;   // runs are kind-homogeneous
        TopoDS_Shape run;
        if (!closedRing && !segChain) {
            run = sweepRun(ptrs.data() + b, e - b, path.wireRadius, wirePolygonSegments);
        }
        if (run.IsNull()) {
            // Per-primitive sweeps for this span (rings revolve exactly).
            run = sweepPiecewise(ptrs.data() + b, e - b, path.wireRadius,
                                 wirePolygonSegments);
        }
        if (!run.IsNull()) builder.Add(compound, run);
    }
    return fuseAllSolids(compound);
}

// ---------------------------------------------------------------------------------------
// Connection replay: MKF's drawn ConnectionReservedSpace rectangles (layer == "") are the
// authoritative routes — the pink (terminal) and blue (link) boxes of the Painter SVG.
// Each rect maps to centreline waypoints in the YZ plane.
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
    const RSpace* stub = nullptr;
    for (const RSpace* s : group) {
        if (rectIsVertical(*s)) stub = s; else run = s;
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
    if (isToroidal) {
        throw std::runtime_error("ConductorBuilder: toroidal cores are not supported yet "
                                 "(real-winding milestone 3)");
    }
    if (bobbinPd.get_column_shape() != MAS::ColumnShape::ROUND) {
        throw std::runtime_error("ConductorBuilder: only ROUND columns are supported yet "
                                 "(rectangular/oblong arrive in real-winding milestone 2)");
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

    auto station = [](const MAS::Turn* t) -> PlanePt {
        const auto& c = t->get_coordinates();
        if (c.size() < 2) {
            throw std::runtime_error("ConductorBuilder: turn '" + t->get_name() +
                                     "' has fewer than 2 coordinates");
        }
        return {c[0], c[1]};
    };

    std::vector<ConductorPath> paths;
    paths.reserve(conductors.size());

    for (size_t ci = 0; ci < conductors.size(); ++ci) {
        const auto& ct = conductors[ci];
        auto [wireW, wireH] =
            TurnBuilder::wireDimensions(wireMap.at(ct.winding), opts.paintCoating);
        double wireRadius = std::min(wireW, wireH) / 2.0;
        double od = 2.0 * wireRadius;

        ConductorPath path;
        path.name = ct.winding + " parallel " + std::to_string(ct.parallel);
        path.wireRadius = wireRadius;

        const auto& turns = ct.turns;
        PlanePt first = station(turns.front());
        PlanePt last = station(turns.back());

        // This conductor's drawn rects: terminals split entrance/exit by which station
        // the group touches; links matched per jump geometrically.
        std::vector<const RSpace*> terminalRects;
        for (const auto& s : drawn) {
            if (s.winding != ct.winding || s.parallel != ct.parallel) continue;
            if (s.isTerminal) terminalRects.push_back(&s);
            // Non-terminal (blue) link rects stay 2D documentation: the conical wrap
            // between the crossings IS the 3D transition.
        }
        // Terminal grouping follows MKF's emission order (verified in
        // get_connection_reserved_spaces): per (winding, parallel) the ENTRANCE lead's
        // rects are pushed first, then the EXIT lead's, each as [vertical stub?,
        // horizontal run] — the run closes its group. Proximity-based grouping is
        // ambiguous when a crossing sits within a wire of the window edge (the entrance
        // edge run then "touches" the exit station too).
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
                path.prims.push_back(std::move(pr));
            }
        };

        // Entrance: MKF's drawn route, walked from the border TO the first station.
        {
            auto wp = terminalWaypoints(entranceGroup, first, path.name + " entrance");
            std::reverse(wp.begin(), wp.end());
            pushPlaneSegs(wp, "entrance lead", 0);
        }

        // Wraps between consecutive CROSSINGS. Under MKF's real-winding model the
        // turnsDescription holds the N+1 window crossings of an N-turn winding (the
        // first entry is the beginning of the first turn), so wrap k sweeps one full
        // 360-degree spiral from crossing k-1 to crossing k: a cylindrical helix within
        // a layer, a conical spiral across a layer transition (radius and height both
        // linear in azimuth). No closed rings exist anywhere and consecutive parallels
        // advance in lockstep, so the whole conductor is one continuous, never
        // self-touching path from entrance lead to exit lead.
        for (size_t i = 0; i + 1 < turns.size(); ++i) {
            PlanePt s = station(turns[i]);
            PlanePt nxt = station(turns[i + 1]);
            Primitive wrap;
            wrap.kind = Primitive::SPIRAL;
            wrap.spiral = {s.x, s.y, kPlaneAz, nxt.x, nxt.y, kPlaneAz + kTwoPi};
            wrap.label = "wrap '" + turns[i]->get_name() + "' -> '" +
                         turns[i + 1]->get_name() + "'";
            wrap.turnOrdinal = i;
            path.prims.push_back(std::move(wrap));
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
