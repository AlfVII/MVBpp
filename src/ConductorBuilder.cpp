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
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Circ.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace mvb {

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kTwoPi = 2.0 * std::numbers::pi;
// Adjacent MKF slots are exactly one wire OD apart, so contact is normal; only
// penetration beyond numeric fuzz is a collision.
constexpr double kContactTol = 1e-7;
// Curved primitives are collision-checked as sampled polylines; the sampling step is
// chosen so each chord sags inward by at most this fraction of the wire radius. The same
// bound is granted back as contact allowance in the gate (a chord can under-report the
// true curve-to-curve distance by exactly this much per side).
constexpr double kMaxSagFraction = 0.02;
// Terminal stubs of neighbouring ranks are spaced this many wire ODs apart in arc length
// (one OD would be touching; half an OD of daylight keeps distinct leads distinct).
constexpr double kStubSpacingODs = 1.5;
// Arc length (in wire ODs) a layer-jump crossover consumes: it must span at least its own
// wire (1 OD) to ride over the wrap start below, plus the same again to come down.
constexpr double kCrossoverArcODs = 2.0;
// Hard cap on the crossover arc so thick-wire/small-core combinations degrade into a
// visible error (via the gate) instead of a quarter-toroid of diagonal wire.
constexpr double kMaxCrossoverArc = kPi / 4.0;
// The bulge ramps outside the plateau re-approach the layer radius; their width is a
// shape choice (any positive width clears, since the wrap is already one OD away
// axially by the plateau edge — see the bulge derivation at the call site).
constexpr double kBulgeRampFraction = 0.3;

// Azimuth convention (OCCT right-handed rotation about +Y, the column axis):
//   pos(r, y, az) = (r cos az, y, -r sin az)
// Azimuth +pi/2 points along -Z — the transition-face convention (design decision B):
// the seam sector, where ALL synthesized connecting geometry lives, is centred there.
gp_Pnt azPoint(double r, double y, double az) {
    return gp_Pnt(r * std::cos(az), y, -r * std::sin(az));
}

// ---------------------------------------------------------------------------------------
// Path model: an ordered list of primitives, each an on-station ARC (exact MKF turn
// position) or a synthesized straight SEGMENT (ramp / stub / edge run).
struct Arc {
    double r = 0, y = 0;
    double azStart = 0, azSweep = 0;   // sweep > 0, direction of increasing azimuth
};
struct Seg {
    gp_Pnt a, b;
};
// Helical segment: linear interpolation in (radius, height, azimuth). Ramps and layer
// jumps must be helical, not straight chords — a chord spanning tens of degrees sags
// radially inward by r*(1-cos(span/2)) and eats the inter-layer clearance.
struct Spiral {
    double r0 = 0, y0 = 0, az0 = 0;
    double r1 = 0, y1 = 0, az1 = 0;
};
struct Primitive {
    enum Kind { ARC, SEG, SPIRAL } kind = ARC;
    Arc arc{};
    Seg seg{};
    Spiral spiral{};
    std::string label;                 // for collision diagnostics
    // Electrical turn ordinal this primitive belongs to (entrance lead = first turn's,
    // exit lead = last turn's). A continuous wire legitimately contacts itself only
    // between CONSECUTIVE turns (spring wraps resting on each other, crossovers riding
    // the wrap below); the gate exempts same-conductor pairs with |ordinal diff| <= 1
    // and checks everything farther apart.
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

// Sample step so the polyline's chord sag keeps the capsule check honest.
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

// Concentric arcs: exact cheap prefilter before polyline sampling.
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

// Same-conductor primitives that share an endpoint are path-connected (they MEET there by
// construction — a continuous wire); only their interiors can genuinely collide, and any
// real interior overlap is still caught by every non-touching pair along the path.
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
            // Grant the sampling sag back as contact allowance (each polyline can sit up
            // to kMaxSagFraction*wireRadius inside its true curve), or exact-contact
            // pairs — bulge plateau over a crossed turn, layers exactly one OD apart —
            // false-positive by microns.
            double sagAllowance = kMaxSagFraction * (A.wireRadius + B.wireRadius);
            double minGap = A.wireRadius + B.wireRadius - sagAllowance - kContactTol;
            for (size_t i = 0; i < A.prims.size(); ++i) {
                size_t jStart = (ci == cj) ? i + 1 : 0;
                for (size_t j = jStart; j < B.prims.size(); ++j) {
                    const auto& pa = A.prims[i];
                    const auto& pb = B.prims[j];
                    // Same conductor: consecutive turns of a spring legitimately touch
                    // (wraps resting on each other, crossovers riding the wrap below);
                    // farther-apart turns must genuinely clear.
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
TopoDS_Face wireProfile(const gp_Pnt& center, const gp_Dir& normal, double radius,
                        int segments) {
    gp_Ax2 plane(center, normal);
    if (segments <= 0) {
        gp_Circ circ(plane, radius);
        return BRepBuilderAPI_MakeFace(
                   BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire())
            .Face();
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
    return BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakeWire(poly.Wire()).Wire()).Face();
}

TopoDS_Shape emitConductor(const ConductorPath& path, int wirePolygonSegments) {
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    auto addSphere = [&](const gp_Pnt& p) {
        builder.Add(compound, BRepPrimAPI_MakeSphere(p, path.wireRadius).Shape());
    };

    for (const auto& pr : path.prims) {
        if (pr.kind == Primitive::ARC) {
            const Arc& a = pr.arc;
            if (a.azSweep < 1e-9) continue;
            // Profile at azStart, revolved about the column axis (+Y through the origin);
            // positive sweep advances the azimuth per azPoint()'s convention.
            gp_Pnt c = azPoint(a.r, a.y, a.azStart);
            gp_Dir tangent(-std::sin(a.azStart), 0.0, -std::cos(a.azStart));
            TopoDS_Face prof = wireProfile(c, tangent, path.wireRadius, wirePolygonSegments);
            BRepPrimAPI_MakeRevol rev(prof, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)),
                                      a.azSweep);
            if (!rev.IsDone() || rev.Shape().IsNull()) {
                throw std::runtime_error("ConductorBuilder: arc revolve failed for " +
                                         path.name + " [" + pr.label + "]");
            }
            builder.Add(compound, rev.Shape());
        } else if (pr.kind == Primitive::SEG) {
            const Seg& s = pr.seg;
            double len = s.a.Distance(s.b);
            if (len < 1e-12) continue;
            gp_Dir dir(gp_XYZ(s.b.XYZ() - s.a.XYZ()));
            builder.Add(compound,
                        BRepPrimAPI_MakeCylinder(gp_Ax2(s.a, dir), path.wireRadius, len).Shape());
            // Sphere joints keep the wire envelope continuous across kinks.
            addSphere(s.a);
            addSphere(s.b);
        } else {
            // Helical ramp/jump: short cylinder chain along the sampled spiral.
            auto pts = samplePrim(pr, path.wireRadius);
            for (size_t j = 0; j + 1 < pts.size(); ++j) {
                double len = pts[j].Distance(pts[j + 1]);
                if (len < 1e-12) continue;
                gp_Dir dir(gp_XYZ(pts[j + 1].XYZ() - pts[j].XYZ()));
                builder.Add(compound,
                            BRepPrimAPI_MakeCylinder(gp_Ax2(pts[j], dir), path.wireRadius, len)
                                .Shape());
                addSphere(pts[j]);
            }
            addSphere(pts.back());
        }
    }
    return compound;
}

struct WindowBox {
    double rInner = 0, rOuter = 0;   // radial extent of the winding window
    double yBottom = 0, yTop = 0;    // axial extent
};

// ---------------------------------------------------------------------------------------
// The planner, templated like buildTurnsImpl so both the MAS and the OpenMagnetics typed
// coils work (their getters return different variant types).
template <typename CoilT, typename WireT>
std::vector<NamedShape> buildAllImpl(const CoilT& coil,
                                     const MAS::CoreBobbinProcessedDescription& bobbinPd,
                                     bool isToroidal,
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

    // Winding window bounds (rectangular cross-section for concentric cores).
    const auto& wws = bobbinPd.get_winding_windows();
    if (wws.empty() || !wws[0].get_width() || !wws[0].get_height()) {
        throw std::runtime_error("ConductorBuilder: bobbin winding window has no width/height");
    }
    double colRadius = bobbinPd.get_column_width().value_or(0.0);
    if (colRadius <= 0.0) {
        throw std::runtime_error("ConductorBuilder: bobbin column width (tube outer radius) "
                                 "is missing or non-positive");
    }
    WindowBox window;
    {
        const auto& ww = wws[0];
        double wwWidth = ww.get_width().value();
        double wwHeight = ww.get_height().value();
        double wwCx = colRadius + wwWidth / 2.0;
        double wwCy = 0.0;
        const auto& wwCoords = ww.get_coordinates();
        if (wwCoords && wwCoords->size() >= 2) {
            wwCx = (*wwCoords)[0];
            wwCy = (*wwCoords)[1];
        }
        window.rInner = wwCx - wwWidth / 2.0;
        window.rOuter = wwCx + wwWidth / 2.0;
        window.yBottom = wwCy - wwHeight / 2.0;
        window.yTop = wwCy + wwHeight / 2.0;
    }

    // Wire per winding (sliced to the MAS base — geometry only needs dimensions).
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
    // order, MKF's contract). Conductor identity order: functionalDescription order, then
    // parallel index.
    struct ConductorTurns {
        std::string winding;
        int64_t parallel;
        std::vector<const MAS::Turn*> turns;
    };
    std::vector<ConductorTurns> conductors;
    std::map<std::string, std::vector<size_t>> byWinding;
    for (const auto& winding : coil.get_functional_description()) {
        for (int64_t k = 0; k < winding.get_number_parallels(); ++k) {
            conductors.push_back({winding.get_name(), k, {}});
            byWinding[winding.get_name()].push_back(conductors.size() - 1);
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

    // --- plan every conductor path ------------------------------------------------------
    const double seamAz = kPi / 2.0;   // -Z transition face (design decision B)
    std::vector<ConductorPath> paths;
    paths.reserve(conductors.size());

    for (const auto& windingEntry : byWinding) {
        const auto& idxs = windingEntry.second;
        const size_t K = idxs.size();

        auto wIt = wireMap.find(windingEntry.first);
        if (wIt == wireMap.end()) {
            throw std::runtime_error("ConductorBuilder: winding '" + windingEntry.first +
                                     "' has no wire in coil.functionalDescription");
        }
        auto [wireW, wireH] = TurnBuilder::wireDimensions(wIt->second, opts.paintCoating);
        double wireRadius = std::min(wireW, wireH) / 2.0;
        double od = 2.0 * wireRadius;

        double rMin = std::numeric_limits<double>::max();
        for (size_t idx : idxs)
            for (const MAS::Turn* t : conductors[idx].turns)
                rMin = std::min(rMin, t->get_coordinates().at(0));

        // Seam sector layout: [exit-stub band | shared crossover WINDOW | entrance-stub
        // band]. Every station's on-station arc has its gap over the ENTIRE window, so all
        // ramps/jumps live inside it as stacked near-parallel chords. Two stacked chords
        // with vertical offset p (the slot pitch) and run w have perpendicular clearance
        // p*w/sqrt(w^2 + (K*p)^2); solving for >= one OD gives the minimum window width.
        // Tight-packed layers (p ~ OD) admit no straight-chord crossover at any width —
        // MKF's real-mode SPREAD alignment avoids that; we throw loudly if it reaches us.
        // Spring model: EVERY move to the next station — intra-layer advance AND layer
        // jump — is a one-wrap SPIRAL whose pitch/radius change is spread progressively
        // over the full 360 degrees, starting exactly at the turn's MKF station (real
        // multi-start spring geometry; parallels are parallel starts that never converge).
        // No crossover window is needed; the seam sector holds only the terminal stub
        // bands: [exit stubs | entrance stubs]. Consecutive wraps touch near the seam by
        // ~pitch*delta/2pi — the physical spring contact — as ADJACENT path primitives
        // (exempt); everything else must genuinely clear or the gate throws.
        const double s = kStubSpacingODs * od / rMin;
        const double delta = 2.0 * (static_cast<double>(K) + 1.0) * s;
        if (delta > 2.0 * kPi / 3.0) {
            throw std::runtime_error(
                "ConductorBuilder: seam sector for winding '" + windingEntry.first +
                "' would span " + std::to_string(delta * 180.0 / kPi) +
                " deg (> 120 deg) — wire too thick relative to the innermost turn radius");
        }
        double sectorStart = seamAz - delta / 2.0;   // conductor arrives here (azimuth grows)
        double sectorEnd = seamAz + delta / 2.0;     // turns start/resume here

        // Lead planning pre-pass. A stub crossing another parallel's slot heights must do
        // so where that parallel's coverage has already ended (exit side) or not yet begun
        // (entrance side), so stubs are ordered by NEARNESS to the exit/entrance edge: the
        // conductor closest to its edge gets the azimuth closest to the sector boundary.
        struct LeadPlan {
            double r0 = 0, y0 = 0, rL = 0, yL = 0;
            bool inBottom = true, outBottom = true;
            size_t inRank = 0, outRank = 0;
        };
        std::vector<LeadPlan> leads(K);
        for (size_t k = 0; k < K; ++k) {
            const auto& ts = conductors[idxs[k]].turns;
            auto& lp = leads[k];
            lp.r0 = ts.front()->get_coordinates().at(0);
            lp.y0 = ts.front()->get_coordinates().at(1);
            lp.rL = ts.back()->get_coordinates().at(0);
            lp.yL = ts.back()->get_coordinates().at(1);
            auto goesBottom = [&](double radius, double yEnd, const char* what) {
                double lo = std::numeric_limits<double>::max();
                double hi = std::numeric_limits<double>::lowest();
                for (const MAS::Turn* t : ts) {
                    if (std::abs(t->get_coordinates().at(0) - radius) <= od / 2.0) {
                        lo = std::min(lo, t->get_coordinates().at(1));
                        hi = std::max(hi, t->get_coordinates().at(1));
                    }
                }
                if (yEnd <= lo + od / 4.0) return true;    // end turn is its lowest slot
                if (yEnd >= hi - od / 4.0) return false;   // end turn is its highest slot
                throw std::runtime_error(
                    std::string("ConductorBuilder: ") + what + " of " +
                    conductors[idxs[k]].winding + " parallel " +
                    std::to_string(conductors[idxs[k]].parallel) + " ends mid-layer (y=" +
                    std::to_string(yEnd) + " between " + std::to_string(lo) + " and " +
                    std::to_string(hi) + ") — no collision-free vertical stub direction exists");
            };
            lp.inBottom = goesBottom(lp.r0, lp.y0, "entrance");
            lp.outBottom = goesBottom(lp.rL, lp.yL, "exit");
        }
        for (size_t k = 0; k < K; ++k) {
            size_t inRank = 0, outRank = 0;
            for (size_t j = 0; j < K; ++j) {
                if (j == k) continue;
                // Nearer to the edge = rank 0. Bottom edge: smaller y is nearer.
                if (leads[j].inBottom == leads[k].inBottom) {
                    bool jNearer = leads[k].inBottom ? (leads[j].y0 < leads[k].y0)
                                                     : (leads[j].y0 > leads[k].y0);
                    if (jNearer) ++inRank;
                }
                if (leads[j].outBottom == leads[k].outBottom) {
                    bool jNearer = leads[k].outBottom ? (leads[j].yL < leads[k].yL)
                                                      : (leads[j].yL > leads[k].yL);
                    if (jNearer) ++outRank;
                }
            }
            leads[k].inRank = inRank;
            leads[k].outRank = outRank;
        }
        // Exit band directly BEFORE the entrance band (both near the sector end): the exit
        // spiral's previous wrap then sits nearly a full pitch above the entrance runs —
        // the widest vertical clearance the spring geometry allows.
        auto entranceStubAz = [&](size_t rank) {
            return sectorEnd - (static_cast<double>(rank) + 1.0) * s;
        };
        auto exitStubAz = [&](size_t rank) {
            return sectorEnd - (static_cast<double>(K + rank) + 1.0) * s;
        };

        for (size_t k = 0; k < K; ++k) {
            const auto& ct = conductors[idxs[k]];
            ConductorPath path;
            path.name = ct.winding + " parallel " + std::to_string(ct.parallel);
            path.wireRadius = wireRadius;

            const auto& turns = ct.turns;
            auto station = [&](const MAS::Turn* t) -> std::pair<double, double> {
                const auto& c = t->get_coordinates();
                if (c.size() < 2) {
                    throw std::runtime_error("ConductorBuilder: turn '" + t->get_name() +
                                             "' has fewer than 2 coordinates");
                }
                return {c[0], c[1]};
            };

            auto [r0, y0] = station(turns.front());
            auto [rL, yL] = station(turns.back());

            // Entrance: radial edge run at the window edge the lead heads toward, then a
            // vertical stub to the first station, then a connector arc to the sector end.
            // MKF's real-winding blocking freed exactly these edge slots.
            {
                bool bottomIn = leads[k].inBottom;
                double yEdge = bottomIn ? window.yBottom + wireRadius : window.yTop - wireRadius;
                // MKF's rule: an end already inside the edge band leaves straight at its own
                // axial level (the blocked edge slots of the crossed layers cover it); only
                // ends deeper in the window take a vertical stub to the edge first.
                double yRun = (std::abs(y0 - yEdge) < od) ? y0 : yEdge;
                double az = entranceStubAz(leads[k].inRank);
                double rBorder = window.rOuter - wireRadius;
                path.prims.push_back({Primitive::SEG, {},
                                      {azPoint(rBorder, yRun, az), azPoint(r0, yRun, az)}, {},
                                      "entrance edge run"});
                if (std::abs(y0 - yRun) > 1e-9) {
                    path.prims.push_back({Primitive::SEG, {},
                                          {azPoint(r0, yRun, az), azPoint(r0, y0, az)}, {},
                                          "entrance stub"});
                }
                if (sectorEnd - az > 1e-9) {
                    path.prims.push_back({Primitive::ARC, {r0, y0, az, sectorEnd - az}, {}, {},
                                          "entrance connector"});
                }
            }

            // Wrap phase: every turn passes its exact MKF station at azimuth `phase`.
            // Intra-layer spirals keep the phase; each layer jump advances it by the
            // crossover bump's arc, exactly like a real winding where the crossover
            // consumes a little azimuth.
            double phase = sectorEnd;
            for (size_t i = 0; i < turns.size(); ++i) {
                auto [r, y] = station(turns[i]);
                if (i + 1 < turns.size()) {
                    auto [rN, yN] = station(turns[i + 1]);
                    bool layerJump = std::abs(rN - r) > od / 2.0;
                    if (!layerJump) {
                        // Intra-layer advance: one-wrap spiral from this turn's exact MKF
                        // station to the next slot (real spring geometry). Where the wrap
                        // crosses another parallel's LAST (flat) turn at the same radius,
                        // it bulges radially outward by one wire OD over a short arc — the
                        // classic multifilar crossover riding over the finished turn.
                        struct Cross { double az; };
                        std::vector<Cross> crossings;
                        for (size_t j = 0; j < K; ++j) {
                            if (j == k) continue;
                            const auto& lj = leads[j];
                            if (std::abs(lj.rL - r) > od / 2.0) continue;
                            // A "crossing" nearer than half a wire to either endpoint is
                            // the junction itself (source/destination slot), not a wire
                            // to ride over.
                            double lo = std::min(y, yN) + od / 2.0;
                            double hi = std::max(y, yN) - od / 2.0;
                            if (lj.yL <= lo || lj.yL >= hi) continue;
                            double frac = (lj.yL - y) / (yN - y);
                            crossings.push_back({phase + kTwoPi * frac});
                        }
                        std::sort(crossings.begin(), crossings.end(),
                                  [](const Cross& a, const Cross& b) { return a.az < b.az; });
                        // Trapezoid bulge profile: full +OD plateau while the wrap is
                        // within one OD of the crossed turn's height (half-width
                        // od/slope), plus short entry/exit ramps — anything narrower
                        // leaves a sub-OD pinch at the bulge foot.
                        double slope = std::abs(yN - y) / kTwoPi;
                        // Plateau half-width: the wrap is within one wire envelope of the
                        // crossed turn while |dy| < od + sampling sag; outside it the
                        // axial gap alone clears, so any ramp width works.
                        double clearDy = od + kMaxSagFraction * 2.0 * wireRadius;
                        double wPlateau = (slope > 1e-12) ? clearDy / slope : kTwoPi;
                        double wRamp = kBulgeRampFraction * wPlateau;
                        auto heightAt = [&](double az) {
                            return y + (yN - y) * (az - phase) / kTwoPi;
                        };
                        double azCur = phase;
                        double rCur = r;
                        size_t part = 0;
                        auto pushSpiral = [&](double azTo, double rTo, const char* what) {
                            azTo = std::min(azTo, phase + kTwoPi);
                            if (azTo - azCur < 1e-9) { rCur = rTo; return; }
                            Primitive pr;
                            pr.kind = Primitive::SPIRAL;
                            pr.spiral = {rCur, heightAt(azCur), azCur,
                                         rTo, heightAt(azTo), azTo};
                            pr.label = "turn '" + turns[i]->get_name() + "' (" + what +
                                       " " + std::to_string(part++) + ")";
                            pr.turnOrdinal = i;
                            path.prims.push_back(std::move(pr));
                            azCur = azTo;
                            rCur = rTo;
                        };
                        for (const auto& cx : crossings) {
                            if (cx.az - wPlateau - wRamp < azCur - 1e-9 ||
                                cx.az + wPlateau + wRamp > phase + kTwoPi + 1e-9) {
                                throw std::runtime_error(
                                    "ConductorBuilder: multifilar crossover bulge for turn '" +
                                    turns[i]->get_name() +
                                    "' does not fit within one wrap (pitch too small "
                                    "relative to the wire OD) — see MKF ABT #187");
                            }
                            pushSpiral(cx.az - wPlateau - wRamp, r, "helical");
                            pushSpiral(cx.az - wPlateau, r + od, "crossover ramp");
                            pushSpiral(cx.az + wPlateau, r + od, "crossover plateau");
                            pushSpiral(cx.az + wPlateau + wRamp, r, "crossover ramp");
                        }
                        pushSpiral(phase + kTwoPi, rN, "helical");
                    } else {
                        // Layer jump: full flat wrap at the exact station, then an
                        // azimuth-local crossover riding over the wrap's own start
                        // (adjacent primitives — the physical crossover contact): first
                        // ALL the radius at constant height, then the height climb at the
                        // full destination radius — one layer clear of every slot it
                        // crosses. The crossover advances the wrap phase.
                        double bumpAng = std::min(kCrossoverArcODs * od / rMin, kMaxCrossoverArc);
                        path.prims.push_back({Primitive::ARC, {r, y, phase, kTwoPi}, {}, {},
                                              "turn '" + turns[i]->get_name() + "'", i});
                        Primitive bumpR;
                        bumpR.kind = Primitive::SPIRAL;
                        bumpR.spiral = {r, y, phase, rN, y, phase + bumpAng / 2.0};
                        bumpR.label = "turn '" + turns[i]->get_name() + "' (crossover radial)";
                        bumpR.turnOrdinal = i;
                        path.prims.push_back(std::move(bumpR));
                        Primitive bumpY;
                        bumpY.kind = Primitive::SPIRAL;
                        bumpY.spiral = {rN, y, phase + bumpAng / 2.0, rN, yN, phase + bumpAng};
                        bumpY.label = "turn '" + turns[i]->get_name() + "' (crossover climb)";
                        bumpY.turnOrdinal = i;
                        path.prims.push_back(std::move(bumpY));
                        phase += bumpAng;
                    }
                } else {
                    // Last turn: one full flat wrap at the exact station. The exit then
                    // continues the spring from the wrap's end (see exit block below).
                    path.prims.push_back({Primitive::ARC, {r, y, phase, kTwoPi}, {}, {},
                                          "turn '" + turns[i]->get_name() + "'", i});
                }
            }

            // Exit: MKF's own convention for an end that crosses nothing radially — leave
            // straight OUTWARD at the last turn's own axial level, exactly where the
            // final wrap ends (parallels' runs stack one slot pitch apart vertically).
            // MKF's terminal blocking reserves the crossed-layer slots when the end is
            // not outermost; the gate verifies.
            {
                double rBorder = window.rOuter - wireRadius;
                path.prims.push_back({Primitive::SEG, {},
                                      {azPoint(rL, yL, phase), azPoint(rBorder, yL, phase)},
                                      {}, "exit run", turns.size() - 1});
            }

            // Window-bounds check on synthesized points (the arcs sit on exact MKF
            // positions — those are MKF's responsibility, not ours to second-guess).
            for (const auto& pr : path.prims) {
                if (pr.kind != Primitive::SEG) continue;
                for (const gp_Pnt& p : {pr.seg.a, pr.seg.b}) {
                    double rr = std::hypot(p.X(), p.Z());
                    if (rr < window.rInner - wireRadius - kContactTol ||
                        rr > window.rOuter + kContactTol ||
                        p.Y() < window.yBottom - kContactTol ||
                        p.Y() > window.yTop + kContactTol) {
                        std::ostringstream s;
                        s << "ConductorBuilder: " << path.name << " [" << pr.label
                          << "] leaves the winding window (r=" << rr << ", y=" << p.Y()
                          << "; window r=[" << window.rInner << "," << window.rOuter
                          << "] y=[" << window.yBottom << "," << window.yTop << "])";
                        throw std::runtime_error(s.str());
                    }
                }
            }

            paths.push_back(std::move(path));
        }
    }

    checkCollisions(paths);

    std::vector<NamedShape> out;
    out.reserve(paths.size());
    for (const auto& p : paths) {
        TopoDS_Shape solid = emitConductor(p, opts.wirePolygonSegments);
        if (TurnBuilder::fuseTurnParts()) {
            solid = TurnBuilder::fuseSolids(solid);
        }
        out.push_back({solid, p.name});
    }
    return out;
}

} // namespace

std::vector<NamedShape> ConductorBuilder::buildAll(
    const MAS::Coil& coil, const MAS::CoreBobbinProcessedDescription& bobbin, bool isToroidal,
    const Options& opts) {
    return buildAllImpl<MAS::Coil, MAS::Wire>(coil, bobbin, isToroidal, opts);
}

std::vector<NamedShape> ConductorBuilder::buildAll(
    const OpenMagnetics::Coil& coil, const MAS::CoreBobbinProcessedDescription& bobbin,
    bool isToroidal, const Options& opts) {
    return buildAllImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(coil, bobbin, isToroidal, opts);
}

} // namespace mvb
