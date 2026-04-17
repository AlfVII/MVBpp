#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GC_MakeCircle.hxx>
#include <gp_Circ.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Ax1.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <Standard_Failure.hxx>
#include <StdFail_NotDone.hxx>
#include <iostream>
#include <cmath>

namespace mvb {

static double get_wire_diameter(const MAS::Wire& wire) {
    auto cd = wire.get_conducting_diameter();
    if (cd) return cd->get_nominal().value_or(0.0);
    auto od = wire.get_outer_diameter();
    if (od) return od->get_nominal().value_or(0.0);
    return 0.001;
}

static std::pair<double, double> get_wire_dimensions(const MAS::Wire& wire, const MAS::Turn& turn) {
    auto td = turn.get_dimensions();
    if (td && td->size() >= 2) {
        return {(*td)[0], (*td)[1]};
    }
    auto ow = wire.get_outer_width();
    auto oh = wire.get_outer_height();
    if (ow && oh) {
        return {ow->get_nominal().value_or(0.001), oh->get_nominal().value_or(0.001)};
    }
    auto cw = wire.get_conducting_width();
    auto ch = wire.get_conducting_height();
    if (cw && ch) {
        return {cw->get_nominal().value_or(0.001), ch->get_nominal().value_or(0.001)};
    }
    double d = get_wire_diameter(wire);
    return {d, d};
}

static bool is_rectangular_wire(const MAS::Wire& wire) {
    return wire.get_type() == MAS::WireType::RECTANGULAR || wire.get_type() == MAS::WireType::PLANAR;
}

static TopoDS_Face build_circle_profile(const gp_Ax2& plane, double radius, int segments = 16) {
    gp_Pnt origin = plane.Location();
    gp_Dir dx = plane.XDirection();
    gp_Dir dy = plane.YDirection();

    if (segments <= 0) {
        gp_Circ circ(plane, radius);
        BRepBuilderAPI_MakeEdge edge(circ);
        BRepBuilderAPI_MakeWire wire(edge);
        BRepBuilderAPI_MakeFace face(wire.Wire());
        return face.Face();
    }

    BRepBuilderAPI_MakePolygon poly;
    double offset = M_PI / segments;
    for (int i = 0; i < segments; ++i) {
        double angle = 2.0 * M_PI * i / segments + offset;
        gp_Vec vx(dx);
        gp_Vec vy(dy);
        vx.Scale(radius * std::cos(angle));
        vy.Scale(radius * std::sin(angle));
        gp_Pnt p = origin.Translated(vx).Translated(vy);
        poly.Add(p);
    }
    poly.Close();
    BRepBuilderAPI_MakeFace face(BRepBuilderAPI_MakeWire(poly.Wire()).Wire());
    return face.Face();
}

static TopoDS_Face build_rect_profile(const gp_Ax2& plane, double width, double height) {
    gp_Pnt origin = plane.Location();
    gp_Dir dx = plane.XDirection();
    gp_Dir dy = plane.YDirection();
    double hw = width / 2.0;
    double hh = height / 2.0;

    gp_Pnt p1 = origin.Translated(gp_Vec(dx).Scaled(-hw)).Translated(gp_Vec(dy).Scaled(-hh));
    gp_Pnt p2 = origin.Translated(gp_Vec(dx).Scaled(-hw)).Translated(gp_Vec(dy).Scaled( hh));
    gp_Pnt p3 = origin.Translated(gp_Vec(dx).Scaled( hw)).Translated(gp_Vec(dy).Scaled( hh));
    gp_Pnt p4 = origin.Translated(gp_Vec(dx).Scaled( hw)).Translated(gp_Vec(dy).Scaled(-hh));

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(p1);
    poly.Add(p2);
    poly.Add(p3);
    poly.Add(p4);
    poly.Close();
    BRepBuilderAPI_MakeFace face(poly.Wire());
    return face.Face();
}

static TopoDS_Shape build_concentric_round_column_turn(double turn_radius, double wire_radius,
                                                         double height_pos, bool rect_wire,
                                                         double wire_width, double wire_height) {
    if (rect_wire) {
        // Rectangular wire on round column: revolve a rectangle face 360° around Y axis.
        // This matches MVB.js behavior and avoids BRepOffsetAPI_MakePipe degeneracies.
        double hw = wire_width / 2.0;
        double hh = wire_height / 2.0;
        gp_Pnt p1(turn_radius - hw, height_pos - hh, 0.0);
        gp_Pnt p2(turn_radius + hw, height_pos - hh, 0.0);
        gp_Pnt p3(turn_radius + hw, height_pos + hh, 0.0);
        gp_Pnt p4(turn_radius - hw, height_pos + hh, 0.0);

        BRepBuilderAPI_MakePolygon poly;
        poly.Add(p1);
        poly.Add(p2);
        poly.Add(p3);
        poly.Add(p4);
        poly.Close();

        BRepBuilderAPI_MakeFace faceMaker(poly.Wire());
        if (!faceMaker.IsDone()) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeFace failed for rectangle profile\n";
            return TopoDS_Shape();
        }
        TopoDS_Face profile = faceMaker.Face();

        gp_Ax1 rev_axis(gp_Pnt(0.0, height_pos, 0.0), gp_Dir(0, 1, 0));
        BRepPrimAPI_MakeRevol revol(profile, rev_axis, 2.0 * M_PI - 1e-6);
        if (!revol.IsDone()) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol failed for rect wire turn\n";
            return TopoDS_Shape();
        }
        try {
            return revol.Shape();
        } catch (const Standard_Failure& e) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol threw Standard_Failure: " << e.GetMessageString() << "\n";
            return TopoDS_Shape();
        } catch (...) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol threw unknown exception\n";
            return TopoDS_Shape();
        }
    }

    // Round wire on round column: revolve a polygonal circular profile
    // around the Y axis, using DEFAULT_WIRE_POLYGON_SEGMENTS. When the
    // constant is <=0 we fall back to an exact torus (matches MVB Python's
    // WIRE_POLYGON_SEGMENTS<=0 branch) — the faceted form is preferable
    // for downstream meshers (gmsh) that don't handle the torus
    // parameterisation's seam cleanly.
    const int segments = DEFAULT_WIRE_POLYGON_SEGMENTS;
    if (segments <= 0) {
        gp_Pnt torus_center(0.0, height_pos, 0.0);
        gp_Ax2 torus_axis(torus_center, gp_Dir(0, 1, 0), gp_Dir(1, 0, 0));
        try {
            return BRepPrimAPI_MakeTorus(torus_axis, turn_radius, wire_radius).Shape();
        } catch (const Standard_Failure& e) {
            std::cerr << "ERROR build_concentric_round_column_turn: MakeTorus threw Standard_Failure: " << e.GetMessageString() << "\n";
            return TopoDS_Shape();
        }
    }

    gp_Ax2 profile_plane(gp_Pnt(turn_radius, height_pos, 0.0),
                          gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    TopoDS_Face profile = build_circle_profile(profile_plane, wire_radius, segments);
    gp_Ax1 rev_axis(gp_Pnt(0.0, height_pos, 0.0), gp_Dir(0, 1, 0));
    BRepPrimAPI_MakeRevol revol(profile, rev_axis, 2.0 * M_PI - 1e-6);
    if (!revol.IsDone()) {
        std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol failed for round wire turn\n";
        return TopoDS_Shape();
    }
    try {
        return revol.Shape();
    } catch (const Standard_Failure& e) {
        std::cerr << "ERROR build_concentric_round_column_turn: MakeRevol threw Standard_Failure: " << e.GetMessageString() << "\n";
        return TopoDS_Shape();
    }
}

static TopoDS_Shape build_concentric_rect_column_turn(double radial_pos, double wire_radius,
                                                       double height_pos,
                                                       double half_col_width, double half_col_depth,
                                                       bool rect_wire,
                                                       double wire_width, double wire_height) {
    double turn_turn_radius = radial_pos - half_col_width;
    double min_bend = std::max(wire_width, wire_height) / 2.0 * 1.02;
    if (turn_turn_radius < min_bend) {
        turn_turn_radius = min_bend;
    }

    double wire_x_pos = half_col_width + turn_turn_radius;
    double wire_z_pos = half_col_depth + turn_turn_radius;
    double y = height_pos;

    gp_Pnt pts[8] = {
        gp_Pnt(+half_col_width, y, +wire_z_pos),
        gp_Pnt(-half_col_width, y, +wire_z_pos),
        gp_Pnt(-wire_x_pos,     y, +half_col_depth),
        gp_Pnt(-wire_x_pos,     y, -half_col_depth),
        gp_Pnt(-half_col_width, y, -wire_z_pos),
        gp_Pnt(+half_col_width, y, -wire_z_pos),
        gp_Pnt(+wire_x_pos,     y, -half_col_depth),
        gp_Pnt(+wire_x_pos,     y, +half_col_depth),
    };

    std::pair<gp_Pnt, gp_Dir> corners[4] = {
        {gp_Pnt(-half_col_width, y, +half_col_depth), gp_Dir(0, 0, 1)},
        {gp_Pnt(-half_col_width, y, -half_col_depth), gp_Dir(-1, 0, 0)},
        {gp_Pnt(+half_col_width, y, -half_col_depth), gp_Dir(0, 0, -1)},
        {gp_Pnt(+half_col_width, y, +half_col_depth), gp_Dir(1, 0, 0)},
    };

    // For round wire, build the turn manually from straight cylinders and corner quarter-tori.
    // This avoids BRepOffsetAPI_MakePipe failures on closed paths in OCCT 7.9.
    if (!rect_wire) {
        BRep_Builder builder;
        TopoDS_Compound compound;
        builder.MakeCompound(compound);

        // 4 straight segments: extrude a circle along each straight edge
        int straight_pairs[4][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}};
        for (int i = 0; i < 4; ++i) {
            gp_Pnt p1 = pts[straight_pairs[i][0]];
            gp_Pnt p2 = pts[straight_pairs[i][1]];
            gp_Vec vec(p1, p2);
            double len = vec.Magnitude();
            if (len < 1e-12) continue;

            // Cylinder along the segment, centered at p1, axis = vec
            gp_Ax2 cylAx2(p1, gp_Dir(vec));
            TopoDS_Shape cyl = BRepPrimAPI_MakeCylinder(cylAx2, wire_radius, len).Shape();
            builder.Add(compound, cyl);
        }

        // 4 corner quarter-tori: revolve a circle around the Y axis at each corner center
        int corner_incoming[4] = {0, 2, 4, 6}; // start index of straight segment before each corner
        for (int i = 0; i < 4; ++i) {
            auto [c_center, c_xref] = corners[i];
            gp_Pnt circ_center = c_center.Translated(gp_Vec(c_xref).Scaled(turn_turn_radius));
            gp_Dir incoming_dir = gp_Dir(gp_Vec(pts[corner_incoming[i]], pts[corner_incoming[i] + 1]));
            gp_Ax2 circ_ax2(circ_center, incoming_dir, gp_Dir(0, 1, 0));
            gp_Circ circ(circ_ax2, wire_radius);
            TopoDS_Edge circ_edge = BRepBuilderAPI_MakeEdge(circ).Edge();
            TopoDS_Wire circ_wire = BRepBuilderAPI_MakeWire(circ_edge).Wire();
            TopoDS_Face circ_face = BRepBuilderAPI_MakeFace(circ_wire).Face();

            gp_Ax1 rev_axis(c_center, gp_Dir(0, 1, 0));
            // All corners sweep -90° around +Y (clockwise when viewed from +Y)
            BRepPrimAPI_MakeRevol revol(circ_face, rev_axis, -M_PI / 2.0);
            if (revol.IsDone()) {
                builder.Add(compound, revol.Shape());
            }
        }

        // Fuse all solids into one
        TopoDS_Shape result;
        for (TopExp_Explorer exp(compound, TopAbs_SOLID); exp.More(); exp.Next()) {
            TopoDS_Shape solid = exp.Current();
            if (result.IsNull()) {
                result = solid;
            } else {
                BRepAlgoAPI_Fuse fuse(result, solid);
                result = fuse.Shape();
            }
        }
        return result.IsNull() ? compound : result;
    }

    // Rectangular wire: build the turn from 4 axis-aligned boxes + 4 corner
    // quarter-swept rectangles. Mirrors the MVB.js approach in
    // replicadBuilder.js::_makeQuarterSweptRectangle and avoids BRepOffsetAPI_MakePipe
    // entirely (it segfaults on closed planar paths in OCCT 7.9).
    double hw_wire = wire_width / 2.0;
    double hh_wire = wire_height / 2.0;

    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    // Straight segment 0: along -X at z=+wire_z_pos
    {
        gp_Pnt corner(-half_col_width, y - hh_wire, wire_z_pos - hw_wire);
        TopoDS_Shape box = BRepPrimAPI_MakeBox(corner, 2.0 * half_col_width, wire_height, wire_width).Shape();
        builder.Add(compound, box);
    }
    // Straight segment 1: along -Z at x=-wire_x_pos
    {
        gp_Pnt corner(-wire_x_pos - hw_wire, y - hh_wire, -half_col_depth);
        TopoDS_Shape box = BRepPrimAPI_MakeBox(corner, wire_width, wire_height, 2.0 * half_col_depth).Shape();
        builder.Add(compound, box);
    }
    // Straight segment 2: along +X at z=-wire_z_pos
    {
        gp_Pnt corner(-half_col_width, y - hh_wire, -wire_z_pos - hw_wire);
        TopoDS_Shape box = BRepPrimAPI_MakeBox(corner, 2.0 * half_col_width, wire_height, wire_width).Shape();
        builder.Add(compound, box);
    }
    // Straight segment 3: along +Z at x=+wire_x_pos
    {
        gp_Pnt corner(wire_x_pos - hw_wire, y - hh_wire, -half_col_depth);
        TopoDS_Shape box = BRepPrimAPI_MakeBox(corner, wire_width, wire_height, 2.0 * half_col_depth).Shape();
        builder.Add(compound, box);
    }

    // 4 corner quarter-swept rectangles: revolve a rect face -90° around Y
    // centered at each corner. The profile sits at circ_center = c_center + c_xref * bend_radius,
    // oriented with normal = incoming_dir, local X = c_xref (radial), local Y = +Y (column axis).
    int corner_incoming[4] = {0, 2, 4, 6};
    for (int i = 0; i < 4; ++i) {
        auto [c_center, c_xref] = corners[i];
        gp_Pnt circ_center = c_center.Translated(gp_Vec(c_xref).Scaled(turn_turn_radius));
        gp_Dir incoming_dir = gp_Dir(gp_Vec(pts[corner_incoming[i]], pts[corner_incoming[i] + 1]));

        gp_Ax2 profPlane(circ_center, incoming_dir, c_xref);
        TopoDS_Face rect_face = build_rect_profile(profPlane, wire_width, wire_height);

        gp_Ax1 rev_axis(c_center, gp_Dir(0, 1, 0));
        BRepPrimAPI_MakeRevol revol(rect_face, rev_axis, -M_PI / 2.0);
        if (revol.IsDone()) {
            builder.Add(compound, revol.Shape());
        } else {
            std::cerr << "ERROR build_concentric_rect_column_turn: quarter-swept rect "
                      << i << " MakeRevol failed\n";
        }
    }

    // Fuse all solids into one
    TopoDS_Shape result;
    for (TopExp_Explorer exp(compound, TopAbs_SOLID); exp.More(); exp.Next()) {
        TopoDS_Shape solid = exp.Current();
        if (result.IsNull()) {
            result = solid;
        } else {
            BRepAlgoAPI_Fuse fuse(result, solid);
            if (fuse.IsDone()) {
                result = fuse.Shape();
            }
        }
    }
    return result.IsNull() ? compound : result;
}

// Helpers for toroidal turn construction (mirrors MVB.js _createToroidalTurn).
static TopoDS_Shape rotate_about_axis(const TopoDS_Shape& s, const gp_Pnt& pt,
                                       const gp_Dir& dir, double rad) {
    if (s.IsNull() || std::abs(rad) < 1e-12) return s;
    gp_Trsf t;
    t.SetRotation(gp_Ax1(pt, dir), rad);
    return BRepBuilderAPI_Transform(s, t).Shape();
}

static TopoDS_Shape make_quarter_torus_at(double major_R, double minor_r,
                                            const gp_Pnt& center, double startAngleDeg) {
    double a = startAngleDeg * M_PI / 180.0;
    gp_Ax2 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1),
              gp_Dir(std::cos(a), std::sin(a), 0.0));
    try {
        // Match MVB.js _makeQuarterSweptRectangle pattern but with a circular
        // wire cross-section: disk profile in the XZ plane at (major_R, 0, 0)
        // (so the revolution axis Z lies IN the face's plane), revolve around
        // Z by 90°, then rotate by startAngle, then translate to center.
        gp_Ax2 circ_ax(gp_Pnt(major_R, 0.0, 0.0), gp_Dir(0, 1, 0));
        gp_Circ circ(circ_ax, minor_r);
        TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circ).Edge();
        TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge).Wire();
        TopoDS_Face face = BRepBuilderAPI_MakeFace(wire, Standard_True).Face();
        gp_Ax1 rev_axis(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
        BRepPrimAPI_MakeRevol mk(face, rev_axis, M_PI / 2.0);
        if (!mk.IsDone()) return TopoDS_Shape();
        TopoDS_Shape shape = mk.Shape();
        if (std::abs(a) > 1e-12) {
            gp_Trsf rot;
            rot.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), a);
            shape = BRepBuilderAPI_Transform(shape, rot).Shape();
        }
        return translate_shape(shape, center.X(), center.Y(), center.Z());
    } catch (const Standard_Failure& e) {
        std::cerr << "ERROR make_quarter_torus_at: " << e.GetMessageString() << "\n";
        return TopoDS_Shape();
    }
}

static TopoDS_Shape build_toroidal_turn(const MAS::Turn& turn, const MAS::Wire& wire,
                                         const MAS::CoreBobbinProcessedDescription& bobbin) {
    bool rect_wire = is_rectangular_wire(wire);
    auto [wire_w, wire_h] = get_wire_dimensions(wire, turn);
    double wire_radius = std::min(wire_w, wire_h) / 2.0;

    const auto& coords = turn.get_coordinates();
    if (coords.size() < 2) return TopoDS_Shape();
    double cx = coords[0], cy = coords[1];
    double innerRadial = std::sqrt(cx * cx + cy * cy);
    double innerAngleDeg = std::atan2(cy, cx) * 180.0 / M_PI;

    double outerRadial = innerRadial;
    double outerAngleDeg = innerAngleDeg;
    auto addCoords = turn.get_additional_coordinates();
    double wwRadialHeight = 0.0;
    const auto& wws = bobbin.get_winding_windows();
    if (!wws.empty()) wwRadialHeight = wws[0].get_radial_height().value_or(0.0);

    if (addCoords && !addCoords->empty() && (*addCoords)[0].size() >= 2) {
        double ax = (*addCoords)[0][0], ay = (*addCoords)[0][1];
        outerRadial = std::sqrt(ax * ax + ay * ay);
        outerAngleDeg = std::atan2(ay, ax) * 180.0 / M_PI;
    } else {
        outerRadial = innerRadial + wwRadialHeight;
    }

    double turnAngleDeg = turn.get_rotation().value_or(0.0);
    double angleDiffRad = (outerAngleDeg - innerAngleDeg) * M_PI / 180.0;
    double turnRotationRad = (turnAngleDeg - 180.0) * M_PI / 180.0;

    double halfDepth = bobbin.get_column_depth();
    double bendRadius = rect_wire ? std::max(wire_w, wire_h) / 2.0 : wire_radius;
    double baseClearance = wire_radius;
    double layerClearance = wwRadialHeight - innerRadial;
    double radialHeight = halfDepth + baseClearance + layerClearance;

    double tubeLength = std::max(1e-7, radialHeight - bendRadius);
    double radialDistance = outerRadial - innerRadial;
    double radialLength = std::max(1e-7, radialDistance - 2.0 * bendRadius);
    double innerX = -innerRadial;
    double outerX = innerX - radialDistance;

    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    // Round wire only — rectangular toroidal wires not yet supported.
    if (rect_wire) {
        std::cerr << "WARN build_toroidal_turn: rectangular wire toroidal not implemented, using round\n";
    }

    auto build_half = [&](double ySign) {
        // 1. Inner tube along Y from origin
        {
            TopoDS_Shape c = BRepPrimAPI_MakeCylinder(wire_radius, tubeLength).Shape();
            c = rotate_about_axis(c, gp_Pnt(0,0,0), gp_Dir(1,0,0), -M_PI/2.0 * ySign);
            c = translate_shape(c, innerX, 0.0, 0.0);
            if (!c.IsNull()) builder.Add(compound, c);
        }
        // 2. Inner corner
        {
            double sa = (ySign > 0) ? 0.0 : 270.0;
            gp_Pnt center(innerX - bendRadius, tubeLength * ySign, 0.0);
            TopoDS_Shape s = make_quarter_torus_at(bendRadius, wire_radius, center, sa);
            s = rotate_about_axis(s, gp_Pnt(innerX, 0, 0), gp_Dir(0,1,0), angleDiffRad / 2.0);
            if (!s.IsNull()) builder.Add(compound, s);
        }
        // 3. Radial segment along X
        {
            TopoDS_Shape c = BRepPrimAPI_MakeCylinder(wire_radius, radialLength).Shape();
            c = rotate_about_axis(c, gp_Pnt(0,0,0), gp_Dir(0,1,0), -M_PI/2.0);
            c = rotate_about_axis(c, gp_Pnt(0,0,0), gp_Dir(0,1,0), angleDiffRad / 2.0);
            c = translate_shape(c, innerX - bendRadius, radialHeight * ySign, 0.0);
            c = rotate_about_axis(c, gp_Pnt(innerX, 0, 0), gp_Dir(0,1,0), angleDiffRad / 2.0);
            if (!c.IsNull()) builder.Add(compound, c);
        }
        // 4. Outer corner
        {
            double sa = (ySign > 0) ? 90.0 : 180.0;
            gp_Pnt center(outerX + bendRadius, tubeLength * ySign, 0.0);
            TopoDS_Shape s = make_quarter_torus_at(bendRadius, wire_radius, center, sa);
            s = rotate_about_axis(s, gp_Pnt(innerX, 0, 0), gp_Dir(0,1,0), angleDiffRad);
            if (!s.IsNull()) builder.Add(compound, s);
        }
        // 5. Outer tube
        {
            TopoDS_Shape c = BRepPrimAPI_MakeCylinder(wire_radius, tubeLength).Shape();
            c = rotate_about_axis(c, gp_Pnt(0,0,0), gp_Dir(1,0,0), -M_PI/2.0 * ySign);
            c = translate_shape(c, outerX, 0.0, 0.0);
            c = rotate_about_axis(c, gp_Pnt(innerX, 0, 0), gp_Dir(0,1,0), angleDiffRad);
            if (!c.IsNull()) builder.Add(compound, c);
        }
    };

    build_half(+1.0);
    build_half(-1.0);

    // Match MVB.js: keep all pieces as a Compound (no fuse). Disjoint
    // solids do not survive BRepAlgoAPI_Fuse cleanly across multiple calls.
    TopoDS_Shape result = compound;
    result = rotate_about_axis(result, gp_Pnt(0,0,0), gp_Dir(0,1,0), turnRotationRad);
    return result;
}

static TopoDS_Shape build_concentric_oblong_turn(double radial_pos, double wire_radius,
                                                  double height_pos,
                                                  double half_col_width, double half_col_depth,
                                                  bool rect_wire,
                                                  double wire_width, double wire_height) {
    double straight_half = half_col_depth - half_col_width;
    if (straight_half <= 0.0) {
        // Effectively round column
        return build_concentric_round_column_turn(radial_pos, wire_radius, height_pos,
                                                   rect_wire, wire_width, wire_height);
    }

    double tube_z_length = 2.0 * straight_half;
    double prof_r = rect_wire ? std::min(wire_width, wire_height) / 2.0 : wire_radius;

    // Helper to create an exact cylinder along Z at the same placement as before.
    auto make_z_cylinder = [&](double cx, double cz) -> TopoDS_Shape {
        gp_Ax2 ax(gp_Pnt(cx, height_pos, cz), gp_Dir(0, 0, 1));
        return BRepPrimAPI_MakeCylinder(ax, prof_r, tube_z_length).Shape();
    };

    TopoDS_Shape tube_px = make_z_cylinder(radial_pos, 0.0);
    TopoDS_Shape tube_nx = make_z_cylinder(-radial_pos, 0.0);

    // Half-torus arcs at ±Z ends
    auto make_half_torus = [&](double cz, double start_angle_deg) -> TopoDS_Shape {
        gp_Pnt circ_center(radial_pos, height_pos, cz);
        double angle_rad = start_angle_deg * M_PI / 180.0;
        gp_Dir normal(std::cos(angle_rad), 0.0, std::sin(angle_rad));
        gp_Ax2 circ_axis(circ_center, normal);

        Handle(Geom_Circle) circle = GC_MakeCircle(circ_axis, prof_r).Value();
        TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circle).Edge();
        TopoDS_Wire circ_wire = BRepBuilderAPI_MakeWire(edge).Wire();
        TopoDS_Face circ_face = BRepBuilderAPI_MakeFace(circ_wire).Face();

        gp_Ax1 rev_axis(gp_Pnt(0, height_pos, cz), gp_Dir(0, 0, 1));
        BRepPrimAPI_MakeRevol revol(circ_face, rev_axis, M_PI);
        try {
            return revol.Shape();
        } catch (const Standard_Failure& e) {
            std::cerr << "ERROR build_concentric_oblong_turn: MakeRevol threw Standard_Failure: " << e.GetMessageString() << "\n";
            return TopoDS_Shape();
        } catch (...) {
            std::cerr << "ERROR build_concentric_oblong_turn: MakeRevol threw unknown exception\n";
            return TopoDS_Shape();
        }
    };

    TopoDS_Shape torus_pz = make_half_torus(straight_half, 90.0);
    TopoDS_Shape torus_nz = make_half_torus(-straight_half, -90.0);

    // Fuse all pieces
    BRepAlgoAPI_Fuse f1(tube_px, tube_nx);
    if (!f1.IsDone()) return tube_px;
    TopoDS_Shape fused = f1.Shape();

    BRepAlgoAPI_Fuse f2(fused, torus_pz);
    if (f2.IsDone()) fused = f2.Shape();

    BRepAlgoAPI_Fuse f3(fused, torus_nz);
    if (f3.IsDone()) fused = f3.Shape();

    return fused;
}

TopoDS_Shape TurnBuilder::buildTurn(const MAS::Turn& turn,
                                    const MAS::Wire& wire,
                                    const MAS::CoreBobbinProcessedDescription& bobbin,
                                    bool isToroidal) {
    const auto& coords = turn.get_coordinates();
    double radial_pos = coords.size() > 0 ? coords[0] : 0.0;
    double height_pos = coords.size() > 1 ? coords[1] : 0.0;

    bool rect_wire = is_rectangular_wire(wire);
    auto [wire_w, wire_h] = get_wire_dimensions(wire, turn);
    double wire_radius = std::min(wire_w, wire_h) / 2.0;

    if (isToroidal) {
        return build_toroidal_turn(turn, wire, bobbin);
    }

    double half_col_width = bobbin.get_column_width().value_or(0.0);
    double half_col_depth = bobbin.get_column_depth();

    if (bobbin.get_column_shape() == MAS::ColumnShape::ROUND) {
        return build_concentric_round_column_turn(radial_pos, wire_radius, height_pos,
                                                   rect_wire, wire_w, wire_h);
    }

    if (bobbin.get_column_shape() == MAS::ColumnShape::OBLONG) {
        return build_concentric_oblong_turn(radial_pos, wire_radius, height_pos,
                                             half_col_width, half_col_depth,
                                             rect_wire, wire_w, wire_h);
    }

    // Default: rectangular column
    return build_concentric_rect_column_turn(radial_pos, wire_radius, height_pos,
                                              half_col_width, half_col_depth,
                                              rect_wire, wire_w, wire_h);
}

} // namespace mvb
