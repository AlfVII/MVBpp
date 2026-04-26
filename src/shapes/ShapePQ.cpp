#include "mvb/shapes/ShapePQ.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <gp_Pnt.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <cmath>

namespace mvb {
namespace shapes {

static TopoDS_Shape buildLateralPiece(double a, double c, double e, double sin_g, double cos_g,
                                       double j, double l, double height) {
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(j / 4.0, l / 4.0, 0.0));
    poly.Add(gp_Pnt(j / 4.0, -l / 4.0, 0.0));
    poly.Add(gp_Pnt(j / 2.0, -l / 2.0, 0.0));
    // When e*cos_g > c the arc extremum lies outside the flat-edge boundary,
    // producing a self-intersecting polygon.  Only include it when it's inside.
    if (e * cos_g < c) {
        poly.Add(gp_Pnt(e * sin_g, -e * cos_g, 0.0));
    }
    poly.Add(gp_Pnt(e * sin_g, -c, 0.0));
    poly.Add(gp_Pnt(a, -c, 0.0));
    poly.Add(gp_Pnt(a, c, 0.0));
    poly.Add(gp_Pnt(e * sin_g, c, 0.0));
    if (e * cos_g < c) {
        poly.Add(gp_Pnt(e * sin_g, e * cos_g, 0.0));
    }
    poly.Add(gp_Pnt(j / 2.0, l / 2.0, 0.0));
    poly.Close();

    BRepBuilderAPI_MakeWire wire(poly.Wire());
    BRepBuilderAPI_MakeFace face(wire.Wire());
    gp_Vec vec(0, 0, height);
    return BRepPrimAPI_MakePrism(face.Face(), vec).Shape();
}

TopoDS_Shape ShapePQ::buildPiece(const MAS::CoreShape& shapeData) const {
    auto dimsOpt = shapeData.get_dimensions();
    if (!dimsOpt) return TopoDS_Shape();
    auto dims = flatten_dimensions(*dimsOpt);

    double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0, f = 0.0, g = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second / 2.0;
    it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second / 2.0;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second / 2.0;
    it = dims.find("F"); if (it != dims.end()) f = it->second / 2.0;
    it = dims.find("G"); if (it != dims.end()) g = it->second;

    if (b == 0.0) return TopoDS_Shape();

    double l = (f * 2.0) + (c * 2.0 - f * 2.0) / 3.0;
    auto itL = dims.find("L");
    if (itL != dims.end()) l = itL->second;

    double j = f;
    auto itJ = dims.find("J");
    if (itJ != dims.end()) j = itJ->second;

    double g_angle = 0.0;
    if (e > 0.0) {
        if (g > 0.0 && g <= e * 2.0) {
            g_angle = std::asin(g / (e * 2.0));
        } else {
            g_angle = std::asin((e + f) / (e * 2.0));
        }
    }
    double sin_g = std::sin(g_angle);
    double cos_g = std::cos(g_angle);

    // Central structure: two cylinders stacked.
    // 1) Connector disk (radius E/2, height B−D): solid yoke base that joins the
    //    column to the lateral wings — this is the "solid cylinder joining the two
    //    parts" of the PQ (column side + wing side).
    // 2) Real column (radius F/2, height D): the cylinder that protrudes into the
    //    winding space.  The winding-window cut will remove the E/2→F/2 annulus at
    //    this height, leaving only this column.
    double yoke_h = b - d;
    TopoDS_Shape central;
    {
        TopoDS_Shape connector = build_polygon_cylinder(
            (yoke_h > 1e-9 ? yoke_h : b), f, m_corePolygonSegments);
        TopoDS_Shape column    = build_polygon_cylinder(d, f, m_corePolygonSegments);
        // Translate column to sit on top of the connector disk.
        if (yoke_h > 1e-9) {
            gp_Trsf shift; shift.SetTranslation(gp_Vec(0, 0, yoke_h));
            column = BRepBuilderAPI_Transform(column, shift).Shape();
        }
        BRepAlgoAPI_Fuse fconn(connector, column);
        central = (yoke_h > 1e-9 && fconn.IsDone()) ? fconn.Shape() : connector;
    }

    // Right lateral piece
    TopoDS_Shape rightLateral = buildLateralPiece(a, c, e, sin_g, cos_g, j, l, b);

    // Left lateral piece (mirror across YZ plane)
    gp_Trsf mirror;
    mirror.SetMirror(gp_Ax2(gp_Pnt(0, 0, b / 2.0), gp_Dir(1, 0, 0)));
    TopoDS_Shape leftLateral = BRepBuilderAPI_Transform(rightLateral, mirror).Shape();

    // Fuse all three
    BRepAlgoAPI_Fuse f1(central, rightLateral);
    if (!f1.IsDone()) return central;
    TopoDS_Shape fused = f1.Shape();

    BRepAlgoAPI_Fuse f2(fused, leftLateral);
    if (f2.IsDone()) fused = f2.Shape();

    // Apply winding window cut
    TopoDS_Shape window = buildWindingWindow(dims);
    if (!window.IsNull()) {
        BRepAlgoAPI_Cut cutter(fused, window);
        if (cutter.IsDone()) fused = cutter.Shape();
    }

    // Apply extras (translate by -B)
    fused = applyExtras(dims, fused);

    // Intrinsic rotation
    fused = rotate_shape(fused, -std::numbers::pi / 2.0, 0.0, 0.0);
    return fused;
}

} // namespace shapes
} // namespace mvb
