#include "mvb/shapes/ShapeDrum.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <algorithm>

namespace mvb {
namespace shapes {

// Build frame (pre intrinsic -90degX rotation): extruded along Z from 0 to B, bottom
// flange at the base. The groove is cut as an annular band, exactly the ShapeP idiom.

TopoDS_Face ShapeDrum::buildProfile(const std::map<std::string, double>& dims) const {
    double flange = dims.count("A") ? dims.at("A") : 0.0;
    if (dims.count("A2") && dims.at("A2") > 0.0) {
        flange = std::max(flange, dims.at("A2"));
    }
    TopoDS_Wire outer = build_polygon_circle(flange / 2.0, m_corePolygonSegments);
    BRepBuilderAPI_MakeFace face(outer);
    return face.Face();
}

TopoDS_Shape ShapeDrum::buildWindingWindow(const std::map<std::string, double>& dims) const {
    double a = dims.count("A") ? dims.at("A") : 0.0;
    double c = dims.count("C") ? dims.at("C") : 0.0;
    double e = dims.count("E") ? dims.at("E") : 0.0;
    double f = dims.count("F") ? dims.at("F") : 0.0;
    if (a == 0.0 || c == 0.0 || e == 0.0) return TopoDS_Shape();

    // Annular groove: from the post surface out past the flange edge, spanning z in [F, F+E].
    TopoDS_Shape outerCyl = build_polygon_cylinder(e, a, m_corePolygonSegments);  // radius a: safely past any flange
    outerCyl = translate_shape(outerCyl, 0.0, 0.0, f);
    TopoDS_Shape post = build_polygon_cylinder(e, c / 2.0, m_corePolygonSegments);
    post = translate_shape(post, 0.0, 0.0, f);

    BRepAlgoAPI_Cut cutter(outerCyl, post);
    return cutter.IsDone() ? cutter.Shape() : TopoDS_Shape();
}

TopoDS_Shape ShapeDrum::applyExtras(const std::map<std::string, double>& dims,
                                    const TopoDS_Shape& piece) const {
    double a = dims.count("A") ? dims.at("A") : 0.0;
    double a2 = (dims.count("A2") && dims.at("A2") > 0.0) ? dims.at("A2") : a;
    double b = dims.count("B") ? dims.at("B") : 0.0;
    double d = dims.count("D") ? dims.at("D") : 0.0;
    double f = dims.count("F") ? dims.at("F") : 0.0;
    double h = dims.count("H") ? dims.at("H") : 0.0;

    TopoDS_Shape result = piece;

    // Asymmetric flanges: the profile used the LARGER disc; trim the other flange's band
    // down to its own diameter. Convention: A = bottom flange, A2 = top flange.
    double larger = std::max(a, a2);
    if (a2 < larger || a < larger) {
        double smaller = std::min(a, a2);
        double bandBase = (a2 < a) ? (b - d) : 0.0;   // trim the top band if A2 is smaller, else the bottom
        double bandHeight = (a2 < a) ? d : f;
        TopoDS_Shape trimOuter = build_polygon_cylinder(bandHeight, larger, m_corePolygonSegments);
        trimOuter = translate_shape(trimOuter, 0.0, 0.0, bandBase);
        TopoDS_Shape trimInner = build_polygon_cylinder(bandHeight, smaller / 2.0, m_corePolygonSegments);
        trimInner = translate_shape(trimInner, 0.0, 0.0, bandBase);
        BRepAlgoAPI_Cut ring(trimOuter, trimInner);
        if (ring.IsDone()) {
            BRepAlgoAPI_Cut trimmed(result, ring.Shape());
            if (trimmed.IsDone()) result = trimmed.Shape();
        }
    }

    // Bore through the whole piece.
    if (h > 0.0) {
        TopoDS_Shape bore =
            build_polygon_cylinder(b, h / 2.0, m_corePolygonSegments, /*circumscribed=*/true);
        BRepAlgoAPI_Cut cutter(result, bore);
        if (cutter.IsDone()) result = cutter.Shape();
    }

    // Centre the single piece on the origin (open-shape cores are emitted as one CLOSED
    // element at {0,0,z}); after the intrinsic rotation the drum axis lands on Y.
    return translate_shape(result, 0.0, 0.0, -b / 2.0);
}

} // namespace shapes
} // namespace mvb
