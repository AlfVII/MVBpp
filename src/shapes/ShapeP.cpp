#include "mvb/shapes/ShapeP.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

TopoDS_Face ShapeP::buildProfile(const std::map<std::string, double>& dims) const {
    double a = 0.0;
    auto it = dims.find("A");
    if (it != dims.end()) a = it->second / 2.0;

    TopoDS_Wire outer = build_polygon_circle(a, DEFAULT_CORE_POLYGON_SEGMENTS);
    BRepBuilderAPI_MakeFace face(outer);
    return face.Face();
}

TopoDS_Shape ShapeP::buildWindingWindow(const std::map<std::string, double>& dims) const {
    double b = 0.0, d = 0.0, e = 0.0, f = 0.0;
    auto it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;

    if (b == 0.0 || d == 0.0 || e == 0.0) return TopoDS_Shape();

    double zCenter = b - d / 2.0;

    // build_polygon_cylinder creates a base-at-zero prism; shift down by d/2 to center it at zCenter
    TopoDS_Shape outerCyl = build_polygon_cylinder(d, e / 2.0, DEFAULT_CORE_POLYGON_SEGMENTS);
    outerCyl = translate_shape(outerCyl, 0.0, 0.0, zCenter - d / 2.0);

    TopoDS_Shape innerCyl = build_polygon_cylinder(d, f / 2.0, DEFAULT_CORE_POLYGON_SEGMENTS);
    innerCyl = translate_shape(innerCyl, 0.0, 0.0, zCenter - d / 2.0);

    BRepAlgoAPI_Cut cutter(outerCyl, innerCyl);
    return cutter.IsDone() ? cutter.Shape() : TopoDS_Shape();
}

TopoDS_Shape ShapeP::applyExtras(const std::map<std::string, double>& dims,
                                 const TopoDS_Shape& piece) const {
    double b = 0.0;
    auto it = dims.find("B");
    if (it != dims.end()) b = it->second;
    return translate_shape(piece, 0.0, 0.0, -b);
}

// Round-column gap cut: cylinder of radius F/2, height = gap length,
// centered at machining.coordinates (matches MVB.js ERShape.applyMachining).
TopoDS_Shape ShapeP::applyMachining(const TopoDS_Shape& piece,
                                    const MAS::Machining& machining,
                                    const std::map<std::string, double>& dims) const {
    const auto& coords = machining.get_coordinates();
    if (coords.size() < 2) return piece;
    double gapLength = machining.get_length();
    if (gapLength <= 0.0) return piece;

    double f = 0.0;
    auto it = dims.find("F"); if (it != dims.end()) f = it->second;
    if (f <= 0.0) return piece;

    double xCoord = coords[0];
    double yCoord = coords[1];
    if (std::abs(xCoord) > 1e-12) {
        // Side-column gaps for P-family aren't physically meaningful (no
        // outer column to cut). Skip rather than attempting a degenerate cut.
        return piece;
    }

    // Centre-column gap: cylinder along Y (the column axis after the
    // intrinsic -90°X rotation), radius F/2, height gapLength, centred at
    // Y = yCoord. Using default BRepPrimAPI_MakeCylinder orients along Z
    // and cuts a transverse slot through the column — same mistake Python
    // avoids by its coord system and ShapeEr handles via an explicit axis.
    gp_Ax2 cylAxis(gp_Pnt(0.0, yCoord - gapLength / 2.0, 0.0), gp_Dir(0, 1, 0));
    TopoDS_Shape tool = BRepPrimAPI_MakeCylinder(cylAxis, f / 2.0, gapLength).Shape();
    BRepAlgoAPI_Cut cutter(piece, tool);
    return cutter.IsDone() ? cutter.Shape() : piece;
}

} // namespace shapes
} // namespace mvb
