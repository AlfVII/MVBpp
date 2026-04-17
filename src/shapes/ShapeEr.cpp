#include "mvb/shapes/ShapeEr.h"
#include "mvb/Utils.h"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Trsf.hxx>
#include <gp_Pnt.hxx>
#include <cmath>

namespace mvb {
namespace shapes {

TopoDS_Shape ShapeEr::buildWindingWindow(const std::map<std::string, double>& dims) const {
    double b = 0.0, c = 0.0, d = 0.0, e = 0.0, f = 0.0, g = 0.0;
    auto it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;
    it = dims.find("G"); if (it != dims.end()) g = it->second;

    if (b == 0.0 || d == 0.0 || e == 0.0) return TopoDS_Shape();

    double zCenter = b - d / 2.0;

    TopoDS_Shape outerCyl = build_polygon_cylinder(d, e / 2.0, 0);
    outerCyl = translate_shape(outerCyl, 0.0, 0.0, zCenter - d / 2.0);

    TopoDS_Shape innerCyl = build_polygon_cylinder(d, f / 2.0, 0);
    innerCyl = translate_shape(innerCyl, 0.0, 0.0, zCenter - d / 2.0);

    BRepAlgoAPI_Cut ringCut(outerCyl, innerCyl);
    if (!ringCut.IsDone()) return TopoDS_Shape();
    TopoDS_Shape ww = ringCut.Shape();

    // ER/EQ variant: when G > F and C > F, add a rectangular slot of
    // (G × C × D) centred on the column, minus the central column cylinder,
    // to the winding window cut. Matches MVB Python Er.get_negative_winding_window.
    if (g > f && c > f) {
        gp_Pnt boxCorner(-g / 2.0, -c / 2.0, zCenter - d / 2.0);
        TopoDS_Shape cube = BRepPrimAPI_MakeBox(boxCorner, g, c, d).Shape();
        BRepAlgoAPI_Cut cubeMinus(cube, innerCyl);
        if (cubeMinus.IsDone()) {
            BRepAlgoAPI_Fuse fuser(ww, cubeMinus.Shape());
            if (fuser.IsDone()) ww = fuser.Shape();
        }
    }
    return ww;
}

TopoDS_Shape ShapeEr::applyMachining(const TopoDS_Shape& piece,
                                     const MAS::Machining& machining,
                                     const std::map<std::string, double>& dims) const {
    const std::vector<double>& coords = machining.get_coordinates();
    if (coords.size() < 2) return piece;

    double gapLength = machining.get_length();
    double xCoord = coords[0];
    double yCoord = coords[1];

    if (std::abs(xCoord) < 1e-12) {
        double f = 0.0;
        auto it = dims.find("F");
        if (it != dims.end()) f = it->second;
        if (f <= 0.0) return piece;

        // Center column: cylinder along Y axis, CENTERED on yCoord to match
        // MVB Python apply_machining (which uses workplane XZ polygon_cylinder
        // → centered on workplane origin → translate → centered on y_coord).
        gp_Ax2 cylAxis(gp_Pnt(0, yCoord - gapLength / 2.0, 0), gp_Dir(0, 1, 0));
        TopoDS_Shape tool = BRepPrimAPI_MakeCylinder(cylAxis, f / 2.0, gapLength).Shape();
        if (tool.IsNull()) return piece;

        BRepAlgoAPI_Cut cutter(piece, tool);
        return cutter.IsDone() ? cutter.Shape() : piece;
    }

    // Side columns: use generic rectangular tool from base class
    return ShapeE::applyMachining(piece, machining, dims);
}

} // namespace shapes
} // namespace mvb
