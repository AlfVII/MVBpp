#include "mvb/shapes/ShapeBuilder.h"
#include "mvb/Utils.h"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <cmath>

namespace mvb {
namespace shapes {

TopoDS_Shape ShapeBuilder::buildPiece(const MAS::CoreShape& shapeData) const {
    auto dimsOpt = shapeData.get_dimensions();
    if (!dimsOpt) return TopoDS_Shape();
    auto dims = flatten_dimensions(*dimsOpt);

    // Build base profile and extrude
    TopoDS_Face profile = buildProfile(dims);
    double height = dims.count("B") ? dims.at("B") : 0.0;
    if (isToroidal()) {
        height = dims.count("C") ? dims.at("C") : height;
    }
    TopoDS_Shape piece = extrude(profile, height);

    // Cut winding window
    TopoDS_Shape window = buildWindingWindow(dims);
    if (!window.IsNull()) {
        piece = BRepAlgoAPI_Cut(piece, window).Shape();
    }

    // Apply extras
    piece = applyExtras(dims, piece);

    // Intrinsic rotation: concentric cores are built in XY-extruded-Z then rotated -90° around X
    if (!isToroidal()) {
        piece = rotate_shape(piece, -std::numbers::pi / 2.0, 0.0, 0.0);
    }

    return piece;
}

TopoDS_Shape ShapeBuilder::applyMachining(const TopoDS_Shape& piece,
                                          const MAS::Machining& machining,
                                          const std::map<std::string, double>& dims) const {
    // Generic rectangular machining tool oriented along Y (column axis)
    const std::vector<double>& coords = machining.get_coordinates();
    if (coords.size() < 2) return piece;

    double gapLength = machining.get_length();
    double xCoord = coords[0];
    double yCoord = coords[1];

    TopoDS_Shape tool;
    double f = dims.count("F") ? dims.at("F") : 0.0;
    double c = dims.count("C") ? dims.at("C") : 0.0;
    double a = dims.count("A") ? dims.at("A") : 0.0;
    // For EFD cores the central column depth is F2 (not C). Use F2 when
    // available so the gap tool matches MVB Python's EFD.apply_machining.
    double centerColumnDepth = dims.count("F2") ? dims.at("F2") : c;

    if (std::abs(xCoord) < 1e-12) {
        // Center column: F (width) × gap × centerColumnDepth
        tool = makeBox(f, gapLength, centerColumnDepth);
    } else {
        // Side column: A/2 × gap × C (or A for EFD-family, matching MVB Python)
        double sideDepth = dims.count("F2") ? a : c;
        double sideX = (xCoord < 0) ? -a / 4.0 : a / 4.0;
        tool = makeBox(a / 2.0, gapLength, sideDepth);
        tool = translate_shape(tool, sideX, 0.0, 0.0);

        // Avoid cutting into center column (with matching depth)
        TopoDS_Shape centerTool = makeBox(f * 1.001, gapLength, centerColumnDepth * 1.001);
        BRepAlgoAPI_Cut cutter(tool, centerTool);
        if (cutter.IsDone()) {
            tool = cutter.Shape();
        }
    }

    tool = translate_shape(tool, 0.0, yCoord, 0.0);
    BRepAlgoAPI_Cut cutter(piece, tool);
    return cutter.IsDone() ? cutter.Shape() : piece;
}

TopoDS_Shape ShapeBuilder::buildWindingWindow(const std::map<std::string, double>&) const {
    return TopoDS_Shape();
}

TopoDS_Shape ShapeBuilder::applyExtras(const std::map<std::string, double>&,
                                        const TopoDS_Shape& piece) const {
    return piece;
}

TopoDS_Shape ShapeBuilder::extrude(const TopoDS_Face& face, double height) {
    gp_Vec vec(0, 0, height);
    return BRepPrimAPI_MakePrism(face, vec).Shape();
}

TopoDS_Shape ShapeBuilder::makeBox(double x, double y, double z) {
    gp_Pnt corner(-x / 2.0, -y / 2.0, -z / 2.0);
    return BRepPrimAPI_MakeBox(corner, x, y, z).Shape();
}

TopoDS_Shape ShapeBuilder::makeCylinder(double height, double radius, int segments) {
    return build_polygon_cylinder(height, radius, segments);
}

} // namespace shapes
} // namespace mvb
