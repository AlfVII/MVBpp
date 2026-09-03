#include "mvb/shapes/ShapeDrumSemishielded.h"
#include "mvb/Utils.h"
#include <numbers>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakePrism.hxx>

namespace mvb {
namespace shapes {

TopoDS_Shape buildSemishieldedShell(const std::map<std::string, double>& dims,
                                    const TopoDS_Shape& drum,
                                    int /*corePolygonSegments*/) {
    double j = dims.count("J") ? dims.at("J") : 0.0;
    double k = dims.count("K") ? dims.at("K") : 0.0;
    double l = dims.count("L") ? dims.at("L") : 0.0;
    double b = dims.count("B") ? dims.at("B") : 0.0;
    if (j <= 0.0 || k <= 0.0 || l <= 0.0) return TopoDS_Shape();

    // Same build frame as ShapeDrum: extruded along Z from 0 to B, so the envelope is centred
    // on the drum's own height rather than on the origin.
    double base = (b - l) / 2.0;

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-j / 2.0, -k / 2.0, base));
    poly.Add(gp_Pnt( j / 2.0, -k / 2.0, base));
    poly.Add(gp_Pnt( j / 2.0,  k / 2.0, base));
    poly.Add(gp_Pnt(-j / 2.0,  k / 2.0, base));
    poly.Close();
    TopoDS_Face face = BRepBuilderAPI_MakeFace(poly.Wire()).Face();
    TopoDS_Shape envelope = BRepPrimAPI_MakePrism(face, gp_Vec(0.0, 0.0, l)).Shape();

    // The drum leaves ShapeBuilder::buildPiece already rotated -90deg about X (the intrinsic
    // rotation that puts a concentric core's axis on Y). The envelope is built in the raw
    // extrusion frame, so it must take the SAME rotation or the shell lands 90deg out of
    // alignment with the drum it encases — which is exactly what the first build produced:
    // drum extent Y=1.65 (its height) against shell extent Z=1.8 (its height).
    if (drum.IsNull()) {
        return rotate_shape(envelope, -std::numbers::pi / 2.0, 0.0, 0.0);
    }

    // The shell is what the epoxy actually occupies: the envelope with the drum carved out.
    // Subtracting rather than returning the solid block matters for rendering — a translucent
    // block would show the drum through a second layer of epoxy that is not there — and it is
    // the honest volume if anyone ever meshes it.
    // Rotate BEFORE subtracting: the drum is already in the rotated frame, so the cut has to
    // happen there or it subtracts a solid that is not where the drum actually is.
    TopoDS_Shape orientedEnvelope = rotate_shape(envelope, -std::numbers::pi / 2.0, 0.0, 0.0);
    BRepAlgoAPI_Cut shell(orientedEnvelope, drum);
    return shell.IsDone() ? shell.Shape() : orientedEnvelope;
}

} // namespace shapes
} // namespace mvb
