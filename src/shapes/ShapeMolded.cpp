#include "mvb/shapes/ShapeMolded.h"
#include "mvb/Utils.h"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>

namespace mvb {
namespace shapes {

// Build frame (pre intrinsic -90degX rotation): the block is extruded along Z from 0 to B by
// ShapeBuilder::buildPiece, so Z is the coil axis here and the footprint lies in XY.
TopoDS_Face ShapeMolded::buildProfile(const std::map<std::string, double>& dims) const {
    double a = dims.count("A") ? dims.at("A") / 2.0 : 0.0;
    double c = dims.count("C") ? dims.at("C") / 2.0 : 0.0;
    if (a <= 0.0 || c <= 0.0) return TopoDS_Face();

    // A rectangular body footprint: A wide, C deep, centred on the axis. Moulded inductors
    // are plain blocks — the shaping is all internal.
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a, -c, 0.0));
    poly.Add(gp_Pnt( a, -c, 0.0));
    poly.Add(gp_Pnt( a,  c, 0.0));
    poly.Add(gp_Pnt(-a,  c, 0.0));
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

TopoDS_Shape ShapeMolded::buildWindingWindow(const std::map<std::string, double>& dims) const {
    double b = dims.count("B") ? dims.at("B") : 0.0;
    double d = dims.count("D") ? dims.at("D") : 0.0;   // cavity height
    double e = dims.count("E") ? dims.at("E") : 0.0;   // cavity outer diameter
    double f = dims.count("F") ? dims.at("F") : 0.0;   // cavity inner diameter (the post)
    if (b <= 0.0 || d <= 0.0 || e <= 0.0) return TopoDS_Shape();

    // The cavity is the annulus the winding occupies: outer diameter E, inner diameter F,
    // height D, centred on the body so the composite closes over and under the coil. F may
    // legitimately be 0 (no post), in which case the cavity is a plain cylinder.
    //
    // Geometry validity (E > F, D < B, E inside A x C) is enforced upstream by MKF's
    // CorePieceMolded, so an out-of-range record cannot reach here from a validated MAS.
    // Drawing the bare block would be a visible defect rather than silent physics, so an
    // incomplete record degrades to the solid body instead of aborting the render.
    if (e <= f || d >= b) {
        return TopoDS_Shape();
    }

    double base = (b - d) / 2.0;
    // Circumscribe the cavity: faceting must never leave a sliver of composite inside the
    // volume the winding occupies, or the turns intersect the core in the render.
    TopoDS_Shape cavity = build_polygon_cylinder(d, e / 2.0, m_corePolygonSegments,
                                                 /*circumscribed=*/true);
    cavity = translate_shape(cavity, 0.0, 0.0, base);
    if (f <= 0.0) {
        return cavity;
    }

    // The post is inscribed for the mirrored reason: it is composite that must not eat into
    // the winding volume.
    TopoDS_Shape post = build_polygon_cylinder(d, f / 2.0, m_corePolygonSegments);
    post = translate_shape(post, 0.0, 0.0, base);

    BRepAlgoAPI_Cut annulus(cavity, post);
    return annulus.IsDone() ? annulus.Shape() : TopoDS_Shape();
}

TopoDS_Shape ShapeMolded::applyExtras(const std::map<std::string, double>& dims,
                                      const TopoDS_Shape& piece) const {
    // Centre the body on the origin, the same closing move ShapeDrum makes.
    //
    // ShapeBuilder::buildPiece extrudes the profile from z = 0 to B, so a piece that does not
    // translate is left sitting entirely on one side of the plane. MKF places the winding
    // window at y = 0 (CorePieceMolded sets the window coordinates with no vertical offset),
    // so an uncentred body renders half a body-height away from its own coil — the block
    // floating above the winding rather than moulded around it.
    //
    // Every other family here ends with this translation; ShapeMolded originally had no
    // applyExtras at all, which is how it was missed.
    double b = dims.count("B") ? dims.at("B") : 0.0;
    if (b <= 0.0) {
        return piece;
    }
    return translate_shape(piece, 0.0, 0.0, -b / 2.0);
}

} // namespace shapes
} // namespace mvb
