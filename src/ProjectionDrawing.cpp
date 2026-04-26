#include "mvb/ProjectionDrawing.h"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pnt.hxx>

namespace mvb {

namespace {

// Fuse many shapes into a single compound (no boolean fuse — keep them
// separate so edges stay intact). Returns a null shape for an empty input.
TopoDS_Shape makeCompound(const std::vector<TopoDS_Shape>& shapes) {
    if (shapes.empty()) return TopoDS_Shape();
    if (shapes.size() == 1) return shapes.front();

    BRep_Builder b;
    TopoDS_Compound compound;
    b.MakeCompound(compound);
    for (const auto& s : shapes) {
        if (!s.IsNull()) b.Add(compound, s);
    }
    return compound;
}

// Build a huge axis-aligned half-space box on the "far" side of the plane
// — i.e. the half past the section plane, as viewed from the plane normal.
// Intersecting this with a shape keeps only the portion BEHIND the section,
// which is what the user sees in a cutaway technical drawing.
TopoDS_Shape farHalfBox(const TopoDS_Shape& shape,
                        SectionPlane plane,
                        double section_offset) {
    Bnd_Box bbox;
    BRepBndLib::Add(shape, bbox);
    if (bbox.IsVoid()) return TopoDS_Shape();

    double xmin, ymin, zmin, xmax, ymax, zmax;
    bbox.Get(xmin, ymin, zmin, xmax, ymax, zmax);

    // Inflate slightly so the cutter is safely larger than the shape.
    const double pad = std::max({xmax - xmin, ymax - ymin, zmax - zmin, 1.0}) * 0.1 + 1.0;
    xmin -= pad; ymin -= pad; zmin -= pad;
    xmax += pad; ymax += pad; zmax += pad;

    // "Far" half = region where the plane-normal coordinate is >= offset,
    // because the viewer stands at -infinity along that normal.
    switch (plane) {
        case SectionPlane::XY: {                  // normal = +Z, viewer at +Z far
            // far half = z <= section_offset (past the plane away from viewer)
            return BRepPrimAPI_MakeBox(
                gp_Pnt(xmin, ymin, zmin),
                gp_Pnt(xmax, ymax, section_offset)).Shape();
        }
        case SectionPlane::XZ: {                  // normal = +Y
            return BRepPrimAPI_MakeBox(
                gp_Pnt(xmin, ymin, zmin),
                gp_Pnt(xmax, section_offset, zmax)).Shape();
        }
        case SectionPlane::YZ: {                  // normal = +X
            return BRepPrimAPI_MakeBox(
                gp_Pnt(xmin, ymin, zmin),
                gp_Pnt(section_offset, ymax, zmax)).Shape();
        }
    }
    return TopoDS_Shape();
}

} // namespace

std::string ProjectionDrawing::drawProjection(
    const TopoDS_Shape& shape,
    SectionPlane plane,
    double width_px,
    double margin_px,
    double stroke_width,
    const std::string& stroke_color)
{
    // SectionBuilder::edgesToSvg samples every TopAbs_EDGE of the input
    // and projects via the same project() mapping we want here, so we can
    // re-use it directly.
    return SectionBuilder::edgesToSvg(shape, plane, width_px, margin_px,
                                       stroke_width, stroke_color);
}

std::string ProjectionDrawing::drawProjection(
    const std::vector<TopoDS_Shape>& shapes,
    SectionPlane plane,
    double width_px,
    double margin_px,
    double stroke_width,
    const std::string& stroke_color)
{
    return drawProjection(makeCompound(shapes), plane, width_px, margin_px,
                           stroke_width, stroke_color);
}

std::string ProjectionDrawing::drawCutawaySection(
    const TopoDS_Shape& shape,
    SectionPlane plane,
    double section_offset,
    double width_px,
    double margin_px,
    double stroke_width,
    const std::string& stroke_color)
{
    if (shape.IsNull()) {
        throw std::runtime_error("drawCutawaySection: empty shape");
    }

    TopoDS_Shape cutter = farHalfBox(shape, plane, section_offset);
    if (cutter.IsNull()) {
        // Fallback: just draw the full projection.
        return drawProjection(shape, plane, width_px, margin_px,
                               stroke_width, stroke_color);
    }

    BRepAlgoAPI_Common cut(shape, cutter);
    cut.Build();
    TopoDS_Shape far = cut.IsDone() ? cut.Shape() : shape;

    return drawProjection(far, plane, width_px, margin_px,
                           stroke_width, stroke_color);
}

std::string ProjectionDrawing::drawCutawaySection(
    const std::vector<TopoDS_Shape>& shapes,
    SectionPlane plane,
    double section_offset,
    double width_px,
    double margin_px,
    double stroke_width,
    const std::string& stroke_color)
{
    // Cut each shape individually so boolean COMMON stays cheap, then
    // assemble a single compound for projection.
    std::vector<TopoDS_Shape> cut;
    cut.reserve(shapes.size());
    for (const auto& s : shapes) {
        if (s.IsNull()) continue;
        TopoDS_Shape cutter = farHalfBox(s, plane, section_offset);
        if (cutter.IsNull()) { cut.push_back(s); continue; }
        BRepAlgoAPI_Common c(s, cutter);
        c.Build();
        cut.push_back(c.IsDone() ? c.Shape() : s);
    }
    return drawProjection(cut, plane, width_px, margin_px,
                           stroke_width, stroke_color);
}

} // namespace mvb
