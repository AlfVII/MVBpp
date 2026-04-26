#pragma once

#include "mvb/SectionBuilder.h"
#include <TopoDS_Shape.hxx>
#include <string>
#include <vector>

namespace mvb {

// 2D projections of 3D geometry rendered as SVG. Two modes:
//
//  - drawProjection: project every edge of the input shape(s) onto `plane`
//    and render the result as a wire-frame SVG. Shows the full outline
//    "as seen looking from the plane normal" — all edges, no hidden-line
//    removal (matches cadquery getSVG conventions).
//
//  - drawCutawaySection: cut the shape(s) with a half-space at the given
//    signed offset along the plane normal, keep the half that is FURTHER
//    from the viewer (i.e. past the section plane), and project that onto
//    `plane`. This is the classic "section view" where the near half of
//    the solid has been sliced away, so the section cut is visible AND
//    everything beyond the section plane is visible behind it.
class ProjectionDrawing {
public:
    // Single-shape projection.
    static std::string drawProjection(
        const TopoDS_Shape& shape,
        SectionPlane plane,
        double width_px    = 800.0,
        double margin_px   = 40.0,
        double stroke_width = 0.5,
        const std::string& stroke_color = "#000000");

    // Multi-shape projection (e.g. assembly: core + bobbin + turns).
    static std::string drawProjection(
        const std::vector<TopoDS_Shape>& shapes,
        SectionPlane plane,
        double width_px    = 800.0,
        double margin_px   = 40.0,
        double stroke_width = 0.5,
        const std::string& stroke_color = "#000000");

    // Cutaway section: cuts away the near half (at plane + offset along +normal),
    // projects the far half onto `plane`.
    static std::string drawCutawaySection(
        const TopoDS_Shape& shape,
        SectionPlane plane,
        double section_offset = 0.0,
        double width_px    = 800.0,
        double margin_px   = 40.0,
        double stroke_width = 0.5,
        const std::string& stroke_color = "#000000");

    static std::string drawCutawaySection(
        const std::vector<TopoDS_Shape>& shapes,
        SectionPlane plane,
        double section_offset = 0.0,
        double width_px    = 800.0,
        double margin_px   = 40.0,
        double stroke_width = 0.5,
        const std::string& stroke_color = "#000000");
};

} // namespace mvb
