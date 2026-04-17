#include "mvb/SectionBuilder.h"

#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRep_Builder.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <GCPnts_TangentialDeflection.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pln.hxx>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace mvb {

namespace {

gp_Pln planeFor(SectionPlane plane) {
    switch (plane) {
        case SectionPlane::XY: return gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
        case SectionPlane::XZ: return gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0));
        case SectionPlane::YZ: return gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0));
    }
    return gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
}

// Project a 3D point onto the chosen plane's 2D coordinates.
// Convention: SVG +x = right, +y = down, so we flip the vertical axis.
std::pair<double, double> project(const gp_Pnt& p, SectionPlane plane) {
    switch (plane) {
        case SectionPlane::XY: return {p.X(), p.Y()};
        case SectionPlane::XZ: return {p.X(), p.Z()};
        case SectionPlane::YZ: return {p.Y(), p.Z()};
    }
    return {p.X(), p.Y()};
}

// Sample one edge into a polyline using OCCT's tangential deflection
// adaptive sampler — handles lines, circles, b-splines, etc.
std::vector<gp_Pnt> sampleEdge(const TopoDS_Edge& edge,
                                double angularDeflection = 0.05,
                                double curvatureDeflection = 1e-4) {
    std::vector<gp_Pnt> pts;
    BRepAdaptor_Curve curve(edge);
    GCPnts_TangentialDeflection sampler(curve, angularDeflection,
                                          curvatureDeflection);
    pts.reserve(sampler.NbPoints());
    for (int i = 1; i <= sampler.NbPoints(); ++i) {
        pts.push_back(sampler.Value(i));
    }
    return pts;
}

} // namespace

TopoDS_Shape SectionBuilder::sectionCore(const TopoDS_Shape& solid,
                                          SectionPlane plane) {
    if (solid.IsNull()) return TopoDS_Shape();
    BRepAlgoAPI_Section section(solid, planeFor(plane), Standard_False);
    section.ComputePCurveOn1(Standard_True);
    section.Approximation(Standard_True);
    section.Build();
    if (!section.IsDone()) return TopoDS_Shape();
    return section.Shape();
}

std::string SectionBuilder::edgesToSvg(const TopoDS_Shape& edges,
                                        SectionPlane plane,
                                        double width_px,
                                        double margin_px,
                                        double stroke_width,
                                        const std::string& stroke_color) {
    if (edges.IsNull()) {
        throw std::runtime_error("edgesToSvg: empty section");
    }

    // Sample every edge to a polyline and collect bbox in the projected plane.
    std::vector<std::vector<std::pair<double, double>>> polylines;
    double xmin = +1e30, xmax = -1e30, ymin = +1e30, ymax = -1e30;
    bool any = false;

    for (TopExp_Explorer exp(edges, TopAbs_EDGE); exp.More(); exp.Next()) {
        TopoDS_Edge edge = TopoDS::Edge(exp.Current());
        std::vector<gp_Pnt> pts;
        try { pts = sampleEdge(edge); }
        catch (...) { continue; }
        if (pts.size() < 2) continue;

        std::vector<std::pair<double, double>> poly;
        poly.reserve(pts.size());
        for (const auto& p : pts) {
            auto [u, v] = project(p, plane);
            poly.emplace_back(u, v);
            if (u < xmin) xmin = u; if (u > xmax) xmax = u;
            if (v < ymin) ymin = v; if (v > ymax) ymax = v;
            any = true;
        }
        polylines.push_back(std::move(poly));
    }

    if (!any) {
        throw std::runtime_error("edgesToSvg: section has no sampled edges");
    }

    // Scale-to-fit, with the SVG y-axis flipped (model +y is up, SVG +y is down).
    double model_w = xmax - xmin;
    double model_h = ymax - ymin;
    if (model_w <= 0 || model_h <= 0) {
        throw std::runtime_error("edgesToSvg: degenerate section bbox");
    }
    double inner_w = width_px - 2.0 * margin_px;
    double scale = inner_w / model_w;
    double height_px = model_h * scale + 2.0 * margin_px;

    auto toSvg = [&](double u, double v) -> std::pair<double, double> {
        double x = (u - xmin) * scale + margin_px;
        double y = height_px - ((v - ymin) * scale + margin_px);
        return {x, y};
    };

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(3);
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "width=\"" << width_px << "\" height=\"" << height_px << "\" "
        << "viewBox=\"0 0 " << width_px << " " << height_px << "\">\n";
    svg << "  <g fill=\"none\" stroke=\"" << stroke_color
        << "\" stroke-width=\"" << stroke_width
        << "\" stroke-linecap=\"round\" stroke-linejoin=\"round\">\n";

    for (const auto& poly : polylines) {
        if (poly.size() < 2) continue;
        svg << "    <path d=\"";
        for (size_t i = 0; i < poly.size(); ++i) {
            auto [x, y] = toSvg(poly[i].first, poly[i].second);
            svg << (i == 0 ? "M" : "L") << x << "," << y << " ";
        }
        svg << "\"/>\n";
    }

    svg << "  </g>\n";
    svg << "</svg>\n";
    return svg.str();
}

void SectionBuilder::writeSectionSvg(const TopoDS_Shape& solid,
                                      SectionPlane plane,
                                      const std::string& outputPath) {
    auto edges = sectionCore(solid, plane);
    auto svg = edgesToSvg(edges, plane);
    std::ofstream f(outputPath);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open " + outputPath + " for writing");
    }
    f << svg;
}

} // namespace mvb
