#include "mvb/Symmetry.h"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <cmath>

namespace mvb {

const char* to_string(SymmetryPlane p) {
    switch (p) {
        case SymmetryPlane::None: return "None";
        case SymmetryPlane::XY:   return "XY";
        case SymmetryPlane::YZ:   return "YZ";
        case SymmetryPlane::XZ:   return "XZ";
    }
    return "?";
}

void plane_axes(SymmetryPlane p, int& i_normal, int& i_in_a, int& i_in_b) {
    switch (p) {
        case SymmetryPlane::XY: i_normal = 2; i_in_a = 0; i_in_b = 1; return;
        case SymmetryPlane::YZ: i_normal = 0; i_in_a = 1; i_in_b = 2; return;
        case SymmetryPlane::XZ: i_normal = 1; i_in_a = 0; i_in_b = 2; return;
        case SymmetryPlane::None:
        default:                i_normal = -1; i_in_a = -1; i_in_b = -1;
    }
}

namespace {

double shape_volume(const TopoDS_Shape& s) {
    if (s.IsNull()) return 0.0;
    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props);
    return std::max(props.Mass(), 0.0);
}

// Build a solid axis-aligned box occupying one half of the bounding box.
TopoDS_Shape make_halfspace(SymmetryPlane plane, SymmetryHalf half,
                            const ShapeBBox& b) {
    const double dx = b.xmax - b.xmin;
    const double dy = b.ymax - b.ymin;
    const double dz = b.zmax - b.zmin;
    const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
    const double pad = std::max(diag * 5.0, 1e-3);
    const double cx = 0.5 * (b.xmin + b.xmax);
    const double cy = 0.5 * (b.ymin + b.ymax);
    const double cz = 0.5 * (b.zmin + b.zmax);

    double lx_min = cx - pad, lx_max = cx + pad;
    double ly_min = cy - pad, ly_max = cy + pad;
    double lz_min = cz - pad, lz_max = cz + pad;

    switch (plane) {
        case SymmetryPlane::YZ:
            if (half == SymmetryHalf::Positive) lx_min = 0.0;
            else                                 lx_max = 0.0;
            break;
        case SymmetryPlane::XZ:
            if (half == SymmetryHalf::Positive) ly_min = 0.0;
            else                                 ly_max = 0.0;
            break;
        case SymmetryPlane::XY:
            if (half == SymmetryHalf::Positive) lz_min = 0.0;
            else                                 lz_max = 0.0;
            break;
        case SymmetryPlane::None:
            break;
    }

    BRepPrimAPI_MakeBox mk(gp_Pnt(lx_min, ly_min, lz_min),
                           gp_Pnt(lx_max, ly_max, lz_max));
    return mk.Solid();
}

double volume_on_side(const std::vector<NamedShape>& shapes,
                      const TopoDS_Shape& halfspace) {
    double v = 0.0;
    for (const auto& ns : shapes) {
        if (ns.shape.IsNull()) continue;
        BRepAlgoAPI_Common op(ns.shape, halfspace);
        op.Build();
        if (!op.IsDone()) continue;
        v += shape_volume(op.Shape());
    }
    return v;
}

// Enumerate all solid sub-shapes of a (possibly compound) result. Returns
// each as a separate TopoDS_Shape so consumers can name them individually.
std::vector<TopoDS_Shape> flatten_solids(const TopoDS_Shape& s) {
    std::vector<TopoDS_Shape> out;
    if (s.IsNull()) return out;
    if (s.ShapeType() == TopAbs_SOLID) {
        out.push_back(s);
        return out;
    }
    for (TopExp_Explorer ex(s, TopAbs_SOLID); ex.More(); ex.Next()) {
        out.push_back(ex.Current());
    }
    return out;
}

} // namespace

ShapeBBox aggregate_bbox(const std::vector<NamedShape>& shapes) {
    Bnd_Box box;
    for (const auto& ns : shapes) {
        if (!ns.shape.IsNull()) BRepBndLib::Add(ns.shape, box);
    }
    ShapeBBox out{0, 0, 0, 0, 0, 0};
    if (!box.IsVoid()) {
        box.Get(out.xmin, out.ymin, out.zmin, out.xmax, out.ymax, out.zmax);
    }
    return out;
}

SymmetryResult analyze_symmetry(const std::vector<NamedShape>& shapes,
                                double volume_tolerance) {
    SymmetryResult out;

    double V_total = 0.0;
    for (const auto& ns : shapes) V_total += shape_volume(ns.shape);
    if (V_total <= 0.0) return out;

    const ShapeBBox bb = aggregate_bbox(shapes);

    const std::array<SymmetryPlane, 3> planes = {
        SymmetryPlane::XY, SymmetryPlane::YZ, SymmetryPlane::XZ};

    for (SymmetryPlane p : planes) {
        SymmetryCandidate c;
        c.plane = p;
        TopoDS_Shape pos = make_halfspace(p, SymmetryHalf::Positive, bb);
        TopoDS_Shape neg = make_halfspace(p, SymmetryHalf::Negative, bb);
        c.volume_plus  = volume_on_side(shapes, pos);
        c.volume_minus = volume_on_side(shapes, neg);
        c.volume_asymmetry =
            std::abs(c.volume_plus - c.volume_minus) /
            std::max(V_total, 1e-30);
        c.valid = c.volume_asymmetry < volume_tolerance;
        out.candidates.push_back(c);
        if (c.valid) out.valid_planes.push_back(p);
    }
    return out;
}

std::vector<NamedShape> cut_half_space(const NamedShape& in,
                                       SymmetryPlane plane,
                                       SymmetryHalf  half,
                                       const ShapeBBox& bbox) {
    std::vector<NamedShape> out;
    if (plane == SymmetryPlane::None || in.shape.IsNull()) {
        out.push_back(in);
        return out;
    }
    TopoDS_Shape halfspace = make_halfspace(plane, half, bbox);
    BRepAlgoAPI_Common op(in.shape, halfspace);
    op.Build();
    if (!op.IsDone()) return out;
    const TopoDS_Shape result = op.Shape();
    const auto solids = flatten_solids(result);
    if (solids.empty()) return out;  // whole shape was on the other side

    if (solids.size() == 1) {
        out.emplace_back(solids.front(), in.name);
    } else {
        for (std::size_t i = 0; i < solids.size(); ++i) {
            out.emplace_back(solids[i],
                             in.name + "_" + std::to_string(i));
        }
    }
    return out;
}

std::vector<NamedShape> cut_half_space(const std::vector<NamedShape>& in,
                                       SymmetryPlane plane,
                                       SymmetryHalf  half,
                                       const ShapeBBox& bbox) {
    std::vector<NamedShape> out;
    out.reserve(in.size());
    for (const auto& ns : in) {
        auto r = cut_half_space(ns, plane, half, bbox);
        out.insert(out.end(),
                   std::make_move_iterator(r.begin()),
                   std::make_move_iterator(r.end()));
    }
    return out;
}

std::vector<NamedShape> cut_to_region(
    const std::vector<NamedShape>& in,
    const std::vector<std::pair<SymmetryPlane, SymmetryHalf>>& cuts,
    const ShapeBBox& bbox) {
    std::vector<NamedShape> current = in;
    for (const auto& [plane, half] : cuts) {
        current = cut_half_space(current, plane, half, bbox);
        if (current.empty()) break;
    }
    return current;
}

} // namespace mvb
