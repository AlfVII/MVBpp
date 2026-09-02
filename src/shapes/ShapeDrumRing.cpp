#include "mvb/shapes/ShapeDrumRing.h"
#include "mvb/Utils.h"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>

namespace mvb {
namespace shapes {

// Build frame: identical to ShapeDrum (extruded along Z, then centred by the base class), so
// the ring is added in the already-centred frame — z from -L/2 to +L/2 keeps it concentric
// and vertically centred on the drum, which is how the matched pairs are assembled.
TopoDS_Shape ShapeDrumRing::applyExtras(const std::map<std::string, double>& dims,
                                        const TopoDS_Shape& piece) const {
    // Drum first: asymmetric-flange trim, bore, and the centring translation.
    TopoDS_Shape result = ShapeDrum::applyExtras(dims, piece);

    double j = dims.count("J") ? dims.at("J") : 0.0;
    double k = dims.count("K") ? dims.at("K") : 0.0;
    double l = dims.count("L") ? dims.at("L") : 0.0;
    // Geometry validity (K > A, letters present) is enforced upstream by MKF's
    // CorePieceDrumRing; drawing the bare drum here would be a visible defect, not silent
    // physics, so an incomplete record degrades to the drum rather than aborting the render.
    if (j <= 0.0 || k <= 0.0 || l <= 0.0 || k >= j) {
        return result;
    }

    TopoDS_Shape ringOuter = build_polygon_cylinder(l, j / 2.0, m_corePolygonSegments);
    ringOuter = translate_shape(ringOuter, 0.0, 0.0, -l / 2.0);
    // The ring bore faces the drum winding: circumscribe so faceting never bites in.
    TopoDS_Shape ringBore =
        build_polygon_cylinder(l, k / 2.0, m_corePolygonSegments, /*circumscribed=*/true);
    ringBore = translate_shape(ringBore, 0.0, 0.0, -l / 2.0);

    BRepAlgoAPI_Cut ring(ringOuter, ringBore);
    if (!ring.IsDone()) {
        return result;
    }

    // Disjoint solids by construction (the clearance gap): the fuse yields the two-solid
    // assembly, which is what a shielded drum physically is.
    BRepAlgoAPI_Fuse assembly(result, ring.Shape());
    return assembly.IsDone() ? assembly.Shape() : result;
}

} // namespace shapes
} // namespace mvb
