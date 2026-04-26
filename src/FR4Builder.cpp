#include "mvb/FR4Builder.h"

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>

#include <algorithm>
#include <cmath>
#include <optional>

namespace mvb {

namespace {

// Bobbin processed description fields we need. The type is on MAS::Coil
// as a variant<Bobbin, string>; callers typically hand us an enriched
// coil where the struct is populated.
std::optional<MAS::CoreBobbinProcessedDescription>
getBobbinPd(const MAS::Coil& coil) {
    const auto& var = coil.get_bobbin();
    const MAS::Bobbin* b = std::get_if<MAS::Bobbin>(&var);
    if (!b) return std::nullopt;
    return b->get_processed_description();
}

} // namespace

TopoDS_Shape FR4Builder::buildFR4Board(
    const MAS::Coil& coil,
    double borderToWireDistance,
    double coreToLayerDistance)
{
    const auto& groupsOpt = coil.get_groups_description();
    if (!groupsOpt || groupsOpt->empty()) return TopoDS_Shape();

    const auto& group0 = (*groupsOpt).front();
    if (group0.get_type() != MAS::WiringTechnology::PRINTED) return TopoDS_Shape();

    const auto& groupCoords = group0.get_coordinates();
    const auto& groupDims   = group0.get_dimensions();
    if (groupCoords.size() < 2 || groupDims.size() < 2) return TopoDS_Shape();

    const double groupX      = groupCoords[0];                // radial
    const double groupZ      = (groupCoords.size() > 1) ? groupCoords[1] : 0.0; // height
    const double groupWidth  = groupDims[0];
    const double groupHeight = groupDims[1];

    // Widest layer across all groups — sets the Y (depth) extent of the board.
    double maxLayerWidth = 0.0;
    for (const auto& g : *groupsOpt) {
        const auto& d = g.get_dimensions();
        if (!d.empty()) maxLayerWidth = std::max(maxLayerWidth, d[0]);
    }

    // Bobbin column metrics.
    double columnWidth  = 0.002;  // 2 mm default
    double columnDepth  = 0.002;
    MAS::ColumnShape columnShape = MAS::ColumnShape::RECTANGULAR;

    auto bobbinPd = getBobbinPd(coil);
    if (bobbinPd) {
        if (auto cw = bobbinPd->get_column_width()) columnWidth = *cw;
        columnDepth = bobbinPd->get_column_depth();
        columnShape = bobbinPd->get_column_shape();
    }

    // Board thickness: group height minus a 4 µm clearance so the FR4 plane
    // doesn't Z-fight with turn planes — with a MIN_FR4_THICKNESS floor.
    const double clearance    = 4e-6;
    const double fr4Thickness = std::max(groupHeight - clearance,
                                          MIN_FR4_THICKNESS - clearance);

    // Radial extents in replicad logic: groupX ± groupWidth/2.
    const double innerRadius     = groupX - groupWidth / 2.0;
    const double outerRadius     = groupX + groupWidth / 2.0 + borderToWireDistance;
    const double holeInnerRadius = innerRadius - coreToLayerDistance;

    // Y-depth of the board: use bobbin column depth + widest layer width.
    double groupDepth = columnDepth + maxLayerWidth;
    if (groupDepth <= 0.0) groupDepth = outerRadius;  // defensive fallback

    TopoDS_Shape board;

    if (columnShape == MAS::ColumnShape::ROUND) {
        // Annulus: outer disc minus inner hole cylinder, oriented along Z,
        // centred at groupZ.
        gp_Ax2 axOuter(gp_Pnt(0.0, 0.0, groupZ - fr4Thickness / 2.0), gp_Dir(0, 0, 1));
        TopoDS_Shape outerDisc = BRepPrimAPI_MakeCylinder(axOuter, outerRadius, fr4Thickness).Shape();

        // Hole extends slightly beyond so the cut is clean.
        const double holeExtra = 2e-3;
        gp_Ax2 axHole(gp_Pnt(0.0, 0.0, groupZ - fr4Thickness / 2.0 - holeExtra / 2.0), gp_Dir(0, 0, 1));
        TopoDS_Shape hole = BRepPrimAPI_MakeCylinder(axHole, holeInnerRadius,
                                                       fr4Thickness + holeExtra).Shape();
        BRepAlgoAPI_Cut cutter(outerDisc, hole);
        board = cutter.IsDone() ? cutter.Shape() : outerDisc;
    }
    else if (columnShape == MAS::ColumnShape::OBLONG) {
        // Rectangular plate, stadium-shaped hole.
        const double boardW = outerRadius * 2.0;
        const double boardD = groupDepth * 2.0;
        gp_Pnt boardCorner(-boardW / 2.0, -boardD / 2.0, groupZ - fr4Thickness / 2.0);
        TopoDS_Shape plate = BRepPrimAPI_MakeBox(boardCorner, boardW, boardD, fr4Thickness).Shape();

        const double fullColDepth = columnDepth * 2.0;
        const double fullColWidth = columnWidth * 2.0;
        const double holeRadius   = columnWidth + coreToLayerDistance;
        const double straightLen  = std::max(fullColDepth - fullColWidth, 0.0);
        const double holeExtra    = 2e-3;
        const double holeThickness = fr4Thickness + holeExtra;

        // Central rectangle
        gp_Pnt rectCorner(-holeRadius, -straightLen / 2.0,
                          groupZ - holeThickness / 2.0);
        TopoDS_Shape rect = BRepPrimAPI_MakeBox(rectCorner, 2.0 * holeRadius,
                                                  straightLen, holeThickness).Shape();

        // Two semicircle caps
        gp_Ax2 axTop(gp_Pnt(0.0, straightLen / 2.0, groupZ - holeThickness / 2.0),
                     gp_Dir(0, 0, 1));
        TopoDS_Shape capTop = BRepPrimAPI_MakeCylinder(axTop, holeRadius, holeThickness).Shape();
        gp_Ax2 axBot(gp_Pnt(0.0, -straightLen / 2.0, groupZ - holeThickness / 2.0),
                     gp_Dir(0, 0, 1));
        TopoDS_Shape capBot = BRepPrimAPI_MakeCylinder(axBot, holeRadius, holeThickness).Shape();

        TopoDS_Shape hole = rect;
        { BRepAlgoAPI_Fuse f(hole, capTop); if (f.IsDone()) hole = f.Shape(); }
        { BRepAlgoAPI_Fuse f(hole, capBot); if (f.IsDone()) hole = f.Shape(); }

        BRepAlgoAPI_Cut cutter(plate, hole);
        board = cutter.IsDone() ? cutter.Shape() : plate;
    }
    else {
        // Rectangular column → rectangular board with rectangular hole.
        const double boardW = outerRadius * 2.0;
        const double boardD = groupDepth * 2.0;
        gp_Pnt boardCorner(-boardW / 2.0, -boardD / 2.0, groupZ - fr4Thickness / 2.0);
        TopoDS_Shape plate = BRepPrimAPI_MakeBox(boardCorner, boardW, boardD, fr4Thickness).Shape();

        const double holeW = (innerRadius - coreToLayerDistance) * 2.0;
        const double holeD = (columnDepth + coreToLayerDistance) * 2.0;
        const double holeExtra    = 2e-3;
        const double holeThickness = fr4Thickness + holeExtra;

        if (holeW > 0.0 && holeD > 0.0) {
            gp_Pnt holeCorner(-holeW / 2.0, -holeD / 2.0,
                              groupZ - holeThickness / 2.0);
            TopoDS_Shape hole = BRepPrimAPI_MakeBox(holeCorner, holeW, holeD, holeThickness).Shape();
            BRepAlgoAPI_Cut cutter(plate, hole);
            board = cutter.IsDone() ? cutter.Shape() : plate;
        } else {
            board = plate;
        }
    }

    return board;
}

} // namespace mvb
