#include "mvb/BobbinBuilder.h"
#include "mvb/Utils.h"
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <cmath>

namespace mvb {

TopoDS_Shape BobbinBuilder::buildBobbin(const MAS::CoreBobbinProcessedDescription& bobbin, double flangeThickness, bool axisIsY) {
    if (flangeThickness < 0.0 || std::isnan(flangeThickness)) {
        throw std::invalid_argument("BobbinBuilder: flangeThickness is invalid");
    }
    bool hasFlanges = flangeThickness > 0.0;

    double colWidth = bobbin.get_column_width().value_or(0.0);
    double colDepth = bobbin.get_column_depth();
    double wallThickness = bobbin.get_wall_thickness();
    if (wallThickness < 0.0 || std::isnan(wallThickness)) {
        throw std::invalid_argument("BobbinBuilder: wallThickness is invalid");
    }

    if (colWidth <= 0.0 || colDepth <= 0.0) {
        throw std::invalid_argument("BobbinBuilder: column_width and column_depth must be > 0");
    }

    const auto& wwList = bobbin.get_winding_windows();
    double wwWidth = 0.0;
    double wwHeight = 0.0;
    double wwDepth = 0.0;
    if (!wwList.empty()) {
        wwWidth = wwList[0].get_width().value_or(0.0);
        wwHeight = wwList[0].get_height().value_or(0.0);
        wwDepth = wwList[0].get_radial_height().value_or(wwWidth);
    }
    if (wwWidth <= 0.0) wwWidth = colWidth;
    if (wwHeight <= 0.0) wwHeight = colDepth;
    if (wwDepth <= 0.0) wwDepth = wwWidth;

    double height = wwHeight;
    if (height <= 0.0) return TopoDS_Shape();

    TopoDS_Shape bobbinShape;

    if (bobbin.get_column_shape() == MAS::ColumnShape::ROUND) {
        double outerR = colWidth;
        double innerR = colWidth - bobbin.get_column_thickness();
        if (innerR <= 0.0) innerR = outerR - wallThickness;

        TopoDS_Shape outerCyl = BRepPrimAPI_MakeCylinder(outerR, height).Shape();
        outerCyl = translate_shape(outerCyl, 0.0, 0.0, -height / 2.0);

        // If innerR >= outerR (e.g. planar bobbin with columnThickness=0),
        // skip the hollow cut and use a solid body — matches MVB Python
        // behaviour where an identical-radius cyl cut returns the uncut
        // cylinder rather than an empty shape.
        TopoDS_Shape body;
        if (innerR + 1e-9 >= outerR) {
            body = outerCyl;
        } else {
            double innerHeight = height * 1.1;
            TopoDS_Shape innerCyl = BRepPrimAPI_MakeCylinder(innerR, innerHeight).Shape();
            innerCyl = translate_shape(innerCyl, 0.0, 0.0, -innerHeight / 2.0);
            body = BRepAlgoAPI_Cut(outerCyl, innerCyl).Shape();
        }

        if (hasFlanges) {
            // Round-column bobbins get round (disc) yokes. Flange extends
            // exactly to the bobbin's winding-window outer edge — no
            // extension, purely derived from MKF-processed data.
            double flangeOuterR = outerR + wwWidth;
            double holeHeight = flangeThickness * 1.2;
            double holeMargin = (holeHeight - flangeThickness) / 2.0;

            TopoDS_Shape topFlangeDisc = BRepPrimAPI_MakeCylinder(flangeOuterR, flangeThickness).Shape();
            topFlangeDisc = translate_shape(topFlangeDisc, 0.0, 0.0, height / 2.0);
            TopoDS_Shape topHole = BRepPrimAPI_MakeCylinder(innerR, holeHeight).Shape();
            topHole = translate_shape(topHole, 0.0, 0.0, height / 2.0 - holeMargin);
            TopoDS_Shape topFlangeNet = BRepAlgoAPI_Cut(topFlangeDisc, topHole).Shape();

            TopoDS_Shape bottomFlangeDisc = BRepPrimAPI_MakeCylinder(flangeOuterR, flangeThickness).Shape();
            bottomFlangeDisc = translate_shape(bottomFlangeDisc, 0.0, 0.0, -(height / 2.0 + flangeThickness));
            TopoDS_Shape bottomHole = BRepPrimAPI_MakeCylinder(innerR, holeHeight).Shape();
            bottomHole = translate_shape(bottomHole, 0.0, 0.0, -(height / 2.0 + flangeThickness) - holeMargin);
            TopoDS_Shape bottomFlangeNet = BRepAlgoAPI_Cut(bottomFlangeDisc, bottomHole).Shape();

            TopoDS_Shape fused = BRepAlgoAPI_Fuse(body, topFlangeNet).Shape();
            bobbinShape = BRepAlgoAPI_Fuse(fused, bottomFlangeNet).Shape();
        } else {
            bobbinShape = body;
        }
    } else {
        // Rectangular (or irregular = EFD) column.
        // MKF: column_width = core_col_half + column_thickness, so
        //   - bobbin tube OUTER (radial) = 2 × column_width
        //   - hole (where the core column sits) = 2 × (column_width − column_thickness)
        // wall_thickness is for the WW outer flange walls (axial), NOT the
        // column tube wall. The column tube wall thickness IS column_thickness.
        double columnThickness = bobbin.get_column_thickness();
        double outerWidth = 2.0 * colWidth;
        double outerDepth = 2.0 * colDepth;
        double holeWidth = 2.0 * std::max(0.0, colWidth - columnThickness);
        double holeDepth = 2.0 * std::max(0.0, colDepth - columnThickness);
        if (holeWidth <= 0.0) holeWidth = outerWidth * 0.8;
        if (holeDepth <= 0.0) holeDepth = outerDepth * 0.8;
        double eps = wallThickness > 0 ? wallThickness * 0.1 : 1e-6;

        gp_Pnt oCorner(-outerWidth / 2.0, -outerDepth / 2.0, -height / 2.0);
        TopoDS_Shape outer = BRepPrimAPI_MakeBox(oCorner, outerWidth, outerDepth, height).Shape();
        gp_Pnt cCorner(-holeWidth / 2.0, -holeDepth / 2.0, -height / 2.0 - eps);
        TopoDS_Shape central = BRepPrimAPI_MakeBox(cCorner, holeWidth, holeDepth, height + 2.0 * eps).Shape();
        TopoDS_Shape body = BRepAlgoAPI_Cut(outer, central).Shape();

        if (hasFlanges) {
            // Flange extends BEYOND the bobbin tube by the same amount on
            // each side in BOTH X and Z directions, because any wire wrapping
            // around the column needs the same yoke surface to rest on.
            // Extension on each side = ww_width (the radial winding window).
            double flangeWidth = 2.0 * (colWidth + wwWidth);
            double flangeDepth = 2.0 * (colDepth + wwWidth);
            double holeEps = flangeThickness * 0.1;

            gp_Pnt tfCorner(-flangeWidth / 2.0, -flangeDepth / 2.0, height / 2.0);
            TopoDS_Shape topFlange = BRepPrimAPI_MakeBox(tfCorner, flangeWidth, flangeDepth, flangeThickness).Shape();
            gp_Pnt thCorner(-holeWidth / 2.0, -holeDepth / 2.0, height / 2.0 - holeEps);
            TopoDS_Shape topHole = BRepPrimAPI_MakeBox(thCorner, holeWidth, holeDepth, flangeThickness + 2.0 * holeEps).Shape();
            TopoDS_Shape topFlangeNet = BRepAlgoAPI_Cut(topFlange, topHole).Shape();

            gp_Pnt bfCorner(-flangeWidth / 2.0, -flangeDepth / 2.0, -(height / 2.0 + flangeThickness));
            TopoDS_Shape bottomFlange = BRepPrimAPI_MakeBox(bfCorner, flangeWidth, flangeDepth, flangeThickness).Shape();
            gp_Pnt bhCorner(-holeWidth / 2.0, -holeDepth / 2.0, -(height / 2.0 + flangeThickness) - holeEps);
            TopoDS_Shape bottomHole = BRepPrimAPI_MakeBox(bhCorner, holeWidth, holeDepth, flangeThickness + 2.0 * holeEps).Shape();
            TopoDS_Shape bottomFlangeNet = BRepAlgoAPI_Cut(bottomFlange, bottomHole).Shape();

            TopoDS_Shape fused = BRepAlgoAPI_Fuse(body, topFlangeNet).Shape();
            bobbinShape = BRepAlgoAPI_Fuse(fused, bottomFlangeNet).Shape();
        } else {
            bobbinShape = body;
        }
    }

    if (axisIsY && !bobbinShape.IsNull()) {
        gp_Trsf rot;
        rot.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0)), -M_PI / 2.0);
        bobbinShape = BRepBuilderAPI_Transform(bobbinShape, rot).Shape();
    }

    return bobbinShape;
}

} // namespace mvb
