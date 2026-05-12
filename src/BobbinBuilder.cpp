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
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <cmath>
#include <numbers>

namespace mvb {

static void checkIsDone(const BRepAlgoAPI_Cut& op, const char* ctx) {
    if (!op.IsDone()) throw std::runtime_error(std::string("BobbinBuilder: Cut failed: ") + ctx);
}
static void checkIsDone(const BRepAlgoAPI_Fuse& op, const char* ctx) {
    if (!op.IsDone()) throw std::runtime_error(std::string("BobbinBuilder: Fuse failed: ") + ctx);
}

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

        TopoDS_Shape body;
        if (innerR + 1e-9 >= outerR) {
            body = outerCyl;
        } else {
            double innerHeight = height * 1.1;
            TopoDS_Shape innerCyl = BRepPrimAPI_MakeCylinder(innerR, innerHeight).Shape();
            innerCyl = translate_shape(innerCyl, 0.0, 0.0, -innerHeight / 2.0);
            BRepAlgoAPI_Cut cutOp(outerCyl, innerCyl);
            checkIsDone(cutOp, "round body hollow");
            body = cutOp.Shape();
        }

        if (hasFlanges) {
            double flangeOuterR = outerR + wwWidth;
            double holeHeight = flangeThickness * 1.2;
            double holeMargin = (holeHeight - flangeThickness) / 2.0;

            TopoDS_Shape topFlangeDisc = BRepPrimAPI_MakeCylinder(flangeOuterR, flangeThickness).Shape();
            topFlangeDisc = translate_shape(topFlangeDisc, 0.0, 0.0, height / 2.0);
            TopoDS_Shape topHole = BRepPrimAPI_MakeCylinder(innerR, holeHeight).Shape();
            topHole = translate_shape(topHole, 0.0, 0.0, height / 2.0 - holeMargin);
            BRepAlgoAPI_Cut topCut(topFlangeDisc, topHole);
            checkIsDone(topCut, "round top flange");

            TopoDS_Shape bottomFlangeDisc = BRepPrimAPI_MakeCylinder(flangeOuterR, flangeThickness).Shape();
            bottomFlangeDisc = translate_shape(bottomFlangeDisc, 0.0, 0.0, -(height / 2.0 + flangeThickness));
            TopoDS_Shape bottomHole = BRepPrimAPI_MakeCylinder(innerR, holeHeight).Shape();
            bottomHole = translate_shape(bottomHole, 0.0, 0.0, -(height / 2.0 + flangeThickness) - holeMargin);
            BRepAlgoAPI_Cut bottomCut(bottomFlangeDisc, bottomHole);
            checkIsDone(bottomCut, "round bottom flange");

            BRepAlgoAPI_Fuse fuseTop(body, topCut.Shape());
            checkIsDone(fuseTop, "round top fuse");
            BRepAlgoAPI_Fuse fuseBottom(fuseTop.Shape(), bottomCut.Shape());
            checkIsDone(fuseBottom, "round bottom fuse");
            bobbinShape = fuseBottom.Shape();
        } else {
            bobbinShape = body;
        }
    } else {
        // Rectangular (or irregular = EFD) column.
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
        BRepAlgoAPI_Cut bodyCut(outer, central);
        checkIsDone(bodyCut, "rect body hollow");

        if (hasFlanges) {
            double flangeWidth = 2.0 * (colWidth + wwWidth);
            double flangeDepth = 2.0 * (colDepth + wwWidth);
            double holeEps = flangeThickness * 0.1;

            gp_Pnt tfCorner(-flangeWidth / 2.0, -flangeDepth / 2.0, height / 2.0);
            TopoDS_Shape topFlange = BRepPrimAPI_MakeBox(tfCorner, flangeWidth, flangeDepth, flangeThickness).Shape();
            gp_Pnt thCorner(-holeWidth / 2.0, -holeDepth / 2.0, height / 2.0 - holeEps);
            TopoDS_Shape topHole = BRepPrimAPI_MakeBox(thCorner, holeWidth, holeDepth, flangeThickness + 2.0 * holeEps).Shape();
            BRepAlgoAPI_Cut topCut(topFlange, topHole);
            checkIsDone(topCut, "rect top flange");

            gp_Pnt bfCorner(-flangeWidth / 2.0, -flangeDepth / 2.0, -(height / 2.0 + flangeThickness));
            TopoDS_Shape bottomFlange = BRepPrimAPI_MakeBox(bfCorner, flangeWidth, flangeDepth, flangeThickness).Shape();
            gp_Pnt bhCorner(-holeWidth / 2.0, -holeDepth / 2.0, -(height / 2.0 + flangeThickness) - holeEps);
            TopoDS_Shape bottomHole = BRepPrimAPI_MakeBox(bhCorner, holeWidth, holeDepth, flangeThickness + 2.0 * holeEps).Shape();
            BRepAlgoAPI_Cut bottomCut(bottomFlange, bottomHole);
            checkIsDone(bottomCut, "rect bottom flange");

            BRepAlgoAPI_Fuse fuseTop(bodyCut.Shape(), topCut.Shape());
            checkIsDone(fuseTop, "rect top fuse");
            BRepAlgoAPI_Fuse fuseBottom(fuseTop.Shape(), bottomCut.Shape());
            checkIsDone(fuseBottom, "rect bottom fuse");
            bobbinShape = fuseBottom.Shape();
        } else {
            bobbinShape = bodyCut.Shape();
        }
    }

    if (axisIsY && !bobbinShape.IsNull()) {
        gp_Trsf rot;
        rot.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0)), -std::numbers::pi / 2.0);
        bobbinShape = BRepBuilderAPI_Transform(bobbinShape, rot).Shape();
    }

    return bobbinShape;
}

} // namespace mvb
