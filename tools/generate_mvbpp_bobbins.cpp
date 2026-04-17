#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <STEPControl_Writer.hxx>
#include <Interface_Static.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <iostream>
#include <cmath>
#include <filesystem>

static TopoDS_Shape translate_shape(const TopoDS_Shape& shape, double x, double y, double z) {
    gp_Trsf trsf;
    trsf.SetTranslation(gp_Vec(x, y, z));
    return BRepBuilderAPI_Transform(shape, trsf).Shape();
}

static double volume(const TopoDS_Shape& s) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props);
    return props.Mass();
}

static bool exportSTEP(const TopoDS_Shape& shape, const std::string& path) {
    STEPControl_Writer writer;
    Interface_Static::SetCVal("xstep.cascade.unit", "M");
    Interface_Static::SetCVal("write.step.unit", "M");
    if (writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone) return false;
    return writer.Write(path.c_str()) == IFSelect_RetDone;
}

static void printBBox(const std::string& label, const TopoDS_Shape& shape) {
    Bnd_Box box;
    BRepBndLib::Add(shape, box);
    double xmin, ymin, zmin, xmax, ymax, zmax;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    int solids = 0;
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    std::cout << label << " bbox (mm):\n";
    std::cout << "  x: [" << xmin * 1000 << ", " << xmax * 1000 << "]\n";
    std::cout << "  y: [" << ymin * 1000 << ", " << ymax * 1000 << "]\n";
    std::cout << "  z: [" << zmin * 1000 << ", " << zmax * 1000 << "]\n";
    std::cout << "  volume: " << props.Mass() * 1e9 << " mm^3";
    std::cout << "  solids: " << solids << "\n\n";
}

static TopoDS_Shape buildRectangularBobbin() {
    // Same params as scripts/generate_mvb_bobbins.py
    double col_width  = 0.005;
    double col_depth  = 0.004;
    double ww_width   = 0.003;
    double ww_depth   = 0.003;
    double ww_height  = 0.010;
    double wall_thick = 0.001;
    double flange_thick = 0.001;
    double flange_ext   = 0.002;

    double outer_width = ww_width + 2.0 * wall_thick;
    double outer_depth = ww_depth + 2.0 * wall_thick;
    double height = ww_height;

    double hole_width = ww_width * 0.8;
    double hole_depth = ww_depth * 0.8;

    // Body = outer box - inner box - central hole
    gp_Pnt oCorner(-outer_width / 2.0, -outer_depth / 2.0, -height / 2.0);
    TopoDS_Shape outer = BRepPrimAPI_MakeBox(oCorner, outer_width, outer_depth, height).Shape();

    gp_Pnt iCorner(-ww_width / 2.0, -ww_depth / 2.0, -height / 2.0 - 0.0005);
    TopoDS_Shape inner = BRepPrimAPI_MakeBox(iCorner, ww_width, ww_depth, height * 1.001).Shape();

    gp_Pnt cCorner(-hole_width / 2.0, -hole_depth / 2.0, -height / 2.0 - 0.001);
    TopoDS_Shape central = BRepPrimAPI_MakeBox(cCorner, hole_width, hole_depth, height * 1.002).Shape();

    TopoDS_Shape body = BRepAlgoAPI_Cut(outer, inner).Shape();
    body = BRepAlgoAPI_Cut(body, central).Shape();

    // Top flange
    double flange_width = outer_width + flange_ext * 2.0;
    double flange_depth = outer_depth + flange_ext * 2.0;
    gp_Pnt tfCorner(-flange_width / 2.0, -flange_depth / 2.0, height / 2.0);
    TopoDS_Shape topFlange = BRepPrimAPI_MakeBox(tfCorner, flange_width, flange_depth, flange_thick).Shape();
    gp_Pnt thCorner(-hole_width / 2.0, -hole_depth / 2.0, height / 2.0 - 0.0005);
    TopoDS_Shape topHole = BRepPrimAPI_MakeBox(thCorner, hole_width, hole_depth, flange_thick * 1.1).Shape();
    topFlange = BRepAlgoAPI_Cut(topFlange, topHole).Shape();

    // Bottom flange
    gp_Pnt bfCorner(-flange_width / 2.0, -flange_depth / 2.0, -(height / 2.0 + flange_thick));
    TopoDS_Shape bottomFlange = BRepPrimAPI_MakeBox(bfCorner, flange_width, flange_depth, flange_thick).Shape();
    gp_Pnt bhCorner(-hole_width / 2.0, -hole_depth / 2.0, -(height / 2.0 + flange_thick) - 0.0005);
    TopoDS_Shape bottomHole = BRepPrimAPI_MakeBox(bhCorner, hole_width, hole_depth, flange_thick * 1.1).Shape();
    bottomFlange = BRepAlgoAPI_Cut(bottomFlange, bottomHole).Shape();

    TopoDS_Shape bobbin = BRepAlgoAPI_Fuse(body, topFlange).Shape();
    bobbin = BRepAlgoAPI_Fuse(bobbin, bottomFlange).Shape();
    return bobbin;
}

static TopoDS_Shape buildRoundBobbin() {
    double col_width  = 0.005;  // outer_r in python
    double col_depth  = 0.004;
    double ww_width   = 0.003;
    double ww_height  = 0.010;
    double wall_thick = 0.001;
    double flange_thick = 0.001;
    double flange_ext   = 0.002;

    double outer_r = col_width;
    double inner_r = col_depth - wall_thick;
    double height = ww_height;

    // Body = outer cylinder - inner cylinder
    TopoDS_Shape outerCyl = BRepPrimAPI_MakeCylinder(outer_r, height).Shape();
    outerCyl = translate_shape(outerCyl, 0.0, 0.0, -height / 2.0);

    TopoDS_Shape innerCyl = BRepPrimAPI_MakeCylinder(inner_r, height * 1.001).Shape();
    innerCyl = translate_shape(innerCyl, 0.0, 0.0, -height / 2.0);

    TopoDS_Shape body = BRepAlgoAPI_Cut(outerCyl, innerCyl).Shape();

    // Flanges
    double flangeOuterX = outer_r + ww_width + flange_ext;
    double flangeHalfY = col_depth / 2.0;

    // Top flange
    gp_Pnt tfCorner(-flangeOuterX, -flangeHalfY, height / 2.0);
    TopoDS_Shape topFlange = BRepPrimAPI_MakeBox(tfCorner, flangeOuterX * 2.0, flangeHalfY * 2.0, flange_thick).Shape();
    // Match CadQuery cylinder() convention: centered at origin
    TopoDS_Shape topHole = BRepPrimAPI_MakeCylinder(inner_r, flange_thick * 1.1).Shape();
    topHole = translate_shape(topHole, 0.0, 0.0, height / 2.0 - flange_thick * 1.1 / 2.0);
    TopoDS_Shape topFlangeNet = BRepAlgoAPI_Cut(topFlange, topHole).Shape();

    // Bottom flange
    gp_Pnt bfCorner(-flangeOuterX, -flangeHalfY, -(height / 2.0 + flange_thick));
    TopoDS_Shape bottomFlange = BRepPrimAPI_MakeBox(bfCorner, flangeOuterX * 2.0, flangeHalfY * 2.0, flange_thick).Shape();
    TopoDS_Shape bottomHole = BRepPrimAPI_MakeCylinder(inner_r, flange_thick * 1.1).Shape();
    bottomHole = translate_shape(bottomHole, 0.0, 0.0, -(height / 2.0 + flange_thick * 1.1) - flange_thick * 1.1 / 2.0);
    TopoDS_Shape bottomFlangeNet = BRepAlgoAPI_Cut(bottomFlange, bottomHole).Shape();

    std::cout << "  [DEBUG] body vol: " << volume(body)*1e9 << " mm^3\n";
    std::cout << "  [DEBUG] topFlange gross: " << volume(topFlange)*1e9 << " mm^3\n";
    std::cout << "  [DEBUG] topHole vol: " << volume(topHole)*1e9 << " mm^3\n";
    std::cout << "  [DEBUG] topFlange net: " << volume(topFlangeNet)*1e9 << " mm^3\n";
    std::cout << "  [DEBUG] bottomFlange gross: " << volume(bottomFlange)*1e9 << " mm^3\n";
    std::cout << "  [DEBUG] bottomHole vol: " << volume(bottomHole)*1e9 << " mm^3\n";
    std::cout << "  [DEBUG] bottomFlange net: " << volume(bottomFlangeNet)*1e9 << " mm^3\n";
    std::cout << "  [DEBUG] sum parts: " << (volume(body)+volume(topFlangeNet)+volume(bottomFlangeNet))*1e9 << " mm^3\n";

    TopoDS_Shape fused1 = BRepAlgoAPI_Fuse(body, topFlangeNet).Shape();
    std::cout << "  [DEBUG] body+top fuse vol: " << volume(fused1)*1e9 << " mm^3\n";
    TopoDS_Shape bobbin = BRepAlgoAPI_Fuse(fused1, bottomFlangeNet).Shape();
    std::cout << "  [DEBUG] final fuse vol: " << volume(bobbin)*1e9 << " mm^3\n";
    return bobbin;
}

int main() {
    std::string outDir = "/home/alf/OpenMagnetics/MVB++/bobbin_comparison";
    std::filesystem::create_directories(outDir);

    TopoDS_Shape rect = buildRectangularBobbin();
    std::string rectPath = outDir + "/mvbpp_rectangular_bobbin.step";
    exportSTEP(rect, rectPath);
    printBBox("MVB++ Rectangular", rect);

    TopoDS_Shape round = buildRoundBobbin();
    std::string roundPath = outDir + "/mvbpp_round_bobbin.step";
    exportSTEP(round, roundPath);
    printBBox("MVB++ Round", round);

    return 0;
}
