// ABT #970 step 1: the FOIL MESHING SPIKE.
//
// Before any sheet-wrap primitive is written, one question decides whether the primitive is
// worth writing at all: does a foil winding MESH? A foil turn is a band as tall as the winding
// window and as thin as the foil, so its aspect ratio is 100:1 to 1000:1 -- the regime where a
// tetrahedral mesher either explodes in element count or refuses the boundary outright. This
// tool builds the bands ALONE (no core, no terminals, no primitive in the product path) at the
// real dimensions MKF's layout now produces, so the mesh question can be answered on its own.
//
//   mvbpp_foilspike <out.step> [turns] [thickness_mm] [height_mm] [innerRadius_mm] [gap_mm]
//
// Defaults are two_switch_forward_transformer_complete's Secondary as MKF lays it out after the
// ABT #881 fix: 8 bands, 0.200 mm thick, 25.935 mm tall, innermost at r = 10.238 mm.
#include "mvb/StepExporter.h"
#include "mvb/NamedShape.h"

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include <TopoDS.hxx>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("usage: %s <out.step> [turns] [thickness_mm] [height_mm] [innerR_mm] [gap_mm]\n",
                    argv[0]);
        return 2;
    }
    const std::string out = argv[1];
    const int turns       = argc > 2 ? std::atoi(argv[2]) : 8;
    const double thick    = argc > 3 ? std::atof(argv[3]) : 0.200;
    const double height   = argc > 4 ? std::atof(argv[4]) : 25.935;
    const double innerR   = argc > 5 ? std::atof(argv[5]) : 10.238;
    // The radial film between one band and the next. MKF's layout puts the bands one thickness
    // apart centre to centre, i.e. touching; a real foil winding has an interlayer film, and
    // touching solids are what the mesher has to fragment. Default to the 25 um MAS insulation.
    const double gap      = argc > 6 ? std::atof(argv[6]) : 0.025;
    // FACETING. A closed band's lateral faces are full 360-degree cylinders -- periodic surfaces,
    // which gmsh dispatches to its fragile mesher off the SURFACE rather than the face
    // (OCCFace.cpp:142), and which is what broke 16 of the corpus's 23 mesh failures until the
    // round wire was faceted. segments > 0 builds each band as a prism between two N-gons
    // instead, so every lateral face is planar. 0 keeps the exact cylinders.
    const int segments    = argc > 7 ? std::atoi(argv[7]) : 0;

    // MVB++ BUILDS IN METRES and the STEP exporter declares mm, so the geometry must be metres
    // here; passing mm straight in wrote a part 1000x too big (bbox diagonal 42.7 m) and the
    // mesher refused it for reasons that had nothing to do with foil.
    const double mm = 1e-3;
    std::vector<mvb::NamedShape> shapes;
    for (int i = 0; i < turns; ++i) {
        const double r0 = (innerR + double(i) * (thick + gap)) * mm;
        const double r1 = r0 + thick * mm;
        const double z0 = -0.5 * height * mm;
        const gp_Ax2 axis(gp_Pnt(0, 0, z0), gp_Dir(0, 0, 1));
        TopoDS_Shape band;
        if (segments <= 0) {
            // A band is the difference of two coaxial cylinders: one closed shell, four faces.
            TopoDS_Shape outer = BRepPrimAPI_MakeCylinder(axis, r1, height * mm).Shape();
            TopoDS_Shape inner = BRepPrimAPI_MakeCylinder(axis, r0, height * mm).Shape();
            BRepAlgoAPI_Cut cut(outer, inner);
            if (!cut.IsDone()) {
                std::fprintf(stderr, "foilspike: band %d failed to build\n", i);
                return 1;
            }
            band = cut.Shape();
        }
        else {
            // Two coaxial N-gons, both with their vertices ON the true circles, extruded: the
            // annular face becomes a planar face with a hole and every lateral face is a plane.
            auto ring = [&](double radius) {
                BRepBuilderAPI_MakePolygon poly;
                for (int k = 0; k < segments; ++k) {
                    const double a = 2.0 * M_PI * double(k) / double(segments);
                    poly.Add(gp_Pnt(radius * std::cos(a), radius * std::sin(a), z0));
                }
                poly.Close();
                return poly.Wire();
            };
            BRepBuilderAPI_MakeFace face(ring(r1), Standard_True);
            face.Add(TopoDS::Wire(ring(r0).Reversed()));
            band = BRepPrimAPI_MakePrism(face.Face(), gp_Vec(0, 0, height * mm)).Shape();
        }
        mvb::NamedShape ns;
        ns.shape = band;
        ns.name = "Foil band " + std::to_string(i);
        shapes.push_back(ns);
    }
    if (!mvb::exportSTEP(shapes, out)) {
        std::fprintf(stderr, "foilspike: could not write %s\n", out.c_str());
        return 1;
    }
    std::printf("foilspike: %d band(s), %.4f mm thick, %.4f mm tall, r %.4f..%.4f mm, gap %.4f mm,"
                " %s, aspect %.0f:1 -> %s\n",
                turns, thick, height, innerR,
                innerR + double(turns - 1) * (thick + gap) + thick, gap,
                segments > 0 ? ("faceted " + std::to_string(segments)).c_str() : "exact cylinders",
                height / thick, out.c_str());
    return 0;
}
