#include "mvb/StepExporter.h"

#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Bnd_Box.hxx>
#include <gp_Trsf.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <StlAPI_Writer.hxx>
#include <TDF_Label.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDataStd_Name.hxx>
#include <TDocStd_Application.hxx>
#include <TDocStd_Document.hxx>
#include <TCollection_ExtendedString.hxx>
#include <TopoDS_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <random>

namespace mvb {

namespace {

std::string label_name(const TDF_Label& label) {
    Handle(TDataStd_Name) attr;
    if (label.FindAttribute(TDataStd_Name::GetID(), attr)) {
        const TCollection_ExtendedString& n = attr->Get();
        std::string out;
        out.reserve(n.Length());
        for (int i = 1; i <= n.Length(); ++i) {
            out.push_back(static_cast<char>(n.Value(i)));
        }
        return out;
    }
    return {};
}

} // namespace

bool exportSTEP(const std::vector<NamedShape>& shapes,
                const std::string& filepath) {
    if (shapes.empty()) return false;

    Handle(TDocStd_Document) doc =
        new TDocStd_Document(TCollection_ExtendedString("BinXCAF"));
    Handle(XCAFDoc_ShapeTool) shapeTool =
        XCAFDoc_DocumentTool::ShapeTool(doc->Main());

    // Model coordinates are in SI (metres) — the whole MVB++/MKF/MAS convention
    // (see exportSTLToBytes below). OCCT's length unit is millimetres and it
    // faithfully labels whatever numbers it is given: writing metre-valued
    // coordinates with STEPCAFControl_Writer's default emits
    // SI_UNIT(.MILLI.,.METRE.) over numbers that are actually metres, so every
    // STEP opens 1000x too small (ABT #317). Setting write.step.unit=M does NOT
    // fix it — OCCT then divides the coordinates by 1000 to match the metre
    // header, preserving the same wrong physical size (verified empirically).
    // The only correct fix is to convert the geometry metres->millimetres here,
    // once, so the emitted numbers match the default mm header. This is the
    // single source of truth for the conversion: callers pass metres and MUST
    // NOT pre-scale.
    gp_Trsf metresToMillimetres;
    metresToMillimetres.SetScale(gp_Pnt(0, 0, 0), 1000.0);

    for (const auto& ns : shapes) {
        if (ns.shape.IsNull()) continue;
        const TopoDS_Shape shapeMm =
            BRepBuilderAPI_Transform(ns.shape, metresToMillimetres).Shape();
        // AddShape(isAssembly=false): register the shape as a free top-level
        // component so STEPCAFControl_Writer emits it as a discrete product
        // (with our name attached) rather than burying it in a compound.
        TDF_Label lab = shapeTool->AddShape(shapeMm, Standard_False);
        if (!ns.name.empty()) {
            TDataStd_Name::Set(lab, TCollection_ExtendedString(ns.name.c_str()));
        }
    }

    STEPCAFControl_Writer writer;
    if (!writer.Transfer(doc)) {
        std::cerr << "mvb::exportSTEP: STEPCAFControl_Writer::Transfer failed\n";
        return false;
    }
    if (writer.Write(filepath.c_str()) != IFSelect_RetDone) {
        std::cerr << "mvb::exportSTEP: write failed to " << filepath << "\n";
        return false;
    }
    return true;
}

bool exportSTEP(const std::vector<TopoDS_Shape>& shapes,
                const std::vector<std::string>& names,
                const std::string& filepath) {
    std::vector<NamedShape> ns;
    ns.reserve(shapes.size());
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        const std::string n =
            (i < names.size()) ? names[i] : std::string{};
        ns.emplace_back(shapes[i], n);
    }
    return exportSTEP(ns, filepath);
}

std::vector<NamedShape> importSTEP(const std::string& filepath) {
    std::vector<NamedShape> out;

    Handle(TDocStd_Document) doc =
        new TDocStd_Document(TCollection_ExtendedString("BinXCAF"));

    STEPCAFControl_Reader reader;
    if (reader.ReadFile(filepath.c_str()) != IFSelect_RetDone) return out;
    if (!reader.Transfer(doc)) return out;

    Handle(XCAFDoc_ShapeTool) shapeTool =
        XCAFDoc_DocumentTool::ShapeTool(doc->Main());

    TDF_LabelSequence freeShapes;
    shapeTool->GetFreeShapes(freeShapes);
    for (Standard_Integer i = 1; i <= freeShapes.Length(); ++i) {
        const TDF_Label lab = freeShapes.Value(i);
        TopoDS_Shape s = shapeTool->GetShape(lab);
        if (s.IsNull()) continue;
        out.emplace_back(s, label_name(lab));
    }
    return out;
}

bool exportSTL(const TopoDS_Shape& compound, const std::string& filepath) {
    if (compound.IsNull()) return false;
    // Own the metres->millimetres conversion here, like exportSTEP and
    // exportSTLToBytes (ABT #317), so every exporter is uniform: callers hand over
    // native SI-metre geometry and MUST pass scale=1.0 (NOT the old 1000 — passing
    // 1000 now double-scales, 1e6x too big). The 0.001 deflection was always in the
    // caller's mm space, so it is unchanged.
    gp_Trsf metresToMillimetres;
    metresToMillimetres.SetScale(gp_Pnt(0, 0, 0), 1000.0);
    const TopoDS_Shape mm = BRepBuilderAPI_Transform(compound, metresToMillimetres).Shape();
    BRepMesh_IncrementalMesh mesh(mm, 0.001);
    StlAPI_Writer writer;
    writer.Write(mm, filepath.c_str());
    return true;
}

std::string exportSTLToBytes(const std::vector<TopoDS_Shape>& shapes,
                              double toleranceMm,
                              double angularTolerance,
                              bool binary)
{
    if (shapes.empty()) return {};

    // Model coordinates are in SI (metres); scale metres -> millimetres so the STL
    // opens at correct physical size in the CAD/slicer tools that consume it (they
    // assume mm — the frontend even passes its mesh tolerance as `tolMm`). Same
    // convention as exportSTEP (ABT #317): a raw-metre STL opened 1000x too small.
    gp_Trsf metresToMillimetres;
    metresToMillimetres.SetScale(gp_Pnt(0, 0, 0), 1000.0);

    // Combine into a single compound — StlAPI_Writer walks sub-solids.
    TopoDS_Compound compound;
    TopoDS_Builder b;
    b.MakeCompound(compound);
    bool any = false;
    for (const auto& s : shapes) {
        if (s.IsNull()) continue;
        b.Add(compound, BRepBuilderAPI_Transform(s, metresToMillimetres).Shape());
        any = true;
    }
    if (!any) return {};

    // Geometry is now in millimetres (scaled above). Compute a scale-appropriate
    // absolute linear deflection from the overall bounding-box diagonal. Using
    // relative deflection on tiny per-face boxes (e.g. individual wire turns) creates
    // sub-micron meshes and OOMs in WASM. Absolute deflection ~0.5 % of overall size,
    // floored at 1.0 mm (was 0.001 m before the mm scaling — physically identical),
    // keeps mesh counts bounded while preserving recognisable geometry.
    Bnd_Box bbox;
    BRepBndLib::Add(compound, bbox);
    Standard_Real xMin, yMin, zMin, xMax, yMax, zMax;
    bbox.Get(xMin, yMin, zMin, xMax, yMax, zMax);
    const Standard_Real dx = xMax - xMin;
    const Standard_Real dy = yMax - yMin;
    const Standard_Real dz = zMax - zMin;
    const Standard_Real diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
    const Standard_Real linDeflection = std::max(diagonal * 0.005, 1.0);
    BRepMesh_IncrementalMesh mesh(compound, linDeflection, Standard_False,
                                   angularTolerance, Standard_False);
    mesh.Perform();

    // StlAPI_Writer has no in-memory API; write to a temp file and slurp.
    std::filesystem::path tmpDir = std::filesystem::temp_directory_path();
    // Use random suffix to avoid race conditions with concurrent calls.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dist;
    std::string tmpName = "mvbpp_stl_" + std::to_string(dist(gen)) + ".stl";
    std::filesystem::path tmpPath = tmpDir / tmpName;
    StlAPI_Writer writer;
    writer.ASCIIMode() = !binary;
    if (!writer.Write(compound, tmpPath.string().c_str())) return {};

    std::ifstream f(tmpPath, std::ios::binary);
    if (!f) { std::filesystem::remove(tmpPath); return {}; }
    std::string data((std::istreambuf_iterator<char>(f)),
                      std::istreambuf_iterator<char>());
    std::filesystem::remove(tmpPath);
    return data;
}

} // namespace mvb
