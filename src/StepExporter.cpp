#include "mvb/StepExporter.h"
#include <STEPCAFControl_Writer.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDataStd_Name.hxx>
#include <TDF_Label.hxx>
#include <TopoDS_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <StlAPI_Writer.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <iostream>

namespace mvb {

bool exportSTEP(const std::vector<TopoDS_Shape>& shapes,
                const std::vector<std::string>& names,
                const std::string& filepath) {
    if (shapes.empty()) {
        return false;
    }

    // Combine all shapes into a single compound (like Python's cq.Compound.makeCompound)
    TopoDS_Compound compound;
    TopoDS_Builder builder;
    builder.MakeCompound(compound);
    for (const auto& shape : shapes) {
        if (!shape.IsNull()) {
            builder.Add(compound, shape);
        }
    }

    Handle(TDocStd_Document) doc = new TDocStd_Document("BinXCAF");
    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());

    TDF_Label label = shapeTool->AddShape(compound, Standard_False);
    TDataStd_Name::Set(label, "Assembly");

    STEPCAFControl_Writer writer;
    writer.Perform(doc, filepath.c_str());
    return true;
}

bool exportSTL(const TopoDS_Shape& compound,
               const std::string& filepath) {
    if (compound.IsNull()) {
        return false;
    }
    BRepMesh_IncrementalMesh mesh(compound, 0.001);
    StlAPI_Writer writer;
    writer.Write(compound, filepath.c_str());
    return true;
}

} // namespace mvb
