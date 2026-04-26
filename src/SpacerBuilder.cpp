#include "mvb/SpacerBuilder.h"
#include "mvb/Utils.h"

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pnt.hxx>

namespace mvb {

std::vector<TopoDS_Shape> SpacerBuilder::buildSpacers(
    const std::vector<MAS::CoreGeometricalDescriptionElement>& geometricalDescription)
{
    std::vector<TopoDS_Shape> out;
    for (const auto& part : geometricalDescription) {
        if (part.get_type() != MAS::CoreGeometricalDescriptionElementType::SPACER) continue;

        // MAS convention for spacer: dimensions = [width(X), thickness(Z), depth(Y)].
        // The spacer is placed with its CENTRE at coordinates (x, y, z)
        // in the MVB++ column-along-Y frame (same frame buildCore translates
        // half-sets into).
        auto dimsOpt = part.get_dimensions();
        const auto& coords = part.get_coordinates();
        if (!dimsOpt || dimsOpt->size() < 3) continue;
        const auto& dims = *dimsOpt;

        const double width    = dims[0]; // X
        const double thickness = dims[1]; // Y-axis thickness (column axis)
        const double depth    = dims[2]; // Z

        if (width <= 0.0 || thickness <= 0.0 || depth <= 0.0) continue;

        const double cx = coords.size() > 0 ? coords[0] : 0.0;
        const double cy = coords.size() > 1 ? coords[1] : 0.0;
        const double cz = coords.size() > 2 ? coords[2] : 0.0;

        gp_Pnt corner(cx - width / 2.0, cy - thickness / 2.0, cz - depth / 2.0);
        TopoDS_Shape box = BRepPrimAPI_MakeBox(corner, width, thickness, depth).Shape();
        if (!box.IsNull()) out.push_back(box);
    }
    return out;
}

TopoDS_Shape SpacerBuilder::buildSpacersCompound(
    const std::vector<MAS::CoreGeometricalDescriptionElement>& geometricalDescription)
{
    auto spacers = buildSpacers(geometricalDescription);
    if (spacers.empty()) return TopoDS_Shape();
    if (spacers.size() == 1) return spacers.front();

    BRep_Builder b;
    TopoDS_Compound compound;
    b.MakeCompound(compound);
    for (const auto& s : spacers) b.Add(compound, s);
    return compound;
}

} // namespace mvb
