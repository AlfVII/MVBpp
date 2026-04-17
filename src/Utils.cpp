#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <variant>
#include <stdexcept>
#include <numbers>
#include <map>

// MKF headers for safe autocomplete
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"

namespace mvb {

double flatten_dimension(const MAS::Dimension& dim) {
    if (std::holds_alternative<MAS::DimensionWithTolerance>(dim)) {
        const auto& dwt = std::get<MAS::DimensionWithTolerance>(dim);
        if (dwt.get_nominal()) {
            return *dwt.get_nominal();
        }
        if (dwt.get_maximum().has_value() && dwt.get_minimum().has_value()) {
            return (*dwt.get_maximum() + *dwt.get_minimum()) / 2.0;
        }
        if (dwt.get_maximum().has_value()) {
            return *dwt.get_maximum();
        }
        if (dwt.get_minimum().has_value()) {
            return *dwt.get_minimum();
        }
        throw std::runtime_error("DimensionWithTolerance missing nominal value");
    }
    return std::get<double>(dim);
}

std::map<std::string, double> flatten_dimensions(const std::map<std::string, MAS::Dimension>& dims) {
    std::map<std::string, double> result;
    for (const auto& [key, dim] : dims) {
        result[key] = flatten_dimension(dim);
    }
    return result;
}

TopoDS_Wire build_polygon_circle(double radius, int segments) {
    if (segments <= 0) {
        gp_Circ circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), radius);
        BRepBuilderAPI_MakeEdge edge(circ);
        BRepBuilderAPI_MakeWire wire(edge);
        return wire.Wire();
    }
    BRepBuilderAPI_MakePolygon poly;
    for (int i = 0; i < segments; ++i) {
        double angle = 2.0 * std::numbers::pi * i / segments;
        poly.Add(gp_Pnt(radius * std::cos(angle), radius * std::sin(angle), 0.0));
    }
    poly.Close();
    return BRepBuilderAPI_MakeWire(poly.Wire()).Wire();
}

TopoDS_Shape build_polygon_cylinder(double height, double radius, int segments) {
    if (segments <= 0) {
        return BRepPrimAPI_MakeCylinder(radius, height).Shape();
    }
    TopoDS_Wire wire = build_polygon_circle(radius, segments);
    BRepBuilderAPI_MakeFace face(wire);
    gp_Vec vec(0, 0, height);
    return BRepPrimAPI_MakePrism(face.Face(), vec).Shape();
}

TopoDS_Shape rotate_shape(const TopoDS_Shape& shape, double rx, double ry, double rz) {
    TopoDS_Shape result = shape;
    if (rx != 0.0) {
        gp_Trsf trsf;
        trsf.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0)), rx);
        result = BRepBuilderAPI_Transform(result, trsf).Shape();
    }
    if (ry != 0.0) {
        gp_Trsf trsf;
        trsf.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,1,0)), ry);
        result = BRepBuilderAPI_Transform(result, trsf).Shape();
    }
    if (rz != 0.0) {
        gp_Trsf trsf;
        trsf.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,0,1)), rz);
        result = BRepBuilderAPI_Transform(result, trsf).Shape();
    }
    return result;
}

gp_Trsf translation_trsf(double x, double y, double z) {
    gp_Trsf trsf;
    trsf.SetTranslation(gp_Vec(x, y, z));
    return trsf;
}

TopoDS_Shape translate_shape(const TopoDS_Shape& shape, double x, double y, double z) {
    return BRepBuilderAPI_Transform(shape, translation_trsf(x, y, z)).Shape();
}

void patch_dimension_nominals(nlohmann::json& j) {
    if (j.is_object()) {
        // If this looks like a dimension object with min/max but no nominal, add one
        bool hasMin = j.contains("minimum") && !j["minimum"].is_null();
        bool hasMax = j.contains("maximum") && !j["maximum"].is_null();
        if (hasMin || hasMax) {
            if (!j.contains("nominal") || j["nominal"].is_null()) {
                if (hasMin && hasMax) {
                    j["nominal"] = (j["minimum"].get<double>() + j["maximum"].get<double>()) / 2.0;
                } else if (hasMax) {
                    j["nominal"] = j["maximum"].get<double>();
                } else {
                    j["nominal"] = j["minimum"].get<double>();
                }
            }
        }
        for (auto& [key, val] : j.items()) {
            patch_dimension_nominals(val);
        }
    } else if (j.is_array()) {
        for (auto& elem : j) {
            patch_dimension_nominals(elem);
        }
    }
}

std::string core_shape_family_to_string(MAS::CoreShapeFamily family) {
    switch (family) {
        case MAS::CoreShapeFamily::C: return "c";
        case MAS::CoreShapeFamily::DRUM: return "drum";
        case MAS::CoreShapeFamily::E: return "e";
        case MAS::CoreShapeFamily::EC: return "ec";
        case MAS::CoreShapeFamily::EFD: return "efd";
        case MAS::CoreShapeFamily::EI: return "ei";
        case MAS::CoreShapeFamily::EL: return "el";
        case MAS::CoreShapeFamily::ELP: return "elp";
        case MAS::CoreShapeFamily::EP: return "ep";
        case MAS::CoreShapeFamily::EPX: return "epx";
        case MAS::CoreShapeFamily::EQ: return "eq";
        case MAS::CoreShapeFamily::ER: return "er";
        case MAS::CoreShapeFamily::ETD: return "etd";
        case MAS::CoreShapeFamily::H: return "h";
        case MAS::CoreShapeFamily::LP: return "lp";
        case MAS::CoreShapeFamily::P: return "p";
        case MAS::CoreShapeFamily::PM: return "pm";
        case MAS::CoreShapeFamily::PQ: return "pq";
        case MAS::CoreShapeFamily::RM: return "rm";
        case MAS::CoreShapeFamily::ROD: return "rod";
        case MAS::CoreShapeFamily::T: return "t";
        case MAS::CoreShapeFamily::U: return "u";
        case MAS::CoreShapeFamily::UT: return "ut";
    }
    return "unknown";
}

OpenMagnetics::Magnetic magnetic_autocomplete_safe(const nlohmann::json& magneticJson) {
    using json = nlohmann::json;

    json coreJson = magneticJson.contains("core") ? magneticJson.at("core") : json::object();
    json coilJson = magneticJson.contains("coil") ? magneticJson.at("coil") : json::object();

    // Construct OpenMagnetics::Core and Coil directly from JSON.
    // The key fix: pass false to Coil constructor to skip wind(),
    // which crashes for raw MAS files with "Basic" bobbins.
    OpenMagnetics::Core core(coreJson);
    OpenMagnetics::Coil coil(coilJson, false);

    OpenMagnetics::Magnetic om;
    om.set_core(core);
    om.set_coil(coil);

    if (magneticJson.contains("distributorsInfo")) {
        om.set_distributors_info(magneticJson.at("distributorsInfo").get<std::optional<std::vector<MAS::DistributorInfo>>>());
    }
    if (magneticJson.contains("manufacturerInfo")) {
        om.set_manufacturer_info(magneticJson.at("manufacturerInfo").get<std::optional<MAS::MagneticManufacturerInfo>>());
    }
    if (magneticJson.contains("rotation")) {
        om.set_rotation(magneticJson.at("rotation").get<std::optional<std::vector<double>>>());
    }

    return OpenMagnetics::magnetic_autocomplete(om, json{});
}

OpenMagnetics::Magnetic magnetic_autocomplete_safe(const MAS::Magnetic& magnetic) {
    json j;
    to_json(j, magnetic);
    return magnetic_autocomplete_safe(j);
}

bool is_shape_usable(const TopoDS_Shape& shape) {
    if (shape.IsNull()) return false;
    try {
        Bnd_Box bb;
        BRepBndLib::Add(shape, bb);
        return !bb.IsVoid();
    } catch (...) {
        return false;
    }
}

TopoDS_Shape cut_bobbin(const TopoDS_Shape& bobbin, const std::vector<TopoDS_Shape>& cutters) {
    if (std::getenv("MVB_NO_BOBBIN_CUT")) return bobbin;
    if (bobbin.IsNull()) return bobbin;

    // Scale into mm for the cut so coordinate magnitudes are well above
    // OCCT's 1e-7 numerical tolerance floor, then scale back to meters.
    const double S = 1000.0;
    gp_Trsf up;   up.SetScale(gp_Pnt(0, 0, 0), S);
    gp_Trsf down; down.SetScale(gp_Pnt(0, 0, 0), 1.0 / S);

    TopTools_ListOfShape tools;
    for (const auto& tool : cutters) {
        if (tool.IsNull()) continue;
        tools.Append(BRepBuilderAPI_Transform(tool, up).Shape());
    }
    if (tools.IsEmpty()) return bobbin;

    TopoDS_Shape bobbinScaled = BRepBuilderAPI_Transform(bobbin, up).Shape();
    TopTools_ListOfShape args;
    args.Append(bobbinScaled);

    BRepAlgoAPI_Cut cutter;
    cutter.SetArguments(args);
    cutter.SetTools(tools);
    cutter.Build();
    if (!cutter.IsDone()) return bobbin;
    TopoDS_Shape candidate = cutter.Shape();
    if (candidate.IsNull()) return bobbin;

    for (TopExp_Explorer exp(candidate, TopAbs_SOLID); exp.More(); exp.Next()) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(exp.Current(), props);
        if (props.Mass() < 0) return bobbin;
    }

    return BRepBuilderAPI_Transform(candidate, down).Shape();
}

} // namespace mvb
