#include "mvb/MagneticBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Symmetry.h"
#include "mvb/Utils.h"
#include "mvb/shapes/ShapeBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/BobbinBuilder.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include <nlohmann/json.hpp>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <filesystem>
#include <stdexcept>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <numbers>

namespace mvb {

using json = nlohmann::json;

static bool isCoreToroidal(const MAS::MagneticCore& core) {
    auto geo = core.get_geometrical_description();
    if (!geo) return false;
    for (const auto& piece : *geo) {
        if (piece.get_type() == MAS::CoreGeometricalDescriptionElementType::TOROIDAL) {
            return true;
        }
    }
    return false;
}

template<typename BobbinT>
static MAS::CoreBobbinProcessedDescription getBobbinProcessedT(const std::variant<BobbinT, std::string>& bobbinVar) {
    const BobbinT* bobbin = std::get_if<BobbinT>(&bobbinVar);
    if (bobbin) {
        auto pd = bobbin->get_processed_description();
        if (pd) return *pd;
    }
    return MAS::CoreBobbinProcessedDescription();
}

static MAS::CoreBobbinProcessedDescription getBobbinProcessed(const MAS::Coil& coil) {
    return getBobbinProcessedT(coil.get_bobbin());
}

static MAS::CoreBobbinProcessedDescription getBobbinProcessed(const OpenMagnetics::Coil& coil) {
    return getBobbinProcessedT(coil.get_bobbin());
}

template<typename BobbinT>
static std::string getBobbinNameT(const std::variant<BobbinT, std::string>& bobbinVar, const std::string& fallback) {
    if (const BobbinT* b = std::get_if<BobbinT>(&bobbinVar)) {
        return b->get_name().value_or(fallback);
    }
    if (const std::string* s = std::get_if<std::string>(&bobbinVar)) {
        return s->empty() ? fallback : *s;
    }
    return fallback;
}

static void patchBobbinDimensions(MAS::CoreBobbinProcessedDescription& bobbinPd, const MAS::MagneticCore& core) {
    double colWidth = bobbinPd.get_column_width().value_or(0.0);

    // MKF sometimes leaves column_depth uninitialized and column_width empty.
    // Fall back to the core's central column dimensions, using MKF's convention:
    // column_width = core_column_width / 2 + wall_thickness
    // (this is the outer radius of the bobbin body tube).
    if (colWidth <= 0.0) {
        auto corePd = core.get_processed_description();
        if (corePd && !corePd->get_columns().empty()) {
            const auto& centralCol = corePd->get_columns()[0];
            double wallThickness = bobbinPd.get_wall_thickness();
            if (wallThickness <= 0.0 || std::isnan(wallThickness)) {
                wallThickness = 0.0005;
            }
            if (centralCol.get_width() > 0.0) {
                bobbinPd.set_column_width(centralCol.get_width() / 2.0 + wallThickness);
            }
            if (centralCol.get_depth() > 0.0) {
                bobbinPd.set_column_depth(centralCol.get_depth() / 2.0 + wallThickness);
            }
        }
    }
}

static double getDimValue(const MAS::Dimension& dim) {
    if (const double* v = std::get_if<double>(&dim)) return *v;
    if (const MAS::DimensionWithTolerance* v = std::get_if<MAS::DimensionWithTolerance>(&dim)) {
        auto n = v->get_nominal();
        if (n) return *n;
    }
    return 0.0;
}


static MAS::Wire defaultWire() {
    MAS::Wire wire;
    wire.set_type(MAS::WireType::ROUND);
    MAS::DimensionWithTolerance dim;
    dim.set_nominal(0.001);
    wire.set_outer_diameter(std::optional<MAS::DimensionWithTolerance>(dim));
    wire.set_conducting_diameter(std::optional<MAS::DimensionWithTolerance>(dim));
    return wire;
}

static std::string drawMagneticCommon(const std::vector<NamedShape>& named,
                                      const std::string& outputPath,
                                      const std::string& format,
                                      double scale) {
    std::vector<TopoDS_Shape> allShapes;
    std::vector<std::string>  allNames;
    allShapes.reserve(named.size());
    allNames.reserve(named.size());
    for (const auto& ns : named) {
        allShapes.push_back(ns.shape);
        allNames.push_back(ns.name);
    }

    if (scale != 1.0) {
        gp_Trsf trsf;
        trsf.SetScale(gp_Pnt(0, 0, 0), scale);
        for (auto& s : allShapes)
            s = BRepBuilderAPI_Transform(s, trsf).Shape();
    }

    std::filesystem::path out = outputPath;
    if (format == "stl") {
        out /= "magnetic.stl";
        exportSTL(allShapes.empty() ? TopoDS_Shape() : allShapes[0], out.string());
        return out.string();
    }
    out /= "magnetic.step";
    exportSTEP(allShapes, allNames, out.string());
    return out.string();
}

std::string MagneticBuilder::drawMagnetic(const MAS::Magnetic& magnetic,
                                          const std::string& outputPath,
                                          const std::string& format,
                                          bool includeBobbin,
                                          double scale,
                                          int symmetryPlanes,
                                          int wirePolygonSegments,
                                          int corePolygonSegments) const {
    auto named = buildAllNamed(magnetic, includeBobbin, symmetryPlanes,
                               wirePolygonSegments, corePolygonSegments);
    return drawMagneticCommon(named, outputPath, format, scale);
}

std::string MagneticBuilder::drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                                          const std::string& outputPath,
                                          const std::string& format,
                                          bool includeBobbin,
                                          double scale,
                                          int symmetryPlanes,
                                          int wirePolygonSegments,
                                          int corePolygonSegments) const {
    auto named = buildAllNamed(magnetic, includeBobbin, symmetryPlanes,
                               wirePolygonSegments, corePolygonSegments);
    return drawMagneticCommon(named, outputPath, format, scale);
}

std::vector<TopoDS_Shape> MagneticBuilder::buildCore(const MAS::MagneticCore& core,
                                                     int corePolygonSegments) const {
    std::vector<TopoDS_Shape> result;
    auto geoOpt = core.get_geometrical_description();
    if (!geoOpt) return result;

    for (const auto& piece : *geoOpt) {
        if (piece.get_type() != MAS::CoreGeometricalDescriptionElementType::HALF_SET
            && piece.get_type() != MAS::CoreGeometricalDescriptionElementType::TOROIDAL) {
            continue;
        }

        auto shapeOpt = piece.get_shape();
        if (!shapeOpt) continue;
        const MAS::CoreShape* shapeData = std::get_if<MAS::CoreShape>(&*shapeOpt);
        if (!shapeData) continue;

        auto builder = shapes::createShapeBuilder(shapeData->get_family(), "", corePolygonSegments);
        if (!builder) continue;

        TopoDS_Shape shape = builder->buildPiece(*shapeData);
        if (shape.IsNull()) continue;

        auto dimsOpt = shapeData->get_dimensions();
        auto dims = dimsOpt ? flatten_dimensions(*dimsOpt) : std::map<std::string, double>{};

        // Apply rotation from geometrical description FIRST; Python MVB applies
        // machining in the post-rotation frame so the gap coordinates refer to
        // the already-flipped piece. Applying machining before rotation leaves
        // the gap tool positioned outside the flipped piece and the cut is a
        // silent no-op (observed on stacked E cores — example 18).
        auto rotOpt = piece.get_rotation();
        if (rotOpt && rotOpt->size() >= 3) {
            shape = rotate_shape(shape, (*rotOpt)[0], (*rotOpt)[1], (*rotOpt)[2]);
        }

        // Apply machining after rotation. Wrap each cut so that an
        // OCCT/std failure on one column doesn't abort the whole core —
        // the un-gapped shape is a strictly-better fallback than no shape.
        auto machiningOpt = piece.get_machining();
        if (machiningOpt) {
            for (const auto& mach : *machiningOpt) {
                try {
                    shape = builder->applyMachining(shape, mach, dims);
                } catch (...) {
                    // Keep the previous shape; gap cut on this column
                    // would have produced a degenerate solid.
                }
            }
        }

        // Apply translation from geometrical description (after machining).
        // Toroidal pieces are rotated {π/2, π/2, 0} upstream so their hole axis
        // points along world Y. MKF still emits the stack offset in coords[2]
        // (MAS Z, the pre-rotation axial) for historical parity with
        // two-piece sets. Swap Y↔Z here so stacked rings displace along the
        // post-rotation axial direction and visually stack along the hole,
        // not perpendicular to it.
        auto coords = piece.get_coordinates();
        if (coords.size() >= 3) {
            if (piece.get_type() == MAS::CoreGeometricalDescriptionElementType::TOROIDAL) {
                shape = translate_shape(shape, coords[0], coords[2], coords[1]);
            } else {
                shape = translate_shape(shape, coords[0], coords[1], coords[2]);
            }
        }

        // Drop phantom sub-solids left by OCCT boolean fuse/cut artifacts
        // (observed on PQ) — but keep every solid that is a meaningful
        // fraction of the largest, so distributed-gap cores (e.g. two
        // central-column subtractive gaps on PQ5050) retain the floating
        // middle chunk between the two gaps.
        {
            std::vector<std::pair<double, TopoDS_Shape>> solids;
            double largestVol = 0.0;
            for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next()) {
                GProp_GProps props;
                BRepGProp::VolumeProperties(e.Current(), props);
                double v = std::abs(props.Mass());
                solids.emplace_back(v, e.Current());
                if (v > largestVol) largestVol = v;
            }
            if (solids.size() > 1) {
                // 1% of the largest volume is well below any real core
                // chunk but well above typical boolean residue.
                const double threshold = largestVol * 0.01;
                TopoDS_Compound compound;
                BRep_Builder builder;
                builder.MakeCompound(compound);
                int kept = 0;
                for (const auto& [v, s] : solids) {
                    if (v >= threshold) { builder.Add(compound, s); ++kept; }
                }
                if (kept >= 1) shape = compound;
            }
        }
        result.push_back(shape);
    }
    return result;
}

TopoDS_Shape MagneticBuilder::buildBobbin(const MAS::Coil& coil, const MAS::MagneticCore& core) const {
    auto bobbinPd = getBobbinProcessed(coil);
    patchBobbinDimensions(bobbinPd, core);
    if (bobbinPd.get_column_width().value_or(0.0) <= 0.0) return TopoDS_Shape();
    double flangeThickness = bobbinPd.get_wall_thickness();
    if (flangeThickness < 0.0 || std::isnan(flangeThickness)) flangeThickness = 0.0;
    return BobbinBuilder::buildBobbin(bobbinPd, flangeThickness, !isCoreToroidal(core));
}

TopoDS_Shape MagneticBuilder::buildBobbin(const OpenMagnetics::Coil& coil, const MAS::MagneticCore& core) const {
    auto bobbinPd = getBobbinProcessed(coil);
    patchBobbinDimensions(bobbinPd, core);
    if (bobbinPd.get_column_width().value_or(0.0) <= 0.0) return TopoDS_Shape();
    double flangeThickness = bobbinPd.get_wall_thickness();
    if (flangeThickness < 0.0 || std::isnan(flangeThickness)) flangeThickness = 0.0;
    return BobbinBuilder::buildBobbin(bobbinPd, flangeThickness, !isCoreToroidal(core));
}

template<typename CoilT, typename WireT>
static std::vector<TopoDS_Shape> buildTurnsImpl(const CoilT& coil, const MAS::MagneticCore& core, std::vector<std::string>* outNames, int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS) {
    std::vector<TopoDS_Shape> result;
    auto turnsOpt = coil.get_turns_description();
    if (!turnsOpt || turnsOpt->empty()) return result;

    TurnBuilder::clearCache();  // Reset per-call so different magnetics don't share cached shapes
    bool toroidal = isCoreToroidal(core);
    auto bobbinPd = getBobbinProcessed(coil);
    patchBobbinDimensions(bobbinPd, core);

    const auto& funcDesc = coil.get_functional_description();
    std::map<std::string, MAS::Wire> wireMap;
    for (const auto& winding : funcDesc) {
        const auto& wireVar = winding.get_wire();
        if (std::holds_alternative<std::string>(wireVar)) {
            try {
                wireMap[winding.get_name()] = OpenMagnetics::find_wire_by_name(std::get<std::string>(wireVar));
            } catch (...) {
                wireMap[winding.get_name()] = defaultWire();
            }
        } else {
            wireMap[winding.get_name()] = std::get<WireT>(wireVar);
        }
    }

    size_t turnIdx = 0;
    for (const auto& turn : *turnsOpt) {
        MAS::Wire wire = defaultWire();
        auto it = wireMap.find(turn.get_winding());
        if (it != wireMap.end()) {
            wire = it->second;
        }

        TopoDS_Shape t = TurnBuilder::buildTurn(turn, wire, bobbinPd, toroidal, wirePolygonSegments);
        if (!t.IsNull()) {
            result.push_back(t);
            if (outNames) {
                const std::string& n = turn.get_name();
                outNames->push_back(n.empty() ? ("Turn_" + std::to_string(turnIdx)) : n);
            }
        }
        ++turnIdx;
    }
    return result;
}

std::vector<TopoDS_Shape> MagneticBuilder::buildTurns(const MAS::Coil& coil, const MAS::MagneticCore& core,
                                                       int wirePolygonSegments) const {
    return buildTurnsImpl<MAS::Coil, MAS::Wire>(coil, core, nullptr, wirePolygonSegments);
}

std::vector<TopoDS_Shape> MagneticBuilder::buildTurns(const OpenMagnetics::Coil& coil, const MAS::MagneticCore& core,
                                                       int wirePolygonSegments) const {
    return buildTurnsImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(coil, core, nullptr, wirePolygonSegments);
}

std::vector<TopoDS_Shape> MagneticBuilder::buildTurns(const MAS::Coil& coil, const MAS::MagneticCore& core,
                                                       std::vector<std::string>& names,
                                                       int wirePolygonSegments) const {
    return buildTurnsImpl<MAS::Coil, MAS::Wire>(coil, core, &names, wirePolygonSegments);
}

std::vector<TopoDS_Shape> MagneticBuilder::buildTurns(const OpenMagnetics::Coil& coil, const MAS::MagneticCore& core,
                                                       std::vector<std::string>& names,
                                                       int wirePolygonSegments) const {
    return buildTurnsImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(coil, core, &names, wirePolygonSegments);
}

// ---- Named-shape overloads ------------------------------------------------

std::vector<NamedShape> MagneticBuilder::buildCoreNamed(const MAS::MagneticCore& core,
                                                         int corePolygonSegments) const {
    auto shapes = buildCore(core, corePolygonSegments);
    std::vector<NamedShape> out;
    out.reserve(shapes.size());
    const std::string base = core.get_name().value_or("Core");
    const bool single = shapes.size() == 1;
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        out.emplace_back(shapes[i],
                         single ? base : base + "_" + std::to_string(i));
    }
    return out;
}

std::vector<NamedShape> MagneticBuilder::buildTurnsNamed(const MAS::Coil& coil,
                                                         const MAS::MagneticCore& core,
                                                         int wirePolygonSegments) const {
    std::vector<std::string> names;
    auto shapes = buildTurns(coil, core, names, wirePolygonSegments);
    std::vector<NamedShape> out;
    out.reserve(shapes.size());
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        const std::string n = (i < names.size() && !names[i].empty())
                                ? names[i]
                                : "Turn_" + std::to_string(i);
        out.emplace_back(shapes[i], n);
    }
    return out;
}

std::vector<NamedShape> MagneticBuilder::buildTurnsNamed(const OpenMagnetics::Coil& coil,
                                                         const MAS::MagneticCore& core,
                                                         int wirePolygonSegments) const {
    std::vector<std::string> names;
    auto shapes = buildTurns(coil, core, names, wirePolygonSegments);
    std::vector<NamedShape> out;
    out.reserve(shapes.size());
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        const std::string n = (i < names.size() && !names[i].empty())
                                ? names[i]
                                : "Turn_" + std::to_string(i);
        out.emplace_back(shapes[i], n);
    }
    return out;
}

NamedShape MagneticBuilder::buildBobbinNamed(const MAS::Coil& coil,
                                             const MAS::MagneticCore& core) const {
    NamedShape ns;
    ns.shape = buildBobbin(coil, core);
    ns.name = getBobbinNameT<MAS::Bobbin>(coil.get_bobbin(), "Bobbin");
    return ns;
}

NamedShape MagneticBuilder::buildBobbinNamed(const OpenMagnetics::Coil& coil,
                                             const MAS::MagneticCore& core) const {
    NamedShape ns;
    ns.shape = buildBobbin(coil, core);
    ns.name = getBobbinNameT<OpenMagnetics::Bobbin>(coil.get_bobbin(), "Bobbin");
    return ns;
}

template<typename CoilT, typename WireT>
static std::vector<NamedShape> buildTurnsNamedImpl(const CoilT& coil, const MAS::MagneticCore& core,
                                                    int wirePolygonSegments) {
    std::vector<std::string> names;
    auto shapes = buildTurnsImpl<CoilT, WireT>(coil, core, &names, wirePolygonSegments);
    std::vector<NamedShape> out;
    out.reserve(shapes.size());
    for (std::size_t i = 0; i < shapes.size(); ++i) {
        const std::string n = (i < names.size() && !names[i].empty())
                                ? names[i]
                                : "Turn_" + std::to_string(i);
        out.emplace_back(shapes[i], n);
    }
    return out;
}

// Apply up to `numPlanes` symmetry cuts to `shapes`, picking the most valid
// planes from analyze_symmetry. Returns the cut shapes.
static std::vector<NamedShape> apply_symmetry(std::vector<NamedShape> shapes, int numPlanes) {
    if (numPlanes <= 0 || shapes.empty()) return shapes;

    SymmetryResult sym = analyze_symmetry(shapes);
    if (sym.valid_planes.empty()) return shapes;

    const int n = std::min(numPlanes, static_cast<int>(sym.valid_planes.size()));
    ShapeBBox bbox = aggregate_bbox(shapes);

    std::vector<std::pair<SymmetryPlane, SymmetryHalf>> cuts;
    cuts.reserve(n);
    for (int i = 0; i < n; ++i)
        cuts.emplace_back(sym.valid_planes[i], SymmetryHalf::Positive);

    return cut_to_region(shapes, cuts, bbox);
}

std::vector<NamedShape> MagneticBuilder::buildAllNamed(const MAS::Magnetic& magnetic,
                                                        bool includeBobbin,
                                                        int symmetryPlanes,
                                                        int wirePolygonSegments,
                                                        int corePolygonSegments) const {
    OpenMagnetics::Magnetic enriched = magnetic_autocomplete_safe(magnetic);
    return buildAllNamed(enriched, includeBobbin, symmetryPlanes,
                         wirePolygonSegments, corePolygonSegments);
}

std::vector<NamedShape> MagneticBuilder::buildAllNamed(const OpenMagnetics::Magnetic& magnetic,
                                                        bool includeBobbin,
                                                        int symmetryPlanes,
                                                        int wirePolygonSegments,
                                                        int corePolygonSegments) const {
    auto all = buildCoreNamed(magnetic.get_core(), corePolygonSegments);

    if (includeBobbin) {
        auto bobbin = buildBobbinNamed(magnetic.get_coil(), magnetic.get_core());
        if (!bobbin.shape.IsNull()) {
            std::vector<TopoDS_Shape> cutters;
            for (const auto& ns : all) cutters.push_back(ns.shape);
            auto turns = buildTurnsImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(
                magnetic.get_coil(), magnetic.get_core(), nullptr, wirePolygonSegments);
            cutters.insert(cutters.end(), turns.begin(), turns.end());
            bobbin.shape = cut_bobbin(bobbin.shape, cutters);
            if (!bobbin.shape.IsNull()) all.push_back(bobbin);
        }
    }

    auto turns = buildTurnsNamedImpl<OpenMagnetics::Coil, OpenMagnetics::Wire>(
        magnetic.get_coil(), magnetic.get_core(), wirePolygonSegments);
    all.insert(all.end(), turns.begin(), turns.end());

    return apply_symmetry(std::move(all), symmetryPlanes);
}

} // namespace mvb
