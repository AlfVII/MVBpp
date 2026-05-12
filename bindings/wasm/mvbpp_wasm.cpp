#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <emscripten/emscripten.h>

#include "mvb/MagneticBuilder.h"
#include "mvb/Symmetry.h"
#include "mvb/StepExporter.h"
#include "mvb/SectionDrawing.h"
#include "mvb/SectionBuilder.h"
#include "mvb/ProjectionDrawing.h"
#include "mvb/SpacerBuilder.h"
#include "mvb/FR4Builder.h"
#include "mvb/BobbinBuilder.h"
#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Coil.h"
#include <nlohmann/json.hpp>

#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace emscripten;
using json = nlohmann::json;

// ─────────────────────────────────────────────────────────────────────────────
// JS-throw helper + guard template.
//
// All registered embind entry points are wrapped through `guard<&fn>::call` so
// any C++ exception thrown inside is converted to a JS Error with the original
// message preserved. Without this, exceptions cross the WASM/JS boundary as a
// bare numeric pointer that is unreadable from JS unless the runtime was built
// with -sEXPORT_EXCEPTION_HANDLING_HELPERS=1 (which would in turn require
// rebuilding OCCT with -fwasm-exceptions to keep the EH ABI consistent).
// ─────────────────────────────────────────────────────────────────────────────
[[noreturn]] inline void throw_js_error(const std::string& msg) {
    // EM_ASM-thrown JS errors propagate up through WASM and surface in the
    // calling JS context as a regular Error. Control does not return.
    EM_ASM({ throw new Error(UTF8ToString($0)); }, msg.c_str());
    __builtin_unreachable();
}

template <auto Fn> struct guard;

template <typename R, typename... Args, R(*Fn)(Args...)>
struct guard<Fn> {
    static R call(Args... args) {
        try {
            return Fn(args...);
        } catch (const std::exception& e) {
            throw_js_error(std::string("[mvbpp] ") + e.what());
        } catch (...) {
            throw_js_error("[mvbpp] unknown C++ exception");
        }
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────
namespace {

val makeUint8Array(const std::string& data) {
    val memoryView = val(typed_memory_view(data.size(),
                                            reinterpret_cast<const uint8_t*>(data.data())));
    val out = val::global("Uint8Array").new_(data.size());
    out.call<void>("set", memoryView);
    return out;
}

TopoDS_Shape applyUniformScale(const TopoDS_Shape& s, double scale) {
    if (s.IsNull() || scale == 1.0) return s;
    gp_Trsf t;
    t.SetScale(gp_Pnt(0, 0, 0), scale);
    return BRepBuilderAPI_Transform(s, t).Shape();
}

std::vector<TopoDS_Shape> applyUniformScale(const std::vector<TopoDS_Shape>& shapes,
                                             double scale) {
    if (scale == 1.0) return shapes;
    std::vector<TopoDS_Shape> out;
    out.reserve(shapes.size());
    for (const auto& s : shapes) out.push_back(applyUniformScale(s, scale));
    return out;
}

OpenMagnetics::Magnetic parseEnriched(const std::string& json_str) {
    nlohmann::json j;
    try {
        j = json::parse(json_str);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("parseEnriched: json::parse failed: ") + e.what());
    }
    // Skip expensive autocomplete if geometricalDescription is already present
    // (pre-enriched by the JS MKF worker before calling MVB++).
    bool alreadyEnriched =
        j.contains("core") &&
        j.at("core").contains("geometricalDescription") &&
        !j.at("core").at("geometricalDescription").is_null();
    if (alreadyEnriched) {
        using json = nlohmann::json;
        json coreJson = j.at("core");
        json coilJson = j.contains("coil") ? j.at("coil") : json::object();
        OpenMagnetics::Core core;
        try {
            core = OpenMagnetics::Core(coreJson);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("parseEnriched[fast]: Core ctor failed: ") + e.what());
        }
        // Validate that the core has usable geometrical data
        if (core.get_geometrical_description().has_value() &&
            !core.get_geometrical_description()->empty()) {
            OpenMagnetics::Coil coil;
            try {
                coil = OpenMagnetics::Coil(coilJson, false);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("parseEnriched[fast]: Coil ctor failed: ") + e.what());
            }
            OpenMagnetics::Magnetic om;
            om.set_core(core);
            om.set_coil(coil);
            return om;
        }
        throw std::runtime_error("parseEnriched[fast]: core had geometricalDescription but it was empty/unusable");
    }
    try {
        return mvb::magnetic_autocomplete_safe(j);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("parseEnriched[autocomplete]: ") + e.what());
    }
}

mvb::SectionPlane parsePlane(const std::string& plane) {
    if (plane == "XY" || plane == "xy") return mvb::SectionPlane::XY;
    if (plane == "XZ" || plane == "xz") return mvb::SectionPlane::XZ;
    if (plane == "YZ" || plane == "yz") return mvb::SectionPlane::YZ;
    throw std::runtime_error("mvbpp_wasm: unknown plane '" + plane + "'");
}

// ── Core (just the iron) ────────────────────────────────────────────────────
std::vector<TopoDS_Shape> buildCoreShapes(const std::string& magnetic_json,
                                           double scale,
                                           int corePolygonSegments) {
    auto mag = parseEnriched(magnetic_json);
    mvb::MagneticBuilder builder;
    auto shapes = builder.buildCore(mag.get_core(), corePolygonSegments);
    return applyUniformScale(shapes, scale);
}

// ── Full assembly named shapes ──────────────────────────────────────────────
std::vector<mvb::NamedShape> buildAssemblyNamed(const std::string& magnetic_json,
                                                 bool includeBobbin,
                                                 int symmetryPlanes,
                                                 int wirePolygonSegments,
                                                 int corePolygonSegments,
                                                 double scale) {
    auto mag = parseEnriched(magnetic_json);
    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(mag, includeBobbin, symmetryPlanes,
                                        wirePolygonSegments, corePolygonSegments);
    if (scale != 1.0) {
        for (auto& ns : named) ns.shape = applyUniformScale(ns.shape, scale);
    }
    return named;
}

std::vector<TopoDS_Shape> toShapes(const std::vector<mvb::NamedShape>& named) {
    std::vector<TopoDS_Shape> out;
    out.reserve(named.size());
    for (const auto& ns : named) out.push_back(ns.shape);
    return out;
}

// ── Assembly component filter ──────────────────────────────────────────────
// Builds the full assembly, then returns only the requested components.
// `components` is a bitmask: 1=CORE, 2=BOBBIN, 4=TURNS. 0 or 7 = all.
std::vector<TopoDS_Shape> buildAssemblyFiltered(const std::string& magnetic_json,
                                                 int components,
                                                 int symmetryPlanes,
                                                 int wirePolygonSegments,
                                                 int corePolygonSegments,
                                                 double scale) {
    if (components == 0) components = 0b111;
    const bool wantCore   = components & 0b001;
    const bool wantBobbin = components & 0b010;
    const bool wantTurns  = components & 0b100;

    auto mag = parseEnriched(magnetic_json);
    mvb::MagneticBuilder builder;

    std::vector<TopoDS_Shape> out;
    if (wantCore) {
        auto core = builder.buildCoreNamed(mag.get_core(), corePolygonSegments);
        for (const auto& ns : core) out.push_back(ns.shape);
    }
    if (wantBobbin) {
        auto b = builder.buildBobbinNamed(mag.get_coil(), mag.get_core());
        if (!b.shape.IsNull()) out.push_back(b.shape);
    }
    if (wantTurns) {
        auto turns = builder.buildTurnsNamed(mag.get_coil(), mag.get_core(),
                                              wirePolygonSegments);
        for (const auto& ns : turns) out.push_back(ns.shape);
    }

    // Symmetry cut applied after filtering so the cut respects the full bbox.
    if (symmetryPlanes > 0) {
        std::vector<mvb::NamedShape> ns;
        ns.reserve(out.size());
        for (const auto& s : out) ns.emplace_back(s, std::string{});
        auto sym = mvb::analyze_symmetry(ns);
        if (!sym.valid_planes.empty()) {
            const int n = std::min(symmetryPlanes, (int)sym.valid_planes.size());
            std::vector<std::pair<mvb::SymmetryPlane, mvb::SymmetryHalf>> cuts;
            for (int i = 0; i < n; ++i) cuts.emplace_back(sym.valid_planes[i], mvb::SymmetryHalf::Positive);
            auto bbox = mvb::aggregate_bbox(ns);
            auto cut = mvb::cut_to_region(ns, cuts, bbox);
            out.clear();
            out.reserve(cut.size());
            for (const auto& cn : cut) out.push_back(cn.shape);
        }
    }

    return applyUniformScale(out, scale);
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Exposed functions
// ─────────────────────────────────────────────────────────────────────────────
namespace {

// ── STEP assembly export (to Uint8Array) ────────────────────────────────────
val buildMagneticSTEP(const std::string& magnetic_json,
                      bool includeBobbin,
                      double scale,
                      int symmetryPlanes,
                      int wirePolygonSegments,
                      int corePolygonSegments)
{
    auto named = buildAssemblyNamed(magnetic_json, includeBobbin, symmetryPlanes,
                                     wirePolygonSegments, corePolygonSegments, scale);
    // Use random suffix to avoid race conditions with concurrent calls.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dist;
    std::string tmpName = "mvbpp_" + std::to_string(dist(gen)) + ".step";
    auto tmp = std::filesystem::temp_directory_path() / tmpName;
    mvb::exportSTEP(named, tmp.string());
    std::ifstream f(tmp, std::ios::binary);
    if (!f) { std::filesystem::remove(tmp); throw std::runtime_error("buildMagneticSTEP: failed to read temp file"); }
    std::string data((std::istreambuf_iterator<char>(f)),
                      std::istreambuf_iterator<char>());
    std::filesystem::remove(tmp);
    return makeUint8Array(data);
}

// ── STL assembly export ─────────────────────────────────────────────────────
val buildMagneticSTL(const std::string& magnetic_json,
                     bool includeBobbin,
                     double scale,
                     int symmetryPlanes,
                     int wirePolygonSegments,
                     int corePolygonSegments,
                     double stlToleranceMm,
                     double stlAngularTolerance,
                     bool stlBinary)
{
    auto named = buildAssemblyNamed(magnetic_json, includeBobbin, symmetryPlanes,
                                     wirePolygonSegments, corePolygonSegments, scale);
    auto shapes = toShapes(named);
    std::string data = mvb::exportSTLToBytes(shapes, stlToleranceMm,
                                               stlAngularTolerance, stlBinary);
    return makeUint8Array(data);
}

// ── Core-only STL (from full magnetic JSON) ─────────────────────────────────
val buildCoreSTL(const std::string& magnetic_json,
                 double scale,
                 int corePolygonSegments,
                 double stlToleranceMm,
                 double stlAngularTolerance,
                 bool stlBinary)
{
    auto shapes = buildCoreShapes(magnetic_json, scale, corePolygonSegments);
    std::string data = mvb::exportSTLToBytes(shapes, stlToleranceMm,
                                               stlAngularTolerance, stlBinary);
    return makeUint8Array(data);
}

// ── Spacers STL (from full magnetic JSON) ───────────────────────────────────
val buildSpacersSTL(const std::string& magnetic_json,
                    double scale,
                    double stlToleranceMm,
                    double stlAngularTolerance,
                    bool stlBinary)
{
    auto mag = parseEnriched(magnetic_json);
    auto geoOpt = mag.get_core().get_geometrical_description();
    if (!geoOpt) return makeUint8Array("");
    auto spacers = mvb::SpacerBuilder::buildSpacers(*geoOpt);
    auto scaled = applyUniformScale(spacers, scale);
    std::string data = mvb::exportSTLToBytes(scaled, stlToleranceMm,
                                               stlAngularTolerance, stlBinary);
    return makeUint8Array(data);
}

// ── Bobbin STL ──────────────────────────────────────────────────────────────
val buildBobbinSTL(const std::string& magnetic_json,
                   double scale,
                   double stlToleranceMm,
                   double stlAngularTolerance,
                   bool stlBinary)
{
    auto mag = parseEnriched(magnetic_json);
    mvb::MagneticBuilder builder;
    TopoDS_Shape bobbin = builder.buildBobbin(mag.get_coil(), mag.get_core());
    std::vector<TopoDS_Shape> shapes;
    if (!bobbin.IsNull()) shapes.push_back(applyUniformScale(bobbin, scale));
    std::string data = mvb::exportSTLToBytes(shapes, stlToleranceMm,
                                               stlAngularTolerance, stlBinary);
    return makeUint8Array(data);
}

// ── Turns STL ───────────────────────────────────────────────────────────────
val buildTurnsSTL(const std::string& magnetic_json,
                  double scale,
                  int wirePolygonSegments,
                  double stlToleranceMm,
                  double stlAngularTolerance,
                  bool stlBinary)
{
    auto mag = parseEnriched(magnetic_json);
    mvb::MagneticBuilder builder;
    auto turns = builder.buildTurns(mag.get_coil(), mag.get_core(), wirePolygonSegments);
    auto scaled = applyUniformScale(turns, scale);
    std::string data = mvb::exportSTLToBytes(scaled, stlToleranceMm,
                                               stlAngularTolerance, stlBinary);
    return makeUint8Array(data);
}

// ── FR4 board STL ───────────────────────────────────────────────────────────
val buildFR4BoardSTL(const std::string& magnetic_json,
                     double scale,
                     double borderToWireDistance,
                     double coreToLayerDistance,
                     double stlToleranceMm,
                     double stlAngularTolerance,
                     bool stlBinary)
{
    auto mag = parseEnriched(magnetic_json);
    TopoDS_Shape board = mvb::FR4Builder::buildFR4Board(
        mag.get_coil(), borderToWireDistance, coreToLayerDistance);
    std::vector<TopoDS_Shape> shapes;
    if (!board.IsNull()) shapes.push_back(applyUniformScale(board, scale));
    std::string data = mvb::exportSTLToBytes(shapes, stlToleranceMm,
                                               stlAngularTolerance, stlBinary);
    return makeUint8Array(data);
}

// ── Symmetry planes ─────────────────────────────────────────────────────────
std::vector<std::string> getSymmetryPlanes(const std::string& magnetic_json) {
    auto mag = parseEnriched(magnetic_json);
    mvb::MagneticBuilder builder;
    auto shapes = builder.buildAllNamed(mag, /*includeBobbin=*/false, /*symmetryPlanes=*/0);
    auto sym = mvb::analyze_symmetry(shapes);
    std::vector<std::string> out;
    for (auto p : sym.valid_planes) out.emplace_back(mvb::to_string(p));
    return out;
}

// ── Supported families ──────────────────────────────────────────────────────
std::vector<std::string> getSupportedFamilies() {
    return mvb::get_supported_families();
}

// ── Dimensioned 2D drawings (full + piece) ──────────────────────────────────
std::string drawDimensionedFrontView(const std::string& magnetic_json,
                                      double width_px,
                                      double label_font_px,
                                      const std::string& projection_color,
                                      const std::string& dimension_color) {
    auto mag = parseEnriched(magnetic_json);
    return mvb::SectionDrawing::drawDimensionedFrontView(
        mag, width_px, label_font_px, projection_color, dimension_color);
}

std::string drawDimensionedTopView(const std::string& magnetic_json,
                                    double width_px,
                                    double label_font_px,
                                    const std::string& projection_color,
                                    const std::string& dimension_color) {
    auto mag = parseEnriched(magnetic_json);
    return mvb::SectionDrawing::drawDimensionedTopView(
        mag, width_px, label_font_px, projection_color, dimension_color);
}

// Gap-only front view: shows gap dimensions without core shape labels.
std::string drawCoreGappingTechnicalDrawing(const std::string& magnetic_json,
                                             double width_px,
                                             double label_font_px,
                                             const std::string& projection_color,
                                             const std::string& dimension_color) {
    auto mag = parseEnriched(magnetic_json);
    return mvb::SectionDrawing::drawCoreGappingTechnicalDrawing(
        mag, width_px, label_font_px, projection_color, dimension_color);
}

// ── 2D projections & cutaway sections ───────────────────────────────────────
std::string drawCoreProjection(const std::string& magnetic_json,
                                const std::string& plane,
                                int corePolygonSegments,
                                double width_px,
                                double stroke_width,
                                const std::string& stroke_color) {
    auto shapes = buildCoreShapes(magnetic_json, /*scale=*/1.0, corePolygonSegments);
    return mvb::ProjectionDrawing::drawProjection(shapes, parsePlane(plane),
                                                    width_px, /*margin_px=*/40.0,
                                                    stroke_width, stroke_color);
}

std::string drawCoreCrossSection(const std::string& magnetic_json,
                                  const std::string& plane,
                                  double section_offset,
                                  int corePolygonSegments,
                                  double width_px,
                                  double stroke_width,
                                  const std::string& stroke_color) {
    auto shapes = buildCoreShapes(magnetic_json, /*scale=*/1.0, corePolygonSegments);
    return mvb::ProjectionDrawing::drawCutawaySection(shapes, parsePlane(plane),
                                                       section_offset, width_px,
                                                       /*margin_px=*/40.0,
                                                       stroke_width, stroke_color);
}

// Single-piece views (0-indexed against geometricalDescription's HALF_SET/TOROIDAL list).
TopoDS_Shape pickPiece(const std::string& magnetic_json, int pieceIndex,
                        int corePolygonSegments) {
    auto shapes = buildCoreShapes(magnetic_json, 1.0, corePolygonSegments);
    if (pieceIndex < 0 || pieceIndex >= (int)shapes.size()) {
        throw std::runtime_error("pickPiece: index out of range");
    }
    return shapes[pieceIndex];
}

std::string drawPieceProjection(const std::string& magnetic_json,
                                 int pieceIndex,
                                 const std::string& plane,
                                 int corePolygonSegments,
                                 double width_px,
                                 double stroke_width,
                                 const std::string& stroke_color) {
    auto shape = pickPiece(magnetic_json, pieceIndex, corePolygonSegments);
    return mvb::ProjectionDrawing::drawProjection(shape, parsePlane(plane),
                                                    width_px, /*margin_px=*/40.0,
                                                    stroke_width, stroke_color);
}

std::string drawPieceCrossSection(const std::string& magnetic_json,
                                   int pieceIndex,
                                   const std::string& plane,
                                   double section_offset,
                                   int corePolygonSegments,
                                   double width_px,
                                   double stroke_width,
                                   const std::string& stroke_color) {
    auto shape = pickPiece(magnetic_json, pieceIndex, corePolygonSegments);
    return mvb::ProjectionDrawing::drawCutawaySection(shape, parsePlane(plane),
                                                       section_offset, width_px,
                                                       /*margin_px=*/40.0,
                                                       stroke_width, stroke_color);
}

// Full assembly projection (component bitmask: 1=core, 2=bobbin, 4=turns).
std::string drawAssemblyProjection(const std::string& magnetic_json,
                                    const std::string& plane,
                                    int components,
                                    int symmetryPlanes,
                                    int wirePolygonSegments,
                                    int corePolygonSegments,
                                    double width_px,
                                    double stroke_width,
                                    const std::string& stroke_color) {
    auto shapes = buildAssemblyFiltered(magnetic_json, components, symmetryPlanes,
                                         wirePolygonSegments, corePolygonSegments,
                                         /*scale=*/1.0);
    return mvb::ProjectionDrawing::drawProjection(shapes, parsePlane(plane),
                                                    width_px, /*margin_px=*/40.0,
                                                    stroke_width, stroke_color);
}

std::string drawAssemblyCrossSection(const std::string& magnetic_json,
                                      const std::string& plane,
                                      double section_offset,
                                      int components,
                                      int symmetryPlanes,
                                      int wirePolygonSegments,
                                      int corePolygonSegments,
                                      double width_px,
                                      double stroke_width,
                                      const std::string& stroke_color) {
    auto shapes = buildAssemblyFiltered(magnetic_json, components, symmetryPlanes,
                                         wirePolygonSegments, corePolygonSegments,
                                         /*scale=*/1.0);
    return mvb::ProjectionDrawing::drawCutawaySection(shapes, parsePlane(plane),
                                                       section_offset, width_px,
                                                       /*margin_px=*/40.0,
                                                       stroke_width, stroke_color);
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Embind registration
// ─────────────────────────────────────────────────────────────────────────────
EMSCRIPTEN_BINDINGS(mvbpp) {

    register_vector<std::string>("VectorString");

    // 3D export (STL + STEP) — everything that mvbWorker.js used to
    // delegate to replicad+MVB.js is covered here.
    function("buildMagneticSTL",  &guard<&buildMagneticSTL>::call);
    function("buildMagneticSTEP", &guard<&buildMagneticSTEP>::call);
    function("buildCoreSTL",      &guard<&buildCoreSTL>::call);
    function("buildSpacersSTL",   &guard<&buildSpacersSTL>::call);
    function("buildBobbinSTL",    &guard<&buildBobbinSTL>::call);
    function("buildTurnsSTL",     &guard<&buildTurnsSTL>::call);
    function("buildFR4BoardSTL",  &guard<&buildFR4BoardSTL>::call);

    // Metadata
    function("getSymmetryPlanes",    &guard<&getSymmetryPlanes>::call);
    function("getSupportedFamilies", &guard<&getSupportedFamilies>::call);

    // 2D drawings — dimensioned (Python MVB technical-drawing surface)
    function("drawDimensionedFrontView",        &guard<&drawDimensionedFrontView>::call);
    function("drawDimensionedTopView",          &guard<&drawDimensionedTopView>::call);
    function("drawCoreGappingTechnicalDrawing", &guard<&drawCoreGappingTechnicalDrawing>::call);

    // 2D projections — wire-frame outline seen looking along the plane normal
    function("drawCoreProjection",      &guard<&drawCoreProjection>::call);
    function("drawPieceProjection",     &guard<&drawPieceProjection>::call);
    function("drawAssemblyProjection",  &guard<&drawAssemblyProjection>::call);

    // 2D cutaway sections — project the half past the section plane, so the
    // drawing shows the cut surface AND everything beyond it, as seen
    // looking from the plane normal.
    function("drawCoreCrossSection",     &guard<&drawCoreCrossSection>::call);
    function("drawPieceCrossSection",    &guard<&drawPieceCrossSection>::call);
    function("drawAssemblyCrossSection", &guard<&drawAssemblyCrossSection>::call);
}
