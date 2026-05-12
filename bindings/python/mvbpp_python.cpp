#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "mvb/Symmetry.h"
#include "mvb/StepExporter.h"
#include "mvb/SectionDrawing.h"
#include "constructive_models/Magnetic.h"
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <random>

namespace py = pybind11;

namespace {

py::bytes export_to_bytes(const std::string& json_str,
                          const std::string& format,
                          bool include_bobbin,
                          double scale,
                          int symmetry_planes,
                          int wire_polygon_segments,
                          int core_polygon_segments)
{
    auto j = nlohmann::json::parse(json_str);
    MAS::Magnetic magnetic = j.get<MAS::Magnetic>();
    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(magnetic, include_bobbin, symmetry_planes,
                                       wire_polygon_segments, core_polygon_segments);

    if (scale != 1.0) {
        gp_Trsf trsf;
        trsf.SetScale(gp_Pnt(0, 0, 0), scale);
        for (auto& ns : named)
            ns.shape = BRepBuilderAPI_Transform(ns.shape, trsf).Shape();
    }

    std::string data;
    if (format == "stl") {
        std::vector<TopoDS_Shape> shapes;
        for (const auto& ns : named) shapes.push_back(ns.shape);
        data = mvb::exportSTLToBytes(shapes);
    } else {
        // STEP has no in-memory API; write to a temp file and slurp.
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint64_t> dist;
        auto tmp = std::filesystem::temp_directory_path() /
                   ("mvbpp_" + std::to_string(dist(gen)) + ".step");
        struct TmpCleanup {
            std::filesystem::path p;
            ~TmpCleanup() { std::filesystem::remove(p); }
        } cleanup{tmp};

        if (!mvb::exportSTEP(named, tmp.string()))
            throw std::runtime_error("mvbpp: exportSTEP failed");

        std::ifstream f(tmp, std::ios::binary);
        if (!f) throw std::runtime_error("mvbpp: failed to read temp STEP");
        data.assign((std::istreambuf_iterator<char>(f)),
                     std::istreambuf_iterator<char>());
    }

    if (data.empty())
        throw std::runtime_error("mvbpp: export produced empty output");
    return py::bytes(data);
}

std::string draw_magnetic_impl(const std::string& json_str,
                               const std::string& output_path,
                               const std::string& format,
                               bool include_bobbin,
                               double scale,
                               int symmetry_planes,
                               int wire_polygon_segments,
                               int core_polygon_segments)
{
    auto data = export_to_bytes(json_str, format, include_bobbin, scale,
                                symmetry_planes, wire_polygon_segments,
                                core_polygon_segments);

    std::filesystem::path out(output_path);
    std::filesystem::create_directories(out.parent_path());

    std::ofstream f(out, std::ios::binary);
    if (!f)
        throw std::runtime_error("mvbpp: cannot write to " + output_path);
    std::string s(data);
    f.write(s.data(), s.size());
    return output_path;
}

py::bytes draw_magnetic_to_bytes_impl(const std::string& json_str,
                                       const std::string& format,
                                       bool include_bobbin,
                                       double scale,
                                       int symmetry_planes,
                                       int wire_polygon_segments,
                                       int core_polygon_segments)
{
    return export_to_bytes(json_str, format, include_bobbin, scale,
                           symmetry_planes, wire_polygon_segments,
                           core_polygon_segments);
}

std::vector<std::string> get_symmetry_planes_impl(const std::string& json_str,
                                                   int wire_polygon_segments,
                                                   int core_polygon_segments)
{
    auto j = nlohmann::json::parse(json_str);
    MAS::Magnetic magnetic = j.get<MAS::Magnetic>();
    mvb::MagneticBuilder builder;
    auto shapes = builder.buildAllNamed(magnetic,
                                         /*includeBobbin=*/false,
                                         /*symmetryPlanes=*/0,
                                         wire_polygon_segments,
                                         core_polygon_segments);
    auto result = mvb::analyze_symmetry(shapes);
    std::vector<std::string> names;
    names.reserve(result.valid_planes.size());
    for (auto p : result.valid_planes)
        names.emplace_back(mvb::to_string(p));
    return names;
}

OpenMagnetics::Magnetic parse_enriched(const std::string& json_str) {
    auto j = nlohmann::json::parse(json_str);
    return mvb::magnetic_autocomplete_safe(j);
}

} // namespace

PYBIND11_MODULE(mvbpp, m)
{
    m.doc() = R"doc(
MVB++ Python bindings — magnetics 3-D geometry builder.

Quick start
-----------
    import mvbpp, json

    mag_json = json.dumps(your_mas_magnetic_dict)

    # Full 3D STEP as bytes (default)
    data = mvbpp.draw_magnetic_to_bytes(mag_json)

    # Half domain (one symmetry cut)
    data = mvbpp.draw_magnetic_to_bytes(mag_json, symmetry_planes=1)

    # Write directly to a file
    mvbpp.draw_magnetic(mag_json, "/path/to/output.step", symmetry_planes=1, scale=1000.0)

    # Inspect which planes are geometrically symmetric before cutting
    planes = mvbpp.get_symmetry_planes(mag_json)   # e.g. ['XY', 'YZ']

    # 2D dimensioned drawings (SVG string)
    svg = mvbpp.draw_dimensioned_front_view(mag_json)
    svg = mvbpp.draw_dimensioned_top_view(mag_json)
)doc";

    m.def("draw_magnetic",
          &draw_magnetic_impl,
          py::arg("json_str"),
          py::arg("output_path"),
          py::arg("format")                = std::string("step"),
          py::arg("include_bobbin")        = true,
          py::arg("scale")                 = 1.0,
          py::arg("symmetry_planes")       = 0,
          py::arg("wire_polygon_segments") = mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
          py::arg("core_polygon_segments") = mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
          "Build geometry and write a STEP/STL file. Returns the output path.");

    m.def("draw_magnetic_to_bytes",
          &draw_magnetic_to_bytes_impl,
          py::arg("json_str"),
          py::arg("format")                = std::string("step"),
          py::arg("include_bobbin")        = true,
          py::arg("scale")                 = 1.0,
          py::arg("symmetry_planes")       = 0,
          py::arg("wire_polygon_segments") = mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
          py::arg("core_polygon_segments") = mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
          "Build geometry and return STEP/STL as bytes (no file I/O).");

    m.def("get_symmetry_planes",
          &get_symmetry_planes_impl,
          py::arg("json_str"),
          py::arg("wire_polygon_segments") = mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
          py::arg("core_polygon_segments") = mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
          "Detect symmetry planes. Returns list of plane names: ['XY', 'YZ', 'XZ'].");

    m.def("draw_dimensioned_front_view",
          [](const std::string& json_str, double width_px, double label_font_px,
             const std::string& projection_color, const std::string& dimension_color) {
              return mvb::SectionDrawing::drawDimensionedFrontView(
                  parse_enriched(json_str), width_px, label_font_px,
                  projection_color, dimension_color);
          },
          py::arg("json_str"), py::arg("width_px") = 800.0,
          py::arg("label_font_px") = 12.0,
          py::arg("projection_color") = std::string("#000000"),
          py::arg("dimension_color") = std::string("#1976d2"),
          "Render a dimensioned front-elevation SVG.");

    m.def("draw_dimensioned_top_view",
          [](const std::string& json_str, double width_px, double label_font_px,
             const std::string& projection_color, const std::string& dimension_color) {
              return mvb::SectionDrawing::drawDimensionedTopView(
                  parse_enriched(json_str), width_px, label_font_px,
                  projection_color, dimension_color);
          },
          py::arg("json_str"), py::arg("width_px") = 800.0,
          py::arg("label_font_px") = 12.0,
          py::arg("projection_color") = std::string("#000000"),
          py::arg("dimension_color") = std::string("#1976d2"),
          "Render a dimensioned top-view SVG.");

    m.def("write_dimensioned_front_view",
          [](const std::string& json_str, const std::string& output_path) {
              mvb::SectionDrawing::writeDimensionedFrontView(
                  parse_enriched(json_str), output_path);
          },
          py::arg("json_str"), py::arg("output_path"),
          "Write dimensioned front view to disk.");

    m.def("write_dimensioned_top_view",
          [](const std::string& json_str, const std::string& output_path) {
              mvb::SectionDrawing::writeDimensionedTopView(
                  parse_enriched(json_str), output_path);
          },
          py::arg("json_str"), py::arg("output_path"),
          "Write dimensioned top view to disk.");
}
