#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "mvb/Symmetry.h"
#include "mvb/StepExporter.h"
#include "mvb/SectionDrawing.h"
#include "constructive_models/Magnetic.h"
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {

std::string draw_magnetic_impl(const std::string& json_str,
                               const std::string& output_path,
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
    return builder.drawMagnetic(magnetic, output_path, format, include_bobbin,
                                 scale, symmetry_planes,
                                 wire_polygon_segments, core_polygon_segments);
}

py::bytes draw_magnetic_to_bytes_impl(const std::string& json_str,
                                      const std::string& format,
                                      bool include_bobbin,
                                      double scale,
                                      int symmetry_planes,
                                      int wire_polygon_segments,
                                      int core_polygon_segments)
{
    auto tmp = std::filesystem::temp_directory_path() / "mvbpp_tmp.step";
    draw_magnetic_impl(json_str, tmp.string(), format, include_bobbin, scale,
                        symmetry_planes, wire_polygon_segments, core_polygon_segments);

    std::ifstream f(tmp, std::ios::binary);
    if (!f)
        throw std::runtime_error("mvbpp: failed to open temp STEP file");
    std::string data((std::istreambuf_iterator<char>(f)),
                      std::istreambuf_iterator<char>());
    std::filesystem::remove(tmp);
    return py::bytes(data);
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

    # Quarter domain (two symmetry cuts)
    data = mvbpp.draw_magnetic_to_bytes(mag_json, symmetry_planes=2)

    # Write directly to a file
    mvbpp.draw_magnetic(mag_json, "/path/to/output", symmetry_planes=1, scale=1000.0)

    # Inspect which planes are geometrically symmetric before cutting
    planes = mvbpp.get_symmetry_planes(mag_json)   # e.g. ['XY', 'YZ']

    # 2D dimensioned drawings (SVG string)
    svg = mvbpp.draw_dimensioned_front_view(mag_json)
    svg = mvbpp.draw_dimensioned_top_view(mag_json)

Every configuration knob — including polygon segment counts for cores
and wire cross-sections — is passed as an explicit keyword argument.
There is no config object.
)doc";

    // draw_magnetic(json_str, output_path, *, format, include_bobbin, scale,
    //               symmetry_planes, wire_polygon_segments, core_polygon_segments) -> str
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
          "Build geometry from a MAS JSON string and write a STEP/STL file.\n"
          "Returns the output path.\n\n"
          "Parameters\n"
          "----------\n"
          "format                : 'step' (default) or 'stl'\n"
          "include_bobbin        : include bobbin geometry (default True)\n"
          "scale                 : uniform scale factor (default 1.0; use 1000 for mm)\n"
          "symmetry_planes       : 0=full, 1=half, 2=quarter domain\n"
          "wire_polygon_segments : <=0 = exact torus, >0 = polygon segments for wire cross-section\n"
          "core_polygon_segments : polygon segments for core cylinders/circles");

    // draw_magnetic_to_bytes(...) -> bytes
    m.def("draw_magnetic_to_bytes",
          &draw_magnetic_to_bytes_impl,
          py::arg("json_str"),
          py::arg("format")                = std::string("step"),
          py::arg("include_bobbin")        = true,
          py::arg("scale")                 = 1.0,
          py::arg("symmetry_planes")       = 0,
          py::arg("wire_polygon_segments") = mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
          py::arg("core_polygon_segments") = mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
          "Build geometry from a MAS JSON string.\n"
          "Returns the STEP file contents as bytes (no file I/O required).");

    // get_symmetry_planes(json_str) -> list[str]
    m.def("get_symmetry_planes",
          &get_symmetry_planes_impl,
          py::arg("json_str"),
          py::arg("wire_polygon_segments") = mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
          py::arg("core_polygon_segments") = mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
          "Detect symmetry planes of a magnetic from its MAS JSON string.\n"
          "Returns a list of plane names: subsets of ['XY', 'YZ', 'XZ'].");

    // ── 2D dimensioned drawings (SVG) ────────────────────────────────────────
    m.def("draw_dimensioned_front_view",
          [](const std::string& json_str,
             double width_px,
             double label_font_px,
             const std::string& projection_color,
             const std::string& dimension_color) {
              auto magnetic = parse_enriched(json_str);
              return mvb::SectionDrawing::drawDimensionedFrontView(
                  magnetic, width_px, label_font_px,
                  projection_color, dimension_color);
          },
          py::arg("json_str"),
          py::arg("width_px")         = 800.0,
          py::arg("label_font_px")    = 12.0,
          py::arg("projection_color") = std::string("#000000"),
          py::arg("dimension_color")  = std::string("#1976d2"),
          "Render a dimensioned front-elevation cross-section SVG.\n"
          "Returns the SVG XML as a string.\n\n"
          "Annotations: core B (or C for toroids), window height D, every gap length,\n"
          "and the chunk sizes between gaps.");

    m.def("draw_dimensioned_top_view",
          [](const std::string& json_str,
             double width_px,
             double label_font_px,
             const std::string& projection_color,
             const std::string& dimension_color) {
              auto magnetic = parse_enriched(json_str);
              return mvb::SectionDrawing::drawDimensionedTopView(
                  magnetic, width_px, label_font_px,
                  projection_color, dimension_color);
          },
          py::arg("json_str"),
          py::arg("width_px")         = 800.0,
          py::arg("label_font_px")    = 12.0,
          py::arg("projection_color") = std::string("#000000"),
          py::arg("dimension_color")  = std::string("#1976d2"),
          "Render a dimensioned top-view SVG (X horizontal, Z vertical).\n"
          "Returns the SVG XML as a string.\n\n"
          "Annotates the full per-family dimension set (A, E, F, C, G, H, J, L, K, F2).");

    m.def("write_dimensioned_front_view",
          [](const std::string& json_str, const std::string& output_path) {
              auto magnetic = parse_enriched(json_str);
              mvb::SectionDrawing::writeDimensionedFrontView(magnetic, output_path);
          },
          py::arg("json_str"),
          py::arg("output_path"),
          "Build the dimensioned front view and write it to disk.");

    m.def("write_dimensioned_top_view",
          [](const std::string& json_str, const std::string& output_path) {
              auto magnetic = parse_enriched(json_str);
              mvb::SectionDrawing::writeDimensionedTopView(magnetic, output_path);
          },
          py::arg("json_str"),
          py::arg("output_path"),
          "Build the dimensioned top view and write it to disk.");
}
