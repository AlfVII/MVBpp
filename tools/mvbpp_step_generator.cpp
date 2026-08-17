#include "mvb/MagneticBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "Magnetic.h"
#include "Mas.h"
#include "Utils.h"
#include "support/Painter.h"
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numbers>
#include <filesystem>

namespace fs = std::filesystem;
using json = nlohmann::json;

static void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog << " [options] <input.json>\n"
              << "Options:\n"
              << "  -o, --output <path>   Output STEP file path\n"
              << "  -d, --output-dir <dir> Output directory (batch mode)\n"
              << "  --no-mkf              Skip MKF enrichment\n"
              << "  --real                Real winding: continuous conductor per (winding, parallel)\n"
              << "  --fem                 FEM geometry: one-piece / conformal conductors (slow); "
                 "default is the fast drawing compound\n"
              << "  --segments <N>        Wire+core polygon segments (0 = exact analytic curves)\n"
              << "  -h, --help            Show this help\n";
}


// ABT #685 (Alf, 2026-08-17): "always produce step and 2D SVG together." The 2D view is the
// faster diagnosis for most 3D faults -- a wrap collision is usually a layout problem visible
// in the cross-section -- so the STEP never ships without it. Painted from the SAME enriched
// magnetic the 3D was built from, so the two cannot describe different coils.
//
// Failures here are reported and swallowed: a projection the Painter does not implement for
// this core (the whole-magnetic views are two-piece-set + rectangular-window only) must not
// cost the caller the STEP it asked for.
static void paintProjections(const MAS::Magnetic& rawMagnetic, bool useRealWinding,
                             const fs::path& stepPath) {
    using OpenMagnetics::Painter;
    using OpenMagnetics::PainterProjection;
    OpenMagnetics::Magnetic magnetic;
    try {
        auto enriched = useRealWinding
            ? mvb::magnetic_autocomplete_safe(nlohmann::json(rawMagnetic), true)
            : mvb::magnetic_autocomplete_safe(nlohmann::json(rawMagnetic));
        magnetic = enriched;
    } catch (const std::exception& e) {
        std::cerr << "2D painter: could not enrich for painting: " << e.what() << "\n";
        return;
    }
    const std::vector<std::pair<PainterProjection, std::string>> views{
        {PainterProjection::XY, ""},          // the winding-window cross-section
        {PainterProjection::XZ, ".top"},      // looking down the column: lanes and bumps
    };
    for (const auto& [projection, suffix] : views) {
        fs::path svg = stepPath;
        svg.replace_extension("");
        svg += suffix + std::string(".svg");
        try {
            fs::remove(svg);
            Painter painter(svg);
            painter.paint_magnetic(magnetic, projection);
            painter.export_svg();
            if (fs::exists(svg)) {
                std::cout << "Generated: " << svg << "\n";
            }
        } catch (const std::exception& e) {
            std::cerr << "2D painter (" << (suffix.empty() ? "XY" : "XZ") << "): " << e.what()
                      << "\n";
        }
    }
}

static bool processFile(const fs::path& inputPath, const fs::path& outputPath, bool useMkf,
                        bool useRealWinding, int segments, bool copperFootprint = false,
                        bool femReady = false) {
    try {
        // Read JSON
        std::ifstream f(inputPath);
        if (!f.is_open()) {
            std::cerr << "Error: Cannot open " << inputPath << "\n";
            return false;
        }
        
        json j;
        f >> j;
        
        // Patch dimensions
        mvb::patch_dimension_nominals(j);
        
        // Get magnetic data
        MAS::Magnetic magnetic;
        if (j.contains("magnetic")) {
            magnetic = j.at("magnetic").get<MAS::Magnetic>();
        } else {
            magnetic = j.get<MAS::Magnetic>();
        }
        
        // Generate STEP
        mvb::MagneticBuilder builder;
        const std::string format       = "step";
        const bool includeBobbin       = true;
        const double scale             = 1.0;     // exporters emit mm natively (ABT #317); 1000 double-scaled 1e6x
        const int symmetryPlanes       = 0;

        std::string result;
        if (useRealWinding) {
            // Real winding requires the MKF wind (turn blocking on); no fallback.
            auto enriched = mvb::magnetic_autocomplete_safe(magnetic, true);
            // MVB_DUMP_TURNS=<path>: dump the ENRICHED turnsDescription (exactly what MKF's wind/
            // blocking placed and what ConductorBuilder consumes) BEFORE drawing, so turn placement
            // can be analysed even when the 3D build throws on a collision. This is the "MKF painter"
            // view: turn centres + per-winding wire, in window cross-section coordinates.
            if (const char* dump = std::getenv("MVB_DUMP_TURNS")) {
                nlohmann::json enr; OpenMagnetics::to_json(enr, enriched);
                nlohmann::json out;
                out["turnsDescription"] = enr["coil"].contains("turnsDescription")
                                              ? enr["coil"]["turnsDescription"] : nlohmann::json::array();
                out["windings"] = nlohmann::json::array();
                OpenMagnetics::Coil coil = enriched.get_coil();
                size_t wi = 0;
                for (const auto& w : enr["coil"]["functionalDescription"]) {
                    const auto wire = coil.resolve_wire(wi++);
                    nlohmann::json ww; ww["name"] = w["name"];
                    if (wire.get_conducting_diameter())
                        ww["conductingDiameter"] = OpenMagnetics::resolve_dimensional_values(wire.get_conducting_diameter().value());
                    if (wire.get_outer_diameter())
                        ww["outerDiameter"] = OpenMagnetics::resolve_dimensional_values(wire.get_outer_diameter().value());
                    if (wire.get_conducting_width())
                        ww["conductingWidth"] = OpenMagnetics::resolve_dimensional_values(wire.get_conducting_width().value());
                    if (wire.get_conducting_height())
                        ww["conductingHeight"] = OpenMagnetics::resolve_dimensional_values(wire.get_conducting_height().value());
                    out["windings"].push_back(ww);
                }
                std::ofstream df(dump); df << out.dump(1);
                std::cerr << "[dump-turns] wrote " << dump << "\n";
            }
            mvb::DrawConfig cfg{format, includeBobbin, scale, symmetryPlanes};
            cfg.useRealWindingGeometry = true;
            cfg.paintCoating = !copperFootprint;  // --copper => bare CONDUCTING footprint (FEM)
            cfg.femReady = femReady;              // --fem => slow one-piece/conformal meshable geometry
            if (segments >= 0) { cfg.wirePolygonSegments = segments; cfg.corePolygonSegments = segments; }
            result = builder.drawMagnetic(enriched, outputPath.parent_path().string(), cfg);
        } else if (useMkf) {
            try {
                // Use MVB++'s safe MKF wrapper to avoid Coil::wind() crashes on raw MAS files
                auto enriched = mvb::magnetic_autocomplete_safe(magnetic);
                result = builder.drawMagnetic(enriched, outputPath.parent_path().string(),
                                               format, includeBobbin, scale, symmetryPlanes);
            } catch (const std::exception& e) {
                std::cerr << "MKF enrichment failed: " << e.what() << "\n";
                std::cerr << "Falling back to raw magnetic...\n";
                result = builder.drawMagnetic(magnetic, outputPath.parent_path().string(),
                                               format, includeBobbin, scale, symmetryPlanes);
            }
        } else {
            result = builder.drawMagnetic(magnetic, outputPath.parent_path().string(),
                                           format, includeBobbin, scale, symmetryPlanes);
        }
        
        if (!result.empty()) {
            // drawMagnetic writes to <parent>/magnetic.step; rename to requested output path
            fs::path generated = outputPath.parent_path() / "magnetic.step";
            if (fs::exists(generated) && generated != outputPath) {
                fs::rename(generated, outputPath);
            }
            std::cout << "Generated: " << outputPath << "\n";
            paintProjections(magnetic, useRealWinding, outputPath);
            return true;
        } else {
            std::cerr << "Failed to generate STEP for " << inputPath << "\n";
            return false;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error processing " << inputPath << ": " << e.what() << "\n";
        return false;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }
    
    fs::path inputPath;
    fs::path outputPath;
    fs::path outputDir;
    bool useMkf = true;
    bool useRealWinding = false;
    bool copperFootprint = false;
    bool femReady = false;
    int segments = -1;  // -1 = builder default; 0 = exact analytic curves

    // Parse arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-o" || arg == "--output") {
            if (++i < argc) outputPath = argv[i];
        } else if (arg == "-d" || arg == "--output-dir") {
            if (++i < argc) outputDir = argv[i];
        } else if (arg == "--no-mkf") {
            useMkf = false;
        } else if (arg == "--real") {
            useRealWinding = true;
        } else if (arg == "--copper") {
            copperFootprint = true;
        } else if (arg == "--fem") {
            femReady = true;
        } else if (arg == "--segments") {
            if (++i < argc) segments = std::stoi(argv[i]);
        } else if (arg[0] != '-') {
            inputPath = arg;
        }
    }
    
    if (inputPath.empty()) {
        std::cerr << "Error: No input file specified\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Determine output path
    if (outputPath.empty()) {
        fs::path dir = outputDir.empty() ? fs::current_path() : outputDir;
        fs::create_directories(dir);
        outputPath = dir / (inputPath.stem().string() + "_mvbpp.step");
    }
    
    bool success = processFile(inputPath, outputPath, useMkf, useRealWinding, segments,
                               copperFootprint, femReady);
    return success ? 0 : 1;
}
