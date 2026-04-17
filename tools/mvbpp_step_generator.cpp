#include "mvb/MagneticBuilder.h"
#include "mvb/StepExporter.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "Magnetic.h"
#include "Mas.h"
#include "Utils.h"
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
              << "  -h, --help            Show this help\n";
}

static bool processFile(const fs::path& inputPath, const fs::path& outputPath, bool useMkf) {
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
        mvb::DrawConfig config;
        config.includeBobbin = true;
        config.scale = 1000.0;  // Scale to mm like Python MVB
        
        std::string result;
        if (useMkf) {
            try {
                // Use MVB++'s safe MKF wrapper to avoid Coil::wind() crashes on raw MAS files
                auto enriched = mvb::magnetic_autocomplete_safe(magnetic);
                result = builder.drawMagnetic(enriched, outputPath.parent_path().string(), config);
            } catch (const std::exception& e) {
                std::cerr << "MKF enrichment failed: " << e.what() << "\n";
                std::cerr << "Falling back to raw magnetic...\n";
                result = builder.drawMagnetic(magnetic, outputPath.parent_path().string(), config);
            }
        } else {
            result = builder.drawMagnetic(magnetic, outputPath.parent_path().string(), config);
        }
        
        if (!result.empty()) {
            // drawMagnetic writes to <parent>/magnetic.step; rename to requested output path
            fs::path generated = outputPath.parent_path() / "magnetic.step";
            if (fs::exists(generated) && generated != outputPath) {
                fs::rename(generated, outputPath);
            }
            std::cout << "Generated: " << outputPath << "\n";
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
    
    bool success = processFile(inputPath, outputPath, useMkf);
    return success ? 0 : 1;
}
