// mvbpp_paint2d -- ABT #869: LOOK at the layout, not just its filling factors.
//
// The designs of the #869 class are refused before anything is drawn, so the ordinary
// step-generator path (which paints alongside the STEP) never gets to show them. This tool paints
// MKF's own 2D winding-window cross-section for a design that does NOT fit, in two states:
//
//   <out>.refused.svg    the layout the fit gate actually rejects: real winding on,
//                        coilWindEvenIfNotFit on, nothing relaxed. Turns are placed but the
//                        packer was free to overstuff a layer, so the picture shows one column
//                        of turns running out of the window rather than the layers the design
//                        really needs.
//   <out>.overflow.svg   WHY IT DOES NOT FIT, for the user to act on (Alf, 2026-08-24: "when the
//                        coil won't fit, I want the user to have a visual image of why, which are
//                        the turns overflowing, so they can fix it"). coilWindEvenIfNotFit +
//                        coilAllowHorizontalOverflow, and coating squish deliberately OFF.
//
//                        That combination is the honest one: with horizontal overflow allowed the
//                        SECTION may run past the window edge, but every LAYER still has to be a
//                        real layer (filling factor <= 1.0, no squish allowance). So try_rewind
//                        keeps adding layers until each one holds only what it can, and the
//                        copper the window cannot take ends up DRAWN OUTSIDE IT -- exactly the
//                        turns the user has to remove, in the place they would sit.
//
// Painted from the SAME enriched magnetic in each case, so what is drawn is what MKF laid out.
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "Magnetic.h"
#include "Utils.h"
#include "constructive_models/Coil.h"
#include "support/Painter.h"
#include "support/Settings.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;
namespace fs = std::filesystem;

namespace {

// Enrich with the given fit relaxations and paint the winding-window cross-section.
// Deliberately calls OpenMagnetics::magnetic_autocomplete directly rather than MVB++'s
// magnetic_autocomplete_safe: the safe wrapper REFUSES a layout whose blocking was declined
// (that is its job), and refusing is precisely the case being drawn here.
void paintOne(const json& magneticJson, const fs::path& svgPath, bool windEvenIfNotFit,
              bool squish, bool horizontalOverflowAllowed, const std::string& label) {
    OpenMagnetics::SettingsGuard<bool> realWinding(
        OpenMagnetics::Settings::GetInstance(),
        &OpenMagnetics::Settings::get_coil_use_real_winding_geometry,
        &OpenMagnetics::Settings::set_coil_use_real_winding_geometry, true);
    OpenMagnetics::SettingsGuard<bool> evenIfNotFit(
        OpenMagnetics::Settings::GetInstance(),
        &OpenMagnetics::Settings::get_coil_wind_even_if_not_fit,
        &OpenMagnetics::Settings::set_coil_wind_even_if_not_fit, windEvenIfNotFit);
    OpenMagnetics::SettingsGuard<bool> coatingSquish(
        OpenMagnetics::Settings::GetInstance(),
        &OpenMagnetics::Settings::get_coil_allow_coating_squish,
        &OpenMagnetics::Settings::set_coil_allow_coating_squish, squish);
    OpenMagnetics::SettingsGuard<bool> horizontalOverflow(
        OpenMagnetics::Settings::GetInstance(),
        &OpenMagnetics::Settings::get_coil_allow_horizontal_overflow,
        &OpenMagnetics::Settings::set_coil_allow_horizontal_overflow, horizontalOverflowAllowed);

    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coil(magneticJson.at("coil"), /*windInConstructor=*/false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    OpenMagnetics::Magnetic enriched;
    try {
        enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
    } catch (const std::exception& e) {
        std::cerr << "  " << label << ": autocomplete threw: " << e.what() << "\n";
        return;
    }
    auto enrichedCoil = enriched.get_mutable_coil();
    const auto turns = enrichedCoil.get_turns_description();
    std::cout << "  " << label << ": turns=" << (turns ? turns->size() : 0)
              << " fits=" << enrichedCoil.are_sections_and_layers_fitting()
              << " copperInside=" << enrichedCoil.are_turns_inside_winding_window()
              << " blockingApplied=" << enrichedCoil.is_real_winding_blocking_applied();
    if (enrichedCoil.get_layers_description()) {
        const auto layers = enrichedCoil.get_layers_description().value();
        size_t conduction = 0;
        for (const auto& layer : layers)
            if (layer.get_type() == MAS::ElectricalType::CONDUCTION) ++conduction;
        std::cout << " conductionLayers=" << conduction;
    }
    std::cout << "\n";

    try {
        fs::remove(svgPath);
        OpenMagnetics::Painter painter(svgPath);
        painter.paint_magnetic(enriched, OpenMagnetics::PainterProjection::XY);
        painter.export_svg();
        if (fs::exists(svgPath)) std::cout << "  wrote " << svgPath.string() << "\n";
    } catch (const std::exception& e) {
        std::cerr << "  " << label << ": painter: " << e.what() << "\n";
    }
}

}   // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: mvbpp_paint2d <input.json> <output-stem>\n";
        return 2;
    }
    std::ifstream f(argv[1]);
    if (!f.is_open()) {
        std::cerr << "cannot open " << argv[1] << "\n";
        return 2;
    }
    json j;
    f >> j;
    mvb::patch_dimension_nominals(j);
    const json magneticJson = j.contains("magnetic") ? j.at("magnetic") : j;

    const fs::path stem = argv[2];
    std::cout << fs::path(argv[1]).stem().string() << ":\n";
    paintOne(magneticJson, fs::path(stem.string() + ".refused.svg"), /*windEvenIfNotFit=*/true,
             /*squish=*/false, /*horizontalOverflowAllowed=*/false, "refused");
    paintOne(magneticJson, fs::path(stem.string() + ".overflow.svg"), /*windEvenIfNotFit=*/true,
             /*squish=*/false, /*horizontalOverflowAllowed=*/true, "overflow");
    return 0;
}
