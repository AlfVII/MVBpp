// mvbpp_fitprobe -- ABT #869 diagnosis. Winds a MAS with real winding geometry (coating squish
// and horizontal overflow forced ON so the layout is PRODUCED rather than refused) and prints
// the window, the sections and the layers MKF actually allocated, plus the copper each one has
// to hold. That is the evidence needed to tell "the design data is genuinely over-full" from
// "MKF proportioned the window badly and the design would fit".
#include "mvb/Utils.h"
#include "MAS.hpp"
#include "Magnetic.h"
#include "Utils.h"
#include "constructive_models/Coil.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: mvbpp_fitprobe <input.json>\n";
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
    json magneticJson = j.contains("magnetic") ? j.at("magnetic") : j;

    OpenMagnetics::Magnetic enriched;
    try {
        enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);
    } catch (const std::exception& e) {
        std::cerr << "autocomplete threw: " << e.what() << "\n";
        return 1;
    }
    auto coil = enriched.get_mutable_coil();
    auto bobbin = coil.resolve_bobbin();
    std::cout << std::fixed << std::setprecision(3);
    if (bobbin.get_processed_description()) {
        const auto processed = bobbin.get_processed_description().value();
        const auto windows = processed.get_winding_windows();
        for (const auto& window : windows) {
            if (!window.get_coordinates()) continue;
            std::cout << "window centre (" << (*window.get_coordinates())[0] * 1e3 << ","
                      << (*window.get_coordinates())[1] * 1e3 << ") mm  size "
                      << (window.get_width() ? *window.get_width() * 1e3 : 0) << " x "
                      << (window.get_height() ? *window.get_height() * 1e3 : 0) << " mm\n";
        }
    }
    auto wires = coil.get_wires();
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size();
         ++windingIndex) {
        const auto& winding = coil.get_functional_description()[windingIndex];
        std::cout << "winding " << winding.get_name() << ": turns " << winding.get_number_turns()
                  << " x " << winding.get_number_parallels() << " parallels, wire OD "
                  << wires[windingIndex].get_maximum_outer_width() * 1e3 << " x "
                  << wires[windingIndex].get_maximum_outer_height() * 1e3 << " mm\n";
    }
    if (coil.get_groups_description()) {
        const auto groups = coil.get_groups_description().value();
        for (const auto& group : groups) {
            std::cout << "  group " << std::setw(28) << std::left << group.get_name() << std::right
                      << " x " << group.get_coordinates()[0] * 1e3 << " y "
                      << group.get_coordinates()[1] * 1e3 << "  dims "
                      << group.get_dimensions()[0] * 1e3 << " x " << group.get_dimensions()[1] * 1e3
                      << "\n";
        }
    }
    if (coil.get_sections_description()) {
        // MAS getters return BY VALUE: `for (... : coil.get_sections_description().value())`
        // binds into an optional temporary that dies at the end of the range-init expression.
        // Copy first (the same trap the MKF window checks document).
        const auto sections = coil.get_sections_description().value();
        for (const auto& section : sections) {
            std::cout << "  section " << std::setw(28) << std::left << section.get_name()
                      << std::right << " x " << section.get_coordinates()[0] * 1e3 << " y "
                      << section.get_coordinates()[1] * 1e3 << "  dims "
                      << section.get_dimensions()[0] * 1e3 << " x "
                      << section.get_dimensions()[1] * 1e3 << "  ff "
                      << (section.get_filling_factor() ? *section.get_filling_factor() : -1)
                      << "  layers " << (section.get_layers_orientation() == MAS::WindingOrientation::OVERLAPPING
                                             ? "overlapping" : "contiguous")
                      << "\n";
        }
    }
    if (coil.get_layers_description()) {
        const auto layers = coil.get_layers_description().value();
        for (const auto& layer : layers) {
            std::cout << "    layer " << std::setw(34) << std::left << layer.get_name()
                      << std::right << " x " << layer.get_coordinates()[0] * 1e3 << " y "
                      << layer.get_coordinates()[1] * 1e3 << "  dims "
                      << layer.get_dimensions()[0] * 1e3 << " x " << layer.get_dimensions()[1] * 1e3
                      << "  ff " << (layer.get_filling_factor() ? *layer.get_filling_factor() : -1)
                      << "  turns " << layer.get_turns_alignment().has_value() << "\n";
        }
    }
    if (coil.get_turns_description()) {
        std::map<std::string, size_t> perLayer;
        double minimumX = 1e9, maximumX = -1e9, minimumY = 1e9, maximumY = -1e9;
        const auto turns = coil.get_turns_description().value();
        for (const auto& turn : turns) {
            perLayer[turn.get_layer() ? *turn.get_layer() : std::string("(none)")]++;
            minimumX = std::min(minimumX, turn.get_coordinates()[0]);
            maximumX = std::max(maximumX, turn.get_coordinates()[0]);
            minimumY = std::min(minimumY, turn.get_coordinates()[1]);
            maximumY = std::max(maximumY, turn.get_coordinates()[1]);
        }
        std::cout << "  turns " << turns.size() << " spanning x ["
                  << minimumX * 1e3 << "," << maximumX * 1e3 << "] y [" << minimumY * 1e3 << ","
                  << maximumY * 1e3 << "] mm\n";
        for (const auto& [layerName, count] : perLayer) {
            std::cout << "    " << std::setw(34) << std::left << layerName << std::right << " holds "
                      << count << " turns\n";
        }
    }
    std::cout << "  fits=" << coil.are_sections_and_layers_fitting()
              << " copperInside=" << coil.are_turns_inside_winding_window()
              << " blockingApplied=" << coil.is_real_winding_blocking_applied() << "\n";
    return 0;
}
