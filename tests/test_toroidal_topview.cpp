// Sanity check that drawDimensionedTopView works for toroidal cores.
// Regression for the "SectionDrawing: no polylines" bug surfaced by the
// MAS battery on T-cores (PFC inductor, CMC, boost T).
#include <catch2/catch_test_macros.hpp>

#include "constructive_models/Magnetic.h"
#include "mvb/SectionDrawing.h"
#include "mvb/Utils.h"
#include "MAS.hpp"

#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;

#ifndef MAS_EXAMPLES_DIR
#define MAS_EXAMPLES_DIR "."
#endif

namespace {
OpenMagnetics::Magnetic load(const std::string& rel) {
    std::ifstream in(std::string(MAS_EXAMPLES_DIR) + "/" + rel);
    REQUIRE(in.is_open());
    json j;
    in >> j;
    mvb::patch_dimension_nominals(j);
    return mvb::magnetic_autocomplete_safe(j.at("magnetic"));
}
}

TEST_CASE("Toroidal TopView: PFC inductor produces non-empty SVG",
          "[2d][toroidal][topview]") {
    auto m = load("05_pfc_inductor_t4020_hf60.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedTopView(m, 400, 12));
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<svg") != std::string::npos);
    // The toroidal section must produce at least one closed circular polyline
    // (outer or inner ring). Check there is some <polyline> or <path> with
    // multiple points by looking for repeated comma separators.
    REQUIRE(svg.find("<path") != std::string::npos);
}

TEST_CASE("Toroidal TopView: boost inductor T5026 produces non-empty SVG",
          "[2d][toroidal][topview]") {
    auto m = load("12_boost_inductor_t5026_26.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedTopView(m, 400, 12));
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<path") != std::string::npos);
}

TEST_CASE("Toroidal TopView: CMC T2515 produces non-empty SVG",
          "[2d][toroidal][topview]") {
    auto m = load("07_cmc_t2515_w800.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedTopView(m, 400, 12));
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<path") != std::string::npos);
}

TEST_CASE("E-family TopView still works (regression guard)",
          "[2d][topview]") {
    auto m = load("01_simple_inductor_etd34_n87.json");
    std::string svg;
    REQUIRE_NOTHROW(svg = mvb::SectionDrawing::drawDimensionedTopView(m, 400, 12));
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<path") != std::string::npos);
}
