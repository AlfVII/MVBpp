#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Symmetry.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include <nlohmann/json.hpp>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <fstream>
#include <string>

using Catch::Matchers::WithinRel;
using json = nlohmann::json;

#ifndef MAS_EXAMPLES_DIR
#define MAS_EXAMPLES_DIR "."
#endif

// ── helpers ──────────────────────────────────────────────────────────────────

static double total_volume(const std::vector<mvb::NamedShape>& shapes) {
    double v = 0.0;
    for (const auto& ns : shapes) {
        if (ns.shape.IsNull()) continue;
        GProp_GProps props;
        BRepGProp::VolumeProperties(ns.shape, props);
        v += props.Mass();
    }
    return v;
}

static MAS::Magnetic load_mas_example(const std::string& filename) {
    std::string path = std::string(MAS_EXAMPLES_DIR) + "/" + filename;
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("Cannot open " + path);
    json j;
    f >> j;
    mvb::patch_dimension_nominals(j);
    return j.at("magnetic").get<MAS::Magnetic>();
}

// ── symmetry detection ────────────────────────────────────────────────────────

TEST_CASE("ETD34 inductor has at least two symmetry planes", "[symmetry]") {
    mvb::MagneticBuilder builder;

    auto shapes = builder.buildAllNamed(load_mas_example("01_simple_inductor_etd34_n87.json"),
                                         /*includeBobbin=*/true, /*symmetryPlanes=*/0);
    REQUIRE_FALSE(shapes.empty());

    auto result = mvb::analyze_symmetry(shapes);
    INFO("Valid planes found: " << result.valid_planes.size());
    CHECK(result.valid_planes.size() >= 2);
}

// ── volume ratios via symmetryPlanes parameter ────────────────────────────────

TEST_CASE("ETD34: symmetryPlanes=1 gives half total volume", "[symmetry]") {
    mvb::MagneticBuilder builder;

    double V_full = total_volume(builder.buildAllNamed(
        load_mas_example("01_simple_inductor_etd34_n87.json"), true, 0));
    REQUIRE(V_full > 0.0);

    double V_half = total_volume(builder.buildAllNamed(
        load_mas_example("01_simple_inductor_etd34_n87.json"), true, 1));
    REQUIRE(V_half > 0.0);

    CHECK_THAT(V_half / V_full, WithinRel(0.5, 0.01));
}

TEST_CASE("ETD34: symmetryPlanes=2 gives quarter total volume", "[symmetry]") {
    mvb::MagneticBuilder builder;

    double V_full = total_volume(builder.buildAllNamed(
        load_mas_example("01_simple_inductor_etd34_n87.json"), true, 0));
    REQUIRE(V_full > 0.0);

    double V_quarter = total_volume(builder.buildAllNamed(
        load_mas_example("01_simple_inductor_etd34_n87.json"), true, 2));
    REQUIRE(V_quarter > 0.0);

    CHECK_THAT(V_quarter / V_full, WithinRel(0.25, 0.01));
}

TEST_CASE("PQ3230 inductor: symmetryPlanes=1 gives half volume", "[symmetry]") {
    mvb::MagneticBuilder builder;

    double V_full = total_volume(builder.buildAllNamed(
        load_mas_example("03_buck_inductor_pq3230_n95.json"), true, 0));
    REQUIRE(V_full > 0.0);

    double V_half = total_volume(builder.buildAllNamed(
        load_mas_example("03_buck_inductor_pq3230_n95.json"), true, 1));
    REQUIRE(V_half > 0.0);

    CHECK_THAT(V_half / V_full, WithinRel(0.5, 0.01));
}

TEST_CASE("Cut shapes are all non-null and non-zero volume", "[symmetry]") {
    mvb::MagneticBuilder builder;

    auto shapes = builder.buildAllNamed(
        load_mas_example("01_simple_inductor_etd34_n87.json"), true, 1);
    REQUIRE_FALSE(shapes.empty());
    for (const auto& ns : shapes) {
        CHECK_FALSE(ns.shape.IsNull());
        GProp_GProps props;
        BRepGProp::VolumeProperties(ns.shape, props);
        CHECK(props.Mass() > 0.0);
    }
}

TEST_CASE("symmetryPlanes=0 returns full assembly", "[symmetry]") {
    mvb::MagneticBuilder builder;

    auto shapes = builder.buildAllNamed(
        load_mas_example("01_simple_inductor_etd34_n87.json"), true, 0);
    // At minimum: 2 core halves + 1 bobbin
    CHECK(shapes.size() >= 3);
}

// ── direct API: cut_to_region ─────────────────────────────────────────────────

TEST_CASE("Explicit cut_to_region gives half volume", "[symmetry]") {
    mvb::MagneticBuilder builder;

    auto full = builder.buildAllNamed(
        load_mas_example("01_simple_inductor_etd34_n87.json"), true, 0);
    double V_full = total_volume(full);
    REQUIRE(V_full > 0.0);

    auto sym = mvb::analyze_symmetry(full);
    REQUIRE_FALSE(sym.valid_planes.empty());

    auto bbox = mvb::aggregate_bbox(full);
    std::vector<std::pair<mvb::SymmetryPlane, mvb::SymmetryHalf>> cuts = {
        {sym.valid_planes[0], mvb::SymmetryHalf::Positive}
    };
    auto cut = mvb::cut_to_region(full, cuts, bbox);
    double V_cut = total_volume(cut);

    CHECK_THAT(V_cut / V_full, WithinRel(0.5, 0.01));
}

TEST_CASE("Explicit cut_to_region two planes gives quarter volume", "[symmetry]") {
    mvb::MagneticBuilder builder;

    auto full = builder.buildAllNamed(
        load_mas_example("01_simple_inductor_etd34_n87.json"), true, 0);
    double V_full = total_volume(full);
    REQUIRE(V_full > 0.0);

    auto sym = mvb::analyze_symmetry(full);
    REQUIRE(sym.valid_planes.size() >= 2);

    auto bbox = mvb::aggregate_bbox(full);
    std::vector<std::pair<mvb::SymmetryPlane, mvb::SymmetryHalf>> cuts = {
        {sym.valid_planes[0], mvb::SymmetryHalf::Positive},
        {sym.valid_planes[1], mvb::SymmetryHalf::Positive},
    };
    auto cut = mvb::cut_to_region(full, cuts, bbox);
    double V_cut = total_volume(cut);

    CHECK_THAT(V_cut / V_full, WithinRel(0.25, 0.01));
}
