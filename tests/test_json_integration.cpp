#include <catch2/catch_test_macros.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs.hxx>

using json = nlohmann::json;

TEST_CASE("Build complete magnetic from concentric_rectangular_column_one_turn.json", "[json][integration][.]" ) {
    std::ifstream f("testData/concentric_rectangular_column_one_turn.json");
    REQUIRE(f.is_open());
    json j;
    try {
        f >> j;
    } catch (const std::exception& e) {
        FAIL("JSON read failed: " << e.what());
    }

    try {
        mvb::patch_dimension_nominals(j);
    } catch (const std::exception& e) {
        FAIL("JSON patch failed: " << e.what());
    }

    MAS::Magnetic magnetic;
    try {
        magnetic = j.at("magnetic").get<MAS::Magnetic>();
    } catch (const std::exception& e) {
        FAIL("MAS parse failed: " << e.what());
    } catch (...) {
        FAIL("MAS parse failed with unknown exception type");
    }

    mvb::MagneticBuilder builder;
    std::vector<TopoDS_Shape> coreShapes;
    TopoDS_Shape bobbinShape;
    std::vector<TopoDS_Shape> turnShapes;
    try {
        coreShapes = builder.buildCore(magnetic.get_core());
        bobbinShape = builder.buildBobbin(magnetic.get_coil(), magnetic.get_core());
        turnShapes = builder.buildTurns(magnetic.get_coil(), magnetic.get_core());
    } catch (const std::exception& e) {
        FAIL("Build failed: " << e.what());
    }

    REQUIRE(!coreShapes.empty());
    REQUIRE(!bobbinShape.IsNull());
    REQUIRE(!turnShapes.empty());

    int solids = 0;
    double totalVolume = 0.0;
    for (const auto& s : coreShapes) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(s, props);
        totalVolume += props.Mass();
        for (TopExp_Explorer exp(s, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    }
    for (const auto& s : {bobbinShape}) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(s, props);
        totalVolume += props.Mass();
        for (TopExp_Explorer exp(s, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    }
    for (const auto& s : turnShapes) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(s, props);
        totalVolume += props.Mass();
        for (TopExp_Explorer exp(s, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    }

    REQUIRE(totalVolume > 0.0);
    REQUIRE(solids >= 3); // core + bobbin + turns
}
