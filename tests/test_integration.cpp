#include <catch2/catch_test_macros.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "MAS.hpp"
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <numbers>
#include <vector>
#include <map>
#include <filesystem>

static MAS::Dimension make_dim(double v) { return MAS::Dimension(v); }

TEST_CASE("E core with bobbin and turns builds complete assembly", "[integration]") {
    MAS::Magnetic magnetic;

    // Core: E 19/8/5
    MAS::CoreShape shape;
    shape.set_family(MAS::CoreShapeFamily::E);
    shape.set_type(MAS::FunctionalDescriptionType::STANDARD);
    std::map<std::string, MAS::Dimension> dims;
    dims["A"] = make_dim(0.019);
    dims["B"] = make_dim(0.008);
    dims["C"] = make_dim(0.005);
    dims["D"] = make_dim(0.005);
    dims["E"] = make_dim(0.012);
    dims["F"] = make_dim(0.006);
    shape.set_dimensions(dims);

    MAS::CoreGeometricalDescriptionElement piece1;
    piece1.set_type(MAS::CoreGeometricalDescriptionElementType::HALF_SET);
    piece1.set_coordinates({0.0, 0.0, 0.0});
    piece1.set_rotation(std::optional<std::vector<double>>(std::vector<double>{0.0, 0.0, 0.0}));
    piece1.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    MAS::CoreGeometricalDescriptionElement piece2;
    piece2.set_type(MAS::CoreGeometricalDescriptionElementType::HALF_SET);
    piece2.set_coordinates({0.0, 0.0, 0.0});
    piece2.set_rotation(std::optional<std::vector<double>>(std::vector<double>{std::numbers::pi, 0.0, 0.0}));
    piece2.set_shape(std::optional<MAS::CoreShapeDataOrNameUnion>(shape));

    MAS::MagneticCore core;
    core.set_geometrical_description(std::optional<std::vector<MAS::CoreGeometricalDescriptionElement>>(std::vector<MAS::CoreGeometricalDescriptionElement>{piece1, piece2}));
    magnetic.set_core(core);

    // Bobbin processed description
    MAS::CoreBobbinProcessedDescription bobbinPd;
    bobbinPd.set_column_width(0.007);
    bobbinPd.set_column_depth(0.004);
    bobbinPd.set_column_thickness(0.006);
    bobbinPd.set_wall_thickness(0.0005);
    bobbinPd.set_column_shape(MAS::ColumnShape::RECTANGULAR);

    MAS::Bobbin bobbin;
    bobbin.set_processed_description(std::optional<MAS::CoreBobbinProcessedDescription>(bobbinPd));

    // Turns
    std::vector<MAS::Turn> turns;
    for (int i = 0; i < 3; ++i) {
        MAS::Turn turn;
        turn.set_name("Turn_" + std::to_string(i));
        turn.set_winding("Primary");
        turn.set_length(0.05);
        turn.set_parallel(1);
        // radial position, height position
        double radial = 0.005;
        double height = -0.002 + i * 0.0012;
        turn.set_coordinates({radial, height});
        turns.push_back(turn);
    }

    MAS::CoilFunctionalDescription coilFunc;
    coilFunc.set_name("Primary");
    coilFunc.set_number_turns(3);
    coilFunc.set_number_parallels(1);
    coilFunc.set_isolation_side(MAS::IsolationSide::PRIMARY);
    coilFunc.set_wire("Round 1.00 - Grade 1");

    MAS::Coil coil;
    coil.set_bobbin(MAS::BobbinDataOrNameUnion(bobbin));
    coil.set_turns_description(std::optional<std::vector<MAS::Turn>>(turns));
    coil.set_functional_description(std::vector<MAS::CoilFunctionalDescription>{coilFunc});
    magnetic.set_coil(coil);

    mvb::MagneticBuilder builder;
    auto coreShapes = builder.buildCore(magnetic.get_core());
    auto bobbinShape = builder.buildBobbin(magnetic.get_coil(), magnetic.get_core());
    auto turnShapes = builder.buildTurns(magnetic.get_coil(), magnetic.get_core());

    REQUIRE(coreShapes.size() == 2);
    REQUIRE(!bobbinShape.IsNull());
    REQUIRE(turnShapes.size() == 3);

    // Combine and sanity check
    std::vector<TopoDS_Shape> all;
    for (const auto& s : coreShapes) all.push_back(s);
    all.push_back(bobbinShape);
    for (const auto& s : turnShapes) all.push_back(s);

    double totalVolume = 0.0;
    int solids = 0;
    for (const auto& s : all) {
        REQUIRE(!s.IsNull());
        GProp_GProps props;
        BRepGProp::VolumeProperties(s, props);
        totalVolume += props.Mass();
        for (TopExp_Explorer exp(s, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
    }

    REQUIRE(totalVolume > 0.0);
    REQUIRE(solids >= 6); // 2 core + 1 bobbin + 3 turns
}
