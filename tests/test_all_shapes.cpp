// Mirrors MVB Python's test_builder.py::test_all_shapes_generated.
// Iterates over MAS/data/core_shapes.ndjson and verifies every shape
// (excluding ui/pqi/ut, same exclusions as MVB Python) builds a non-empty
// solid via the C++ shape factory.
#include <catch2/catch_test_macros.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "mvb/shapes/ShapeBuilder.h"
#include "MAS.hpp"
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

static const std::set<std::string> EXCLUDED_FAMILIES = {"ui", "pqi", "ut"};

TEST_CASE("All MAS core shapes build a non-empty solid",
          "[shapes][all_shapes]") {
    std::ifstream f("/home/alf/OpenMagnetics/MAS2/data/core_shapes.ndjson");
    REQUIRE(f.is_open());

    int total = 0, skipped = 0, failed = 0;
    std::vector<std::string> failures;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json j = json::parse(line, nullptr, false);
        if (j.is_discarded()) continue;

        std::string family = j.value("family", "");
        std::string name = j.value("name", "?");
        if (EXCLUDED_FAMILIES.count(family)) { ++skipped; continue; }

        // Patch null dim subfields (same as the step generator does)
        mvb::patch_dimension_nominals(j);

        MAS::CoreShape shape;
        try {
            shape = j.get<MAS::CoreShape>();
        } catch (const std::exception& e) {
            failures.push_back(name + " (parse: " + e.what() + ")");
            ++failed; ++total;
            continue;
        }

        auto subtype = shape.get_family_subtype().value_or("");
        auto builder = mvb::shapes::createShapeBuilder(shape.get_family(), subtype);
        if (!builder) {
            failures.push_back(name + " (no builder)");
            ++failed; ++total;
            continue;
        }

        TopoDS_Shape piece;
        try {
            piece = builder->buildPiece(shape);
        } catch (const std::exception& e) {
            failures.push_back(name + " (build: " + e.what() + ")");
            ++failed; ++total;
            continue;
        }

        if (piece.IsNull()) {
            failures.push_back(name + " (null)");
            ++failed; ++total;
            continue;
        }

        GProp_GProps props;
        BRepGProp::VolumeProperties(piece, props);
        if (props.Mass() <= 0.0) {
            failures.push_back(name + " (zero volume)");
            ++failed;
        }
        ++total;
    }

    std::cerr << "[all_shapes] total=" << total << " skipped=" << skipped
              << " failed=" << failed << "\n";
    for (const auto& f : failures) std::cerr << "  FAIL: " << f << "\n";

    REQUIRE(total > 900);          // sanity: should hit ~952 shapes
    REQUIRE(failures.empty());
}

// Mirrors MVB Python's test_get_families: every family in core_shapes.ndjson
// (except ui/pqi) must be routable via the factory.
TEST_CASE("All MAS families are supported by the factory",
          "[shapes][get_families]") {
    std::ifstream f("/home/alf/OpenMagnetics/MAS2/data/core_shapes.ndjson");
    REQUIRE(f.is_open());

    std::set<std::string> families;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        json j = json::parse(line, nullptr, false);
        if (j.is_discarded()) continue;
        families.insert(j.value("family", ""));
    }

    std::vector<std::string> unsupported;
    for (const auto& fam : families) {
        if (fam == "ui" || fam == "pqi") continue;  // MVB Python excludes these
        // Look up via a minimal CoreShape just to test factory routing
        json stub = {{"family", fam}, {"name", "stub"},
                     {"type", "standard"}, {"dimensions", json::object()}};
        MAS::CoreShape s;
        try { s = stub.get<MAS::CoreShape>(); }
        catch (...) { unsupported.push_back(fam + " (parse failed)"); continue; }
        auto b = mvb::shapes::createShapeBuilder(s.get_family(), "");
        if (!b) unsupported.push_back(fam);
    }

    std::cerr << "[get_families] families=" << families.size()
              << " unsupported=" << unsupported.size() << "\n";
    for (const auto& u : unsupported) std::cerr << "  UNSUPPORTED: " << u << "\n";

    REQUIRE(unsupported.empty());
}
