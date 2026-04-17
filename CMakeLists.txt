cmake_minimum_required(VERSION 3.18)
project(MVB++ VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

include(FetchContent)
include(ExternalProject)

# ============================================================================
# OpenCASCADE 7.9.3 — build from source, minimal configuration
# ============================================================================
set(OCCT_PREFIX "${CMAKE_BINARY_DIR}/occt-install")

ExternalProject_Add(
    occt_external
    URL https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V7_9_3.tar.gz
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DBUILD_LIBRARY_TYPE=Shared
        -DBUILD_MODULE_Draw=OFF
        -DBUILD_MODULE_DETools=OFF
        -DBUILD_DOC_Overview=OFF
        -DBUILD_PERF_TESTS=OFF
        -DBUILD_TESTING=OFF
        -DBUILD_SAMPLES_MFC=OFF
        -DBUILD_SAMPLES_QT=OFF
        -DBUILD_Inspector=OFF
        -DBUILD_RESOURCES=OFF
        -DUSE_FREETYPE=OFF
        -DUSE_FREEIMAGE=OFF
        -DUSE_TBB=OFF
        -DUSE_OPENVR=OFF
        -DCMAKE_INSTALL_PREFIX=${OCCT_PREFIX}
        -DINSTALL_DIR_LIB=lib
        -DINSTALL_DIR_INCLUDE=include/opencascade
    BUILD_BYPRODUCTS
        ${OCCT_PREFIX}/lib/libTKernel.so
    USES_TERMINAL_BUILD TRUE
)

add_library(occt INTERFACE)
add_dependencies(occt occt_external)
target_include_directories(occt INTERFACE ${OCCT_PREFIX}/include/opencascade)
target_link_directories(occt INTERFACE ${OCCT_PREFIX}/lib)

set(OCCT_ESSENTIAL_LIBS
    TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep
    TKGeomAlgo TKTopAlgo TKPrim TKBO TKBool TKMesh TKShHealing
    TKCDF TKLCAF TKXSBase TKXCAF TKCAF TKVCAF TKBin TKXml TKStd
    TKOffset TKRWMesh TKDESTEP TKDEIGES TKDESTL
)
target_link_libraries(occt INTERFACE ${OCCT_ESSENTIAL_LIBS})

# ============================================================================
# nlohmann/json — reused by both MKF and MVB++
# ============================================================================
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json/
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(json)

# ============================================================================
# Catch2 v3 for tests
# ============================================================================
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.8.0
)
FetchContent_MakeAvailable(Catch2)

# ============================================================================
# MKF — brings MAS.hpp (via quicktype) and physics/autocomplete logic
# ============================================================================
set(MAS_DIRECTORY "${CMAKE_BINARY_DIR}/MAS" CACHE PATH "" FORCE)
set(MAS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../MAS" CACHE PATH "" FORCE)

# Tell MKF not to build its tests or Python interface
set(INCLUDE_MKF_TESTS OFF CACHE BOOL "" FORCE)
set(INCLUDE_PYMKF OFF CACHE BOOL "" FORCE)
set(ENABLE_MATPLOTPP OFF CACHE BOOL "" FORCE)
set(ENABLE_NGSPICE OFF CACHE BOOL "" FORCE)

add_subdirectory(../MKF ${CMAKE_BINARY_DIR}/MKF)

# ============================================================================
# MVB++ library
# ============================================================================
add_library(mvb++ STATIC
    src/MagneticBuilder.cpp
    src/SectionBuilder.cpp
    src/SectionDrawing.cpp
    src/StepExporter.cpp
    src/Utils.cpp
    src/TurnBuilder.cpp
    src/BobbinBuilder.cpp
    src/shapes/ShapeBuilder.cpp
    src/shapes/ShapeE.cpp
    src/shapes/ShapeEp.cpp
    src/shapes/ShapeEpx.cpp
    src/shapes/ShapeLp.cpp
    src/shapes/ShapeEc.cpp
    src/shapes/ShapeEr.cpp
    src/shapes/ShapeP.cpp
    src/shapes/ShapePm.cpp
    src/shapes/ShapePQ.cpp
    src/shapes/ShapeRM.cpp
    src/shapes/ShapeT.cpp
    src/shapes/ShapeU.cpp
    src/shapes/ShapeUr.cpp
    src/shapes/ShapeUt.cpp
    src/shapes/Factory.cpp
)

target_include_directories(mvb++ PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${MAS_DIRECTORY}
    ${json_SOURCE_DIR}/include/nlohmann
    ${json_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/_deps/spline-src/src
    ${CMAKE_BINARY_DIR}/_deps/magic-enum-src/include/magic_enum
    ${CMAKE_BINARY_DIR}/_deps/svg-src/src
)

target_link_libraries(mvb++ PUBLIC
    occt
    MKF
    nlohmann_json::nlohmann_json
)

# ============================================================================
# Tests
# ============================================================================
enable_testing()
add_executable(mvb_tests
    tests/test_main.cpp
    tests/test_core_shapes.cpp
    tests/test_all_shapes.cpp
    tests/test_all_gapped.cpp
    tests/test_step_comparison.cpp
    tests/test_integration.cpp
    tests/test_json_integration.cpp
)

target_compile_definitions(mvb_tests PRIVATE
    REFERENCE_STEPS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/reference_steps"
)

target_link_libraries(mvb_tests PRIVATE mvb++ Catch2::Catch2WithMain)

add_executable(mvbpp_generate_steps tests/generate_mvbpp_steps.cpp)
target_link_libraries(mvbpp_generate_steps PRIVATE mvb++)

add_executable(mvbpp_step_generator tools/mvbpp_step_generator.cpp)
target_link_libraries(mvbpp_step_generator PRIVATE mvb++)

add_executable(mvbpp_section tools/mvbpp_section.cpp)
target_link_libraries(mvbpp_section PRIVATE mvb++)

add_executable(mvbpp_generate_bobbins tools/generate_mvbpp_bobbins.cpp)
target_link_libraries(mvbpp_generate_bobbins PRIVATE occt)

add_test(NAME mvb_tests COMMAND mvb_tests)

# Symlink MVB testData into our build tree for runtime access
add_custom_command(TARGET mvb_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../MVB/tests/testData
        ${CMAKE_CURRENT_BINARY_DIR}/testData
    COMMENT "Linking MVB testData"
)

