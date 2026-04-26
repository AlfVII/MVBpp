cmake_minimum_required(VERSION 3.18)
project(mvbpp_python)

# ── pybind11 ──────────────────────────────────────────────────────────────────
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# ── Python module ─────────────────────────────────────────────────────────────
pybind11_add_module(mvbpp mvbpp_python.cpp)

target_include_directories(mvbpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${json_SOURCE_DIR}/include        # provided by parent if built together
)

target_link_libraries(mvbpp PRIVATE mvb++)

# install alongside the wheel build
install(TARGETS mvbpp DESTINATION .)
