#pragma once

#include <TopoDS_Shape.hxx>
#include <vector>
#include <string>

namespace mvb {

// Export a vector of shapes as a STEP file with per-shape names.
bool exportSTEP(const std::vector<TopoDS_Shape>& shapes,
                const std::vector<std::string>& names,
                const std::string& filepath);

// Export a single compound shape as an STL file.
bool exportSTL(const TopoDS_Shape& compound,
               const std::string& filepath);

} // namespace mvb
