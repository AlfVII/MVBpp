#pragma once

#include "MAS.hpp"
#include "Utils.h"
#include <TopoDS_Shape.hxx>
#include <vector>
#include <string>
#include <map>

// Forward declarations for MKF-enriched magnetic overloads
namespace OpenMagnetics { class Magnetic; class Coil; }

namespace mvb {

struct DrawConfig {
    std::string format = "step";        // "step" or "stl"
    bool separateComponents = true;     // STEP only: core/bobbin/turns as separate roots
    bool includeBobbin = true;
    double machiningOvershoot = 1.0;    // 1.0 = exact MAS dimensions
    double scale = 1.0;                 // Uniform scale applied before export (MVB uses 1000 for mm)
};

class MagneticBuilder {
public:
    std::string drawMagnetic(const MAS::Magnetic& magnetic,
                             const std::string& outputPath,
                             const DrawConfig& config = {}) const;

    // Overload that accepts an already-enriched OpenMagnetics::Magnetic
    // to avoid object-slicing issues with MAS::Magnetic
    std::string drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                             const std::string& outputPath,
                             const DrawConfig& config = {}) const;

    std::vector<TopoDS_Shape> buildCore(const MAS::MagneticCore& core) const;
    TopoDS_Shape buildBobbin(const MAS::Coil& coil, const MAS::MagneticCore& core) const;
    TopoDS_Shape buildBobbin(const OpenMagnetics::Coil& coil, const MAS::MagneticCore& core) const;
    std::vector<TopoDS_Shape> buildTurns(const MAS::Coil& coil, const MAS::MagneticCore& core) const;
    std::vector<TopoDS_Shape> buildTurns(const OpenMagnetics::Coil& coil, const MAS::MagneticCore& core) const;
};

} // namespace mvb
