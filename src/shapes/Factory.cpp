#include "mvb/shapes/ShapeBuilder.h"
#include "mvb/shapes/ShapeC.h"
#include "mvb/shapes/ShapeE.h"
#include "mvb/shapes/ShapeEp.h"
#include "mvb/shapes/ShapeLp.h"
#include "mvb/shapes/ShapeEpx.h"
#include "mvb/shapes/ShapeEr.h"
#include "mvb/shapes/ShapeEc.h"
#include "mvb/shapes/ShapeP.h"
#include "mvb/shapes/ShapePm.h"
#include "mvb/shapes/ShapePQ.h"
#include "mvb/shapes/ShapeRM.h"
#include "mvb/shapes/ShapeT.h"
#include "mvb/shapes/ShapeU.h"
#include "mvb/shapes/ShapeUr.h"
#include "mvb/shapes/ShapeUt.h"
#include <memory>

namespace mvb {
namespace shapes {

std::unique_ptr<ShapeBuilder> createShapeBuilder(MAS::CoreShapeFamily family, const std::string& subtype) {
    switch (family) {
        case MAS::CoreShapeFamily::E:
        case MAS::CoreShapeFamily::EFD:
        case MAS::CoreShapeFamily::EI:
        case MAS::CoreShapeFamily::EL:
        case MAS::CoreShapeFamily::ELP:
        case MAS::CoreShapeFamily::PLANAR_E:
        case MAS::CoreShapeFamily::PLANAR_EL:
            return std::make_unique<ShapeE>();
        case MAS::CoreShapeFamily::ER:
        case MAS::CoreShapeFamily::ETD:
        case MAS::CoreShapeFamily::EQ:
        case MAS::CoreShapeFamily::PLANAR_ER:
            return std::make_unique<ShapeEr>();
        case MAS::CoreShapeFamily::T:
            return std::make_unique<ShapeT>();
        case MAS::CoreShapeFamily::U:
        
            return std::make_unique<ShapeU>();
        case MAS::CoreShapeFamily::UR:
            return std::make_unique<ShapeUr>();
        case MAS::CoreShapeFamily::UT:
            return std::make_unique<ShapeUt>();
        case MAS::CoreShapeFamily::C:
            return std::make_unique<ShapeC>();
        case MAS::CoreShapeFamily::P:
            return std::make_unique<ShapeP>();
        case MAS::CoreShapeFamily::PM:
            return std::make_unique<ShapePm>(subtype);
        case MAS::CoreShapeFamily::EP:
            return std::make_unique<ShapeEp>();
        case MAS::CoreShapeFamily::EPX:
            return std::make_unique<ShapeEpx>();
        case MAS::CoreShapeFamily::LP:
            return std::make_unique<ShapeLp>();
        case MAS::CoreShapeFamily::EC:
            return std::make_unique<ShapeEc>();
        case MAS::CoreShapeFamily::PQ:
        case MAS::CoreShapeFamily::PQI:
            return std::make_unique<ShapePQ>();
        case MAS::CoreShapeFamily::RM:
            return std::make_unique<ShapeRM>();
        default:
            return nullptr;
    }
}

} // namespace shapes
} // namespace mvb
