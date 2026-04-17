#include "mvb/SectionDrawing.h"
#include "mvb/SectionBuilder.h"
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <GCPnts_TangentialDeflection.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <variant>
#include <vector>

namespace mvb {

namespace {

enum class ViewKind { FRONT, TOP };

struct Dim {
    std::string label;
    double from_x, from_y;   // model-space anchor points (mm)
    double to_x,   to_y;
    bool   horizontal;       // arrow runs horizontally vs vertically
    double offset;           // perpendicular offset from the bracket (model mm)
};

// Edge sampler that projects onto the plane coordinates of the view.
// FrontView: (X, Y) from section plane XY. TopView: (X, Z) from XZ.
std::vector<std::pair<double,double>> sampleEdge(const TopoDS_Edge& edge, ViewKind view) {
    std::vector<std::pair<double,double>> out;
    try {
        BRepAdaptor_Curve curve(edge);
        GCPnts_TangentialDeflection sampler(curve, 0.05, 1e-4);
        out.reserve(sampler.NbPoints());
        for (int i = 1; i <= sampler.NbPoints(); ++i) {
            const gp_Pnt& p = sampler.Value(i);
            if (view == ViewKind::FRONT) {
                out.emplace_back(p.X(), p.Y());
            } else {
                out.emplace_back(p.X(), p.Z());
            }
        }
    } catch (...) { /* degenerate edge */ }
    return out;
}

std::string formatGapLabel(double length_m) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(2);
    if (length_m < 0.0001) s << (length_m * 1e6) << " μm";
    else                    s << (length_m * 1000.0) << " mm";
    return s.str();
}

std::string formatLenLabel(double length_m) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(2) << (length_m * 1000.0);
    return s.str();
}

// Returns the lower-case family classification used by the Python port.
std::string familyKey(MAS::CoreShapeFamily f) {
    switch (f) {
        case MAS::CoreShapeFamily::T:  return "t";
        case MAS::CoreShapeFamily::UR: return "ur";
        case MAS::CoreShapeFamily::UT: return "ut";
        case MAS::CoreShapeFamily::EFD: return "efd";
        case MAS::CoreShapeFamily::EP:  return "ep";
        case MAS::CoreShapeFamily::EPX: return "epx";
        case MAS::CoreShapeFamily::PQ:  return "pq";
        case MAS::CoreShapeFamily::P:   return "p";
        case MAS::CoreShapeFamily::PM:  return "pm";
        case MAS::CoreShapeFamily::RM:  return "rm";
        default: return "e";   // E-family default: E/ETD/ER/EQ/EC/LP/etc., and U/C
    }
}

// Returns a CoreShape by value (MAS getters return by value, so pointers into
// temporaries would dangle). Prefers a shape entry with populated dimensions.
std::optional<MAS::CoreShape> resolveCoreShape(const OpenMagnetics::Magnetic& magnetic) {
    std::optional<MAS::CoreShape> fallback;
    auto consider = [&](const MAS::CoreShape& s) -> std::optional<MAS::CoreShape> {
        if (s.get_dimensions() && !s.get_dimensions()->empty()) return s;
        if (!fallback) fallback = s;
        return std::nullopt;
    };

    const auto& funcShape = magnetic.get_core().get_functional_description().get_shape();
    if (const auto* p = std::get_if<MAS::CoreShape>(&funcShape)) {
        if (auto hit = consider(*p)) return hit;
    }
    auto geomOpt = magnetic.get_core().get_geometrical_description();
    if (geomOpt) {
        for (const auto& entry : *geomOpt) {
            auto shapeOpt = entry.get_shape();
            if (!shapeOpt) continue;
            if (const auto* p = std::get_if<MAS::CoreShape>(&*shapeOpt)) {
                if (auto hit = consider(*p)) return hit;
            }
        }
    }
    return fallback;
}

// Returns a map of CoreShape dims in millimetres (e.g. flat["A"] = 34.2).
std::map<std::string, double> flattenMm(const MAS::CoreShape& shape) {
    std::map<std::string, double> out;
    auto dimsOpt = shape.get_dimensions();
    if (!dimsOpt) return out;
    auto flat = flatten_dimensions(*dimsOpt);
    for (const auto& [k, v] : flat) out[k] = v * 1000.0;
    return out;
}

double getMm(const std::map<std::string, double>& d, const char* k, double dflt = 0.0) {
    auto it = d.find(k);
    return it == d.end() ? dflt : it->second;
}

// ==========================================================================
// FrontView dimensions — core-shape B/D (or C for toroidal).
// ==========================================================================
std::vector<Dim> buildCoreShapeFrontDims(const MAS::CoreShape& shape) {
    std::vector<Dim> dims;
    auto flat = flattenMm(shape);
    if (flat.empty()) return dims;

    const std::string fam = familyKey(shape.get_family());
    double A = getMm(flat, "A"), B = getMm(flat, "B"),
           C = getMm(flat, "C"), D = getMm(flat, "D"), E = getMm(flat, "E");

    if (fam == "t") {
        // Toroidal: C = height along the column axis.
        if (C > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "C = " << C;
            dims.push_back({l.str(), 0.0, -C / 2.0, 0.0, +C / 2.0, false,
                             std::max(A * 0.6, 15.0)});
        }
        return dims;
    }

    // E-family and UR/UT: the assembled core has two half-sets — upper at
    // Y ∈ [0, B] and lower at [-B, 0]. B and D describe the upper half, so
    // anchor their arrows on the upper half only.
    if (B > 0 && A > 0) {
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "B = " << B;
        dims.push_back({l.str(), A / 2.0, 0.0, A / 2.0, +B, false, A * 0.12});
    }
    if (D > 0 && E > 0) {
        // Upper winding window opens at the mating plane y=0 and extends up
        // to y=D (back plate fills y ∈ [D, B]).
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "D = " << D;
        dims.push_back({l.str(), E / 2.0, 0.0, E / 2.0, +D, false, A * 0.05});
    }
    return dims;
}

// ==========================================================================
// FrontView — gap length arrows.
// ==========================================================================
std::vector<Dim> buildGapDims(const std::vector<MAS::CoreGap>& gapping,
                               double shapeA_mm) {
    std::vector<Dim> dims;
    const double halfA_mm = shapeA_mm / 2.0;
    const double breathing_mm = std::max(3.0, shapeA_mm * 0.08);

    // Dedupe visually-identical gap arrows. Two gaps with the same (type,
    // length, y) are typically the symmetric left/right residuals on the
    // lateral columns — one label conveys the information for both.
    std::set<std::tuple<int, long long, long long>> seen;
    auto q = [](double v) -> long long { return static_cast<long long>(std::llround(v * 1e6)); };

    int renderedIdx = 0;
    for (const auto& gap : gapping) {
        auto coordsOpt = gap.get_coordinates();
        double gx_mm = 0.0, gy_mm = 0.0;
        if (coordsOpt && coordsOpt->size() >= 2) {
            gx_mm = (*coordsOpt)[0] * 1000.0;
            gy_mm = (*coordsOpt)[1] * 1000.0;
        }
        double length_m = gap.get_length();
        if (length_m <= 0.0) continue;

        auto sig = std::make_tuple(static_cast<int>(gap.get_type()), q(length_m * 1000.0), q(gy_mm));
        if (!seen.insert(sig).second) continue;

        double len_mm = length_m * 1000.0;
        double targetX_mm = halfA_mm + breathing_mm + renderedIdx * (breathing_mm * 1.2);
        double perpOffset_mm = targetX_mm - gx_mm;
        double y_bot = gy_mm - len_mm / 2.0;
        double y_top = gy_mm + len_mm / 2.0;

        dims.push_back({
            formatGapLabel(length_m),
            gx_mm, y_bot, gx_mm, y_top,
            false, perpOffset_mm
        });
        ++renderedIdx;
    }
    return dims;
}

// ==========================================================================
// FrontView — chunk arrows between gaps, plus top and bottom cap arrows.
// Mirrors freecad_builder.py add_dimensions_and_export_view lines 580-633.
// Chunks are shown for columns that have more than one gap.
// ==========================================================================
std::vector<Dim> buildGapChunkDims(const std::vector<MAS::CoreGap>& gapping,
                                   double columnSemiHeight_mm,
                                   double shapeA_mm) {
    std::vector<Dim> dims;
    if (columnSemiHeight_mm <= 0.0) return dims;

    struct ColGap {
        double y_mm;       // centre of gap along column axis
        double len_mm;     // gap length
        double secHalf_mm; // section half-width (perpendicular to axis)
    };
    // Group gaps by (x, z) column key.
    std::map<std::pair<long long, long long>, std::vector<ColGap>> byColumn;
    auto quantize = [](double v) -> long long {
        return static_cast<long long>(std::llround(v * 1e6));
    };
    for (const auto& gap : gapping) {
        auto coordsOpt = gap.get_coordinates();
        if (!coordsOpt || coordsOpt->size() < 3) continue;
        double length_m = gap.get_length();
        if (length_m <= 0.0) continue;
        double gx_mm = (*coordsOpt)[0] * 1000.0;
        double gy_mm = (*coordsOpt)[1] * 1000.0;
        double gz_mm = (*coordsOpt)[2] * 1000.0;
        double secHalf_mm = 0.0;
        auto secOpt = gap.get_section_dimensions();
        if (secOpt && !secOpt->empty()) secHalf_mm = (*secOpt)[0] * 1000.0 / 2.0;
        byColumn[{quantize(gx_mm), quantize(gz_mm)}].push_back(
            {gy_mm, length_m * 1000.0, secHalf_mm});
    }

    const double inner_breathing_mm = std::max(2.0, shapeA_mm * 0.04);

    for (auto& [key, col] : byColumn) {
        if (col.size() < 2) continue;
        std::sort(col.begin(), col.end(),
                  [](const ColGap& a, const ColGap& b) { return a.y_mm < b.y_mm; });
        double gx_mm = key.first / 1e6;

        auto arrowX = [&](double secHalf) {
            // Sit arrow just beyond the gap's section-half, on the outside.
            double dir = (gx_mm >= 0) ? +1.0 : -1.0;
            return gx_mm + dir * (secHalf + inner_breathing_mm);
        };

        // Between-gap chunks.
        for (size_t i = 0; i + 1 < col.size(); ++i) {
            const auto& g = col[i];
            const auto& n = col[i + 1];
            double y_bot_of_chunk = g.y_mm + g.len_mm / 2.0;   // top of this gap
            double y_top_of_chunk = n.y_mm - n.len_mm / 2.0;   // bottom of next gap
            double chunk = y_top_of_chunk - y_bot_of_chunk;
            if (chunk <= 1e-6) continue;
            double x = arrowX(g.secHalf_mm);
            double perpOffset_mm = x - gx_mm;
            dims.push_back({
                formatLenLabel(chunk / 1000.0),
                gx_mm, y_bot_of_chunk, gx_mm, y_top_of_chunk,
                false, perpOffset_mm
            });
        }

        // Top cap — from the topmost gap to the top of the column.
        {
            const auto& top = col.back();
            double y_gap_top = top.y_mm + top.len_mm / 2.0;
            double y_col_top = +columnSemiHeight_mm;
            double chunk = y_col_top - y_gap_top;
            if (chunk > 1e-6) {
                double x = arrowX(top.secHalf_mm);
                double perpOffset_mm = x - gx_mm;
                dims.push_back({
                    formatLenLabel(chunk / 1000.0),
                    gx_mm, y_gap_top, gx_mm, y_col_top,
                    false, perpOffset_mm
                });
            }
        }
        // Bottom cap — from the bottom of the column to the lowermost gap.
        {
            const auto& bot = col.front();
            double y_col_bot = -columnSemiHeight_mm;
            double y_gap_bot = bot.y_mm - bot.len_mm / 2.0;
            double chunk = y_gap_bot - y_col_bot;
            if (chunk > 1e-6) {
                double x = arrowX(bot.secHalf_mm);
                double perpOffset_mm = x - gx_mm;
                dims.push_back({
                    formatLenLabel(chunk / 1000.0),
                    gx_mm, y_col_bot, gx_mm, y_gap_bot,
                    false, perpOffset_mm
                });
            }
        }
    }
    return dims;
}

// ==========================================================================
// TopView dimensions, per family. Coordinate system: X horizontal, Z vertical
// in the drawing (MVB++'s column axis is Y, perpendicular to the top view).
// Python uses (X, Y) in top view; since its column axis is Z, what Python
// calls "y" here corresponds to MVB++'s "z".
// ==========================================================================
std::vector<Dim> buildCoreShapeTopDims(const MAS::CoreShape& shape) {
    std::vector<Dim> dims;
    auto flat = flattenMm(shape);
    if (flat.empty()) return dims;

    const std::string fam = familyKey(shape.get_family());
    double A  = getMm(flat, "A"),  B  = getMm(flat, "B"),
           C  = getMm(flat, "C"),  D  = getMm(flat, "D"),
           E  = getMm(flat, "E"),  F  = getMm(flat, "F"),
           G  = getMm(flat, "G"),  H  = getMm(flat, "H"),
           J  = getMm(flat, "J"),  K  = getMm(flat, "K"),
           L  = getMm(flat, "L"),  F2 = getMm(flat, "F2");
    (void)B; (void)D;   // only used in FrontView

    const double shape_semi_height = (fam != "p") ? C / 2.0 : A / 2.0;
    double h_offset = std::max(A * 0.15, 5.0);   // running offset, away from right
    double v_offset = std::max(A * 0.15, 5.0);   // running offset, away from top
    double increment = std::max(A * 0.10, 4.0);

    auto pushDistY = [&](double xs, double ys, double xe, double ye,
                          const std::string& label, double off,
                          double /*label_align*/ = 0.0) {
        dims.push_back({label, xs, ys, xe, ye, false, off});
    };
    auto pushDistX = [&](double xs, double ys, double xe, double ye,
                          const std::string& label, double off,
                          double /*label_align*/ = 0.0) {
        dims.push_back({label, xs, ys, xe, ye, true, off});
    };

    // ---------- T (toroidal) ----------
    if (fam == "t") {
        // Python: B then A, both DistanceX centred on origin.
        double sh = A / 2.0;
        if (B > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "B = " << B;
            pushDistX(-B / 2.0, 0.0, +B / 2.0, 0.0, l.str(), v_offset + sh);
            v_offset += increment;
        }
        if (A > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "A = " << A;
            pushDistX(-A / 2.0, 0.0, +A / 2.0, 0.0, l.str(), v_offset + sh);
        }
        return dims;
    }

    // ---------- UR ----------
    if (fam == "ur") {
        double F_ur = F > 0 ? F : C;
        if (C > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "C = " << C;
            pushDistY(A / 2.0, -C / 2.0, A / 2.0, +C / 2.0, l.str(), h_offset);
            h_offset += increment;
        }
        if (G > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "G = " << G;
            pushDistX(-A / 2.0 + F_ur / 2.0 - G / 2.0, 0.0,
                      -A / 2.0 + F_ur / 2.0 + G / 2.0, 0.0,
                      l.str(), v_offset + shape_semi_height);
            pushDistX(A / 2.0 - F_ur / 2.0 - G / 2.0, 0.0,
                      A / 2.0 - F_ur / 2.0 + G / 2.0, 0.0,
                      l.str(), v_offset + shape_semi_height);
            v_offset += increment;
        }
        if (F > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "F = " << F;
            pushDistX(-A / 2.0, 0.0, -A / 2.0 + F, 0.0,
                      l.str(), v_offset + shape_semi_height);
        }
        if (H > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "H = " << H;
            pushDistX(A / 2.0 - H, 0.0, A / 2.0, 0.0,
                      l.str(), v_offset + shape_semi_height);
        }
        double left_col  = F > 0 ? F : C;
        double right_col = H > 0 ? H : C;
        double E_ur = (E > 0) ? E : (A - left_col - right_col);
        if (E_ur > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "E ≥ " << E_ur;
            pushDistX(-A / 2.0 + left_col, 0.0, A / 2.0 - right_col, 0.0,
                      l.str(), v_offset + shape_semi_height);
            v_offset += increment;
        }
        std::ostringstream lA; lA << std::fixed << std::setprecision(2) << "A = " << A;
        pushDistX(-A / 2.0, 0.0, A / 2.0, 0.0, lA.str(), v_offset + shape_semi_height);
        return dims;
    }

    // ---------- UT ----------
    if (fam == "ut") {
        if (C > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "C = " << C;
            pushDistY(A / 2.0, -C / 2.0, A / 2.0, +C / 2.0, l.str(), h_offset);
            h_offset += increment;
        }
        if (F > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "F = " << F;
            pushDistX(-A / 2.0, 0.0, -A / 2.0 + F, 0.0,
                      l.str(), v_offset + shape_semi_height);
        }
        double left_col = (F > 0) ? F : 0.0;
        double right_col = A - E - left_col;
        if (E > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "E = " << E;
            pushDistX(-A / 2.0 + left_col, 0.0, A / 2.0 - right_col, 0.0,
                      l.str(), v_offset + shape_semi_height);
            v_offset += increment;
        }
        std::ostringstream lA; lA << std::fixed << std::setprecision(2) << "A = " << A;
        pushDistX(-A / 2.0, 0.0, A / 2.0, 0.0, lA.str(), v_offset + shape_semi_height);
        return dims;
    }

    // ---------- E-family default (E, ETD, ER, EQ, EC, LP, PLANAR_*, EFD, EP,
    //            EPX, PQ, P, PM, RM, U, C) ----------
    // Order matches Python _e_family_dims TopView branch.
    double correction = 0.0;
    double k = 0.0;

    if (L > 0) {
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "L = " << L;
        pushDistY(J / 2.0, -L / 2.0, J / 2.0, +L / 2.0,
                  l.str(), h_offset + A / 2.0 - J / 2.0);
        h_offset += increment;
    }
    if (K > 0) {
        double height_of_dimension = -C / 2.0;
        double k_val = K;
        if (fam == "efd") { height_of_dimension = +C / 2.0; k_val = -K; }
        else if (fam == "ep") { height_of_dimension = 0.0; k_val = -K; }
        if (K < 0) correction = K / 2.0;
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "K = " << K;
        pushDistY(F / 2.0, height_of_dimension + correction,
                  F / 2.0, height_of_dimension + k_val + correction,
                  l.str(), h_offset + A / 2.0 - F / 2.0);
        h_offset += increment;
    }
    if (F2 > 0) {
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "F2 = " << F2;
        if (fam == "efd") {
            pushDistY(0.0, C / 2.0 + correction - K - F2,
                      0.0, C / 2.0 + correction - K,
                      l.str(), h_offset + A / 2.0);
        } else {
            pushDistY(0.0, -F2 / 2.0, 0.0, +F2 / 2.0,
                      l.str(), h_offset + A / 2.0);
        }
        h_offset += increment;
    }
    if (C > 0) {
        double c_corr = correction;
        if (fam == "ep") c_corr = C / 2.0 - K;
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "C = " << C;
        pushDistY(A / 2.0, -C / 2.0 + c_corr, A / 2.0, +C / 2.0 + c_corr,
                  l.str(), h_offset);
        h_offset += increment;
    }
    if (H > 0) {
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "H = " << H;
        pushDistX(-H / 2.0, 0.0, +H / 2.0, 0.0,
                  l.str(), v_offset + shape_semi_height);
        v_offset += increment;
    }
    if (J > 0 && fam == "pq") {
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "J = " << J;
        pushDistX(-J / 2.0, L / 2.0, +J / 2.0, L / 2.0,
                  l.str(), v_offset + shape_semi_height - L / 2.0);
        v_offset += increment;
    }

    if (fam == "ep" || fam == "epx") k = C / 2.0 - K;
    else if (fam == "efd" && K < 0)   k = -C / 2.0 - K * 2.0;
    else                              k = 0.0;

    if (fam != "p") {
        if (F > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "F = " << F;
            pushDistX(-F / 2.0, -k, +F / 2.0, -k,
                      l.str(), v_offset + k + shape_semi_height);
            v_offset += increment;
        }
        if (G > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "G = " << G;
            pushDistX(-G / 2.0, shape_semi_height, +G / 2.0, shape_semi_height,
                      l.str(), v_offset);
            v_offset += increment;
        }
    } else {
        if (G > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "G = " << G;
            pushDistX(-G / 2.0, E / 2.0, +G / 2.0, E / 2.0,
                      l.str(), v_offset + A / 2.0 - E / 2.0);
            v_offset += increment;
        }
        if (F > 0) {
            std::ostringstream l; l << std::fixed << std::setprecision(2) << "F = " << F;
            pushDistX(-F / 2.0, -k, +F / 2.0, -k,
                      l.str(), v_offset + k + shape_semi_height);
            v_offset += increment;
        }
    }
    if (E > 0) {
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "E = " << E;
        pushDistX(-E / 2.0, -k, +E / 2.0, -k,
                  l.str(), v_offset + k + shape_semi_height);
        v_offset += increment;
    }
    if (A > 0) {
        std::ostringstream l; l << std::fixed << std::setprecision(2) << "A = " << A;
        pushDistX(-A / 2.0, 0.0, +A / 2.0, 0.0,
                  l.str(), v_offset + shape_semi_height);
    }
    return dims;
}

// ==========================================================================
// Build the fused core solid for the given enriched magnetic.
// ==========================================================================
TopoDS_Shape buildFusedCore(const OpenMagnetics::Magnetic& magnetic) {
    MagneticBuilder builder;
    auto pieces = builder.buildCore(magnetic.get_core());
    if (pieces.empty()) {
        throw std::runtime_error("SectionDrawing: no core pieces built");
    }
    TopoDS_Shape core = pieces[0];
    for (size_t i = 1; i < pieces.size(); ++i) {
        BRepAlgoAPI_Fuse fuser(core, pieces[i]);
        if (fuser.IsDone()) core = fuser.Shape();
    }
    return core;
}

// ==========================================================================
// SVG emitter shared between FrontView and TopView.
// ==========================================================================
std::string emitSvg(const std::vector<std::vector<std::pair<double,double>>>& polylines,
                    double xmin, double xmax, double ymin, double ymax,
                    const std::vector<Dim>& dims,
                    double width_px, double label_font_px,
                    const std::string& projection_color,
                    const std::string& dimension_color) {
    // Expand bbox for dim lines / labels.
    double dmxmin = xmin, dmxmax = xmax, dmymin = ymin, dmymax = ymax;
    for (const auto& d : dims) {
        double fx = d.from_x, fy = d.from_y, tx = d.to_x, ty = d.to_y;
        if (d.horizontal) { fy += d.offset; ty += d.offset; }
        else              { fx += d.offset; tx += d.offset; }
        dmxmin = std::min({dmxmin, fx, tx});
        dmxmax = std::max({dmxmax, fx, tx});
        dmymin = std::min({dmymin, fy, ty});
        dmymax = std::max({dmymax, fy, ty});
    }
    double model_w = dmxmax - dmxmin;
    double model_h = dmymax - dmymin;
    double pad_x = model_w * 0.05;
    double pad_y = model_h * 0.08;
    double draw_xmin = dmxmin - pad_x, draw_xmax = dmxmax + pad_x;
    double draw_ymin = dmymin - pad_y, draw_ymax = dmymax + pad_y;
    double draw_w = draw_xmax - draw_xmin;
    double draw_h = draw_ymax - draw_ymin;
    if (draw_w <= 0 || draw_h <= 0) {
        throw std::runtime_error("SectionDrawing: degenerate bbox");
    }
    double margin_px = 30.0;
    double inner_w = width_px - 2 * margin_px;
    double scale = inner_w / draw_w;
    double height_px = draw_h * scale + 2 * margin_px;

    auto toSvg = [&](double x, double y) -> std::pair<double,double> {
        double sx = (x - draw_xmin) * scale + margin_px;
        double sy = height_px - ((y - draw_ymin) * scale + margin_px);
        return {sx, sy};
    };

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(2);
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width_px
        << "\" height=\"" << height_px << "\" viewBox=\"0 0 "
        << width_px << " " << height_px << "\">\n";

    // Cross-section outline.
    svg << "  <g fill=\"none\" stroke=\"" << projection_color
        << "\" stroke-width=\"0.8\" stroke-linecap=\"round\" stroke-linejoin=\"round\">\n";
    for (const auto& poly : polylines) {
        if (poly.size() < 2) continue;
        svg << "    <path d=\"";
        for (size_t i = 0; i < poly.size(); ++i) {
            auto [x, y] = toSvg(poly[i].first, poly[i].second);
            svg << (i == 0 ? "M" : "L") << x << "," << y << " ";
        }
        svg << "\"/>\n";
    }
    svg << "  </g>\n";

    // FreeCAD-style arrows + text.
    constexpr double arrow_size = 6.0;
    constexpr double arrow_length = 15.0;

    svg << "  <g font-family=\"osifont, sans-serif\" font-size=\"" << label_font_px
        << "\" fill=\"" << dimension_color << "\" stroke=\"" << dimension_color
        << "\" stroke-width=\"0.8\" stroke-linecap=\"round\" stroke-linejoin=\"bevel\">\n";

    for (const auto& d : dims) {
        double fx = d.from_x, fy = d.from_y, tx = d.to_x, ty = d.to_y;
        if (d.horizontal) { fy += d.offset; ty += d.offset; }
        else              { fx += d.offset; tx += d.offset; }

        auto [sx1, sy1] = toSvg(fx, fy);
        auto [sx2, sy2] = toSvg(tx, ty);
        auto [ax,  ay ] = toSvg(d.from_x, d.from_y);
        auto [bx,  by ] = toSvg(d.to_x,   d.to_y);

        svg << "    <path fill=\"none\" d=\"";
        svg << "M" << sx1 << "," << sy1 << " L" << sx2 << "," << sy2 << " ";
        svg << "M" << ax  << "," << ay  << " L" << sx1 << "," << sy1 << " ";
        svg << "M" << bx  << "," << by  << " L" << sx2 << "," << sy2;
        svg << "\"/>\n";

        double dxp = sx2 - sx1, dyp = sy2 - sy1;
        double lenp = std::sqrt(dxp * dxp + dyp * dyp);
        if (lenp > 1e-6) {
            dxp /= lenp; dyp /= lenp;
            double px = -dyp, py = dxp;
            auto head = [&](double tipx, double tipy, double sign) {
                double bxp = tipx + sign * dxp * arrow_length;
                double byp = tipy + sign * dyp * arrow_length;
                double l1x = bxp + px * arrow_size, l1y = byp + py * arrow_size;
                double l2x = bxp - px * arrow_size, l2y = byp - py * arrow_size;
                svg << "    <path d=\"M" << tipx << "," << tipy
                    << " L" << l1x << "," << l1y
                    << " L" << l2x << "," << l2y << " Z\"/>\n";
            };
            head(sx1, sy1, +1.0);
            head(sx2, sy2, -1.0);
        }

        double midx = (sx1 + sx2) / 2.0;
        double midy = (sy1 + sy2) / 2.0;
        if (d.horizontal) {
            svg << "    <text x=\"" << midx << "\" y=\"" << (midy - 4)
                << "\" text-anchor=\"middle\" stroke=\"none\">" << d.label
                << "</text>\n";
        } else {
            svg << "    <text x=\"" << midx << "\" y=\"" << midy
                << "\" text-anchor=\"middle\" stroke=\"none\" "
                << "transform=\"rotate(-90 " << midx << " " << midy << ")\" "
                << "dy=\"-4\">" << d.label << "</text>\n";
        }
    }
    svg << "  </g>\n</svg>\n";
    return svg.str();
}

// ==========================================================================
// Main rendering paths.
// ==========================================================================
std::string drawView(const OpenMagnetics::Magnetic& magnetic,
                     ViewKind view,
                     double width_px, double label_font_px,
                     const std::string& projection_color,
                     const std::string& dimension_color) {
    TopoDS_Shape core = buildFusedCore(magnetic);

    SectionPlane plane = (view == ViewKind::FRONT)
        ? SectionPlane::XY    // plane z=0 → X horizontal, Y vertical (column axis)
        : SectionPlane::XZ;   // plane y=0 → X horizontal, Z vertical (depth)
    auto edges = SectionBuilder::sectionCore(core, plane);
    if (edges.IsNull()) throw std::runtime_error("SectionDrawing: section returned nothing");

    // Sample edges → polylines in mm; collect bbox.
    std::vector<std::vector<std::pair<double,double>>> polylines;
    double xmin = +1e30, xmax = -1e30, ymin = +1e30, ymax = -1e30;
    for (TopExp_Explorer exp(edges, TopAbs_EDGE); exp.More(); exp.Next()) {
        auto poly = sampleEdge(TopoDS::Edge(exp.Current()), view);
        if (poly.size() < 2) continue;
        std::vector<std::pair<double,double>> mmPoly;
        mmPoly.reserve(poly.size());
        for (auto& p : poly) {
            double x = p.first * 1000.0, y = p.second * 1000.0;
            mmPoly.emplace_back(x, y);
            if (x < xmin) xmin = x; if (x > xmax) xmax = x;
            if (y < ymin) ymin = y; if (y > ymax) ymax = y;
        }
        polylines.push_back(std::move(mmPoly));
    }
    if (polylines.empty()) throw std::runtime_error("SectionDrawing: no polylines");

    // Resolve CoreShape for labels and pick up A, family.
    auto shapeOpt = resolveCoreShape(magnetic);
    double shapeA_mm = xmax - xmin;
    if (shapeOpt) {
        auto flat = flattenMm(*shapeOpt);
        if (auto it = flat.find("A"); it != flat.end()) shapeA_mm = it->second;
    }

    std::vector<Dim> dims;
    if (shapeOpt) {
        auto core_dims = (view == ViewKind::FRONT)
            ? buildCoreShapeFrontDims(*shapeOpt)
            : buildCoreShapeTopDims(*shapeOpt);
        dims.insert(dims.end(), core_dims.begin(), core_dims.end());
    }

    // Gap annotations live on the FrontView only.
    if (view == ViewKind::FRONT) {
        const auto& gapping = magnetic.get_core().get_functional_description().get_gapping();
        auto gap_dims = buildGapDims(gapping, shapeA_mm);
        dims.insert(dims.end(), gap_dims.begin(), gap_dims.end());

        double colSemi_mm = 0.0;
        auto pd = magnetic.get_core().get_processed_description();
        if (pd && !pd->get_columns().empty()) {
            colSemi_mm = pd->get_columns().front().get_height() * 1000.0 / 2.0;
        }
        auto chunk_dims = buildGapChunkDims(gapping, colSemi_mm, shapeA_mm);
        dims.insert(dims.end(), chunk_dims.begin(), chunk_dims.end());
    }

    return emitSvg(polylines, xmin, xmax, ymin, ymax, dims,
                   width_px, label_font_px, projection_color, dimension_color);
}

} // namespace

std::string SectionDrawing::drawDimensionedFrontView(
        const OpenMagnetics::Magnetic& magnetic,
        double width_px, double label_font_px,
        const std::string& projection_color,
        const std::string& dimension_color) {
    return drawView(magnetic, ViewKind::FRONT, width_px, label_font_px,
                    projection_color, dimension_color);
}

std::string SectionDrawing::drawDimensionedTopView(
        const OpenMagnetics::Magnetic& magnetic,
        double width_px, double label_font_px,
        const std::string& projection_color,
        const std::string& dimension_color) {
    return drawView(magnetic, ViewKind::TOP, width_px, label_font_px,
                    projection_color, dimension_color);
}

void SectionDrawing::writeDimensionedFrontView(
        const OpenMagnetics::Magnetic& magnetic,
        const std::string& outputPath) {
    auto svg = drawDimensionedFrontView(magnetic);
    std::ofstream f(outputPath);
    if (!f.is_open()) throw std::runtime_error("Cannot open " + outputPath);
    f << svg;
}

void SectionDrawing::writeDimensionedTopView(
        const OpenMagnetics::Magnetic& magnetic,
        const std::string& outputPath) {
    auto svg = drawDimensionedTopView(magnetic);
    std::ofstream f(outputPath);
    if (!f.is_open()) throw std::runtime_error("Cannot open " + outputPath);
    f << svg;
}

} // namespace mvb
