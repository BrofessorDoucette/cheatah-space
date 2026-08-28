/// @file irbem_t89_test.cpp
/// @brief The suite for `space/irbem/ext_t89.hpp` — Tsyganenko (1989).
///
/// An external field model is a fit to data, so "is it right?" cannot be answered by comparing it to
/// an analytic truth: there isn't one. What CAN be answered, and is answered here, is every question
/// about whether the implementation is the model:
///
///  - **`div B = 0`.** The field is the curl of a vector potential, so its divergence is identically
///    zero — for every tilt, every Kp bin and every point. A central-difference divergence therefore
///    has to fall as `h^2` and it does, over four decades. This one test covers every analytic
///    derivative in the file at once (`dW/dx`, `dW/dy`, `dz_s/dx`, `dz_s/dy`, `dD/dx`, `dD/dy`, the
///    whole `B_z` assembly), and a sign error in any of them makes the residual stop scaling.
///  - **The paper's own consistency conditions.** Tsyganenko states that eq. (20)'s `C_16..C_19` are
///    determined from `C_6..C_15` and `dx` by `div B = 0` — four identities per Kp bin, 24 in all.
///    Checking them is a transcription check on Table 1 that cannot pass by accident.
///  - **Exact symmetries.** At zero tilt the model is a mirror about the GSM equator, and at ANY
///    tilt it is a mirror about the noon-midnight meridian. Both hold BITWISE, not to a tolerance,
///    because the arithmetic on the two sides of the mirror is the same arithmetic.
///  - **The physics the model exists for.** The near-tail `B_z` depression deepens monotonically
///    with Kp, which is the paper's own central result (its Figs. 3 and 4).
///  - **The differential against the IRBEM oracle**, across ALL SEVEN Kp bins, when the oracle is
///    present. That test asserts an envelope taken from a measurement, and the measurement says the
///    two do NOT agree closely — see the header's own account of why. Asserting a tight tolerance
///    there would be asserting something false.
///
/// The oracle test `dlopen`s IRBEM at runtime rather than linking it: this binary must build and
/// pass on a machine that has never heard of IRBEM, and `dlopen` lives in libc on this platform so
/// it costs no link flag. `libdl` is not needed and is not asked for.

#include <gtest/gtest.h>

#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

#include "alloc_counter.hpp"
#include "space/irbem/ext_t89.hpp"

namespace {

namespace ir = cheatah::space::irbem;

using ir::Frame;
using ir::Position;
using ir::Status;
using ir::T89Parameters;
using ir::t89_bin_count;
using ir::t89_bin_is_published;
using ir::t89_coefficient_sets;
using ir::t89_components;
using ir::t89_field;
using ir::t89_field_at;
using ir::t89_field_batch;
using ir::t89_field_host;
using ir::t89_fixed;
using ir::t89_linear_count;
using ir::t89_param_block;
using ir::t89_param_count;
using ir::t89_parameters;
using ir::t89_published_set_count;

/// A sink the optimizer cannot see through, so the allocation test's calls actually happen.
volatile double sink = 0.0;

/// A GSM point, spelled so a test reads as coordinates rather than as a constructor call.
Position<Frame::GSM> at(double x, double y, double z) {
    return Position<Frame::GSM>{cheatah::fixarray::vec3d(x, y, z)};
}

/// The tilts the suite sweeps: zero, and both signs of a realistic seasonal-diurnal excursion.
/// The tilt drives every `sin(psi)` and `tan(psi)` term, so a single tilt leaves half the model
/// untested and zero tilt leaves rather more than half.
constexpr std::array<double, 5> kTilts{-0.55, -0.21, 0.0, 0.21, 0.55};

/// `div B` at one point by central differences, with the eq. (20) block optionally switched off.
///
/// Switching it off matters: the published `C_16..C_19` are rounded to four significant figures, so
/// eq. (20) carries an irreducible ~1e-3 nT/R_E of divergence that is a property of the TABLE and
/// would otherwise mask whether the analytic derivatives elsewhere are exact.
double divergence(const T89Parameters<double>& p, double ps, double x, double y, double z,
                  double h) {
    const double s = std::sin(ps);
    const double c = std::cos(ps);
    const auto b = [&](double a, double d, double e) { return t89_components<double>(p, s, c, a, d, e); };
    return ((b(x + h, y, z)[0] - b(x - h, y, z)[0]) + (b(x, y + h, z)[1] - b(x, y - h, z)[1]) +
            (b(x, y, z + h)[2] - b(x, y, z - h)[2])) /
           (2.0 * h);
}

/// The same parameter set with the eq. (20) expansion zeroed out.
T89Parameters<double> without_eq20(const T89Parameters<double>& p) {
    T89Parameters<double> q = p;
    for (std::size_t k = 5; k < t89_linear_count; ++k) q.c[k] = 0.0;
    return q;
}

/// The largest `|div B|` over the sampled box, at difference step @p h.
double worst_divergence(double h, bool with_eq20) {
    double worst = 0.0;
    for (double ps : kTilts) {
        for (int bin = 1; bin <= static_cast<int>(t89_bin_count); ++bin) {
            const T89Parameters<double> full = t89_parameters<double>(bin);
            const T89Parameters<double> p = with_eq20 ? full : without_eq20(full);
            for (double x = -25.0; x <= 12.0; x += 3.7)
                for (double y = -11.0; y <= 11.0; y += 3.3)
                    for (double z = -9.0; z <= 9.0; z += 3.1) {
                        if (std::sqrt((x * x) + (y * y) + (z * z)) < 2.0) continue;
                        worst = std::max(worst, std::fabs(divergence(p, ps, x, y, z, h)));
                    }
        }
    }
    return worst;
}

/// A deterministic scatter of GSM points over the inner magnetosphere and near tail. A 64-bit LCG
/// so a disagreement is reproducible, and never a lattice, so no component is systematically zero.
std::vector<Position<Frame::GSM>> scatter(std::size_t n) {
    std::vector<Position<Frame::GSM>> out;
    out.reserve(n);
    std::uint64_t s = 0x9E3779B97F4A7C15ULL;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) / 9007199254740992.0;
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double r = 2.5 + (17.5 * next());
        const double th = std::acos(1.0 - (2.0 * next()));
        const double ph = 6.283185307179586 * next();
        out.push_back(at(r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph),
                         r * std::cos(th)));
    }
    return out;
}

// ---- the oracle, opened at runtime and never linked --------------------------------------------

/// `get_field1_` and `coord_trans_vec1_`, as the vendored `matlab/libirbem.h` documents them.
using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

/// The oracle handle, or nulls when IRBEM is not on this machine.
struct Oracle {
    void* handle = nullptr;
    GetField1 get_field = nullptr;
    CoordTransVec1 coord_trans = nullptr;

    /// Whether the oracle can be called.
    [[nodiscard]] bool usable() const { return get_field != nullptr && coord_trans != nullptr; }
};

/// Open the oracle once for the process. `CHEATAH_SPACE_IRBEM_ORACLE` overrides the path.
const Oracle& oracle() {
    static const Oracle* const opened = [] {
        auto* o = new Oracle;
        const char* env = std::getenv("CHEATAH_SPACE_IRBEM_ORACLE");
        const std::string path = env != nullptr ? env : "/tmp/irbem-builds/libirbem-O2.so";
        o->handle = dlopen(path.c_str(), RTLD_NOW);
        if (o->handle == nullptr) return o;
        o->get_field = reinterpret_cast<GetField1>(dlsym(o->handle, "get_field1_"));
        o->coord_trans = reinterpret_cast<CoordTransVec1>(dlsym(o->handle, "coord_trans_vec1_"));
        return o;
    }();
    return *opened;
}

}  // namespace

// ================================================================================================
// The published parameters
// ================================================================================================

TEST(IrbemT89, FixedParametersAreThePublishedOnes) {
    // Tsyganenko (1989) section 3: "it was finally decided to fix the following parameters by the
    // values: L_y = 10, D_x = 13, L_RC = 5, L_T = 6.3, gamma_T = 4, delta = 0.01, gamma_1 = 1."
    EXPECT_EQ(t89_fixed.l_y, 10.0);
    EXPECT_EQ(t89_fixed.d_x, 13.0);
    EXPECT_EQ(t89_fixed.l_rc, 5.0);
    EXPECT_EQ(t89_fixed.l_t, 6.3);
    EXPECT_EQ(t89_fixed.gamma_t, 4.0);
    EXPECT_EQ(t89_fixed.delta, 0.01);
    EXPECT_EQ(t89_fixed.gamma_1, 1.0);
    // Section 4: "the non-linear parameters were fixed at values R_T = 30, x_0c = 4, L_xc^2 = 50,
    // D_yc^2 = 20".
    EXPECT_EQ(t89_fixed.r_t, 30.0);
    EXPECT_EQ(t89_fixed.x_0c, 4.0);
    EXPECT_EQ(t89_fixed.l_xc2, 50.0);
    EXPECT_EQ(t89_fixed.d_yc2, 20.0);
    // The literals inside eqs. (11) and (13), named rather than inlined.
    EXPECT_EQ(t89_fixed.hinge_scale2, 16.0);
    EXPECT_EQ(t89_fixed.h1_offset, 16.0);
    EXPECT_EQ(t89_fixed.h1_scale2, 36.0);
}

TEST(IrbemT89, PublishedCoefficientsAreDivergenceFree) {
    // Tsyganenko (1989) after eq. (20): "the last four coefficients C_16-C_19 are not independent,
    // since they are expressed through the first ones in accordance with equation div B = 0".
    // Taking the divergence of eq. (20) and collecting the four independent monomials —
    // z cos(psi), sin(psi), y^2 sin(psi) and z^2 sin(psi) — gives exactly four identities:
    //
    //     C_6/dx  + C_10 + 2 C_16 = 0
    //     C_7/dx  + C_11 +   C_17 = 0
    //     C_8/dx  + 3 C_12 + C_18 = 0
    //     C_9/dx  + C_13 + 3 C_19 = 0
    //
    // Table 1 is printed to four significant figures, so they close to ~1e-3 and not to zero. That
    // residual is the table's rounding and is the reason the tolerance below is 1.5e-3 rather than
    // an epsilon. The measured worst residual over all seven bins is 1.283e-03, so the cap has
    // about 17% of headroom — and that is enough: perturbing the LAST place of one coefficient
    // (bin 3's C_17, -0.6412 -> -0.6421) takes it to 1.521e-03 and fails this test. Any typo in a
    // more significant digit is off by orders of magnitude.
    double worst = 0.0;
    for (std::size_t b = 0; b < t89_bin_count; ++b) {
        const ir::T89Coefficients& s = t89_coefficient_sets[b];
        const std::array<double, 4> identity{
            (s.c[5] / s.delta_x) + s.c[9] + (2.0 * s.c[15]),
            (s.c[6] / s.delta_x) + s.c[10] + s.c[16],
            (s.c[7] / s.delta_x) + (3.0 * s.c[11]) + s.c[17],
            (s.c[8] / s.delta_x) + s.c[12] + (3.0 * s.c[18]),
        };
        for (std::size_t k = 0; k < identity.size(); ++k) {
            EXPECT_NEAR(identity[k], 0.0, 1.5e-3)
                << "Kp bin " << (b + 1) << ", divergence identity " << k;
            worst = std::max(worst, std::fabs(identity[k]));
        }
    }
    std::printf("[ MEASURED ] worst eq.(20) divergence identity residual: %.3e\n", worst);
}

TEST(IrbemT89, BinSevenRepeatsTheMostDisturbedPublishedSet) {
    EXPECT_EQ(t89_published_set_count, 6U);
    EXPECT_EQ(t89_bin_count, 7U);
    for (int bin = 1; bin <= 6; ++bin) EXPECT_TRUE(t89_bin_is_published(bin));
    EXPECT_FALSE(t89_bin_is_published(7));
    EXPECT_FALSE(t89_bin_is_published(0));
    EXPECT_FALSE(t89_bin_is_published(8));

    // Bins 6 and 7 carry the SAME published Kp >= 5- column, bit for bit. This is the documented
    // gap, asserted so it cannot become an accident: if a seventh set is ever obtained from a
    // publishable source, this test fails and says where to look.
    const ir::T89Coefficients& six = t89_coefficient_sets[5];
    const ir::T89Coefficients& seven = t89_coefficient_sets[6];
    for (std::size_t k = 0; k < t89_linear_count; ++k) EXPECT_EQ(six.c[k], seven.c[k]) << "C_" << (k + 1);
    EXPECT_EQ(six.delta_x, seven.delta_x);
    EXPECT_EQ(six.a_rc, seven.a_rc);
    EXPECT_EQ(six.d_0, seven.d_0);
    EXPECT_EQ(six.gamma_rc, seven.gamma_rc);
    EXPECT_EQ(six.r_c, seven.r_c);
    EXPECT_EQ(six.g, seven.g);
    EXPECT_EQ(six.a_t, seven.a_t);
    EXPECT_EQ(six.d_y, seven.d_y);
    EXPECT_EQ(six.x_0, seven.x_0);
}

TEST(IrbemT89, PublishedCoefficientsAreOrderedByDisturbance) {
    // Table 1's own trend, stated in the paper's section 5.1 and worth pinning because it is what a
    // transposed column would break: the current sheet half-thickness D_0 "shows a rapid monotonic
    // decrease with increasing Kp, from D ~ 2.1 for Kp = 0, 0+ up to D ~ 0.3 for Kp >= 5-", and the
    // ring current's scale radius a_RC "decreases monotonically ... from ~8.2 in very quiet up to
    // ~5.8 in the most disturbed conditions".
    for (std::size_t b = 1; b < t89_published_set_count; ++b) {
        EXPECT_LT(t89_coefficient_sets[b].d_0, t89_coefficient_sets[b - 1].d_0) << "D_0 at bin " << (b + 1);
        EXPECT_LE(t89_coefficient_sets[b].a_rc, t89_coefficient_sets[b - 1].a_rc)
            << "a_RC at bin " << (b + 1);
    }
    EXPECT_NEAR(t89_coefficient_sets[0].d_0, 2.08, 1e-12);
    EXPECT_NEAR(t89_coefficient_sets[5].d_0, 0.3325, 1e-12);
    EXPECT_NEAR(t89_coefficient_sets[0].a_rc, 8.161, 1e-12);
    EXPECT_NEAR(t89_coefficient_sets[5].a_rc, 5.831, 1e-12);
}

TEST(IrbemT89, ParametersRoundTripThroughFloat) {
    for (int bin = 1; bin <= static_cast<int>(t89_bin_count); ++bin) {
        const T89Parameters<double> d = t89_parameters<double>(bin);
        const T89Parameters<float> f = t89_parameters<float>(bin);
        for (std::size_t k = 0; k < t89_linear_count; ++k) {
            EXPECT_EQ(f.c[k], static_cast<float>(d.c[k])) << "bin " << bin << " C_" << (k + 1);
        }
        EXPECT_EQ(f.delta_x, static_cast<float>(d.delta_x));
        EXPECT_EQ(f.a_rc, static_cast<float>(d.a_rc));
        EXPECT_EQ(f.d_0, static_cast<float>(d.d_0));
        EXPECT_EQ(f.gamma_rc, static_cast<float>(d.gamma_rc));
        EXPECT_EQ(f.r_c, static_cast<float>(d.r_c));
        EXPECT_EQ(f.g, static_cast<float>(d.g));
        EXPECT_EQ(f.a_t, static_cast<float>(d.a_t));
        EXPECT_EQ(f.d_y, static_cast<float>(d.d_y));
        EXPECT_EQ(f.x_0, static_cast<float>(d.x_0));
    }
    // The fp64 lane must be the table itself, not a rounded copy of it.
    EXPECT_EQ(t89_parameters<double>(1).c[0], -98.72);
    EXPECT_EQ(t89_parameters<double>(1).delta_x, 24.74);
}

TEST(IrbemT89, ParametersClampAnOutOfRangeBin) {
    // t89_field decides the STATUS; this function must never index out of the table, whatever it
    // is handed. Below range clamps to the quietest set, above range to the most disturbed.
    EXPECT_EQ(t89_parameters<double>(0).c[0], t89_parameters<double>(1).c[0]);
    EXPECT_EQ(t89_parameters<double>(-7).delta_x, t89_parameters<double>(1).delta_x);
    EXPECT_EQ(t89_parameters<double>(8).c[0], t89_parameters<double>(7).c[0]);
    EXPECT_EQ(t89_parameters<double>(1000).d_0, t89_parameters<double>(7).d_0);
    // The float instantiation is a separate function and gets the same guarantee.
    EXPECT_EQ(t89_parameters<float>(0).c[0], t89_parameters<float>(1).c[0]);
    EXPECT_EQ(t89_parameters<float>(99).a_t, t89_parameters<float>(7).a_t);
}

TEST(IrbemT89, TheClosureSheetAxisOriginIsADefinedZero) {
    // Each closure sheet's potential has one genuine singular point: the origin of its own
    // coordinate, GSM (0, 0, -+R_T), where both the radius and the offset vanish together. It sits
    // 30 R_E out of the modelling region and is unreachable by any trace, but the function is total
    // there anyway — that sheet contributes nothing and the OTHER one still does, so the answer is
    // a defined, finite field rather than a NaN that only surfaces downstream.
    for (double z : {-t89_fixed.r_t, t89_fixed.r_t}) {
        const std::array<double, 3> b =
            t89_components<double>(t89_parameters<double>(3), 0.3, std::sqrt(1.0 - 0.09), 0.0, 0.0, z);
        EXPECT_TRUE(std::isfinite(b[0])) << "z = " << z;
        EXPECT_TRUE(std::isfinite(b[1])) << "z = " << z;
        EXPECT_TRUE(std::isfinite(b[2])) << "z = " << z;
    }
    // And the fp32 lane, which is the one the device mirrors.
    std::vector<float> pos{0.0F, 0.0F, static_cast<float>(-t89_fixed.r_t)};
    std::vector<float> out(3, 0.0F);
    ASSERT_TRUE(t89_field_host(pos, out, 0.3F, 0.954F, 3));
    for (float v : out) EXPECT_TRUE(std::isfinite(v));
}

// ================================================================================================
// The mathematics: div B = 0, and the symmetries
// ================================================================================================

TEST(IrbemT89, DivergenceVanishesEverywhere) {
    // The field is the curl of a vector potential, so its divergence is identically zero and a
    // central difference of it must fall as h^2. Four decades of h, with eq. (20) switched off so
    // the published table's own rounding does not set a floor: the ratio between successive steps
    // must be ~100, and is checked to be at least 50 to leave room for the fp64 roundoff term
    // (~eps |B| / h) that eventually takes over.
    const double d2 = worst_divergence(1e-2, false);
    const double d3 = worst_divergence(1e-3, false);
    const double d4 = worst_divergence(1e-4, false);
    std::printf("[ MEASURED ] worst |div B| without eq.(20): h=1e-2 %.3e  h=1e-3 %.3e  h=1e-4 %.3e"
                " nT/Re\n",
                d2, d3, d4);
    EXPECT_GT(d2 / d3, 50.0) << "divergence is not falling as h^2 — an analytic derivative is wrong";
    EXPECT_GT(d3 / d4, 50.0) << "divergence is not falling as h^2 — an analytic derivative is wrong";
    EXPECT_LT(d4, 1e-6);

    // With eq. (20) back in, the residual stops falling: it is the four-significant-figure rounding
    // of C_16..C_19 (see PublishedCoefficientsAreDivergenceFree), not a defect, and it is bounded.
    const double f3 = worst_divergence(1e-3, true);
    const double f4 = worst_divergence(1e-4, true);
    std::printf("[ MEASURED ] worst |div B| WITH eq.(20): h=1e-3 %.3e  h=1e-4 %.3e nT/Re\n", f3, f4);
    EXPECT_LT(f4, 3e-3);
    EXPECT_NEAR(f3, f4, 1e-4) << "the eq.(20) residual should be h-independent — it is the table's "
                                 "rounding, not a truncation error";
}

TEST(IrbemT89, ZeroTiltIsMirrorSymmetricAboutTheEquator) {
    // At psi = 0 the warped sheet is flat (tan(psi) = 0 and sin(psi) = 0, both EXACTLY), so the
    // whole model is even about z = 0 in B_z and odd in B_x and B_y. Bitwise, because the two
    // evaluations perform the same operations on operands that differ only in one sign bit.
    for (int bin = 1; bin <= static_cast<int>(t89_bin_count); ++bin) {
        const T89Parameters<double> p = t89_parameters<double>(bin);
        for (double x : {-18.0, -6.5, -1.5, 3.25, 9.0}) {
            for (double y : {-7.5, -1.25, 0.0, 2.5, 8.0}) {
                for (double z : {0.5, 2.25, 6.0}) {
                    const std::array<double, 3> up = t89_components<double>(p, 0.0, 1.0, x, y, z);
                    const std::array<double, 3> dn = t89_components<double>(p, 0.0, 1.0, x, y, -z);
                    EXPECT_EQ(up[0], -dn[0]);
                    EXPECT_EQ(up[1], -dn[1]);
                    EXPECT_EQ(up[2], dn[2]);
                }
            }
        }
    }
}

TEST(IrbemT89, DawnDuskSymmetryHoldsAtEveryTilt) {
    // The model has no dawn-dusk asymmetry at all: every y dependence is through y^2 or y^4 (the
    // sheet's transverse bend, the thickness, the truncation factors) or is linear in y in B_y.
    // So y -> -y gives (B_x, -B_y, B_z) exactly, at ANY tilt. This is the test that would catch a
    // sign slip in dW/dy, dz_s/dy or dD/dy that the divergence test happens not to see.
    for (double ps : kTilts) {
        const double s = std::sin(ps);
        const double c = std::cos(ps);
        for (int bin = 1; bin <= static_cast<int>(t89_bin_count); ++bin) {
            const T89Parameters<double> p = t89_parameters<double>(bin);
            for (double x : {-16.0, -4.5, 2.0, 8.75}) {
                for (double y : {0.75, 3.5, 9.25}) {
                    for (double z : {-5.5, -1.0, 0.0, 4.25}) {
                        const std::array<double, 3> dusk = t89_components<double>(p, s, c, x, y, z);
                        const std::array<double, 3> dawn = t89_components<double>(p, s, c, x, -y, z);
                        EXPECT_EQ(dusk[0], dawn[0]);
                        EXPECT_EQ(dusk[1], -dawn[1]);
                        EXPECT_EQ(dusk[2], dawn[2]);
                    }
                }
            }
        }
    }
}

TEST(IrbemT89, OnTheNoonMidnightMeridianTheFieldStaysInThatPlane) {
    // A corollary of the dawn-dusk symmetry that is worth its own name because it is the property a
    // field-line trace in the meridian plane relies on: at y = 0 exactly, B_y is exactly zero.
    for (double ps : kTilts) {
        const double s = std::sin(ps);
        const double c = std::cos(ps);
        for (int bin = 1; bin <= static_cast<int>(t89_bin_count); ++bin) {
            const T89Parameters<double> p = t89_parameters<double>(bin);
            for (double x : {-20.0, -6.0, 2.5, 10.0}) {
                for (double z : {-6.0, 0.0, 3.5}) {
                    EXPECT_EQ(t89_components<double>(p, s, c, x, 0.0, z)[1], 0.0);
                }
            }
        }
    }
}

// ================================================================================================
// The physics the model exists for
// ================================================================================================

TEST(IrbemT89, NearTailBzDepressionDeepensWithKp) {
    // The paper's central result (its Fig. 4): the external B_z near the midnight point of the
    // geosynchronous orbit falls monotonically with Kp, from about -20 nT at Kp = 0 to about -60 nT
    // at Kp >= 5-. Measured here at x = -6.6 R_E, zero tilt, so the number is directly comparable.
    double previous = 0.0;
    for (int bin = 1; bin <= static_cast<int>(t89_published_set_count); ++bin) {
        const double bz = t89_field_at(at(-6.6, 0.0, 0.0), 0.0, 1.0, bin).v[2];
        std::printf("[ MEASURED ] Kp bin %d: external Bz at GSM x = -6.6 Re is %+8.3f nT\n", bin, bz);
        EXPECT_LT(bz, 0.0) << "the near-tail external field must DEPRESS Bz";
        if (bin > 1) EXPECT_LT(bz, previous) << "the depression must deepen with Kp, bin " << bin;
        previous = bz;
    }
    // The published figure's endpoints, to the width of the plotted symbols.
    EXPECT_NEAR(t89_field_at(at(-6.6, 0.0, 0.0), 0.0, 1.0, 1).v[2], -24.0, 6.0);
    EXPECT_NEAR(t89_field_at(at(-6.6, 0.0, 0.0), 0.0, 1.0, 6).v[2], -60.0, 12.0);
}

TEST(IrbemT89, TheDaysideFieldIsCompressedAndTheNightsideStretched) {
    // The signature of a magnetosphere rather than a dipole: on the dayside the external field ADDS
    // to the northward internal field (the magnetopause currents compress it), on the nightside it
    // subtracts (the tail current stretches it out).
    for (int bin = 1; bin <= static_cast<int>(t89_bin_count); ++bin) {
        EXPECT_GT(t89_field_at(at(8.0, 0.0, 0.0), 0.0, 1.0, bin).v[2], 0.0) << "dayside, bin " << bin;
        EXPECT_LT(t89_field_at(at(-8.0, 0.0, 0.0), 0.0, 1.0, bin).v[2], 0.0) << "nightside, bin " << bin;
    }
}

TEST(IrbemT89, TheWarpedSheetFollowsTheDipoleEquatorNearEarthAndTheGsmEquatorFarOut) {
    // Eq. (11) is the whole geometric content of the paper's title, and it is checked here through
    // the field rather than through the formula: the tail current sheet is where B_x reverses sign,
    // so tracking that zero along the tail IS tracking the sheet. Near the Earth the sheet lies in
    // the DIPOLE equator (z_SM = 0, i.e. z_GSM = x tan(psi)); far down-tail it flattens towards a
    // plane parallel to the GSM equator. So at a tilted epoch the sheet's GSM height must be a
    // smaller fraction of x tan(psi) at x = -25 than at x = -8.
    const double ps = 0.5;
    const double s = std::sin(ps);
    const double c = std::cos(ps);
    const T89Parameters<double> p = t89_parameters<double>(4);
    const auto sheet_height = [&](double x) {
        // Bisect on B_x, which changes sign across the sheet.
        double lo = -12.0;
        double hi = 12.0;
        const auto bx = [&](double z) { return t89_components<double>(p, s, c, x, 0.0, z)[0]; };
        EXPECT_LT(bx(lo) * bx(hi), 0.0) << "no B_x reversal bracketed at x = " << x;
        for (int i = 0; i < 80; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (bx(lo) * bx(mid) <= 0.0) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        return 0.5 * (lo + hi);
    };
    // The dipole equator is z_SM = 0. The verified rotation is z_SM = x_GSM sin(psi) + z_GSM
    // cos(psi), so in GSM that plane is z_GSM = -x_GSM tan(psi): with the north dipole leaning
    // sunward (psi > 0) the dipole equator rises NORTHWARD down-tail, and the sheet must follow it
    // there and flatten towards z_GSM = 0 further out.
    const double near = sheet_height(-8.0);
    const double far = sheet_height(-25.0);
    const double dipole_near = 8.0 * (s / c);
    const double dipole_far = 25.0 * (s / c);
    std::printf("[ MEASURED ] sheet z_GSM at x=-8: %+.3f (dipole equator %+.3f, GSM equator 0); at "
                "x=-25: %+.3f (dipole equator %+.3f)\n",
                near, dipole_near, far, dipole_far);
    EXPECT_GT(near, 0.0) << "at psi > 0 the near-tail sheet must lie NORTH of the GSM equator";
    EXPECT_GT(far, 0.0);
    EXPECT_LT(std::fabs(far / dipole_far), std::fabs(near / dipole_near))
        << "the sheet must flatten towards the GSM equator down-tail — that is the hinge";
    EXPECT_GT(std::fabs(near / dipole_near), 0.5) << "near the Earth the sheet should track the "
                                                     "dipole equator";
}

// ================================================================================================
// The API surface
// ================================================================================================

TEST(IrbemT89, ReferenceLaneMatchesTheComponentForm) {
    for (int bin = 1; bin <= static_cast<int>(t89_bin_count); ++bin) {
        const T89Parameters<double> p = t89_parameters<double>(bin);
        for (double ps : kTilts) {
            const double s = std::sin(ps);
            const double c = std::cos(ps);
            const std::array<double, 3> raw = t89_components<double>(p, s, c, 4.5, -2.25, 1.5);
            const cheatah::fixarray::vec3d wrapped = t89_field_at(at(4.5, -2.25, 1.5), s, c, bin).v;
            EXPECT_EQ(raw[0], wrapped[0]);
            EXPECT_EQ(raw[1], wrapped[1]);
            EXPECT_EQ(raw[2], wrapped[2]);
        }
    }
}

TEST(IrbemT89, KpBinsSelectDistinctCoefficientSets) {
    // Six distinct published sets must produce six distinct fields. Bin 7 repeats bin 6 by design,
    // so it must produce an IDENTICAL field — which is the assertion that keeps the documented gap
    // from being quietly papered over with an interpolation.
    std::vector<cheatah::fixarray::vec3d> fields;
    for (int bin = 1; bin <= static_cast<int>(t89_bin_count); ++bin) {
        fields.push_back(t89_field_at(at(-6.6, 1.0, 0.5), 0.3, std::sqrt(1.0 - 0.09), bin).v);
    }
    for (std::size_t i = 0; i + 1 < t89_published_set_count; ++i) {
        EXPECT_NE(fields[i][2], fields[i + 1][2]) << "bins " << (i + 1) << " and " << (i + 2);
    }
    EXPECT_EQ(fields[5][0], fields[6][0]);
    EXPECT_EQ(fields[5][1], fields[6][1]);
    EXPECT_EQ(fields[5][2], fields[6][2]);
}

TEST(IrbemT89, OnlyKpDrivesTheModel) {
    // T89 reads Kp and NOTHING else. Dst, Pdyn, By and Bz parameterize T96 and later; sweeping them
    // here must leave the answer bit-identical, and asserting that is the honest way to "sweep the
    // storm drivers" for a model that does not read them. A future edit that starts reading one of
    // those slots fails this test rather than silently changing the model.
    const ir::Epoch epoch{2015.5, 43200.0, 2015, 180};
    ir::RotationTable identity{};
    for (cheatah::fixarray::mat3d& m : identity) m = cheatah::fixarray::mat3d::identity();

    const auto field_for = [&](double dst, double pdyn, double by, double bz) {
        ir::DriverSet d{};
        d[static_cast<std::size_t>(ir::Driver::Kp)] = 35.0;
        d[static_cast<std::size_t>(ir::Driver::Dst)] = dst;
        d[static_cast<std::size_t>(ir::Driver::Pdyn)] = pdyn;
        d[static_cast<std::size_t>(ir::Driver::ByIMF)] = by;
        d[static_cast<std::size_t>(ir::Driver::BzIMF)] = bz;
        const ir::ContextResult built = ir::make_field_context(epoch, 0.3, identity, d);
        EXPECT_TRUE(built.has_value()) << ir::describe(built.error());
        return t89_field(at(5.5, -2.0, 1.75), built.value());
    };

    const ir::Result<ir::FieldVector<Frame::GSM>> quiet = field_for(0.0, 2.0, 0.0, 0.0);
    // Dst 0 .. -400, Pdyn 0.5 .. 40, and southward Bz densely: the storm envelope the library
    // exists for. Every one of these must give the same numbers, because T89 cannot see them.
    for (double dst : {0.0, -50.0, -100.0, -200.0, -300.0, -400.0}) {
        for (double pdyn : {0.5, 2.0, 6.0, 12.0, 25.0, 40.0}) {
            for (double bz : {5.0, 0.0, -2.0, -5.0, -10.0, -15.0, -20.0, -30.0, -40.0}) {
                const ir::Result<ir::FieldVector<Frame::GSM>> stormy =
                    field_for(dst, pdyn, -0.5 * bz, bz);
                EXPECT_EQ(stormy.status, quiet.status);
                EXPECT_EQ(stormy.value.v[0], quiet.value.v[0]);
                EXPECT_EQ(stormy.value.v[1], quiet.value.v[1]);
                EXPECT_EQ(stormy.value.v[2], quiet.value.v[2]);
            }
        }
    }
}

TEST(IrbemT89, KpSweepCoversEveryBinAndReportsTheField) {
    // The sweep that DOES change the answer: Kp 0 through 9 in the third-of-a-unit steps Kp really
    // takes, in IRBEM's Kp x 10 scaling. Every bin must be reached, and the field must be a
    // step function of Kp — constant inside a bin, discontinuous between bins.
    std::array<bool, t89_bin_count + 1> seen{};
    double previous_bz = 0.0;
    int previous_bin = 0;
    for (int k = 0; k <= 90; ++k) {
        const double kp10 = static_cast<double>(k);
        const int bin = ir::t89_kp_bin(kp10);
        ASSERT_GE(bin, 1);
        ASSERT_LE(bin, static_cast<int>(t89_bin_count));
        seen[static_cast<std::size_t>(bin)] = true;
        const ir::Result<ir::FieldVector<Frame::GSM>> r = t89_field(at(-6.6, 0.0, 0.0), 0.0, kp10);
        EXPECT_EQ(r.status, Status::Ok) << "Kp x 10 = " << kp10;
        if (bin == previous_bin) {
            EXPECT_EQ(r.value.v[2], previous_bz) << "the field must be constant inside a Kp bin";
        }
        previous_bin = bin;
        previous_bz = r.value.v[2];
    }
    for (std::size_t b = 1; b <= t89_bin_count; ++b) EXPECT_TRUE(seen[b]) << "Kp bin " << b << " unreachable";
}

TEST(IrbemT89, OutOfRangeKpIsReportedButStillEvaluated) {
    // status.hpp's standing rule: OutOfValidityRange never suppresses the value.
    for (double kp10 : {-1.0, -25.0, 91.0, 200.0}) {
        const ir::Result<ir::FieldVector<Frame::GSM>> r = t89_field(at(-6.6, 0.0, 0.0), 0.2, kp10);
        EXPECT_EQ(r.status, Status::OutOfValidityRange) << "Kp x 10 = " << kp10;
        EXPECT_NE(r.value.v[2], 0.0) << "the value must still be computed";
    }
    // And inside the range there is no caveat.
    EXPECT_EQ(t89_field(at(-6.6, 0.0, 0.0), 0.2, 0.0).status, Status::Ok);
    EXPECT_EQ(t89_field(at(-6.6, 0.0, 0.0), 0.2, 90.0).status, Status::Ok);
}

TEST(IrbemT89, PositionOutsideTheFittedRegionIsReported) {
    // T89's published spatial envelope is r_GEO <= 70 R_E (IRBEM's kext table, restating the
    // paper's 4-70 R_E data range). Beyond it the functional form still evaluates and says so.
    const ir::Result<ir::FieldVector<Frame::GSM>> far = t89_field(at(-80.0, 0.0, 0.0), 0.1, 20.0);
    EXPECT_EQ(far.status, Status::OutOfValidityRange);
    EXPECT_TRUE(std::isfinite(far.value.v[2]));
    EXPECT_EQ(t89_field(at(-60.0, 0.0, 0.0), 0.1, 20.0).status, Status::Ok);
    // Inside the solid Earth is a different failure: there is no answer, not an unreliable one.
    EXPECT_EQ(t89_field(at(0.5, 0.0, 0.0), 0.1, 20.0).status, Status::DomainError);
}

TEST(IrbemT89, NonFiniteInputIsADomainError) {
    const double nan = std::nan("");
    const double inf = std::numeric_limits<double>::infinity();
    for (const Position<Frame::GSM>& p :
         {at(nan, 0.0, 0.0), at(0.0, inf, 0.0), at(4.0, 0.0, -nan)}) {
        const ir::Result<ir::FieldVector<Frame::GSM>> r = t89_field(p, 0.2, 20.0);
        EXPECT_EQ(r.status, Status::DomainError);
        EXPECT_EQ(r.value.v[0], 0.0);
        EXPECT_EQ(r.value.v[1], 0.0);
        EXPECT_EQ(r.value.v[2], 0.0);
    }
    EXPECT_EQ(t89_field(at(5.0, 0.0, 0.0), nan, 20.0).status, Status::DomainError);
    EXPECT_EQ(t89_field(at(5.0, 0.0, 0.0), 0.2, nan).status, Status::DomainError);
}

TEST(IrbemT89, RightAngleTiltIsADomainError) {
    // Eq. (11) carries tan(psi). At |psi| = pi/2 it does not exist, so there is no answer — as
    // opposed to an unreliable one. Physically the tilt never exceeds about 35 degrees, but the
    // check is on the mathematics, not on the physics.
    const double quarter_turn = std::numbers::pi / 2.0;
    EXPECT_EQ(t89_field(at(5.0, 0.0, 0.0), quarter_turn, 20.0).status, Status::DomainError);
    EXPECT_EQ(t89_field(at(5.0, 0.0, 0.0), -quarter_turn, 20.0).status, Status::DomainError);
    EXPECT_EQ(t89_field(at(5.0, 0.0, 0.0), 2.0, 20.0).status, Status::DomainError);
    EXPECT_EQ(t89_field(at(5.0, 0.0, 0.0), std::nextafter(quarter_turn, 0.0), 20.0).status,
              Status::Ok);
}

TEST(IrbemT89, AnOverflowingExtrapolationIsADomainErrorNotANaN) {
    // Eq. (20)'s exp(x/dx) overflows a long way outside the model, and the failure it produces is
    // not a big number -- it is a NaN in B_y, because the assembly multiplies the infinite envelope
    // by a y that is exactly zero. B_y == 0 on the noon-midnight meridian is a property another
    // test in this file asserts, so a NaN there is a broken invariant, not an extrapolation. It
    // must be refused at the boundary rather than discovered a hundred RK4 steps downstream.
    for (double x : {1.0e5, 1.0e6, 1.0e30}) {
        const ir::Result<ir::FieldVector<Frame::GSM>> r = t89_field(at(x, 0.0, 0.0), 0.2, 20.0);
        EXPECT_EQ(r.status, Status::DomainError) << "x = " << x;
        EXPECT_EQ(r.value.v[0], 0.0);
        EXPECT_EQ(r.value.v[1], 0.0);
        EXPECT_EQ(r.value.v[2], 0.0);
    }
    // The three components do not overflow together, so each is checked. B_z alone goes first at
    // (0, 0, 1e120), where the eq. (20) cubic z^3 overflows while the quadratic z^2 in B_x does
    // not; B_y and B_z go together at (0, 1e60, 1e150), where the y z^2 term overflows and B_x's
    // largest term is still only z^2. A check that only looked at B_x would pass both of these.
    for (const Position<Frame::GSM>& p : {at(0.0, 0.0, 1.0e120), at(0.0, 1.0e60, 1.0e150)}) {
        const std::array<double, 3> raw =
            t89_components<double>(t89_parameters<double>(4), 0.3, 0.954, p.v[0], p.v[1], p.v[2]);
        EXPECT_TRUE(std::isfinite(raw[0])) << "B_x must be the component that still fits";
        EXPECT_FALSE(std::isfinite(raw[0]) && std::isfinite(raw[1]) && std::isfinite(raw[2]));
        const ir::Result<ir::FieldVector<Frame::GSM>> r = t89_field(p, 0.3, 20.0);
        EXPECT_EQ(r.status, Status::DomainError);
        EXPECT_EQ(r.value.v[2], 0.0);
    }

    // Far out but still representable stays what it was: an extrapolation, reported and returned.
    const ir::Result<ir::FieldVector<Frame::GSM>> big = t89_field(at(1.0e3, 0.0, 0.0), 0.2, 20.0);
    EXPECT_EQ(big.status, Status::OutOfValidityRange);
    EXPECT_TRUE(std::isfinite(big.value.v[2]));
    EXPECT_NE(big.value.v[2], 0.0);
    // And the anti-sunward direction never overflows at all -- the envelope decays that way.
    EXPECT_EQ(t89_field(at(-1.0e6, 0.0, 0.0), 0.2, 20.0).status, Status::OutOfValidityRange);
}

TEST(IrbemT89, BatchReportsTheSameEnvelopeTheScalarLaneDoes) {
    // One status for N points can only be the worst of them, but it must not be BETTER than the
    // worst: a batch that quietly returns Ok for a point the scalar entry point refuses is a
    // silent hole, and it is the one a drift-shell sweep would fall into.
    const std::vector<Position<Frame::GSM>> good{at(5.0, 1.0, 1.0), at(-8.0, 2.0, -1.0)};
    std::vector<ir::FieldVector<Frame::GSM>> out(good.size());
    EXPECT_EQ(t89_field_batch(good, 0.2, 20.0, out).status, Status::Ok);

    // One point beyond r_GEO = 70 R_E: the whole batch is out of validity, every point computed.
    const std::vector<Position<Frame::GSM>> far{at(5.0, 1.0, 1.0), at(-80.0, 0.0, 0.0)};
    EXPECT_EQ(t89_field_batch(far, 0.2, 20.0, out).status, Status::OutOfValidityRange);
    EXPECT_EQ(t89_field(far[1], 0.2, 20.0).status, Status::OutOfValidityRange);
    for (const ir::FieldVector<Frame::GSM>& b : out) EXPECT_NE(b.v[2], 0.0);

    // One point inside the Earth, and one not finite: a domain error, and every output zeroed,
    // exactly as the scalar lane zeroes its one.
    for (const Position<Frame::GSM>& bad :
         {at(0.5, 0.0, 0.0), at(std::nan(""), 0.0, 0.0),
          at(std::numeric_limits<double>::infinity(), 0.0, 0.0)}) {
        const std::vector<Position<Frame::GSM>> mixed{at(5.0, 1.0, 1.0), bad};
        std::vector<ir::FieldVector<Frame::GSM>> mixed_out(mixed.size(), ir::FieldVector<Frame::GSM>{
                                                                             cheatah::fixarray::vec3d{
                                                                                 1.0, 1.0, 1.0}});
        const ir::Result<bool> r = t89_field_batch(mixed, 0.2, 20.0, mixed_out);
        EXPECT_EQ(r.status, Status::DomainError);
        EXPECT_EQ(t89_field(bad, 0.2, 20.0).status, Status::DomainError);
        for (const ir::FieldVector<Frame::GSM>& b : mixed_out) {
            EXPECT_EQ(b.v[0], 0.0);
            EXPECT_EQ(b.v[1], 0.0);
            EXPECT_EQ(b.v[2], 0.0);
        }
    }

    // Two bad points, the first one first: the fold must stay poisoned across the rest of the scan
    // rather than being re-decided by the last point it looks at.
    const std::vector<Position<Frame::GSM>> two_bad{at(std::nan(""), 0.0, 0.0),
                                                    at(std::numeric_limits<double>::infinity(), 1.0,
                                                       0.0),
                                                    at(5.0, 1.0, 1.0)};
    std::vector<ir::FieldVector<Frame::GSM>> two_out(two_bad.size());
    EXPECT_EQ(t89_field_batch(two_bad, 0.2, 20.0, two_out).status, Status::DomainError);

    // An out-of-range Kp is still reported on an EMPTY batch: nothing to compute is not a reason to
    // drop the caveat the same call would carry with one point in it.
    const std::vector<Position<Frame::GSM>> none;
    std::vector<ir::FieldVector<Frame::GSM>> none_out;
    EXPECT_EQ(t89_field_batch(none, 0.2, 900.0, none_out).status, Status::OutOfValidityRange);
    EXPECT_EQ(t89_field_batch(none, 0.2, 20.0, none_out).status, Status::Ok);
}

TEST(IrbemT89, ContextOverloadAgreesWithTheExplicitOne) {
    const ir::Epoch epoch{2015.5, 43200.0, 2015, 180};
    ir::RotationTable identity{};
    for (cheatah::fixarray::mat3d& m : identity) m = cheatah::fixarray::mat3d::identity();
    ir::DriverSet drivers{};
    drivers[static_cast<std::size_t>(ir::Driver::Kp)] = 47.0;
    const ir::ContextResult built = ir::make_field_context(epoch, -0.42, identity, drivers);
    ASSERT_TRUE(built.has_value()) << ir::describe(built.error());

    const Position<Frame::GSM> p = at(3.75, -1.5, 2.25);
    const ir::Result<ir::FieldVector<Frame::GSM>> viaContext = t89_field(p, built.value());
    const ir::Result<ir::FieldVector<Frame::GSM>> viaScalars = t89_field(p, -0.42, 47.0);
    EXPECT_EQ(viaContext.status, viaScalars.status);
    EXPECT_EQ(viaContext.value.v[0], viaScalars.value.v[0]);
    EXPECT_EQ(viaContext.value.v[1], viaScalars.value.v[1]);
    EXPECT_EQ(viaContext.value.v[2], viaScalars.value.v[2]);
}

// ================================================================================================
// The batch lanes
// ================================================================================================

TEST(IrbemT89, HostFloatLaneTracksTheReferenceLane) {
    // The fp32 host lane exists to mirror the device kernel, so what is measured here is the fp32
    // vs fp64 gap of the SAME expressions — the arithmetic floor the device is then held to.
    const std::vector<Position<Frame::GSM>> pts = scatter(4096);
    std::vector<float> pos(3 * pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    std::vector<float> out(3 * pts.size());
    const double ps = 0.35;
    ASSERT_TRUE(t89_field_host(pos, out, static_cast<float>(std::sin(ps)),
                               static_cast<float>(std::cos(ps)), 4));

    double worst_abs = 0.0;
    double worst_rel = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        // The fp64 reference is evaluated at the SAME fp32-rounded position, so the comparison
        // measures the arithmetic and not the input rounding.
        const std::array<double, 3> ref = t89_components<double>(
            t89_parameters<double>(4), std::sin(ps), std::cos(ps), pos[(3 * i) + 0],
            pos[(3 * i) + 1], pos[(3 * i) + 2]);
        double d2 = 0.0;
        double r2 = 0.0;
        for (int c = 0; c < 3; ++c) {
            const double e = out[(3 * i) + c] - ref[static_cast<std::size_t>(c)];
            d2 += e * e;
            r2 += ref[static_cast<std::size_t>(c)] * ref[static_cast<std::size_t>(c)];
        }
        worst_abs = std::max(worst_abs, std::sqrt(d2));
        worst_rel = std::max(worst_rel, std::sqrt(d2) / (std::sqrt(r2) + 1e-9));
    }
    std::printf("[ MEASURED ] fp32 host lane vs fp64 reference: max |dB| %.3e nT, max relative "
                "%.3e\n",
                worst_abs, worst_rel);
    EXPECT_LT(worst_abs, 1e-2);
}

TEST(IrbemT89, HostFloatLaneRejectsMismatchedSpans) {
    std::vector<float> pos(7, 0.0F);   // not a whole number of points
    std::vector<float> out(7, 0.0F);
    EXPECT_FALSE(t89_field_host(pos, out, 0.0F, 1.0F, 1));
    std::vector<float> pos6(6, 1.0F);
    std::vector<float> out3(3, 0.0F);
    EXPECT_FALSE(t89_field_host(pos6, out3, 0.0F, 1.0F, 1));
    std::vector<float> out6(6, 0.0F);
    EXPECT_TRUE(t89_field_host(pos6, out6, 0.0F, 1.0F, 1));
}

TEST(IrbemT89, ParameterBlockCarriesTheTiltThenTheCoefficients) {
    const std::array<float, t89_param_count> block = t89_param_block(0.25F, 0.75F, 3);
    EXPECT_EQ(block.size(), 30U);
    EXPECT_EQ(block[0], 0.25F);
    EXPECT_EQ(block[1], 0.75F);
    const T89Parameters<float> p = t89_parameters<float>(3);
    for (std::size_t k = 0; k < t89_linear_count; ++k) EXPECT_EQ(block[2 + k], p.c[k]) << "C_" << (k + 1);
    EXPECT_EQ(block[21], p.delta_x);
    EXPECT_EQ(block[22], p.a_rc);
    EXPECT_EQ(block[23], p.d_0);
    EXPECT_EQ(block[24], p.gamma_rc);
    EXPECT_EQ(block[25], p.r_c);
    EXPECT_EQ(block[26], p.g);
    EXPECT_EQ(block[27], p.a_t);
    EXPECT_EQ(block[28], p.d_y);
    EXPECT_EQ(block[29], p.x_0);
}

TEST(IrbemT89, BatchAgreesWithTheReferenceLane) {
    const std::vector<Position<Frame::GSM>> pts = scatter(1000);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    const ir::Result<bool> r = t89_field_batch(pts, 0.28, 35.0, out);
    EXPECT_EQ(r.status, Status::Ok);
    // Whether the device ran is reported rather than assumed. Without a device stack it is false,
    // and the host lane below must then be bit-identical to the reference.
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const ir::FieldVector<Frame::GSM> ref =
            t89_field_at(pts[i], std::sin(0.28), std::cos(0.28), ir::t89_kp_bin(35.0));
        if (r.value) {
            for (int c = 0; c < 3; ++c) {
                EXPECT_NEAR(out[i].v[c], ref.v[c], 1e-2) << "device lane, point " << i;
            }
        } else {
            EXPECT_EQ(out[i].v[0], ref.v[0]);
            EXPECT_EQ(out[i].v[1], ref.v[1]);
            EXPECT_EQ(out[i].v[2], ref.v[2]);
        }
    }
}

TEST(IrbemT89, BatchRejectsMismatchedSpans) {
    const std::vector<Position<Frame::GSM>> pts = scatter(4);
    std::vector<ir::FieldVector<Frame::GSM>> shorter(3);
    EXPECT_EQ(t89_field_batch(pts, 0.1, 20.0, shorter).status, Status::DomainError);
    std::vector<ir::FieldVector<Frame::GSM>> right(4);
    EXPECT_EQ(t89_field_batch(pts, std::nan(""), 20.0, right).status, Status::DomainError);
    EXPECT_EQ(t89_field_batch(pts, 0.1, std::nan(""), right).status, Status::DomainError);
    EXPECT_EQ(t89_field_batch(pts, std::numbers::pi / 2.0, 20.0, right).status, Status::DomainError);
    // An empty batch is a no-op, not an error, and never touches a device.
    const ir::Result<bool> empty = t89_field_batch({}, 0.1, 20.0, {});
    EXPECT_EQ(empty.status, Status::Ok);
    EXPECT_FALSE(empty.value);
    // An out-of-range Kp is carried through the batch the same way it is through the scalar entry.
    EXPECT_EQ(t89_field_batch(pts, 0.1, 300.0, right).status, Status::OutOfValidityRange);
}

TEST(IrbemT89, NothingOnTheHeapInTheHotPath) {
    // The claim is "no heap in a hot path", so it gets a counter. The SECOND call is the one that
    // matters: a routine that allocates a workspace per invocation passes a single-call check.
    const std::vector<Position<Frame::GSM>> pts = scatter(256);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    (void)t89_field_batch(pts, 0.2, 30.0, out);
    (void)t89_field(at(5.0, 1.0, 1.0), 0.2, 30.0);

    const std::size_t before = cheatah_space_test::allocation_count();
    for (int i = 0; i < 64; ++i) {
        sink = sink + t89_field(at(5.0 + 0.01 * i, 1.0, 1.0), 0.2, 30.0).value.v[2];
    }
    (void)t89_field_batch(pts, 0.2, 30.0, out);   // below the crossover: the host lane
    sink = sink + out[0].v[2];
    EXPECT_EQ(before, cheatah_space_test::allocation_count());
}

// ================================================================================================
// The differential against the IRBEM oracle
// ================================================================================================

TEST(IrbemT89, DiffersFromTheIrbemOracleByTheMeasuredEnvelope) {
    const Oracle& o = oracle();
    if (!o.usable()) {
        GTEST_SKIP() << "IRBEM oracle not present (set CHEATAH_SPACE_IRBEM_ORACLE to its .so); "
                        "the oracle is a dev-only black box and is never linked";
    }
    // Three epochs spanning the tilt range. The external field is isolated as kext=4 minus kext=0
    // so the internal IGRF term cancels EXACTLY, and the tilt is read out of the oracle itself
    // (the SM z axis in GSM is (sin psi, 0, cos psi)) so a frame difference cannot masquerade as a
    // model difference.
    struct Epoch {
        int doy;
        double ut;
    };
    const std::array<Epoch, 3> epochs{{{80, 39183.0}, {180, 43200.0}, {355, 7200.0}}};
    const std::array<double, 7> kp10{{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0}};

    // The envelope below is a MEASUREMENT, not an agreement target. space.irbem implements
    // Tsyganenko (1989) as published; IRBEM's kext=4 evaluates the later "T89c" revision, whose
    // parameters were never published and whose source is GPL-3.0 and therefore unreadable here.
    // tools/oracle/t89_diff.cpp measures the gap and shows it is structural rather than a re-fit.
    // These caps sit ~1.4x above the measured RMS so a regression in THIS implementation fails
    // while the known, documented model difference does not.
    const std::array<double, 7> rms_cap{{3.8, 4.5, 6.0, 8.7, 16.0, 21.3, 55.3}};

    for (std::size_t bi = 0; bi < kp10.size(); ++bi) {
        double sum2 = 0.0;
        double sig2 = 0.0;
        std::size_t n = 0;
        for (const Epoch& e : epochs) {
            int iyear = 2015;
            int idoy = e.doy;
            double ut = e.ut;
            int one = 1;
            double ps = 0.0;
            {
                int si = 4;
                int so = 2;
                std::array<double, 3> in{0.0, 0.0, 1.0};
                std::array<double, 3> outv{};
                o.coord_trans(&one, &si, &so, &iyear, &idoy, &ut, in.data(), outv.data());
                ps = std::atan2(outv[0], outv[2]);
            }
            for (double x = -30.0; x <= 14.0; x += 4.0)
                for (double y = -12.0; y <= 12.0; y += 4.0)
                    for (double z = -8.0; z <= 8.0; z += 4.0) {
                        const double r = std::sqrt((x * x) + (y * y) + (z * z));
                        if (r < 3.0 || r > 35.0) continue;
                        std::array<double, 3> gsm{x, y, z};
                        std::array<double, 3> geo{};
                        {
                            int si = 2;
                            int so = 1;
                            o.coord_trans(&one, &si, &so, &iyear, &idoy, &ut, gsm.data(), geo.data());
                        }
                        std::array<int, 5> options{0, 0, 0, 0, 0};
                        int sysaxes = 1;
                        int k0 = 0;
                        int k4 = 4;
                        std::vector<double> mag(25, 0.0);
                        mag[0] = kp10[bi];
                        std::array<double, 3> b0{};
                        std::array<double, 3> b4{};
                        double m0 = 0.0;
                        double m4 = 0.0;
                        double x1 = geo[0];
                        double x2 = geo[1];
                        double x3 = geo[2];
                        o.get_field(&k0, options.data(), &sysaxes, &iyear, &idoy, &ut, &x1, &x2, &x3,
                                    mag.data(), b0.data(), &m0);
                        o.get_field(&k4, options.data(), &sysaxes, &iyear, &idoy, &ut, &x1, &x2, &x3,
                                    mag.data(), b4.data(), &m4);
                        std::array<double, 3> dgeo{b4[0] - b0[0], b4[1] - b0[1], b4[2] - b0[2]};
                        std::array<double, 3> ora{};
                        {
                            int si = 1;
                            int so = 2;
                            o.coord_trans(&one, &si, &so, &iyear, &idoy, &ut, dgeo.data(),
                                          ora.data());
                        }
                        const ir::FieldVector<Frame::GSM> mine = t89_field_at(
                            at(x, y, z), std::sin(ps), std::cos(ps), static_cast<int>(bi) + 1);
                        for (int c = 0; c < 3; ++c) {
                            const double d = mine.v[c] - ora[static_cast<std::size_t>(c)];
                            sum2 += d * d;
                            sig2 += ora[static_cast<std::size_t>(c)] * ora[static_cast<std::size_t>(c)];
                        }
                        ++n;
                    }
        }
        ASSERT_GT(n, 0U);
        const double rms = std::sqrt(sum2 / static_cast<double>(n));
        std::printf("[ MEASURED ] Kp bin %zu: RMS |dB| vs IRBEM kext=4 = %7.3f nT (%.1f%% of the "
                    "external field, %zu points)%s\n",
                    bi + 1, rms, 100.0 * std::sqrt(sum2 / sig2), n,
                    t89_bin_is_published(static_cast<int>(bi) + 1)
                        ? ""
                        : "  <- no published coefficient set for this bin");
        EXPECT_LT(rms, rms_cap[bi]) << "Kp bin " << (bi + 1)
                                    << ": the gap against the oracle has GROWN beyond the "
                                       "documented T89-vs-T89c difference";
        // The two are the same model family, so they must at least agree in sign and order of
        // magnitude everywhere the field is not near a null — a much weaker claim than the RMS cap
        // and one no re-parameterization could break.
        EXPECT_GT(std::sqrt(sig2), 0.0);
        EXPECT_LT(std::sqrt(sum2 / sig2), 0.75) << "Kp bin " << (bi + 1);
    }
}

#if CHEATAH_SPACE_IRBEM_T89_GPU
namespace {

/// Point `CHEATAH_SPACE_IRBEM_SPV_DIR` somewhere for the life of the object, and put it back —
/// value or absence — afterwards. The seam resolves the shader directory on every launch, so this
/// is how the "the shaders were never built" path is reached on a machine that HAS them.
class SpvDirScope {
  public:
    /// @param dir the directory to advertise as the shader directory.
    explicit SpvDirScope(const std::string& dir) {
        if (const char* prev = std::getenv(kVar)) {
            had_ = true;
            prev_ = prev;
        }
        ::setenv(kVar, dir.c_str(), 1);
    }
    SpvDirScope(const SpvDirScope&) = delete;
    SpvDirScope& operator=(const SpvDirScope&) = delete;
    SpvDirScope(SpvDirScope&&) = delete;
    SpvDirScope& operator=(SpvDirScope&&) = delete;
    /// Restore what was there before.
    ~SpvDirScope() {
        if (had_) {
            ::setenv(kVar, prev_.c_str(), 1);
            return;
        }
        ::unsetenv(kVar);
    }

  private:
    /// The variable the seam consults first.
    static constexpr const char* kVar = "CHEATAH_SPACE_IRBEM_SPV_DIR";
    /// Whether it was set before this scope began.
    bool had_ = false;
    /// Its previous value, meaningful only when @ref had_.
    std::string prev_;
};

}  // namespace

TEST(IrbemT89, BatchFallsBackToTheHostWhenTheShaderWasNeverBuilt) {
    // A missing .spv is a deployment problem, not a reason to refuse to compute — so the batch must
    // quietly run the host lane and SAY it did. Reached here by pointing the seam at a directory
    // that has no shaders in it, which is the same thing a half-finished build looks like.
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_t89_f32");
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    {
        const SpvDirScope nowhere(std::filesystem::temp_directory_path().string() +
                                  "/cheatah-space-no-such-shader-dir");
        const ir::Result<bool> r = t89_field_batch(pts, 0.28, 35.0, out);
        EXPECT_EQ(r.status, Status::Ok);
        EXPECT_FALSE(r.value) << "with no compiled shader the batch must run on the host";
    }
    // And the host lane it fell back to is the fp64 reference, bit for bit.
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref =
            t89_field_at(pts[i], std::sin(0.28), std::cos(0.28), ir::t89_kp_bin(35.0));
        ASSERT_EQ(out[i].v[0], ref.v[0]) << "point " << i;
        ASSERT_EQ(out[i].v[1], ref.v[1]) << "point " << i;
        ASSERT_EQ(out[i].v[2], ref.v[2]) << "point " << i;
    }
}

TEST(IrbemT89, BatchUsesTheDeviceWhenOneIsAvailable) {
    // The batch entry point reports whether the device serviced the call, and this asserts it
    // rather than trusting it: a silent fallback to the host is what makes a speedup claim
    // worthless, and it is exactly the failure a missing .spv or a renamed kernel produces.
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_t89_f32");
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    const ir::Result<bool> r = t89_field_batch(pts, 0.28, 35.0, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_TRUE(r.value) << "the batch fell back to the host with a device present";

    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref =
            t89_field_at(pts[i], std::sin(0.28), std::cos(0.28), ir::t89_kp_bin(35.0));
        for (int c = 0; c < 3; ++c) worst = std::max(worst, std::fabs(out[i].v[c] - ref.v[c]));
    }
    std::printf("[ MEASURED ] device batch of %zu vs fp64 reference: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-2);

    // Below the crossover the same call must run on the host, and then it is BIT-identical.
    const std::vector<Position<Frame::GSM>> few = scatter(16);
    std::vector<ir::FieldVector<Frame::GSM>> few_out(few.size());
    const ir::Result<bool> host = t89_field_batch(few, 0.28, 35.0, few_out);
    EXPECT_FALSE(host.value);
}

TEST(IrbemT89, TheDeviceLaneRefusesABadPointBeforeItDispatches) {
    // The envelope fold rides in the staging loop, so the device lane has its own copy of the
    // decision — and its own way to get it wrong. A batch big enough to go to the device, with one
    // point inside the Earth, must come back DomainError with every output zeroed and WITHOUT
    // having reached the GPU (value false), exactly as the host lane below the crossover does.
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_t89_f32");
    std::vector<Position<Frame::GSM>> pts = scatter(n);
    pts[n / 3] = at(0.25, 0.0, 0.0);
    std::vector<ir::FieldVector<Frame::GSM>> out(n, ir::FieldVector<Frame::GSM>{
                                                        cheatah::fixarray::vec3d{7.0, 7.0, 7.0}});
    const ir::Result<bool> r = t89_field_batch(pts, 0.28, 35.0, out);
    EXPECT_EQ(r.status, Status::DomainError);
    EXPECT_FALSE(r.value) << "a refused batch must not have reached the device";
    for (const ir::FieldVector<Frame::GSM>& b : out) {
        ASSERT_EQ(b.v[0], 0.0);
        ASSERT_EQ(b.v[1], 0.0);
        ASSERT_EQ(b.v[2], 0.0);
    }
    // And one point past 70 R_E on the same device-sized batch is a caveat, not a refusal: the
    // device still runs and every answer comes back.
    std::vector<Position<Frame::GSM>> far = scatter(n);
    far[n / 2] = at(-90.0, 0.0, 0.0);
    const ir::Result<bool> f = t89_field_batch(far, 0.28, 35.0, out);
    EXPECT_EQ(f.status, Status::OutOfValidityRange);
    EXPECT_TRUE(f.value) << "an out-of-validity batch must still be computed on the device";
}

TEST(IrbemT89, DeviceKernelAgreesWithTheHostLane) {
    if (!ir::gpu::available()) {
        GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    }
    const std::size_t n = 1 << 16;
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<float> pos(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    const double ps = 0.31;
    const float sp = static_cast<float>(std::sin(ps));
    const float cp = static_cast<float>(std::cos(ps));
    std::vector<float> host(3 * n);
    std::vector<float> device(3 * n);
    ASSERT_TRUE(t89_field_host(pos, host, sp, cp, 4));
    const std::array<float, t89_param_count> block = t89_param_block(sp, cp, 4);
    ir::gpu::dispatch_batch("irbem_t89_f32", pos, device, std::span<const float>(block));

    double worst = 0.0;
    for (std::size_t i = 0; i < 3 * n; ++i) {
        worst = std::max(worst, std::fabs(static_cast<double>(device[i]) - host[i]));
    }
    std::printf("[ MEASURED ] device vs host, %zu points: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-3) << "the device is not evaluating the same expressions";
}
#endif
