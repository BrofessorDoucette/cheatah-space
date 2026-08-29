/// @file irbem_mead_test.cpp
/// @brief The suite for `space/irbem/ext_mead.hpp` — Mead & Fairfield (1975).
///
/// A polynomial model is the one kind of empirical field whose implementation CAN be checked
/// against exact statements rather than tolerances, and this suite leans on that:
///
///  - **The paper's own consistency conditions.** The fit was constrained to be divergence-free,
///    which for this functional form is three linear identities per Kp bin, twelve in all. They
///    hold to the rounding of the printed tables and are a transcription check that cannot pass by
///    accident.
///  - **`div B = 0` by central differences, and its exactness.** Every component is a quadratic,
///    so a second-order stencil is EXACT — not `O(h^2)` but zero truncation error — and the residual
///    is the same number at `h = 0.5` and `h = 5e-4`. With the identities enforced exactly it is
///    roundoff; with the printed table it is the table's rounding plus the one structural term the
///    position-only aberration leaves; and breaking a single coefficient makes it grow.
///  - **Exact symmetries.** Reflecting `z` while reversing the tilt negates `B_x` and `B_y` and
///    preserves `B_z` BITWISE, at every tilt, which is the tilt dependence verified explicitly.
///  - **The aberration.** At zero tilt `B_y` vanishes on the plane `y = -x tan(4 deg)` and nowhere
///    near `y = 0`, which pins both the angle and its sign.
///  - **The differential against the IRBEM oracle**, all four Kp bins, when the oracle is present:
///    parity, asserted at the 1e-6 budget and measured two orders below it.
///
/// The oracle test `dlopen`s IRBEM at runtime rather than linking it, exactly as the T89 suite does.

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
#include <stdexcept>
#include <string>
#include <vector>

#include "alloc_counter.hpp"
#include "irbem_domain_corpus.hpp"
#include "space/irbem/api.hpp"
#include "space/irbem/ext_mead.hpp"
#include "space/irbem/lstar.hpp"

namespace {

namespace ir = cheatah::space::irbem;
namespace corpus = cheatah_space_test;

using ir::Frame;
using ir::mead_aberration_cos;
using ir::mead_aberration_deg;
using ir::mead_aberration_sin;
using ir::mead_bin_count;
using ir::mead_coefficient_sets;
using ir::mead_components;
using ir::mead_deg_per_rad;
using ir::mead_field;
using ir::mead_field_at;
using ir::mead_field_batch;
using ir::mead_field_host;
using ir::mead_kp_bin;
using ir::mead_param_block;
using ir::mead_param_count;
using ir::mead_parameters;
using ir::MeadParameters;
using ir::Position;
using ir::Status;

/// A sink the optimizer cannot see through, so the allocation test's calls actually happen.
volatile double sink = 0.0;

/// A GSM point, spelled so a test reads as coordinates rather than as a constructor call.
Position<Frame::GSM> at(double x, double y, double z) {
    return Position<Frame::GSM>{cheatah::fixarray::vec3d(x, y, z)};
}

/// The tilts the suite sweeps: zero, and both signs of a realistic seasonal-diurnal excursion.
constexpr std::array<double, 5> kTilts{-0.55, -0.21, 0.0, 0.21, 0.55};

/// The three tilt-derived numbers the evaluator takes, from one angle.
struct Tilt {
    double s;
    double c;
    double deg;
};
Tilt tilt_of(double rad) { return {std::sin(rad), std::cos(rad), rad * mead_deg_per_rad}; }

/// `div B` at one point by central differences.
double divergence(const MeadParameters<double>& p, double ps, double x, double y, double z,
                  double h) {
    const Tilt t = tilt_of(ps);
    const auto b = [&](double a, double d, double e) {
        return mead_components<double>(p, t.s, t.c, t.deg, a, d, e);
    };
    return ((b(x + h, y, z)[0] - b(x - h, y, z)[0]) + (b(x, y + h, z)[1] - b(x, y - h, z)[1]) +
            (b(x, y, z + h)[2] - b(x, y, z - h)[2])) /
           (2.0 * h);
}

/// The same coefficient set with the divergence identities enforced EXACTLY, aberration included.
///
/// In the SM frame with the position aberrated by `alpha` and the components not, the divergence
/// of the model is
///     z [c (a2 + b1) + 2 c5] + psi [c (a4 + b2) + c6] + psi x_m [c (2 a5 + b3) + c7]
///     + psi y_m s [2 a6 - b3]
/// with `c = cos(alpha)`, `s = sin(alpha)`. Four coefficients are re-derived to zero all four
/// brackets, in dependency order.
MeadParameters<double> exactly_divergence_free(MeadParameters<double> p) {
    const double c = mead_aberration_cos;
    p.by_txy = 2.0 * p.bx_tyy;
    p.bz_txz = -c * ((2.0 * p.bx_txx) + p.by_txy);
    p.bz_tz = -c * (p.bx_tx + p.by_ty);
    p.bz_zz = -0.5 * c * (p.bx_xz + p.by_yz);
    return p;
}

/// The largest `|div B|` over the sampled box, at difference step @p h, for the sets @p make gives.
template <class Make>
double worst_divergence(double h, Make make) {
    double worst = 0.0;
    for (double ps : kTilts) {
        for (int bin = 1; bin <= static_cast<int>(mead_bin_count); ++bin) {
            const MeadParameters<double> p = make(mead_parameters<double>(bin));
            // Integer induction, coordinates derived per iteration (cert-flp30).
            for (int ix = 0; ix <= 8; ++ix)
                for (int iy = 0; iy <= 6; ++iy)
                    for (int iz = 0; iz <= 5; ++iz) {
                        const double x = -14.0 + (3.5 * ix);
                        const double y = -9.0 + (3.0 * iy);
                        const double z = -8.0 + (3.2 * iz);
                        if (std::sqrt((x * x) + (y * y) + (z * z)) < 2.0) continue;
                        worst = std::max(worst, std::fabs(divergence(p, ps, x, y, z, h)));
                    }
        }
    }
    return worst;
}

/// A deterministic scatter of GSM points over the fitted region. A 64-bit LCG so a disagreement is
/// reproducible, and never a lattice, so no component is systematically zero.
std::vector<Position<Frame::GSM>> scatter(std::size_t n) {
    std::vector<Position<Frame::GSM>> out;
    out.reserve(n);
    std::uint64_t s = 0x9E3779B97F4A7C15ULL;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11U) / 9007199254740992.0;
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double r = 2.5 + (14.0 * next());
        const double th = std::acos(1.0 - (2.0 * next()));
        const double ph = 6.283185307179586 * next();
        out.push_back(at(r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph),
                         r * std::cos(th)));
    }
    return out;
}

/// Kp x 10 values that sit INSIDE each of the four published groups — never on an edge.
constexpr std::array<double, mead_bin_count> kInsideBin{{0.0, 10.0, 23.0, 60.0}};

// ---- the oracle, opened at runtime and never linked --------------------------------------------

using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

/// The oracle handle, or nulls when IRBEM is not on this machine.
struct Oracle {
    void* handle = nullptr;
    GetField1 get_field = nullptr;
    CoordTransVec1 coord_trans = nullptr;
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

/// IGRF at 2015.0, for the total-field tests; asserted rather than dereferenced blind.
ir::Igrf<10> igrf_2015() {
    const auto r = ir::Igrf<10>::at(2015.0);
    if (!r.has_value()) throw std::runtime_error("IGRF at 2015.0 is not available to this suite");
    return r.value();
}

/// The epoch's rotations for the total-field tests.
ir::Rotations epoch_rotations(const ir::Igrf<10>& m) {
    const auto r = ir::api::rotations_at(2015, 180, 43200.0, m);
    EXPECT_EQ(Status::Ok, r.status);
    return r.value;
}

}  // namespace

// ================================================================================================
// The published parameters
// ================================================================================================

TEST(IrbemMead, PublishedCoefficientsAreDivergenceFree) {
    // The paper's least-squares fit was constrained so that div B = 0, and for a field of this
    // form in the frame it was fitted in that is three identities per bin:
    //     a2 + b1 + 2 c5 = 0,   a4 + b2 + c6 = 0,   2 a5 + b3 + c7 = 0.
    // They hold to the rounding of a table printed to four decimals — a transcription check.
    for (std::size_t bin = 0; bin < mead_bin_count; ++bin) {
        const ir::MeadCoefficients& p = mead_coefficient_sets[bin];
        EXPECT_NEAR(p.bx_xz + p.by_yz + (2.0 * p.bz_zz), 0.0, 1.5e-4) << "bin " << (bin + 1);
        EXPECT_NEAR(p.bx_tx + p.by_ty + p.bz_tz, 0.0, 1.5e-4) << "bin " << (bin + 1);
        EXPECT_NEAR((2.0 * p.bx_txx) + p.by_txy + p.bz_txz, 0.0, 1.5e-5) << "bin " << (bin + 1);
    }
    // And a coefficient copied one digit wrong would not: the control that says the tolerances
    // above are tight enough to catch a transcription slip in the fourth decimal.
    ir::MeadCoefficients broken = mead_coefficient_sets[0];
    broken.bz_zz = 0.0805;  // printed: 0.0795
    EXPECT_GT(std::fabs(broken.bx_xz + broken.by_yz + (2.0 * broken.bz_zz)), 1.5e-4);
}

TEST(IrbemMead, PublishedCoefficientsAreOrderedByDisturbance) {
    // The paper's central result: the equatorial depression at the origin, c1, deepens from
    // -9.41 nT in the quietest group to -22.9 nT in the most disturbed, and the tail-like
    // B_x = a1 z shear grows with it. Both monotone across the four columns.
    for (std::size_t bin = 1; bin < mead_bin_count; ++bin) {
        EXPECT_LT(mead_coefficient_sets[bin].bz_1, mead_coefficient_sets[bin - 1].bz_1);
        EXPECT_GT(mead_coefficient_sets[bin].bx_z, mead_coefficient_sets[bin - 1].bx_z);
    }
    EXPECT_EQ(mead_coefficient_sets[0].bz_1, -9.41);
    EXPECT_EQ(mead_coefficient_sets[3].bz_1, -22.90);
}

TEST(IrbemMead, AberrationConstantsAreExactTrigonometry) {
    // std::sin is not constexpr, so the sine and cosine of the aberration angle are literals. This
    // is what keeps them honest: bit-for-bit the library's own value.
    const double rad = mead_aberration_deg * (std::numbers::pi / 180.0);
    EXPECT_EQ(mead_aberration_sin, std::sin(rad));
    EXPECT_EQ(mead_aberration_cos, std::cos(rad));
    EXPECT_EQ(mead_aberration_deg, 4.0);
    EXPECT_EQ(mead_deg_per_rad, 180.0 / std::numbers::pi);
}

TEST(IrbemMead, KpBinsFollowThePublishedGroups) {
    // {0, 0+} | {1-, 1, 1+, 2-} | {2, 2+, 3-} | {>= 3}: in Kp x 10, edges at 4, 20 and 30, with
    // the edge value belonging to the bin ABOVE it — which is where the oracle puts it
    // (tools/oracle/mead_diff.cpp, pass 1).
    EXPECT_EQ(mead_kp_bin(0.0), 1);
    EXPECT_EQ(mead_kp_bin(3.0), 1);   // 0+
    EXPECT_EQ(mead_kp_bin(3.999), 1);
    EXPECT_EQ(mead_kp_bin(4.0), 2);
    EXPECT_EQ(mead_kp_bin(7.0), 2);   // 1-
    EXPECT_EQ(mead_kp_bin(17.0), 2);  // 2-
    EXPECT_EQ(mead_kp_bin(19.999), 2);
    EXPECT_EQ(mead_kp_bin(20.0), 3);  // 2
    EXPECT_EQ(mead_kp_bin(27.0), 3);  // 3-
    EXPECT_EQ(mead_kp_bin(29.999), 3);
    EXPECT_EQ(mead_kp_bin(30.0), 4);  // 3
    EXPECT_EQ(mead_kp_bin(90.0), 4);
    EXPECT_EQ(mead_kp_bin(300.0), 4);
    // The values check_validity refuses still map somewhere definite.
    EXPECT_EQ(mead_kp_bin(-5.0), 1);
    EXPECT_EQ(mead_kp_bin(std::nan("")), 1);
    EXPECT_EQ(mead_kp_bin(std::numeric_limits<double>::infinity()), 4);
    EXPECT_EQ(mead_kp_bin(-std::numeric_limits<double>::infinity()), 1);
}

TEST(IrbemMead, ParametersRoundTripThroughFloat) {
    for (int bin = 1; bin <= static_cast<int>(mead_bin_count); ++bin) {
        const MeadParameters<double> d = mead_parameters<double>(bin);
        const MeadParameters<float> f = mead_parameters<float>(bin);
        EXPECT_EQ(f.bx_z, static_cast<float>(d.bx_z));
        EXPECT_EQ(f.bx_xz, static_cast<float>(d.bx_xz));
        EXPECT_EQ(f.bx_t, static_cast<float>(d.bx_t));
        EXPECT_EQ(f.bx_tx, static_cast<float>(d.bx_tx));
        EXPECT_EQ(f.bx_txx, static_cast<float>(d.bx_txx));
        EXPECT_EQ(f.bx_tyy, static_cast<float>(d.bx_tyy));
        EXPECT_EQ(f.bx_tzz, static_cast<float>(d.bx_tzz));
        EXPECT_EQ(f.by_yz, static_cast<float>(d.by_yz));
        EXPECT_EQ(f.by_ty, static_cast<float>(d.by_ty));
        EXPECT_EQ(f.by_txy, static_cast<float>(d.by_txy));
        EXPECT_EQ(f.bz_1, static_cast<float>(d.bz_1));
        EXPECT_EQ(f.bz_x, static_cast<float>(d.bz_x));
        EXPECT_EQ(f.bz_xx, static_cast<float>(d.bz_xx));
        EXPECT_EQ(f.bz_yy, static_cast<float>(d.bz_yy));
        EXPECT_EQ(f.bz_zz, static_cast<float>(d.bz_zz));
        EXPECT_EQ(f.bz_tz, static_cast<float>(d.bz_tz));
        EXPECT_EQ(f.bz_txz, static_cast<float>(d.bz_txz));
        // And the double set IS the table, not a copy of it that could drift.
        EXPECT_EQ(d.bz_1, mead_coefficient_sets[static_cast<std::size_t>(bin) - 1].bz_1);
    }
}

TEST(IrbemMead, ParametersClampAnOutOfRangeBin) {
    EXPECT_EQ(mead_parameters<double>(0).bz_1, mead_coefficient_sets[0].bz_1);
    EXPECT_EQ(mead_parameters<double>(-7).bz_1, mead_coefficient_sets[0].bz_1);
    EXPECT_EQ(mead_parameters<double>(5).bz_1, mead_coefficient_sets[3].bz_1);
    EXPECT_EQ(mead_parameters<double>(1000).bz_1, mead_coefficient_sets[3].bz_1);
    EXPECT_EQ(mead_parameters<float>(99).bx_z, static_cast<float>(mead_coefficient_sets[3].bx_z));
}

// ================================================================================================
// The physics the evaluator must have
// ================================================================================================

TEST(IrbemMead, DivergenceVanishesEverywhere) {
    // Three statements, each measured, none of which an implementation could fake:
    //
    // 1. EXACTNESS. Every component is a quadratic in the coordinates, so a central difference
    //    has zero truncation error and the stencil returns the analytic divergence at ANY h. The
    //    residual at h = 0.5 and at h = 5e-4 must be the same number to roundoff — the degenerate,
    //    and stronger, form of "falls as h^2".
    // 2. WITH THE IDENTITIES ENFORCED EXACTLY, the residual is roundoff. That verifies every
    //    derivative the evaluator implicitly takes — the frame rotation, the aberration rotation
    //    of the position, all seventeen terms — because any wrong one leaves a residual that
    //    grows with the coordinate.
    // 3. WITH THE PRINTED TABLE, the residual is small and bounded: the tables' four-decimal
    //    rounding plus the one term the position-only aberration leaves, psi y_m sin(4deg)
    //    (2 a6 - b3). Its size is a fact about the published model, reported here as measured.
    const auto published = [](const MeadParameters<double>& p) { return p; };
    const double pub_coarse = worst_divergence(0.5, published);
    const double pub_fine = worst_divergence(5e-4, published);
    const double exact_coarse = worst_divergence(0.5, exactly_divergence_free);
    const double exact_fine = worst_divergence(5e-4, exactly_divergence_free);
    std::printf("[ MEASURED ] worst |div B| nT/Re, printed table: h=0.5 %.3e  h=5e-4 %.3e ; "
                "identities enforced: h=0.5 %.3e  h=5e-4 %.3e\n",
                pub_coarse, pub_fine, exact_coarse, exact_fine);

    // (2) roundoff: |B| ~ 30 nT differenced over 2h = 1e-3 Re amplifies 1e-16 to ~1e-11.
    EXPECT_LT(exact_coarse, 1e-12);
    EXPECT_LT(exact_fine, 1e-10);
    // (1) exactness: the coarse and fine stencils agree on the printed table's divergence.
    EXPECT_NEAR(pub_coarse, pub_fine, 1e-9);
    // (3) bounded by what the table's rounding and the aberration term can produce.
    EXPECT_LT(pub_coarse, 0.1);
    EXPECT_GT(pub_coarse, 1e-3) << "the printed table is NOT exactly divergence-free; a residual "
                                   "this small means the stencil is not sampling the model";

    // Pointwise exactness, at one point, both stencils, to the last few bits.
    const MeadParameters<double> p = exactly_divergence_free(mead_parameters<double>(2));
    for (double h : {1.0, 1e-1, 1e-2, 1e-3}) {
        EXPECT_LT(std::fabs(divergence(p, 0.31, -7.0, 3.0, 2.5, h)), 1e-11) << "h = " << h;
    }

    // THE CONTROL: break one coefficient and the residual must grow with z, or this test proves
    // nothing. bz_zz enters div B as 2 bz_zz z_SM, so a 1e-3 slip at GSM (-7, 3, 8) and tilt
    // 0.31 — z_SM = -7 sin(0.31) + 8 cos(0.31) = 5.48 — is 1.10e-2 nT/Re, and the stencil must
    // return exactly that, since the model is quadratic.
    MeadParameters<double> broken = p;
    broken.bz_zz += 1e-3;
    const double z_sm = (-7.0 * std::sin(0.31)) + (8.0 * std::cos(0.31));
    EXPECT_NEAR(divergence(broken, 0.31, -7.0, 3.0, 8.0, 1e-2), 2e-3 * z_sm, 1e-9);
    EXPECT_GT(std::fabs(divergence(broken, 0.31, -7.0, 3.0, 8.0, 1e-2)), 1e-2);
}

TEST(IrbemMead, ZeroTiltIsMirrorSymmetricAboutTheEquator) {
    // At psi = 0 the tilt terms vanish and what is left is odd in z for B_x and B_y and even for
    // B_z. The arithmetic on the two sides of the mirror is the same arithmetic with a sign, so
    // the symmetry holds BITWISE.
    const MeadParameters<double> p = mead_parameters<double>(3);
    for (const Position<Frame::GSM>& q : scatter(64)) {
        const std::array<double, 3> up = mead_components<double>(p, 0.0, 1.0, 0.0, q.v[0], q.v[1], q.v[2]);
        const std::array<double, 3> dn = mead_components<double>(p, 0.0, 1.0, 0.0, q.v[0], q.v[1], -q.v[2]);
        EXPECT_EQ(up[0], -dn[0]);
        EXPECT_EQ(up[1], -dn[1]);
        EXPECT_EQ(up[2], dn[2]);
    }
}

TEST(IrbemMead, TheTiltDependenceIsExactlyLinear) {
    // Two statements about the tilt, both exact.
    //
    // (a) The polynomials are LINEAR in the tilt angle: at fixed frame rotation, B at psi_2 is the
    //     mean of B at psi_1 and psi_3 when psi_2 is their mean. Measured against the oracle the
    //     angle — not its sine — is what enters (ext_mead.hpp's brief), and this is that structure.
    const MeadParameters<double> p = mead_parameters<double>(4);
    const Tilt frame = tilt_of(0.3);
    for (const Position<Frame::GSM>& q : scatter(32)) {
        const auto b = [&](double deg) {
            return mead_components<double>(p, frame.s, frame.c, deg, q.v[0], q.v[1], q.v[2]);
        };
        const std::array<double, 3> lo = b(-20.0);
        const std::array<double, 3> mid = b(5.0);
        const std::array<double, 3> hi = b(30.0);
        for (int c = 0; c < 3; ++c) {
            EXPECT_NEAR(mid[static_cast<std::size_t>(c)],
                        0.5 * (lo[static_cast<std::size_t>(c)] + hi[static_cast<std::size_t>(c)]), 1e-12);
        }
    }
    // (b) Reflecting z while reversing the tilt is an exact symmetry of the whole model, frame
    //     rotation included: (B_x, B_y, B_z)(x, y, -z; -psi) == (-B_x, -B_y, +B_z)(x, y, z; psi),
    //     bitwise, at three tilts — the explicit tilt check the model's simplicity makes possible.
    for (double ps : {0.15, 0.31, 0.55}) {
        const Tilt plus = tilt_of(ps);
        const Tilt minus = tilt_of(-ps);
        for (const Position<Frame::GSM>& q : scatter(48)) {
            const std::array<double, 3> a =
                mead_components<double>(p, plus.s, plus.c, plus.deg, q.v[0], q.v[1], q.v[2]);
            const std::array<double, 3> m =
                mead_components<double>(p, minus.s, minus.c, minus.deg, q.v[0], q.v[1], -q.v[2]);
            EXPECT_EQ(a[0], -m[0]) << "tilt " << ps;
            EXPECT_EQ(a[1], -m[1]) << "tilt " << ps;
            EXPECT_EQ(a[2], m[2]) << "tilt " << ps;
        }
        // And the tilt terms are not idle: the field at +psi differs from the field at zero tilt.
        const std::array<double, 3> t0 = mead_components<double>(p, 0.0, 1.0, 0.0, -6.6, 0.0, 0.0);
        const std::array<double, 3> t1 =
            mead_components<double>(p, plus.s, plus.c, plus.deg, -6.6, 0.0, 0.0);
        EXPECT_GT(std::fabs(t1[0] - t0[0]), 1.0) << "B_x gains a psi-proportional term of ~0.5 nT/deg";
    }
}

TEST(IrbemMead, TheAberrationRotatesTheNoonMidnightPlaneByFourDegrees) {
    // At zero tilt B_y = b1 y_m z, and y_m = x sin(4deg) + y cos(4deg). So B_y vanishes on the
    // plane y = -x tan(4deg) and is b1 x sin(4deg) z on the GSM meridian y = 0. Both pin the
    // angle AND its sign; an implementation that rotated the other way, or not at all, fails.
    const MeadParameters<double> p = mead_parameters<double>(1);
    const double tan4 = mead_aberration_sin / mead_aberration_cos;
    for (double x : {-10.0, -5.0, 3.0, 8.0}) {
        for (double z : {-3.0, 1.5, 4.0}) {
            const std::array<double, 3> on_plane =
                mead_components<double>(p, 0.0, 1.0, 0.0, x, -x * tan4, z);
            EXPECT_LT(std::fabs(on_plane[1]), 1e-13 * std::fabs(p.by_yz * x * z));
            const std::array<double, 3> on_meridian =
                mead_components<double>(p, 0.0, 1.0, 0.0, x, 0.0, z);
            EXPECT_NEAR(on_meridian[1], p.by_yz * x * mead_aberration_sin * z,
                        1e-12 * std::fabs(p.by_yz * x * z));
            EXPECT_GT(std::fabs(on_meridian[1]), 1e-3) << "the aberration must be visible on y = 0";
        }
    }
}

TEST(IrbemMead, NearTailBzDepressionDeepensWithKp) {
    // The paper's Figure-level result: the external B_z at the origin and in the near tail is
    // negative and deepens monotonically through the four groups.
    double prev = 0.0;
    for (int bin = 1; bin <= static_cast<int>(mead_bin_count); ++bin) {
        const std::array<double, 3> b =
            mead_components<double>(mead_parameters<double>(bin), 0.0, 1.0, 0.0, -6.6, 0.0, 0.0);
        EXPECT_LT(b[2], prev) << "bin " << bin;
        prev = b[2];
    }
    EXPECT_LT(prev, -20.0);
}

TEST(IrbemMead, TheDaysideFieldIsCompressedAndTheNightsideStretched) {
    // Geosynchronous noon versus midnight, quiet bin, zero tilt: the magnetopause currents add a
    // northward B_z on the dayside and the tail depresses it on the nightside.
    const MeadParameters<double> p = mead_parameters<double>(1);
    const std::array<double, 3> noon = mead_components<double>(p, 0.0, 1.0, 0.0, 6.6, 0.0, 0.0);
    const std::array<double, 3> midnight = mead_components<double>(p, 0.0, 1.0, 0.0, -6.6, 0.0, 0.0);
    EXPECT_GT(noon[2], 0.0);
    EXPECT_LT(midnight[2], 0.0);
    EXPECT_GT(noon[2] - midnight[2], 15.0);
}

TEST(IrbemMead, ReferenceLaneMatchesTheComponentForm) {
    const Tilt t = tilt_of(0.27);
    for (const Position<Frame::GSM>& q : scatter(32)) {
        const ir::FieldVector<Frame::GSM> v = mead_field_at(q, t.s, t.c, t.deg, 2);
        const std::array<double, 3> b =
            mead_components<double>(mead_parameters<double>(2), t.s, t.c, t.deg, q.v[0], q.v[1], q.v[2]);
        EXPECT_EQ(v.v[0], b[0]);
        EXPECT_EQ(v.v[1], b[1]);
        EXPECT_EQ(v.v[2], b[2]);
    }
}

TEST(IrbemMead, KpBinsSelectDistinctCoefficientSets) {
    // Kp is a BIN: two values inside one group give bitwise the same field, and the four groups
    // give four different ones. A continuous-in-Kp "smoothness" test would be a BUG here.
    const Position<Frame::GSM> q = at(-6.0, 2.0, 1.0);
    std::array<double, mead_bin_count> bz{};
    for (std::size_t bin = 0; bin < mead_bin_count; ++bin) {
        bz[bin] = mead_field(q, 0.2, kInsideBin[bin]).value.v[2];
    }
    for (std::size_t a = 0; a < mead_bin_count; ++a)
        for (std::size_t b = a + 1; b < mead_bin_count; ++b) EXPECT_NE(bz[a], bz[b]);
    // Inside each group, identical.
    EXPECT_EQ(mead_field(q, 0.2, 0.0).value.v[2], mead_field(q, 0.2, 3.0).value.v[2]);
    EXPECT_EQ(mead_field(q, 0.2, 7.0).value.v[2], mead_field(q, 0.2, 17.0).value.v[2]);
    EXPECT_EQ(mead_field(q, 0.2, 20.0).value.v[2], mead_field(q, 0.2, 27.0).value.v[2]);
    EXPECT_EQ(mead_field(q, 0.2, 30.0).value.v[2], mead_field(q, 0.2, 90.0).value.v[2]);
}

TEST(IrbemMead, OnlyKpDrivesTheModel) {
    // The corpus's four regimes vary Dst, Pdyn, By and Bz along with Kp. Only Kp reaches this
    // model: a context built from each regime's full driver vector gives the field its Kp alone
    // gives, and the storm and extreme regimes — Kp 6 and 8.5, both `>= 3` — give the SAME field,
    // which is a property of the published binning and not a defect.
    const ir::Epoch epoch{2015.5, 43200.0, 2015, 180};
    ir::RotationTable identity{};
    for (cheatah::fixarray::mat3d& m : identity) m = cheatah::fixarray::mat3d::identity();
    const Position<Frame::GSM> q = at(-6.6, 1.0, 0.5);
    std::array<double, 4> bz{};
    for (std::size_t r = 0; r < corpus::regime_drivers.size(); ++r) {
        const corpus::MagInput& m = corpus::regime_drivers[r];
        ir::DriverSet drivers{};
        drivers[static_cast<std::size_t>(ir::Driver::Kp)] = m.kp * 10.0;
        drivers[static_cast<std::size_t>(ir::Driver::Dst)] = m.dst;
        drivers[static_cast<std::size_t>(ir::Driver::Pdyn)] = m.pdyn;
        drivers[static_cast<std::size_t>(ir::Driver::ByIMF)] = m.by_imf;
        drivers[static_cast<std::size_t>(ir::Driver::BzIMF)] = m.bz_imf;
        const ir::ContextResult built = ir::make_field_context(epoch, 0.2, identity, drivers);
        ASSERT_TRUE(built.has_value());
        const ir::Result<ir::FieldVector<Frame::GSM>> viaContext = mead_field(q, built.value());
        EXPECT_EQ(viaContext.status, Status::Ok) << "regime " << r;
        EXPECT_EQ(viaContext.value.v[2], mead_field(q, 0.2, m.kp * 10.0).value.v[2]);
        bz[r] = viaContext.value.v[2];
    }
    EXPECT_NE(bz[0], bz[1]) << "quiet (Kp 1) and moderate (Kp 3.5) are different groups";
    EXPECT_NE(bz[1], bz[0]);
    EXPECT_EQ(bz[1], bz[2]) << "moderate, storm and extreme all fall in the Kp >= 3 group";
    EXPECT_EQ(bz[2], bz[3]);
    // And the storm events of record, likewise: all four are Kp >= 5, the most disturbed group.
    for (const corpus::StormEvent& e : corpus::storm_events) {
        EXPECT_EQ(mead_kp_bin(e.mag.kp * 10.0), 4) << e.name;
    }
}

// ================================================================================================
// The envelope
// ================================================================================================

TEST(IrbemMead, OutOfRangeKpIsReportedButStillEvaluated) {
    for (double kp10 : {-1.0, -25.0, 91.0, 200.0}) {
        const ir::Result<ir::FieldVector<Frame::GSM>> r = mead_field(at(-6.6, 0.0, 0.0), 0.2, kp10);
        EXPECT_EQ(r.status, Status::OutOfValidityRange) << "Kp x 10 = " << kp10;
        EXPECT_NE(r.value.v[2], 0.0) << "the value must still be computed";
    }
    EXPECT_EQ(mead_field(at(-6.6, 0.0, 0.0), 0.2, 0.0).status, Status::Ok);
    EXPECT_EQ(mead_field(at(-6.6, 0.0, 0.0), 0.2, 90.0).status, Status::Ok);
}

TEST(IrbemMead, TheRadialEnvelopeIsCheckedFromBothSides) {
    // The published spatial envelope is r <= 17 R_E. Just inside: Ok, with a value. Just outside:
    // OutOfValidityRange, with the SAME functional form still evaluated. The bound is closed.
    const ir::Result<ir::FieldVector<Frame::GSM>> in = mead_field(at(-16.99, 0.0, 0.0), 0.1, 20.0);
    EXPECT_EQ(in.status, Status::Ok);
    EXPECT_TRUE(std::isfinite(in.value.v[2]));
    EXPECT_EQ(mead_field(at(17.0, 0.0, 0.0), 0.1, 20.0).status, Status::Ok);
    const ir::Result<ir::FieldVector<Frame::GSM>> out = mead_field(at(-17.01, 0.0, 0.0), 0.1, 20.0);
    EXPECT_EQ(out.status, Status::OutOfValidityRange);
    const Tilt t = tilt_of(0.1);
    EXPECT_EQ(out.value.v[2], mead_field_at(at(-17.01, 0.0, 0.0), t.s, t.c, t.deg, 3).v[2]);
    EXPECT_NE(out.value.v[2], 0.0);
    // Off-axis, same radius, same verdict.
    EXPECT_EQ(mead_field(at(10.0, 10.0, 9.0), 0.1, 20.0).status, Status::Ok);           // r = 16.76
    EXPECT_EQ(mead_field(at(10.0, 10.0, 10.0), 0.1, 20.0).status, Status::OutOfValidityRange);  // 17.32
    // Inside the solid Earth is a different failure: there is no answer, not an unreliable one.
    const ir::Result<ir::FieldVector<Frame::GSM>> inside = mead_field(at(0.5, 0.0, 0.0), 0.1, 20.0);
    EXPECT_EQ(inside.status, Status::DomainError);
    EXPECT_EQ(inside.value.v[2], 0.0);
    // Both caveats at once: the FIRST reason — the drivers are checked before the position.
    EXPECT_EQ(mead_field(at(-20.0, 0.0, 0.0), 0.1, 200.0).status, Status::OutOfValidityRange);
}

TEST(IrbemMead, NonFiniteInputIsADomainError) {
    const double nan = std::nan("");
    const double inf = std::numeric_limits<double>::infinity();
    for (const Position<Frame::GSM>& p : {at(nan, 0.0, 0.0), at(0.0, inf, 0.0), at(4.0, 0.0, -nan)}) {
        const ir::Result<ir::FieldVector<Frame::GSM>> r = mead_field(p, 0.2, 20.0);
        EXPECT_EQ(r.status, Status::DomainError);
        EXPECT_EQ(r.value.v[0], 0.0);
        EXPECT_EQ(r.value.v[1], 0.0);
        EXPECT_EQ(r.value.v[2], 0.0);
    }
    EXPECT_EQ(mead_field(at(5.0, 0.0, 0.0), nan, 20.0).status, Status::DomainError);
    EXPECT_EQ(mead_field(at(5.0, 0.0, 0.0), 0.2, nan).status, Status::DomainError);
    EXPECT_EQ(mead_field(at(5.0, 0.0, 0.0), inf, 20.0).status, Status::DomainError);
}

TEST(IrbemMead, ATiltBeyondARightAngleIsADomainError) {
    // Unlike T89 there is no tan(psi): a right-angle tilt itself has a value. Beyond it the number
    // is not a dipole tilt at all — the angle between the dipole axis and z_GSM lies in
    // [-90, 90] degrees by definition, the bound FieldContext enforces — so it is refused, and
    // refused BEFORE the envelope is consulted, on the scalar and the batch lane alike.
    const double quarter_turn = std::numbers::pi / 2.0;
    EXPECT_EQ(mead_field(at(5.0, 0.0, 0.0), quarter_turn, 20.0).status, Status::Ok);
    EXPECT_EQ(mead_field(at(5.0, 0.0, 0.0), -quarter_turn, 20.0).status, Status::Ok);
    const ir::Result<ir::FieldVector<Frame::GSM>> r =
        mead_field(at(5.0, 0.0, 0.0), std::nextafter(quarter_turn, 4.0), 20.0);
    EXPECT_EQ(r.status, Status::DomainError);
    EXPECT_EQ(r.value.v[0], 0.0);
    EXPECT_EQ(r.value.v[2], 0.0);
    EXPECT_EQ(mead_field(at(5.0, 0.0, 0.0), -2.0, 20.0).status, Status::DomainError);
    const std::vector<Position<Frame::GSM>> pts = scatter(8);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    EXPECT_EQ(mead_field_batch(pts, 2.0, 20.0, out).status, Status::DomainError);
    EXPECT_EQ(mead_field_batch(pts, quarter_turn, 20.0, out).status, Status::Ok);
    // Within that bound and a finite radius the quadratics cannot overflow: the most distant
    // finite point there is, at the largest tilt, is merely out of validity — and finite.
    const ir::Result<ir::FieldVector<Frame::GSM>> far =
        mead_field(at(1e150, -1e150, 1e150), quarter_turn, 20.0);
    EXPECT_EQ(far.status, Status::OutOfValidityRange);
    EXPECT_TRUE(std::isfinite(far.value.v[0]));
    EXPECT_TRUE(std::isfinite(far.value.v[1]));
    EXPECT_TRUE(std::isfinite(far.value.v[2]));
}

TEST(IrbemMead, AllLocalTimesAndOffEquatorAgreeAcrossLanes) {
    // Every corpus local time at geosynchronous and at two off-equator heights, every bin, three
    // tilts: the scalar entry, the reference lane and the host batch lane must be the same
    // numbers, and the status must be Ok everywhere inside the envelope.
    std::vector<Position<Frame::GSM>> pts;
    for (double mlt : corpus::local_times) {
        const double phi = (mlt - 12.0) * (std::numbers::pi / 12.0);  // noon on +x
        for (double z : {-2.0, 0.0, 2.0}) pts.push_back(at(6.6 * std::cos(phi), 6.6 * std::sin(phi), z));
    }
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    for (double ps : {-0.45, 0.0, 0.3}) {
        for (std::size_t bin = 0; bin < mead_bin_count; ++bin) {
            const ir::Result<bool> r = mead_field_batch(pts, ps, kInsideBin[bin], out);
            EXPECT_EQ(r.status, Status::Ok);
            for (std::size_t i = 0; i < pts.size(); ++i) {
                const ir::Result<ir::FieldVector<Frame::GSM>> s = mead_field(pts[i], ps, kInsideBin[bin]);
                EXPECT_EQ(s.status, Status::Ok);
                EXPECT_EQ(s.value.v[0], out[i].v[0]);
                EXPECT_EQ(s.value.v[1], out[i].v[1]);
                EXPECT_EQ(s.value.v[2], out[i].v[2]);
                EXPECT_LT(s.value.magnitude(), 120.0) << "an external field of tens of nT at 6.6 Re";
                EXPECT_GT(s.value.magnitude(), 1.0);
            }
        }
    }
}

TEST(IrbemMead, ContextOverloadAgreesWithTheExplicitOne) {
    const ir::Epoch epoch{2015.5, 43200.0, 2015, 180};
    ir::RotationTable identity{};
    for (cheatah::fixarray::mat3d& m : identity) m = cheatah::fixarray::mat3d::identity();
    ir::DriverSet drivers{};
    drivers[static_cast<std::size_t>(ir::Driver::Kp)] = 23.0;
    const ir::ContextResult built = ir::make_field_context(epoch, -0.42, identity, drivers);
    ASSERT_TRUE(built.has_value()) << ir::describe(built.error());

    const Position<Frame::GSM> p = at(3.75, -1.5, 2.25);
    const ir::Result<ir::FieldVector<Frame::GSM>> viaContext = mead_field(p, built.value());
    const ir::Result<ir::FieldVector<Frame::GSM>> viaScalars = mead_field(p, -0.42, 23.0);
    EXPECT_EQ(viaContext.status, viaScalars.status);
    EXPECT_EQ(viaContext.value.v[0], viaScalars.value.v[0]);
    EXPECT_EQ(viaContext.value.v[1], viaScalars.value.v[1]);
    EXPECT_EQ(viaContext.value.v[2], viaScalars.value.v[2]);
}

// ================================================================================================
// The batch lanes
// ================================================================================================

TEST(IrbemMead, HostFloatLaneTracksTheReferenceLane) {
    const std::size_t n = 4096;
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<float> pos(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    const Tilt t = tilt_of(0.31);
    std::vector<float> out(3 * n);
    ASSERT_TRUE(mead_field_host(pos, out, static_cast<float>(t.s), static_cast<float>(t.c),
                                static_cast<float>(t.deg), 4));
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        // The reference evaluated at the ROUNDED position, so the comparison is of arithmetic.
        const ir::FieldVector<Frame::GSM> ref =
            mead_field_at(at(pos[(3 * i) + 0], pos[(3 * i) + 1], pos[(3 * i) + 2]), t.s, t.c, t.deg, 4);
        for (int c = 0; c < 3; ++c) {
            worst = std::max(worst, std::fabs(static_cast<double>(out[(3 * i) + static_cast<std::size_t>(c)]) -
                                              ref.v[c]));
        }
    }
    std::printf("[ MEASURED ] fp32 host lane vs fp64 reference, %zu points: max |dB| = %.3e nT\n", n,
                worst);
    EXPECT_LT(worst, 2e-4);  // fp32 on a ~100 nT field: a few ulps of 1e-5
}

TEST(IrbemMead, HostFloatLaneRejectsMismatchedSpans) {
    std::vector<float> pos(7);  // not a whole number of points
    std::vector<float> out(7, 3.0F);
    EXPECT_FALSE(mead_field_host(pos, out, 0.0F, 1.0F, 0.0F, 1));
    EXPECT_EQ(out[0], 3.0F) << "nothing may be written on refusal";
    std::vector<float> pos6(6);
    std::vector<float> out3(3);
    EXPECT_FALSE(mead_field_host(pos6, out3, 0.0F, 1.0F, 0.0F, 1));
    std::vector<float> out6(6);
    EXPECT_TRUE(mead_field_host(pos6, out6, 0.0F, 1.0F, 0.0F, 1));
}

TEST(IrbemMead, ParameterBlockCarriesTheTiltThenTheCoefficients) {
    const std::array<float, mead_param_count> block = mead_param_block(0.25F, 0.75F, 14.5F, 2);
    const MeadParameters<float> p = mead_parameters<float>(2);
    EXPECT_EQ(block[0], 0.25F);
    EXPECT_EQ(block[1], 0.75F);
    EXPECT_EQ(block[2], 14.5F);
    EXPECT_EQ(block[3], p.bx_z);
    EXPECT_EQ(block[9], p.bx_tzz);
    EXPECT_EQ(block[10], p.by_yz);
    EXPECT_EQ(block[12], p.by_txy);
    EXPECT_EQ(block[13], p.bz_1);
    EXPECT_EQ(block[19], p.bz_txz);
    static_assert(mead_param_count == 20, "the kernel reads twenty floats");
#if CHEATAH_SPACE_IRBEM_MEAD_GPU
    EXPECT_EQ(ir::gpu::kernel_info("irbem_mead_f32").params, mead_param_count);
    EXPECT_EQ(ir::gpu::kernel_info("irbem_mead_f32").bindings, 4U);
#endif
}

TEST(IrbemMead, BatchAgreesWithTheReferenceLane) {
    // Below the crossover (or on a machine with no device) the batch is the fp64 host loop and
    // must be bit-identical to the scalar entry point.
    const std::vector<Position<Frame::GSM>> pts = scatter(64);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    const ir::Result<bool> r = mead_field_batch(pts, 0.28, 23.0, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_FALSE(r.value) << "64 points is below every crossover";
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const ir::Result<ir::FieldVector<Frame::GSM>> s = mead_field(pts[i], 0.28, 23.0);
        EXPECT_EQ(out[i].v[0], s.value.v[0]);
        EXPECT_EQ(out[i].v[1], s.value.v[1]);
        EXPECT_EQ(out[i].v[2], s.value.v[2]);
    }
}

TEST(IrbemMead, BatchRejectsMismatchedSpans) {
    const std::vector<Position<Frame::GSM>> pts = scatter(4);
    std::vector<ir::FieldVector<Frame::GSM>> shorter(3);
    EXPECT_EQ(mead_field_batch(pts, 0.1, 20.0, shorter).status, Status::DomainError);
    std::vector<ir::FieldVector<Frame::GSM>> right(4);
    EXPECT_EQ(mead_field_batch(pts, std::nan(""), 20.0, right).status, Status::DomainError);
    EXPECT_EQ(mead_field_batch(pts, 0.1, std::nan(""), right).status, Status::DomainError);
    const ir::Result<bool> empty = mead_field_batch({}, 0.1, 20.0, {});
    EXPECT_EQ(empty.status, Status::Ok);
    EXPECT_FALSE(empty.value);
    EXPECT_EQ(mead_field_batch(pts, 0.1, 300.0, right).status, Status::OutOfValidityRange);
}

TEST(IrbemMead, BatchReportsTheSameEnvelopeTheScalarLaneDoes) {
    const std::vector<Position<Frame::GSM>> good{at(5.0, 1.0, 1.0), at(-8.0, 2.0, -1.0)};
    std::vector<ir::FieldVector<Frame::GSM>> out(good.size());
    EXPECT_EQ(mead_field_batch(good, 0.2, 20.0, out).status, Status::Ok);

    // One point beyond r = 17 R_E: the whole batch is out of validity, every point computed.
    const std::vector<Position<Frame::GSM>> far{at(5.0, 1.0, 1.0), at(-18.0, 0.0, 0.0)};
    EXPECT_EQ(mead_field_batch(far, 0.2, 20.0, out).status, Status::OutOfValidityRange);
    EXPECT_EQ(mead_field(far[1], 0.2, 20.0).status, Status::OutOfValidityRange);
    for (const ir::FieldVector<Frame::GSM>& b : out) EXPECT_NE(b.v[2], 0.0);

    for (const Position<Frame::GSM>& bad :
         {at(0.5, 0.0, 0.0), at(std::nan(""), 0.0, 0.0), at(std::numeric_limits<double>::infinity(), 0.0, 0.0)}) {
        const std::vector<Position<Frame::GSM>> mixed{at(5.0, 1.0, 1.0), bad};
        std::vector<ir::FieldVector<Frame::GSM>> mixed_out(
            mixed.size(), ir::FieldVector<Frame::GSM>{cheatah::fixarray::vec3d{1.0, 1.0, 1.0}});
        EXPECT_EQ(mead_field_batch(mixed, 0.2, 20.0, mixed_out).status, Status::DomainError);
        EXPECT_EQ(mead_field(bad, 0.2, 20.0).status, Status::DomainError);
        for (const ir::FieldVector<Frame::GSM>& b : mixed_out) {
            EXPECT_EQ(b.v[0], 0.0);
            EXPECT_EQ(b.v[1], 0.0);
            EXPECT_EQ(b.v[2], 0.0);
        }
    }
    // Two bad points, the first one first: the fold must stay poisoned.
    const std::vector<Position<Frame::GSM>> two_bad{
        at(std::nan(""), 0.0, 0.0), at(std::numeric_limits<double>::infinity(), 1.0, 0.0), at(5.0, 1.0, 1.0)};
    std::vector<ir::FieldVector<Frame::GSM>> two_out(two_bad.size());
    EXPECT_EQ(mead_field_batch(two_bad, 0.2, 20.0, two_out).status, Status::DomainError);
    // The fold directly: both ends of the radius band, and the closed bound.
    ir::MeadPositionFold fold;
    fold.add(at(17.0, 0.0, 0.0));
    fold.add(at(3.0, 0.0, 0.0));
    EXPECT_EQ(fold.verdict(), Status::Ok);
    fold.add(at(17.01, 0.0, 0.0));
    EXPECT_EQ(fold.verdict(), Status::OutOfValidityRange);
    // An out-of-range Kp is still reported on an EMPTY batch.
    const std::vector<Position<Frame::GSM>> none;
    std::vector<ir::FieldVector<Frame::GSM>> none_out;
    EXPECT_EQ(mead_field_batch(none, 0.2, 900.0, none_out).status, Status::OutOfValidityRange);
    EXPECT_EQ(mead_field_batch(none, 0.2, 20.0, none_out).status, Status::Ok);
}

TEST(IrbemMead, NothingOnTheHeapInTheHotPath) {
    const std::vector<Position<Frame::GSM>> pts = scatter(256);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    std::vector<float> fpos(std::size_t{3} * 64, 1.5F);
    std::vector<float> fout(std::size_t{3} * 64);
    (void)mead_field_batch(pts, 0.2, 30.0, out);
    (void)mead_field(at(5.0, 1.0, 1.0), 0.2, 30.0);

    const std::size_t before = cheatah_space_test::allocation_count();
    for (int i = 0; i < 64; ++i) {
        sink = sink + mead_field(at(5.0 + (0.01 * i), 1.0, 1.0), 0.2, 30.0).value.v[2];
    }
    (void)mead_field_batch(pts, 0.2, 30.0, out);  // below the crossover: the host lane
    sink = sink + (mead_field_host(fpos, fout, 0.1F, 0.99F, 5.7F, 2) ? fout[5] : 0.0);
    sink = sink + out[0].v[2];
    EXPECT_EQ(before, cheatah_space_test::allocation_count());
}

// ================================================================================================
// The superposition
// ================================================================================================

TEST(IrbemMead, TotalFieldSuperposesInternalAndExternal) {
    const ir::Igrf<10> igrf = igrf_2015();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::TotalFieldMead<10> total(igrf, rot, 23.0);
    static_assert(ir::GeoFieldModel<ir::TotalFieldMead<10>>, "a tracer must be able to follow it");
    static_assert(ir::TotalFieldMead<10>::degree == 10);

    const Position<Frame::GEO> p{cheatah::fixarray::vec3d{6.0, 0.0, 0.0}};
    const ir::FieldVector<Frame::GEO> b_int = igrf.evaluate(p);
    const ir::FieldVector<Frame::GEO> b_tot = total.evaluate(p);
    EXPECT_EQ(23.0, total.kp_times_ten());
    EXPECT_EQ(igrf.g(1, 0), total.g(1, 0));
    EXPECT_EQ(igrf.h(2, 1), total.h(2, 1));
    EXPECT_EQ(igrf.g(1, 0), total.internal().g(1, 0));
    EXPECT_EQ(&rot, &total.rotations());
    // The external field is a real contribution at L = 6 and a perturbation, not a replacement.
    const double rel = std::fabs(b_tot.magnitude() - b_int.magnitude()) / b_int.magnitude();
    EXPECT_GT(rel, 1e-3);
    EXPECT_LT(rel, 0.5);
    // And it is EXACTLY the external model rotated into GEO: rebuild the sum by hand.
    const Position<Frame::GSM> p_gsm = ir::transform<Frame::GSM>(p, rot);
    const double tilt_rad = rot.dipole_tilt_deg * (std::numbers::pi / 180.0);
    const ir::Result<ir::FieldVector<Frame::GSM>> ext = mead_field(p_gsm, tilt_rad, 23.0);
    ASSERT_EQ(ext.status, Status::Ok);
    const ir::FieldVector<Frame::GEO> ext_geo = ir::transform<Frame::GEO>(ext.value, rot);
    for (int c = 0; c < 3; ++c) EXPECT_EQ(b_tot.v[c], b_int.v[c] + ext_geo.v[c]);
    EXPECT_EQ(total.external_status(p), Status::Ok);
}

TEST(IrbemMead, TotalFieldReportsWhenTheExternalModelDeclines) {
    const ir::Igrf<10> igrf = igrf_2015();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::TotalFieldMead<10> total(igrf, rot, 23.0);
    // Beyond 17 Re the external model is out of validity: reported, and still ADDED — the
    // extrapolation is the caller's decision, and the status is how they learn it was made.
    const Position<Frame::GEO> far{cheatah::fixarray::vec3d{-18.0, 0.0, 0.0}};
    EXPECT_EQ(total.external_status(far), Status::OutOfValidityRange);
    EXPECT_NE(total.evaluate(far).magnitude(), igrf.evaluate(far).magnitude());
    // A point the model refuses outright: the internal field alone, never a zero or a NaN.
    const Position<Frame::GEO> inside{cheatah::fixarray::vec3d{0.5, 0.0, 0.0}};
    EXPECT_EQ(total.external_status(inside), Status::DomainError);
    EXPECT_EQ(total.evaluate(inside).magnitude(), igrf.evaluate(inside).magnitude());
    const double nan = std::nan("");
    const Position<Frame::GEO> bad{cheatah::fixarray::vec3d{nan, 0.0, 0.0}};
    EXPECT_EQ(total.external_status(bad), Status::DomainError);
}

TEST(IrbemMead, TotalFieldTracesADriftShellThroughEveryRegime) {
    // The whole reason an external model exists: a field line traced through the TOTAL field, at
    // geosynchronous, in each corpus regime. Every trace closes, the invariants are finite, and
    // the quiet regime's B_min differs from the disturbed ones' — while the three regimes that
    // share the Kp >= 3 group give identical invariants, as the binning says they must.
    const ir::Igrf<10> igrf = igrf_2015();
    const ir::Rotations rot = epoch_rotations(igrf);
    const Position<Frame::GEO> p{cheatah::fixarray::vec3d{6.6, 0.0, 0.0}};
    std::vector<double> bmin;
    for (const corpus::MagInput& m : corpus::regime_drivers) {
        const ir::TotalFieldMead<10> total(igrf, rot, m.kp * 10.0);
        const ir::Result<ir::FieldLine> line = ir::trace_invariant(total, p, 45.0);
        ASSERT_EQ(line.status, Status::Ok) << "Kp " << m.kp;
        EXPECT_TRUE(std::isfinite(line.value.invariant_i));
        EXPECT_GT(line.value.b_min, 0.0);
        bmin.push_back(line.value.b_min);
    }
    const ir::Result<ir::FieldLine> internal_only = ir::trace_invariant(igrf, p, 45.0);
    ASSERT_EQ(internal_only.status, Status::Ok);
    EXPECT_NE(bmin[0], internal_only.value.b_min) << "the external field must reach the trace";
    EXPECT_NE(bmin[0], bmin[1]) << "quiet and moderate are different Kp groups";
    EXPECT_EQ(bmin[1], bmin[2]) << "moderate, storm and extreme share the Kp >= 3 group";
    EXPECT_EQ(bmin[2], bmin[3]);
    // GEO (6.6, 0, 0) at 12:00 UT is local NOON: the magnetopause compression RAISES the
    // equatorial field there. The mirror point, GEO (-6.6, 0, 0), is local midnight, where the
    // tail-like depression lowers it. Both signs, or the frames are wrong.
    EXPECT_GT(bmin[3], internal_only.value.b_min);
    const Position<Frame::GEO> midnight{cheatah::fixarray::vec3d{-6.6, 0.0, 0.0}};
    const ir::TotalFieldMead<10> disturbed(igrf, rot, 60.0);
    const ir::Result<ir::FieldLine> night = ir::trace_invariant(disturbed, midnight, 45.0);
    const ir::Result<ir::FieldLine> night_internal = ir::trace_invariant(igrf, midnight, 45.0);
    ASSERT_EQ(night.status, Status::Ok);
    ASSERT_EQ(night_internal.status, Status::Ok);
    EXPECT_LT(night.value.b_min, night_internal.value.b_min);
}

// ================================================================================================
// The differential against the IRBEM oracle
// ================================================================================================

TEST(IrbemMead, MatchesTheIrbemOracleToParity) {
    const Oracle& o = oracle();
    if (!o.usable()) {
        GTEST_SKIP() << "IRBEM oracle not present (set CHEATAH_SPACE_IRBEM_ORACLE to its .so); "
                        "the oracle is a dev-only black box and is never linked";
    }
    // Three epochs spanning the tilt range (+0.0002, +25.64, -30.42 degrees). The external field is
    // isolated as kext=1 minus kext=0 so the internal IGRF term cancels EXACTLY, and the tilt is
    // read out of the oracle itself so a frame difference cannot masquerade as a model difference.
    struct Epoch {
        int doy;
        double ut;
    };
    const std::array<Epoch, 3> epochs{{{80, 39183.0}, {180, 43200.0}, {355, 7200.0}}};

    // THIS is a parity target, not a measured gap: both the form and the coefficients are
    // published, and tools/oracle/mead_diff.cpp shows the oracle evaluates exactly them. The
    // budget is the internal-field standard, 1e-6 relative; the measurement sits at ~2e-9 RMS.
    for (std::size_t bi = 0; bi < mead_bin_count; ++bi) {
        double sum2 = 0.0;
        double sig2 = 0.0;
        double worst_rel = 0.0;
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
            const Tilt t = tilt_of(ps);
            for (int ix = 0; ix <= 8; ++ix)
                for (int iy = 0; iy <= 6; ++iy)
                    for (int iz = 0; iz <= 4; ++iz) {
                        const double x = -14.0 + (3.5 * ix);
                        const double y = -9.0 + (3.0 * iy);
                        const double z = -8.0 + (4.0 * iz);
                        const double r = std::sqrt((x * x) + (y * y) + (z * z));
                        if (r < 1.5 || r > 17.0) continue;
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
                        int k1 = 1;
                        std::vector<double> mag(25, 0.0);
                        mag[0] = kInsideBin[bi];
                        std::array<double, 3> b0{};
                        std::array<double, 3> b1{};
                        double m0 = 0.0;
                        double m1 = 0.0;
                        double x1 = geo[0];
                        double x2 = geo[1];
                        double x3 = geo[2];
                        o.get_field(&k0, options.data(), &sysaxes, &iyear, &idoy, &ut, &x1, &x2, &x3,
                                    mag.data(), b0.data(), &m0);
                        o.get_field(&k1, options.data(), &sysaxes, &iyear, &idoy, &ut, &x1, &x2, &x3,
                                    mag.data(), b1.data(), &m1);
                        std::array<double, 3> dgeo{b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2]};
                        std::array<double, 3> ora{};
                        {
                            int si = 1;
                            int so = 2;
                            o.coord_trans(&one, &si, &so, &iyear, &idoy, &ut, dgeo.data(), ora.data());
                        }
                        const ir::FieldVector<Frame::GSM> mine =
                            mead_field_at(at(x, y, z), t.s, t.c, t.deg, static_cast<int>(bi) + 1);
                        double d2 = 0.0;
                        double o2 = 0.0;
                        for (int c = 0; c < 3; ++c) {
                            const double d = mine.v[c] - ora[static_cast<std::size_t>(c)];
                            d2 += d * d;
                            o2 += ora[static_cast<std::size_t>(c)] * ora[static_cast<std::size_t>(c)];
                        }
                        sum2 += d2;
                        sig2 += o2;
                        worst_rel = std::max(worst_rel, std::sqrt(d2 / o2));
                        ++n;
                    }
        }
        ASSERT_GT(n, 0U);
        const double rms_rel = std::sqrt(sum2 / sig2);
        std::printf("[ MEASURED ] Kp bin %zu: vs IRBEM kext=1, RMS rel %.3e, worst rel %.3e (%zu points, "
                    "3 tilts)\n",
                    bi + 1, rms_rel, worst_rel, n);
        EXPECT_LT(worst_rel, 1e-6) << "Kp bin " << (bi + 1) << ": outside the parity budget";
        EXPECT_LT(rms_rel, 1e-7) << "Kp bin " << (bi + 1);
    }
}

#if CHEATAH_SPACE_IRBEM_MEAD_GPU
namespace {

/// Point `CHEATAH_SPACE_IRBEM_SPV_DIR` somewhere for the life of the object, and put it back.
class SpvDirScope {
  public:
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
    ~SpvDirScope() {
        if (had_) {
            ::setenv(kVar, prev_.c_str(), 1);
            return;
        }
        ::unsetenv(kVar);
    }

  private:
    static constexpr const char* kVar = "CHEATAH_SPACE_IRBEM_SPV_DIR";
    bool had_ = false;
    std::string prev_;
};

/// The registry routes this kernel to the host at EVERY size — measured, see the
/// `irbem_mead_f32` row — so under Lane::Auto the device branch of mead_field_batch is never
/// taken in production. The seam's `CHEATAH_SPACE_IRBEM_GPU_CROSSOVER` override exists so that a
/// suite can still drive that branch and verify it; this scope sets it for one test.
class CrossoverScope {
  public:
    explicit CrossoverScope(const char* value) {
        if (const char* prev = std::getenv(kVar)) {
            had_ = true;
            prev_ = prev;
        }
        ::setenv(kVar, value, 1);
    }
    CrossoverScope(const CrossoverScope&) = delete;
    CrossoverScope& operator=(const CrossoverScope&) = delete;
    CrossoverScope(CrossoverScope&&) = delete;
    CrossoverScope& operator=(CrossoverScope&&) = delete;
    ~CrossoverScope() {
        if (had_) {
            ::setenv(kVar, prev_.c_str(), 1);
            return;
        }
        ::unsetenv(kVar);
    }

  private:
    static constexpr const char* kVar = "CHEATAH_SPACE_IRBEM_GPU_CROSSOVER";
    bool had_ = false;
    std::string prev_;
};

/// The crossover the device tests force, and a batch size comfortably above it.
constexpr const char* kForcedCrossover = "512";
constexpr std::size_t kDeviceBatch = 4 * 512;

}  // namespace

TEST(IrbemMead, TheRegistryRoutesThisKernelToTheHost) {
    // The measured verdict, pinned: ~50 flops per 36 bytes never pays for the bus, so a batch of
    // any size runs on the host under Lane::Auto. A registry edit that "enables" the device for
    // this kernel without a new measurement fails here.
    EXPECT_EQ(ir::gpu::kernel_info("irbem_mead_f32").crossover_points, ir::gpu::never_faster_on_device);
    EXPECT_FALSE(ir::gpu::prefer_gpu("irbem_mead_f32", std::size_t{1} << 24));
    const std::vector<Position<Frame::GSM>> pts = scatter(kDeviceBatch);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    const ir::Result<bool> r = mead_field_batch(pts, 0.28, 23.0, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_FALSE(r.value) << "a never-faster kernel must not reach the device under Lane::Auto";
}

TEST(IrbemMead, BatchFallsBackToTheHostWhenTheShaderWasNeverBuilt) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const CrossoverScope forced(kForcedCrossover);
    const std::size_t n = kDeviceBatch;
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    {
        const SpvDirScope nowhere(std::filesystem::temp_directory_path().string() +
                                  "/cheatah-space-no-such-shader-dir");
        const ir::Result<bool> r = mead_field_batch(pts, 0.28, 23.0, out);
        EXPECT_EQ(r.status, Status::Ok);
        EXPECT_FALSE(r.value) << "with no compiled shader the batch must run on the host";
    }
    const Tilt t = tilt_of(0.28);
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = mead_field_at(pts[i], t.s, t.c, t.deg, mead_kp_bin(23.0));
        ASSERT_EQ(out[i].v[0], ref.v[0]) << "point " << i;
        ASSERT_EQ(out[i].v[1], ref.v[1]) << "point " << i;
        ASSERT_EQ(out[i].v[2], ref.v[2]) << "point " << i;
    }
}

TEST(IrbemMead, BatchUsesTheDeviceWhenOneIsAvailable) {
    // With the crossover forced below the batch size, the device lane runs and SAYS it did.
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const CrossoverScope forced(kForcedCrossover);
    const std::size_t n = kDeviceBatch;
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    const ir::Result<bool> r = mead_field_batch(pts, 0.28, 23.0, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_TRUE(r.value) << "the batch fell back to the host with a device present";
    const Tilt t = tilt_of(0.28);
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = mead_field_at(pts[i], t.s, t.c, t.deg, mead_kp_bin(23.0));
        for (int c = 0; c < 3; ++c) worst = std::max(worst, std::fabs(out[i].v[c] - ref.v[c]));
    }
    std::printf("[ MEASURED ] device batch of %zu vs fp64 reference: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-3);
    const std::vector<Position<Frame::GSM>> few = scatter(16);
    std::vector<ir::FieldVector<Frame::GSM>> few_out(few.size());
    EXPECT_FALSE(mead_field_batch(few, 0.28, 23.0, few_out).value);
}

TEST(IrbemMead, TheDeviceLaneRefusesABadPointBeforeItDispatches) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const CrossoverScope forced(kForcedCrossover);
    const std::size_t n = kDeviceBatch;
    std::vector<Position<Frame::GSM>> pts = scatter(n);
    pts[n / 3] = at(0.25, 0.0, 0.0);
    std::vector<ir::FieldVector<Frame::GSM>> out(n, ir::FieldVector<Frame::GSM>{cheatah::fixarray::vec3d{7.0, 7.0, 7.0}});
    const ir::Result<bool> r = mead_field_batch(pts, 0.28, 23.0, out);
    EXPECT_EQ(r.status, Status::DomainError);
    EXPECT_FALSE(r.value) << "a refused batch must not have reached the device";
    for (const ir::FieldVector<Frame::GSM>& b : out) {
        ASSERT_EQ(b.v[0], 0.0);
        ASSERT_EQ(b.v[1], 0.0);
        ASSERT_EQ(b.v[2], 0.0);
    }
    std::vector<Position<Frame::GSM>> far = scatter(n);
    far[n / 2] = at(-18.0, 0.0, 0.0);
    const ir::Result<bool> f = mead_field_batch(far, 0.28, 23.0, out);
    EXPECT_EQ(f.status, Status::OutOfValidityRange);
    EXPECT_TRUE(f.value) << "an out-of-validity batch must still be computed on the device";
}

TEST(IrbemMead, DeviceKernelAgreesWithTheHostLane) {
    // Independent of the registry's routing: the kernel is dispatched directly and compared with
    // the fp32 host lane it mirrors expression for expression.
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 1 << 16;
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<float> pos(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    const Tilt t = tilt_of(0.31);
    const auto sp = static_cast<float>(t.s);
    const auto cp = static_cast<float>(t.c);
    const auto td = static_cast<float>(t.deg);
    std::vector<float> host(3 * n);
    std::vector<float> device(3 * n);
    ASSERT_TRUE(mead_field_host(pos, host, sp, cp, td, 4));
    const std::array<float, mead_param_count> block = mead_param_block(sp, cp, td, 4);
    ir::gpu::dispatch_batch("irbem_mead_f32", pos, device, std::span<const float>(block));
    double worst = 0.0;
    for (std::size_t i = 0; i < 3 * n; ++i) {
        worst = std::max(worst, std::fabs(static_cast<double>(device[i]) - host[i]));
    }
    std::printf("[ MEASURED ] device vs host, %zu points: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-3) << "the device is not evaluating the same expressions";
}
#endif
