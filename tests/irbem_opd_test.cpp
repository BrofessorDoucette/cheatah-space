/// @file irbem_opd_test.cpp
/// @brief The suite for `space/irbem/ext_opd.hpp` — the Olson-Pfitzer dynamic model.
///
/// The model's form and coefficients were never published (see the header's brief), so there is no
/// analytic truth AND no published reference to compare against — only IRBEM's distributed code,
/// which this clean room may run as a black box and may not read. What can be established, and is:
///
///  - **`div B = 0`**, model-independent and oracle-free: both pieces are curls of published vector
///    potentials, a similarity transform preserves the property, and a central-difference
///    divergence must fall as `h^2` — across all four activity regimes and five tilts.
///  - **The structure the oracle turns out to have**, which this implementation has by
///    construction: `n` and `V` enter only through `n V^2`; the field is affine in Dst; the Dst
///    gradient does not depend on the pressure. Each is asserted here as an identity on THIS
///    implementation, and `tools/oracle/opd_diff.cpp` measures the same three on the oracle.
///  - **Continuity in every driver.** This is a continuous-driver model. A step-function
///    identity like T89's Kp-bin test would be a BUG here, so the sweeps assert the opposite: no
///    step in a fine sweep is out of proportion to its neighbours, and the Dst dependence is
///    exactly linear.
///  - **Exact reduction to the quiet field** at the reference conditions, which is what makes the
///    "compressed quiet field plus a Dst ring" reading of the formula checkable rather than a
///    story.
///  - **The validity envelope from both sides of every bound**, with the value still returned
///    outside — where the oracle returns a sentinel.
///  - **The differential against the oracle**, asserted as the MEASURED regression envelope and
///    never as agreement: the harness's numbers say the two are ~50% apart in the belts, and a
///    tight tolerance here would be a false claim.
///
/// The oracle test `dlopen`s IRBEM at runtime rather than linking it: this binary builds and passes
/// on a machine that has never heard of IRBEM.

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
#include "space/irbem/api.hpp"
#include "space/irbem/ext_opd.hpp"
#include "space/irbem/lstar.hpp"

namespace {

namespace ir = cheatah::space::irbem;
namespace fx = cheatah::fixarray;

using ir::Frame;
using ir::OpdParameters;
using ir::Position;
using ir::Status;
using ir::opd_components;
using ir::opd_constants;
using ir::opd_dst_star;
using ir::opd_field;
using ir::opd_field_at;
using ir::opd_field_batch;
using ir::opd_field_host;
using ir::opd_param_block;
using ir::opd_param_count;
using ir::opd_parameters;
using ir::opd_pressure_npa;
using ir::opd_quiet_bin;
using ir::opd_reference_pressure_npa;
using ir::opd_ring_coefficient;
using ir::opd_ring_increment;
using ir::opd_ring_unit;
using ir::opd_scale;

/// A sink the optimizer cannot see through, so the allocation test's calls actually happen.
volatile double sink = 0.0;

/// A GSM point, spelled so a test reads as coordinates rather than as a constructor call.
Position<Frame::GSM> at(double x, double y, double z) {
    return Position<Frame::GSM>{fx::vec3d(x, y, z)};
}

/// The tilts the suite sweeps: zero, and both signs of a realistic seasonal-diurnal excursion.
constexpr std::array<double, 5> kTilts{-0.55, -0.21, 0.0, 0.21, 0.55};

/// One driver triple: density (cm^-3), speed (km/s), Dst (nT).
struct Drivers {
    double n;
    double v;
    double dst;
};

/// The corpus's four regimes clipped to the model's documented envelope, plus its two most
/// extreme in-envelope corners. The storm and extreme regimes lie OUTSIDE the envelope in `V` and
/// Dst; they are clipped rather than dropped because the edge of the envelope is where a storm
/// study actually runs this model.
constexpr std::array<Drivers, 6> kDrivers{{
    {5.0, 380.0, -8.0},     // quiet
    {8.0, 450.0, -42.0},    // moderate
    {20.0, 500.0, -100.0},  // storm, clipped
    {45.0, 500.0, -100.0},  // extreme, clipped
    {50.0, 500.0, 20.0},    // most compressed corner, northward-most Dst
    {5.0, 300.0, -100.0},   // least compressed corner, deepest Dst
}};

/// `div B` at one point by central differences.
double divergence(const OpdParameters<double>& p, double ps, double x, double y, double z,
                  double h) {
    const double s = std::sin(ps);
    const double c = std::cos(ps);
    const auto b = [&](double a, double d, double e) { return opd_components<double>(p, s, c, a, d, e); };
    return ((b(x + h, y, z)[0] - b(x - h, y, z)[0]) + (b(x, y + h, z)[1] - b(x, y - h, z)[1]) +
            (b(x, y, z + h)[2] - b(x, y, z - h)[2])) /
           (2.0 * h);
}

/// The same parameters with the quiet field's eq. (20) block zeroed. Tsyganenko's `C_16..C_19`
/// are published to four significant figures, so eq. (20) carries an irreducible ~1e-3 nT/R_E of
/// divergence that is a property of the TABLE; with it off, every analytic derivative in both
/// pieces is exposed.
OpdParameters<double> without_eq20(const OpdParameters<double>& p) {
    OpdParameters<double> q = p;
    for (std::size_t k = 5; k < ir::t89_linear_count; ++k) q.quiet.c[k] = 0.0;
    return q;
}

/// The largest `|div B|` over the sampled box, every regime and tilt, at difference step @p h.
double worst_divergence(double h, bool with_eq20) {
    double worst = 0.0;
    for (double ps : kTilts) {
        for (const Drivers& d : kDrivers) {
            const OpdParameters<double> full = opd_parameters<double>(d.n, d.v, d.dst);
            const OpdParameters<double> p = with_eq20 ? full : without_eq20(full);
            // Integer induction, coordinates derived per iteration (cert-flp30).
            for (int ix = 0; ix <= 10; ++ix)
                for (int iy = 0; iy <= 6; ++iy)
                    for (int iz = 0; iz <= 5; ++iz) {
                        const double x = -25.0 + (3.7 * ix);
                        const double y = -11.0 + (3.3 * iy);
                        const double z = -9.0 + (3.1 * iz);
                        if (std::sqrt((x * x) + (y * y) + (z * z)) < 2.0) continue;
                        worst = std::max(worst, std::fabs(divergence(p, ps, x, y, z, h)));
                    }
        }
    }
    return worst;
}

/// A deterministic scatter of GSM points over the inner magnetosphere and near tail.
std::vector<Position<Frame::GSM>> scatter(std::size_t n, double r_lo = 2.5, double r_hi = 20.0) {
    std::vector<Position<Frame::GSM>> out;
    out.reserve(n);
    std::uint64_t s = 0x9E3779B97F4A7C15ULL;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) / 9007199254740992.0;
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double r = r_lo + ((r_hi - r_lo) * next());
        const double th = std::acos(1.0 - (2.0 * next()));
        const double ph = 6.283185307179586 * next();
        out.push_back(at(r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph),
                         r * std::cos(th)));
    }
    return out;
}

/// `|a - b|` for two field triples.
double dist(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return std::sqrt(((a[0] - b[0]) * (a[0] - b[0])) + ((a[1] - b[1]) * (a[1] - b[1])) +
                     ((a[2] - b[2]) * (a[2] - b[2])));
}

/// `|v|` of a field triple.
double norm(const std::array<double, 3>& v) {
    return std::sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
}

/// `|v|` of a fixarray vector.
double mag(const fx::vec3d& v) { return std::sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2])); }

/// The epoch rotations the total-field tests share.
ir::Rotations epoch_rotations(const ir::Igrf<10>& m) {
    const auto r = ir::api::rotations_at(2015, 180, 43200.0, m);
    EXPECT_EQ(Status::Ok, r.status);
    return r.value;
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

/// Inside the Shue et al. (1997) magnetopause with a 10% margin, at the pressure-balance standoff
/// — the same mask `tools/oracle/opd_diff.cpp` applies, because outside the boundary the oracle's
/// quiet-model exponentials run away and a comparison there measures nothing about either model.
bool inside_magnetopause(double x, double y, double z, double n, double v) {
    const double r = std::sqrt((x * x) + (y * y) + (z * z));
    const double r0 = 11.15 * std::pow(opd_pressure_npa(n, v) / opd_reference_pressure_npa(), -1.0 / 6.0);
    const double mp = r0 * std::pow(2.0 / (1.0 + (x / r)), 0.58);
    return r < 0.9 * mp;
}

}  // namespace

// ================================================================================================
// The published constants and the driver algebra
// ================================================================================================

TEST(IrbemOpd, TheConstantsAreThePublishedOnes) {
    // CODATA 2018 m_p = 1.67262192369e-27 kg, in nPa per (cm^-3 km^2 s^-2).
    EXPECT_EQ(opd_constants.pressure_npa_per_cc_km2s2, 1.67262192369e-6);
    // O'Brien & McPherron (2000): Dst* = Dst - 7.26 sqrt(Pdyn) + 11.
    EXPECT_EQ(opd_constants.dst_pressure_coefficient, 7.26);
    EXPECT_EQ(opd_constants.dst_offset, 11.0);
    // The reference: the envelope's quiet end in density, the canonical nominal speed.
    EXPECT_EQ(opd_constants.reference_density, 5.0);
    EXPECT_EQ(opd_constants.reference_velocity, 400.0);
    // The quiet field is the first published column of Tsyganenko (1989) Table 1, and nothing
    // else — this is what keeps the model free of Kp steps.
    EXPECT_EQ(opd_quiet_bin, 1);
    EXPECT_TRUE(ir::t89_bin_is_published(opd_quiet_bin));
}

TEST(IrbemOpd, DynamicPressureIsTheProtonRamPressure) {
    // n = 5 cm^-3, V = 400 km/s: 1.67262192369e-6 * 5 * 160000 = 1.338097538952 nPa.
    EXPECT_NEAR(opd_pressure_npa(5.0, 400.0), 1.338097538952, 1e-12);
    EXPECT_EQ(opd_reference_pressure_npa(), opd_pressure_npa(5.0, 400.0));
    // Quadratic in V, linear in n: exactly representable ratios.
    EXPECT_EQ(opd_pressure_npa(10.0, 400.0), 2.0 * opd_pressure_npa(5.0, 400.0));
    EXPECT_EQ(opd_pressure_npa(5.0, 800.0), 4.0 * opd_pressure_npa(5.0, 400.0));
    // The envelope's corners, for the record: 0.75 nPa to 20.9 nPa.
    std::printf("[ MEASURED ] Pdyn at (5, 300) = %.3f nPa, at (50, 500) = %.3f nPa\n",
                opd_pressure_npa(5.0, 300.0), opd_pressure_npa(50.0, 500.0));
    EXPECT_NEAR(opd_pressure_npa(5.0, 300.0), 0.7527, 1e-3);
    EXPECT_NEAR(opd_pressure_npa(50.0, 500.0), 20.908, 1e-3);
}

TEST(IrbemOpd, CompressionFollowsTheSixthRootOfPressure) {
    // s = (P / P_0)^(1/6): exactly 1 at the reference, and the sixth root of an exact ratio
    // elsewhere. 64 x the reference pressure is n = 320 at V = 400 (or n = 5 at V = 3200): s = 2.
    EXPECT_EQ(opd_scale(5.0, 400.0), 1.0);
    EXPECT_NEAR(opd_scale(320.0, 400.0), 2.0, 1e-12);
    EXPECT_NEAR(opd_scale(5.0, 3200.0), 2.0, 1e-12);
    // The documented envelope's corners.
    std::printf("[ MEASURED ] s at (5, 300) = %.4f, at (50, 500) = %.4f\n", opd_scale(5.0, 300.0),
                opd_scale(50.0, 500.0));
    EXPECT_NEAR(opd_scale(5.0, 300.0), std::pow(0.5625, 1.0 / 6.0), 1e-12);
    EXPECT_NEAR(opd_scale(50.0, 500.0), std::pow(15.625, 1.0 / 6.0), 1e-12);
    // Monotone in both drivers, which is what "compression" means.
    EXPECT_LT(opd_scale(5.0, 400.0), opd_scale(6.0, 400.0));
    EXPECT_LT(opd_scale(5.0, 400.0), opd_scale(5.0, 401.0));
}

TEST(IrbemOpd, DstIsPressureCorrectedBeforeItDrivesTheRing) {
    // Dst* = Dst - 7.26 sqrt(P) + 11, at P = 4 nPa exactly: -50 - 14.52 + 11 = -53.52.
    EXPECT_NEAR(opd_dst_star(-50.0, 4.0), -53.52, 1e-12);
    EXPECT_NEAR(opd_dst_star(0.0, 0.0), 11.0, 1e-12);
    // The ring increment is the DIFFERENCE from the quiet reference: zero there, the offset gone.
    EXPECT_EQ(opd_ring_increment(5.0, 400.0, 0.0), 0.0);
    EXPECT_NEAR(opd_ring_increment(5.0, 400.0, -100.0), -100.0, 1e-12);
    // A pressure rise alone LOWERS the increment (the magnetopause contribution is subtracted),
    // which is the sign O'Brien & McPherron's correction has.
    EXPECT_LT(opd_ring_increment(50.0, 400.0, 0.0), 0.0);
    EXPECT_NEAR(opd_ring_increment(50.0, 400.0, 0.0),
                -7.26 * (std::sqrt(opd_pressure_npa(50.0, 400.0)) -
                         std::sqrt(opd_reference_pressure_npa())),
                1e-12);
}

TEST(IrbemOpd, TheRingIsNormalisedToDstAtTheCentre) {
    // Dessler-Parker-Sckopke: the ring's field at the Earth's centre IS its Dst* increment, and it
    // points along the DIPOLE axis — the ring lies in the dipole equator. With
    // C = ΔDst* (a + D)^3 / 2 the unit disc's field at the origin, projected onto the dipole axis
    // (sin psi, 0, cos psi) in GSM, times C must be exactly ΔDst* at every tilt, with nothing
    // across the axis. Asserting GSM B_z alone would be wrong by cos(psi) at any tilt but zero.
    const double a = 8.161;
    const double d = 2.08;
    for (double ps : kTilts) {
        const double sn = std::sin(ps);
        const double cs = std::cos(ps);
        const std::array<double, 3> unit = opd_ring_unit<double>(a, d, sn, cs, 0.0, 0.0, 0.0);
        const double along = (unit[0] * sn) + (unit[2] * cs);
        const double across = (unit[0] * cs) - (unit[2] * sn);
        EXPECT_NEAR(along * opd_ring_coefficient(1.0, a, d), 1.0, 1e-12) << "tilt " << ps;
        EXPECT_NEAR(along * opd_ring_coefficient(-73.0, a, d), -73.0, 1e-10) << "tilt " << ps;
        EXPECT_NEAR(across, 0.0, 1e-15) << "tilt " << ps;
        EXPECT_EQ(unit[1], 0.0);
        if (ps == 0.0) EXPECT_EQ(unit[0], 0.0);
    }
    // The shape the normalisation buys: at geosynchronous noon the same unit ring is a third of
    // its central value, decaying to nothing in the far tail — the measured oracle's dB/dDst at
    // (6.6, 0, 0) is 0.30 nT/nT and at (-20, 0, 2) is -0.03; this shape gives 0.33 and -0.01.
    const double c = opd_ring_coefficient(1.0, a, d);
    const double geo = opd_ring_unit<double>(a, d, 0.0, 1.0, 6.6, 0.0, 0.0)[2] * c;
    const double tail = opd_ring_unit<double>(a, d, 0.0, 1.0, -20.0, 0.0, 2.0)[2] * c;
    std::printf("[ MEASURED ] unit ring B_z: centre 1.000, GEO noon %.4f, (-20, 0, 2) %+.4f nT/nT\n", geo,
                tail);
    EXPECT_NEAR(geo, 0.3325, 5e-3);
    EXPECT_LT(std::fabs(tail), 0.05);
    // And the geometry IS the quiet set's, not a second copy of two numbers.
    const OpdParameters<double> p = opd_parameters<double>(10.0, 400.0, 0.0);
    EXPECT_EQ(p.quiet.a_rc, a);
    EXPECT_EQ(p.quiet.d_0, d);
}

TEST(IrbemOpd, ParametersRoundTripThroughFloat) {
    for (const Drivers& d : kDrivers) {
        const OpdParameters<double> dp = opd_parameters<double>(d.n, d.v, d.dst);
        const OpdParameters<float> fp = opd_parameters<float>(d.n, d.v, d.dst);
        // The float block is the double block rounded ONCE, not float arithmetic.
        EXPECT_EQ(fp.scale, static_cast<float>(dp.scale));
        EXPECT_EQ(fp.ring, static_cast<float>(dp.ring));
        for (std::size_t k = 0; k < ir::t89_linear_count; ++k) {
            EXPECT_EQ(fp.quiet.c[k], static_cast<float>(dp.quiet.c[k])) << "C_" << (k + 1);
        }
        EXPECT_EQ(fp.quiet.a_rc, static_cast<float>(dp.quiet.a_rc));
        EXPECT_EQ(fp.quiet.d_0, static_cast<float>(dp.quiet.d_0));
        // And the double block is the published quiet set itself.
        EXPECT_EQ(dp.quiet.c[0], ir::t89_coefficient_sets[0].c[0]);
        EXPECT_EQ(dp.scale, opd_scale(d.n, d.v));
        EXPECT_EQ(dp.ring, opd_ring_coefficient(opd_ring_increment(d.n, d.v, d.dst), dp.quiet.a_rc,
                                                dp.quiet.d_0));
    }
}

// ================================================================================================
// The mathematics: div B = 0, the reduction, the symmetries
// ================================================================================================

TEST(IrbemOpd, DivergenceVanishesEverywhere) {
    // Both pieces are curls and the similarity transform preserves that, so with the quiet
    // field's eq. (20) rounding switched off the central-difference divergence must fall as h^2
    // across every regime and tilt: successive ratios of ~100, checked to be at least 50 to leave
    // room for the fp64 roundoff term that eventually takes over.
    const double d2 = worst_divergence(1e-2, false);
    const double d3 = worst_divergence(1e-3, false);
    const double d4 = worst_divergence(1e-4, false);
    std::printf("[ MEASURED ] worst |div B| without eq.(20): h=1e-2 %.3e  h=1e-3 %.3e  h=1e-4 %.3e"
                " nT/Re\n",
                d2, d3, d4);
    EXPECT_GT(d2 / d3, 50.0) << "divergence is not falling as h^2 — a derivative or the transform is wrong";
    EXPECT_GT(d3 / d4, 50.0) << "divergence is not falling as h^2 — a derivative or the transform is wrong";
    EXPECT_LT(d4, 1e-6);

    // With eq. (20) back in, the residual is the published table's four-figure rounding, scaled
    // by s^4 under the transform (at most 6.25 at the envelope's most compressed corner): bounded,
    // and h-independent, which is the signature of a rounding and not of a truncation error.
    const double f3 = worst_divergence(1e-3, true);
    const double f4 = worst_divergence(1e-4, true);
    std::printf("[ MEASURED ] worst |div B| WITH eq.(20): h=1e-3 %.3e  h=1e-4 %.3e nT/Re\n", f3, f4);
    EXPECT_LT(f4, 2e-2);
    EXPECT_NEAR(f3, f4, 3e-4) << "the eq.(20) residual should be h-independent";
}

TEST(IrbemOpd, ReferenceConditionsReduceToTheQuietField) {
    // At (5, 400, 0) the scale is exactly 1 and the ring increment exactly 0, so the model IS
    // Tsyganenko (1989) Kp bin 1, bit for bit — which is what makes "a compressed quiet field plus
    // a Dst ring" a checkable statement rather than a story.
    const OpdParameters<double> p = opd_parameters<double>(5.0, 400.0, 0.0);
    EXPECT_EQ(p.scale, 1.0);
    EXPECT_EQ(p.ring, 0.0);
    for (double ps : kTilts) {
        const double s = std::sin(ps);
        const double c = std::cos(ps);
        for (const Position<Frame::GSM>& q : scatter(64)) {
            const fx::vec3d mine = opd_field_at(q, s, c, 5.0, 400.0, 0.0).v;
            const fx::vec3d quiet = ir::t89_field_at(q, s, c, opd_quiet_bin).v;
            EXPECT_EQ(mine[0], quiet[0]);
            EXPECT_EQ(mine[1], quiet[1]);
            EXPECT_EQ(mine[2], quiet[2]);
        }
    }
    // And away from the reference it is NOT the quiet field: the drivers do something.
    const fx::vec3d moved = opd_field_at(at(6.6, 0.0, 0.0), 0.0, 1.0, 20.0, 450.0, -60.0).v;
    const fx::vec3d quiet = ir::t89_field_at(at(6.6, 0.0, 0.0), 0.0, 1.0, opd_quiet_bin).v;
    EXPECT_GT(std::fabs(moved[2] - quiet[2]), 1.0);
}

TEST(IrbemOpd, ZeroTiltIsMirrorSymmetricAboutTheEquator) {
    // At psi = 0 both pieces are even about z = 0 in B_z and odd in B_x, B_y — bitwise, because
    // the arithmetic on the two sides differs by one sign bit.
    for (const Drivers& d : kDrivers) {
        const OpdParameters<double> p = opd_parameters<double>(d.n, d.v, d.dst);
        for (double x : {-18.0, -6.5, -1.5, 3.25, 9.0}) {
            for (double y : {-7.5, -1.25, 0.0, 2.5, 8.0}) {
                for (double z : {0.5, 2.25, 6.0}) {
                    const std::array<double, 3> up = opd_components<double>(p, 0.0, 1.0, x, y, z);
                    const std::array<double, 3> dn = opd_components<double>(p, 0.0, 1.0, x, y, -z);
                    EXPECT_EQ(up[0], -dn[0]);
                    EXPECT_EQ(up[1], -dn[1]);
                    EXPECT_EQ(up[2], dn[2]);
                }
            }
        }
    }
}

TEST(IrbemOpd, DawnDuskSymmetryHoldsAtEveryTilt) {
    // Neither the quiet field nor a symmetric ring has a dawn-dusk asymmetry: y -> -y gives
    // (B_x, -B_y, B_z) exactly at ANY tilt, and on the noon-midnight meridian B_y is exactly zero.
    for (double ps : kTilts) {
        const double s = std::sin(ps);
        const double c = std::cos(ps);
        for (const Drivers& d : kDrivers) {
            const OpdParameters<double> p = opd_parameters<double>(d.n, d.v, d.dst);
            for (double x : {-16.0, -4.5, 2.0, 8.75}) {
                for (double y : {0.75, 3.5, 9.25}) {
                    for (double z : {-5.5, -1.0, 0.0, 4.25}) {
                        const std::array<double, 3> dusk = opd_components<double>(p, s, c, x, y, z);
                        const std::array<double, 3> dawn = opd_components<double>(p, s, c, x, -y, z);
                        EXPECT_EQ(dusk[0], dawn[0]);
                        EXPECT_EQ(dusk[1], -dawn[1]);
                        EXPECT_EQ(dusk[2], dawn[2]);
                    }
                    EXPECT_EQ(opd_components<double>(p, s, c, x, 0.0, 1.5)[1], 0.0);
                }
            }
        }
    }
}

// ================================================================================================
// The drivers: what enters, how, and how smoothly
// ================================================================================================

TEST(IrbemOpd, SolarWindEntersOnlyThroughDynamicPressure) {
    // (20, 400), (80, 200) and (5, 800) have the SAME n V^2 = 3.2e6, so the field must be the
    // same to roundoff — the model may see no other combination of the two. The oracle has this
    // property bit for bit (tools/oracle/opd_diff.cpp, section 1a); here the three products round
    // through different intermediates, so the identity is to 1e-12 relative rather than exact.
    for (double ps : {0.0, 0.3}) {
        const double s = std::sin(ps);
        const double c = std::cos(ps);
        for (const Position<Frame::GSM>& q : scatter(48)) {
            const fx::vec3d a = opd_field_at(q, s, c, 20.0, 400.0, -30.0).v;
            const fx::vec3d b = opd_field_at(q, s, c, 80.0, 200.0, -30.0).v;
            const fx::vec3d e = opd_field_at(q, s, c, 5.0, 800.0, -30.0).v;
            const double scale = std::fabs(a[0]) + std::fabs(a[1]) + std::fabs(a[2]);
            for (int k = 0; k < 3; ++k) {
                EXPECT_NEAR(a[k], b[k], 1e-12 * scale);
                EXPECT_NEAR(a[k], e[k], 1e-12 * scale);
            }
        }
    }
    // And a different product is a different field: the drivers are not ignored.
    EXPECT_GT(std::fabs(opd_field_at(at(6.6, 0.0, 0.0), 0.0, 1.0, 20.0, 400.0, 0.0).v[2] -
                        opd_field_at(at(6.6, 0.0, 0.0), 0.0, 1.0, 20.0, 450.0, 0.0).v[2]),
              1.0);
}

TEST(IrbemOpd, FieldIsAffineInDstWithAPressureIndependentGradient) {
    // Dst enters only through the ring coefficient, linearly, and the ring's SHAPE does not read
    // the pressure: so B(Dst) is affine with a gradient that is the same at n = 5 and n = 50. Both
    // are exact properties of the oracle too (opd_diff.cpp, sections 1b and 1c).
    for (const Position<Frame::GSM>& q : scatter(48)) {
        const std::array<double, 3> lo = opd_components<double>(opd_parameters<double>(10.0, 450.0, -100.0), 0.2, std::cos(0.2), q.v[0], q.v[1], q.v[2]);
        const std::array<double, 3> hi = opd_components<double>(opd_parameters<double>(10.0, 450.0, 0.0), 0.2, std::cos(0.2), q.v[0], q.v[1], q.v[2]);
        const std::array<double, 3> mid = opd_components<double>(opd_parameters<double>(10.0, 450.0, -50.0), 0.2, std::cos(0.2), q.v[0], q.v[1], q.v[2]);
        const std::array<double, 3> lin{0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])};
        EXPECT_LT(dist(mid, lin), 1e-11 * (norm(hi) + 1.0));

        std::array<double, 3> g_lo{};
        std::array<double, 3> g_hi{};
        for (int c = 0; c < 3; ++c) {
            const std::size_t k = static_cast<std::size_t>(c);
            g_lo[k] = (opd_components<double>(opd_parameters<double>(5.0, 400.0, -100.0), 0.2, std::cos(0.2), q.v[0], q.v[1], q.v[2])[k] -
                       opd_components<double>(opd_parameters<double>(5.0, 400.0, 0.0), 0.2, std::cos(0.2), q.v[0], q.v[1], q.v[2])[k]) / -100.0;
            g_hi[k] = (opd_components<double>(opd_parameters<double>(50.0, 400.0, -100.0), 0.2, std::cos(0.2), q.v[0], q.v[1], q.v[2])[k] -
                       opd_components<double>(opd_parameters<double>(50.0, 400.0, 0.0), 0.2, std::cos(0.2), q.v[0], q.v[1], q.v[2])[k]) / -100.0;
        }
        EXPECT_LT(dist(g_lo, g_hi), 1e-11);
        // The gradient is the unit ring, exactly: C per nT is (a + D)^3 / 2.
        const std::array<double, 3> unit = opd_ring_unit<double>(8.161, 2.08, 0.2, std::cos(0.2), q.v[0], q.v[1], q.v[2]);
        const double per_nt = opd_ring_coefficient(1.0, 8.161, 2.08);
        for (int c = 0; c < 3; ++c) {
            const std::size_t k = static_cast<std::size_t>(c);
            EXPECT_NEAR(g_lo[k], unit[k] * per_nt, 1e-9);
        }
    }
}

TEST(IrbemOpd, TheDstGradientIsTheRingShapeAtEveryPressure) {
    // The same statement through the public entry point, over the envelope's pressure range and
    // both tilt signs, at the points a storm study looks at.
    for (double ps : {-0.4, 0.0, 0.4}) {
        for (double n : {5.0, 12.0, 30.0, 50.0}) {
            for (double x : {6.6, -6.6, 0.0, 4.0}) {
                const Position<Frame::GSM> q = at(x, x == 0.0 ? 6.6 : 1.0, 0.5);
                const ir::Result<ir::FieldVector<Frame::GSM>> a = opd_field(q, ps, n, 400.0, 0.0);
                const ir::Result<ir::FieldVector<Frame::GSM>> b = opd_field(q, ps, n, 400.0, -100.0);
                ASSERT_EQ(a.status, Status::Ok);
                ASSERT_EQ(b.status, Status::Ok);
                const std::array<double, 3> unit =
                    opd_ring_unit<double>(8.161, 2.08, std::sin(ps), std::cos(ps), q.v[0], q.v[1], q.v[2]);
                const double per_nt = opd_ring_coefficient(1.0, 8.161, 2.08);
                for (int c = 0; c < 3; ++c) {
                    EXPECT_NEAR((b.value.v[c] - a.value.v[c]) / -100.0,
                                unit[static_cast<std::size_t>(c)] * per_nt, 1e-9)
                        << "n " << n << " tilt " << ps << " x " << x;
                }
            }
        }
    }
}

TEST(IrbemOpd, TheFieldIsSmoothInEveryDriver) {
    // This is a continuous-driver model, and a Kp-style step anywhere in it would be a bug. Each
    // driver is swept finely across its whole envelope at four points, and every step of the
    // sweep is held to within 25% of its NEIGHBOURS: a smooth function's adjacent steps differ by
    // O(h) (for the sixth-root response here, by ~4% at the steepest end), a step function's
    // differ by the whole jump. Comparing to the mean step instead would be wrong for a smooth
    // but strongly curved response, which this one is. Dst is exactly linear and gets the tighter
    // check.
    const std::array<Position<Frame::GSM>, 4> pts{at(6.6, 0.0, 0.0), at(-6.6, 0.0, 0.0),
                                                  at(4.0, 3.0, 2.0), at(-12.0, 4.0, -2.0)};
    for (const Position<Frame::GSM>& q : pts) {
        // Density, 5 -> 50 in 180 steps of 0.25.
        {
            std::vector<double> steps;
            fx::vec3d prev = opd_field(q, 0.25, 5.0, 400.0, -30.0).value.v;
            for (int k = 1; k <= 180; ++k) {
                const fx::vec3d cur = opd_field(q, 0.25, 5.0 + (0.25 * k), 400.0, -30.0).value.v;
                steps.push_back(mag(cur - prev));
                prev = cur;
            }
            for (std::size_t k = 1; k + 1 < steps.size(); ++k) {
                EXPECT_LT(steps[k], 1.25 * std::max(steps[k - 1], steps[k + 1]))
                    << "a jump in the density sweep at n = " << (5.0 + (0.25 * static_cast<double>(k + 1)));
            }
            EXPECT_GT(steps[0], 0.0) << "density must move the field";
        }
        // Speed, 300 -> 500 in 200 steps of 1 km/s.
        {
            std::vector<double> steps;
            fx::vec3d prev = opd_field(q, 0.25, 12.0, 300.0, -30.0).value.v;
            for (int k = 1; k <= 200; ++k) {
                const fx::vec3d cur = opd_field(q, 0.25, 12.0, 300.0 + k, -30.0).value.v;
                steps.push_back(mag(cur - prev));
                prev = cur;
            }
            for (std::size_t k = 1; k + 1 < steps.size(); ++k) {
                EXPECT_LT(steps[k], 1.25 * std::max(steps[k - 1], steps[k + 1]))
                    << "a jump in the speed sweep at V = " << (300.0 + static_cast<double>(k + 1));
            }
            EXPECT_GT(steps[0], 0.0) << "speed must move the field";
        }
        // Dst, -100 -> +20 in 120 steps of 1 nT: exactly linear, so every step is the same step.
        {
            fx::vec3d prev = opd_field(q, 0.25, 12.0, 400.0, -100.0).value.v;
            const fx::vec3d first = opd_field(q, 0.25, 12.0, 400.0, -99.0).value.v - prev;
            prev = opd_field(q, 0.25, 12.0, 400.0, -99.0).value.v;
            for (int k = 2; k <= 120; ++k) {
                const fx::vec3d cur = opd_field(q, 0.25, 12.0, 400.0, -100.0 + k).value.v;
                const fx::vec3d step = cur - prev;
                EXPECT_NEAR(step[0], first[0], 1e-9);
                EXPECT_NEAR(step[1], first[1], 1e-9);
                EXPECT_NEAR(step[2], first[2], 1e-9);
                prev = cur;
            }
            EXPECT_GT(mag(first), 0.0) << "Dst must move the field";
        }
    }
}

TEST(IrbemOpd, CompressionScalesTheCentralFieldAsTheCubeOfTheStandoff) {
    // Mead (1964): the field the boundary currents produce at the Earth's centre scales as the
    // inverse CUBE of the standoff distance, `R_mp^-3` — equivalently, with `R_mp ∝ P^(-1/6)`,
    // as the square root of the dynamic pressure, which is the very law Burton et al.'s `b sqrt(P)`
    // term in Dst* encodes. This is what fixes the exponent in `s^3 B_q(s r)`: any other power of
    // `s` passes div-B (a similarity transform is solenoidal at every exponent), passes the
    // reference reduction (s = 1 there) and slips under the oracle envelope (measured: s^2 moves
    // the extreme-regime RMS from 67.7 to 75.5 nT against a 95 nT cap). So it is pinned here, at
    // the centre, with the Dst-driven ring switched off by choosing Dst so that ΔDst* = 0.
    const OpdParameters<double> quiet = opd_parameters<double>(5.0, 400.0, 0.0);
    const std::array<double, 3> b0 = opd_components<double>(quiet, 0.3, std::cos(0.3), 0.0, 0.0, 0.0);
    ASSERT_GT(norm(b0), 1.0) << "the quiet boundary field at the centre is tens of nT";
    for (double n : {10.0, 20.0, 40.0, 50.0}) {
        // The Dst that makes the ring increment vanish: Dst = 7.26 (sqrt(P) - sqrt(P_0)).
        const double dst = opd_constants.dst_pressure_coefficient *
                           (std::sqrt(opd_pressure_npa(n, 400.0)) - std::sqrt(opd_reference_pressure_npa()));
        const OpdParameters<double> p = opd_parameters<double>(n, 400.0, dst);
        ASSERT_NEAR(p.ring, 0.0, 1e-9) << "n " << n;
        const std::array<double, 3> b = opd_components<double>(p, 0.3, std::cos(0.3), 0.0, 0.0, 0.0);
        const double ratio = norm(b) / norm(b0);
        const double s = opd_scale(n, 400.0);
        EXPECT_NEAR(ratio, s * s * s, 1e-9 * (s * s * s)) << "n " << n << ": not Mead's R_mp^-3";
        EXPECT_NEAR(ratio, std::sqrt(n / 5.0), 1e-9 * ratio) << "n " << n << ": not sqrt(P)";
        // And the direction at the centre is unchanged by a similarity transform: the same
        // boundary field, only stronger.
        for (int c = 0; c < 3; ++c) {
            const std::size_t k = static_cast<std::size_t>(c);
            EXPECT_NEAR(b[k], (s * s * s) * b0[k], 1e-9 * norm(b));
        }
    }
}

TEST(IrbemOpd, CompressionRaisesTheDaysideFieldAndDeepensTheRing) {
    // The physics the drivers exist for. Raising the pressure compresses the magnetopause, so the
    // external B_z at geosynchronous NOON must rise monotonically with n V^2 — the oracle goes from
    // 13.2 nT at (5, 400) to 69.2 nT at (50, 400); this model from 8.0 to 61.2. Deepening Dst
    // must LOWER B_z everywhere in the inner magnetosphere, noon and midnight alike.
    double previous = -1e9;
    for (double n : {5.0, 10.0, 20.0, 40.0, 50.0}) {
        const double bz = opd_field(at(6.6, 0.0, 0.0), 0.0, n, 400.0, 0.0).value.v[2];
        std::printf("[ MEASURED ] n = %2.0f cm^-3, V = 400: external Bz at GSM noon 6.6 Re is %+7.3f nT\n", n, bz);
        EXPECT_GT(bz, previous);
        previous = bz;
    }
    EXPECT_GT(opd_field(at(6.6, 0.0, 0.0), 0.0, 5.0, 400.0, 0.0).value.v[2], 0.0) << "dayside compression";
    EXPECT_LT(opd_field(at(-6.6, 0.0, 0.0), 0.0, 5.0, 400.0, 0.0).value.v[2], 0.0) << "nightside stretching";
    for (double x : {6.6, -6.6, 4.0, -4.0}) {
        EXPECT_LT(opd_field(at(x, 0.0, 0.0), 0.0, 10.0, 400.0, -100.0).value.v[2],
                  opd_field(at(x, 0.0, 0.0), 0.0, 10.0, 400.0, 0.0).value.v[2])
            << "a deeper Dst must depress Bz at x = " << x;
    }
}

// ================================================================================================
// The API surface
// ================================================================================================

TEST(IrbemOpd, ReferenceLaneMatchesTheComponentForm) {
    for (const Drivers& d : kDrivers) {
        const OpdParameters<double> p = opd_parameters<double>(d.n, d.v, d.dst);
        for (double ps : kTilts) {
            const double s = std::sin(ps);
            const double c = std::cos(ps);
            const std::array<double, 3> raw = opd_components<double>(p, s, c, 4.5, -2.25, 1.5);
            const fx::vec3d wrapped = opd_field_at(at(4.5, -2.25, 1.5), s, c, d.n, d.v, d.dst).v;
            EXPECT_EQ(raw[0], wrapped[0]);
            EXPECT_EQ(raw[1], wrapped[1]);
            EXPECT_EQ(raw[2], wrapped[2]);
        }
    }
}

TEST(IrbemOpd, ValidityIsReportedFromBothSidesOfEveryBound) {
    // IRBEM's kext table: 5 <= Dsw <= 50, 300 <= Vsw <= 500, -100 <= Dst <= 20, r <= 60 Re. The
    // bounds are closed: sitting ON one is Ok, strictly past it is OutOfValidityRange with the
    // value still computed. The oracle, measured, accepts the bounds and refuses 1e-4 past them.
    const Position<Frame::GSM> q = at(-6.6, 1.0, 0.5);
    const auto check = [&](double n, double v, double dst, Status expected) {
        const ir::Result<ir::FieldVector<Frame::GSM>> r = opd_field(q, 0.2, n, v, dst);
        EXPECT_EQ(r.status, expected) << "n " << n << " V " << v << " Dst " << dst;
        EXPECT_NE(r.value.v[2], 0.0) << "the value must still be computed";
        EXPECT_TRUE(std::isfinite(r.value.v[2]));
    };
    check(5.0, 400.0, 0.0, Status::Ok);
    check(4.999, 400.0, 0.0, Status::OutOfValidityRange);
    check(50.0, 400.0, 0.0, Status::Ok);
    check(50.001, 400.0, 0.0, Status::OutOfValidityRange);
    check(10.0, 300.0, 0.0, Status::Ok);
    check(10.0, 299.9, 0.0, Status::OutOfValidityRange);
    check(10.0, 500.0, 0.0, Status::Ok);
    check(10.0, 500.1, 0.0, Status::OutOfValidityRange);
    check(10.0, 400.0, 20.0, Status::Ok);
    check(10.0, 400.0, 20.1, Status::OutOfValidityRange);
    check(10.0, 400.0, -100.0, Status::Ok);
    check(10.0, 400.0, -100.1, Status::OutOfValidityRange);
    // The corners, all four drivers at once, from inside.
    check(50.0, 500.0, 20.0, Status::Ok);
    check(5.0, 300.0, -100.0, Status::Ok);
    // Radius: 60 Re is the published limit. Down the tail, where the model is still finite.
    EXPECT_EQ(opd_field(at(-59.9, 0.0, 2.0), 0.1, 10.0, 400.0, -20.0).status, Status::Ok);
    const ir::Result<ir::FieldVector<Frame::GSM>> far = opd_field(at(-60.1, 0.0, 2.0), 0.1, 10.0, 400.0, -20.0);
    EXPECT_EQ(far.status, Status::OutOfValidityRange);
    EXPECT_TRUE(std::isfinite(far.value.v[0]));
    EXPECT_NE(far.value.v[0], 0.0);
    // Inside the solid Earth is a different failure: there is no answer, not an unreliable one.
    EXPECT_EQ(opd_field(at(0.5, 0.0, 0.0), 0.1, 10.0, 400.0, -20.0).status, Status::DomainError);
    // A driver caveat and a position caveat together: the driver's is reported first, being the
    // cheaper and more fundamental check.
    EXPECT_EQ(opd_field(at(-60.1, 0.0, 2.0), 0.1, 60.0, 400.0, -20.0).status, Status::OutOfValidityRange);
}

TEST(IrbemOpd, NonFiniteInputIsADomainError) {
    const double nan = std::nan("");
    const double inf = std::numeric_limits<double>::infinity();
    for (const Position<Frame::GSM>& p : {at(nan, 0.0, 0.0), at(0.0, inf, 0.0), at(4.0, 0.0, -nan)}) {
        const ir::Result<ir::FieldVector<Frame::GSM>> r = opd_field(p, 0.2, 10.0, 400.0, -20.0);
        EXPECT_EQ(r.status, Status::DomainError);
        EXPECT_EQ(r.value.v[0], 0.0);
        EXPECT_EQ(r.value.v[1], 0.0);
        EXPECT_EQ(r.value.v[2], 0.0);
    }
    EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), nan, 10.0, 400.0, -20.0).status, Status::DomainError);
    EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), 0.2, nan, 400.0, -20.0).status, Status::DomainError);
    EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), 0.2, 10.0, inf, -20.0).status, Status::DomainError);
    EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), 0.2, 10.0, 400.0, -inf).status, Status::DomainError);
}

TEST(IrbemOpd, AVanishingSolarWindIsADomainError) {
    // No density or no speed means no pressure, no magnetopause and no compression to compute:
    // s would be zero (or, for a negative density, the sixth root of a negative number). These
    // are refused as having no answer, distinct from a small-but-positive value outside the
    // envelope, which is merely an extrapolation.
    for (double n : {0.0, -5.0}) {
        EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), 0.2, n, 400.0, -20.0).status, Status::DomainError);
    }
    for (double v : {0.0, -400.0}) {
        EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), 0.2, 10.0, v, -20.0).status, Status::DomainError);
    }
    EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), 0.2, 1e-3, 400.0, -20.0).status, Status::OutOfValidityRange);
    EXPECT_TRUE(std::isfinite(opd_field(at(5.0, 0.0, 0.0), 0.2, 1e-3, 400.0, -20.0).value.v[2]));
}

TEST(IrbemOpd, RightAngleTiltIsADomainError) {
    // The quiet field carries tan(psi); at |psi| = pi/2 it does not exist.
    const double quarter_turn = std::numbers::pi / 2.0;
    EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), quarter_turn, 10.0, 400.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), -quarter_turn, 10.0, 400.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), 2.0, 10.0, 400.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(opd_field(at(5.0, 0.0, 0.0), std::nextafter(quarter_turn, 0.0), 10.0, 400.0, 0.0).status,
              Status::Ok);
}

TEST(IrbemOpd, AnOverflowingExtrapolationIsADomainErrorNotANaN) {
    // The quiet field's exp(x / dx) overflows far enough sunward, and under the compression it
    // overflows at x / s — sooner at high pressure. The failure is a NaN in B_y (an infinite
    // envelope times an exactly-zero y), which must be refused at the boundary rather than found
    // a hundred RK4 steps downstream.
    for (double x : {1.0e5, 1.0e6, 1.0e30}) {
        for (double n : {5.0, 50.0}) {
            const ir::Result<ir::FieldVector<Frame::GSM>> r = opd_field(at(x, 0.0, 0.0), 0.2, n, 400.0, 0.0);
            EXPECT_EQ(r.status, Status::DomainError) << "x = " << x << " n = " << n;
            EXPECT_EQ(r.value.v[0], 0.0);
            EXPECT_EQ(r.value.v[1], 0.0);
            EXPECT_EQ(r.value.v[2], 0.0);
        }
    }
    // The three components do not overflow together: B_z alone at (0, 0, 1e120), where the quiet
    // field's cubic z^3 overflows while B_x's quadratic does not; B_y and B_z together at
    // (0, 1e60, 1e150). A guard that looked only at B_x would pass both.
    for (const Position<Frame::GSM>& p : {at(0.0, 0.0, 1.0e120), at(0.0, 1.0e60, 1.0e150)}) {
        const std::array<double, 3> raw = opd_components<double>(opd_parameters<double>(10.0, 400.0, -20.0), 0.3, 0.954, p.v[0], p.v[1], p.v[2]);
        EXPECT_TRUE(std::isfinite(raw[0])) << "B_x must be the component that still fits";
        EXPECT_FALSE(std::isfinite(raw[0]) && std::isfinite(raw[1]) && std::isfinite(raw[2]));
        const ir::Result<ir::FieldVector<Frame::GSM>> r = opd_field(p, 0.3, 10.0, 400.0, -20.0);
        EXPECT_EQ(r.status, Status::DomainError);
        EXPECT_EQ(r.value.v[2], 0.0);
    }
    // Far out but still representable stays what it was: an extrapolation, reported and returned.
    const ir::Result<ir::FieldVector<Frame::GSM>> big = opd_field(at(1.0e3, 0.0, 0.0), 0.2, 10.0, 400.0, 0.0);
    EXPECT_EQ(big.status, Status::OutOfValidityRange);
    EXPECT_TRUE(std::isfinite(big.value.v[2]));
    EXPECT_NE(big.value.v[2], 0.0);
    EXPECT_EQ(opd_field(at(-1.0e6, 0.0, 0.0), 0.2, 10.0, 400.0, 0.0).status, Status::OutOfValidityRange);
}

TEST(IrbemOpd, ContextOverloadAgreesWithTheExplicitOne) {
    const ir::Epoch epoch{2015.5, 43200.0, 2015, 180};
    ir::RotationTable identity{};
    for (fx::mat3d& m : identity) m = fx::mat3d::identity();
    ir::DriverSet drivers{};
    drivers[static_cast<std::size_t>(ir::Driver::Dst)] = -47.0;
    drivers[static_cast<std::size_t>(ir::Driver::Dsw)] = 13.0;
    drivers[static_cast<std::size_t>(ir::Driver::Vsw)] = 437.0;
    const ir::ContextResult built = ir::make_field_context(epoch, -0.42, identity, drivers);
    ASSERT_TRUE(built.has_value()) << ir::describe(built.error());

    const Position<Frame::GSM> p = at(3.75, -1.5, 2.25);
    const ir::Result<ir::FieldVector<Frame::GSM>> via_context = opd_field(p, built.value());
    const ir::Result<ir::FieldVector<Frame::GSM>> via_scalars = opd_field(p, -0.42, 13.0, 437.0, -47.0);
    EXPECT_EQ(via_context.status, via_scalars.status);
    EXPECT_EQ(via_context.value.v[0], via_scalars.value.v[0]);
    EXPECT_EQ(via_context.value.v[1], via_scalars.value.v[1]);
    EXPECT_EQ(via_context.value.v[2], via_scalars.value.v[2]);
    // Kp, Pdyn, By and Bz are NOT read: this model has three drivers and sees no others.
    drivers[static_cast<std::size_t>(ir::Driver::Kp)] = 70.0;
    drivers[static_cast<std::size_t>(ir::Driver::Pdyn)] = 30.0;
    drivers[static_cast<std::size_t>(ir::Driver::BzIMF)] = -40.0;
    const ir::ContextResult stormy = ir::make_field_context(epoch, -0.42, identity, drivers);
    ASSERT_TRUE(stormy.has_value());
    EXPECT_EQ(opd_field(p, stormy.value()).value.v[2], via_context.value.v[2]);
}

// ================================================================================================
// The batch lanes
// ================================================================================================

TEST(IrbemOpd, HostFloatLaneTracksTheReferenceLane) {
    const std::vector<Position<Frame::GSM>> pts = scatter(4096);
    std::vector<float> pos(3 * pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    std::vector<float> out(3 * pts.size());
    const double ps = 0.35;
    const OpdParameters<float> pf = opd_parameters<float>(22.0, 470.0, -85.0);
    ASSERT_TRUE(opd_field_host(pos, out, static_cast<float>(std::sin(ps)), static_cast<float>(std::cos(ps)), pf));

    const OpdParameters<double> pd = opd_parameters<double>(22.0, 470.0, -85.0);
    double worst_abs = 0.0;
    double worst_rel = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const std::array<double, 3> ref = opd_components<double>(
            pd, std::sin(ps), std::cos(ps), pos[(3 * i) + 0], pos[(3 * i) + 1], pos[(3 * i) + 2]);
        double d2 = 0.0;
        double r2 = 0.0;
        for (int c = 0; c < 3; ++c) {
            const std::size_t k = static_cast<std::size_t>(c);
            const double e = out[(3 * i) + k] - ref[k];
            d2 += e * e;
            r2 += ref[k] * ref[k];
        }
        worst_abs = std::max(worst_abs, std::sqrt(d2));
        worst_rel = std::max(worst_rel, std::sqrt(d2) / (std::sqrt(r2) + 1e-9));
    }
    std::printf("[ MEASURED ] fp32 host lane vs fp64 reference: max |dB| %.3e nT, max relative %.3e\n",
                worst_abs, worst_rel);
    EXPECT_LT(worst_abs, 2e-2);
}

TEST(IrbemOpd, HostFloatLaneRejectsMismatchedSpans) {
    const OpdParameters<float> p = opd_parameters<float>(10.0, 400.0, 0.0);
    std::vector<float> pos(7, 0.0F);
    std::vector<float> out(7, 0.0F);
    EXPECT_FALSE(opd_field_host(pos, out, 0.0F, 1.0F, p));
    std::vector<float> pos6(6, 1.0F);
    std::vector<float> out3(3, 0.0F);
    EXPECT_FALSE(opd_field_host(pos6, out3, 0.0F, 1.0F, p));
    std::vector<float> out6(6, 0.0F);
    EXPECT_TRUE(opd_field_host(pos6, out6, 0.0F, 1.0F, p));
}

TEST(IrbemOpd, ParameterBlockIsTheQuietBlockThenTheScaleAndTheRing) {
    const OpdParameters<float> p = opd_parameters<float>(30.0, 450.0, -60.0);
    const std::array<float, opd_param_count> block = opd_param_block(0.25F, 0.75F, p);
    EXPECT_EQ(block.size(), 32U);
    const std::array<float, ir::t89_param_count> quiet = ir::t89_param_block(0.25F, 0.75F, opd_quiet_bin);
    for (std::size_t k = 0; k < ir::t89_param_count; ++k) EXPECT_EQ(block[k], quiet[k]) << "slot " << k;
    EXPECT_EQ(block[30], p.scale);
    EXPECT_EQ(block[31], p.ring);
    // The registry's parameter count is this header's, asserted where both are visible.
#if CHEATAH_SPACE_IRBEM_OPD_GPU
    EXPECT_EQ(ir::gpu::kernel_info("irbem_opd_f32").params, opd_param_count);
    EXPECT_EQ(ir::gpu::kernel_info("irbem_opd_f32").bindings, 4U);
#endif
}

TEST(IrbemOpd, BatchReportsTheSameEnvelopeTheScalarLaneDoes) {
    const std::vector<Position<Frame::GSM>> good{at(5.0, 1.0, 1.0), at(-8.0, 2.0, -1.0)};
    std::vector<ir::FieldVector<Frame::GSM>> out(good.size());
    EXPECT_EQ(opd_field_batch(good, 0.2, 10.0, 400.0, -20.0, out).status, Status::Ok);

    // One point beyond 60 Re: the whole batch is out of validity, every point computed.
    const std::vector<Position<Frame::GSM>> far{at(5.0, 1.0, 1.0), at(-70.0, 0.0, 0.0)};
    EXPECT_EQ(opd_field_batch(far, 0.2, 10.0, 400.0, -20.0, out).status, Status::OutOfValidityRange);
    EXPECT_EQ(opd_field(far[1], 0.2, 10.0, 400.0, -20.0).status, Status::OutOfValidityRange);
    for (const ir::FieldVector<Frame::GSM>& b : out) EXPECT_NE(b.v[2], 0.0);

    // One point inside the Earth, and one not finite: a domain error, every output zeroed.
    for (const Position<Frame::GSM>& bad :
         {at(0.5, 0.0, 0.0), at(std::nan(""), 0.0, 0.0), at(std::numeric_limits<double>::infinity(), 0.0, 0.0)}) {
        const std::vector<Position<Frame::GSM>> mixed{at(5.0, 1.0, 1.0), bad};
        std::vector<ir::FieldVector<Frame::GSM>> mixed_out(mixed.size(), ir::FieldVector<Frame::GSM>{fx::vec3d{1.0, 1.0, 1.0}});
        const ir::Result<bool> r = opd_field_batch(mixed, 0.2, 10.0, 400.0, -20.0, mixed_out);
        EXPECT_EQ(r.status, Status::DomainError);
        EXPECT_EQ(opd_field(bad, 0.2, 10.0, 400.0, -20.0).status, Status::DomainError);
        for (const ir::FieldVector<Frame::GSM>& b : mixed_out) {
            EXPECT_EQ(b.v[0], 0.0);
            EXPECT_EQ(b.v[1], 0.0);
            EXPECT_EQ(b.v[2], 0.0);
        }
    }
    // Two bad points, the first one first: the fold stays poisoned.
    const std::vector<Position<Frame::GSM>> two_bad{at(std::nan(""), 0.0, 0.0),
                                                    at(std::numeric_limits<double>::infinity(), 1.0, 0.0),
                                                    at(5.0, 1.0, 1.0)};
    std::vector<ir::FieldVector<Frame::GSM>> two_out(two_bad.size());
    EXPECT_EQ(opd_field_batch(two_bad, 0.2, 10.0, 400.0, -20.0, two_out).status, Status::DomainError);

    // An out-of-range driver is still reported on an EMPTY batch.
    const std::vector<Position<Frame::GSM>> none;
    std::vector<ir::FieldVector<Frame::GSM>> none_out;
    EXPECT_EQ(opd_field_batch(none, 0.2, 60.0, 400.0, -20.0, none_out).status, Status::OutOfValidityRange);
    EXPECT_EQ(opd_field_batch(none, 0.2, 10.0, 400.0, -20.0, none_out).status, Status::Ok);
}

TEST(IrbemOpd, BatchAgreesWithTheReferenceLane) {
    const std::vector<Position<Frame::GSM>> pts = scatter(1000);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    const ir::Result<bool> r = opd_field_batch(pts, 0.28, 18.0, 430.0, -55.0, out);
    EXPECT_EQ(r.status, Status::Ok);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const ir::FieldVector<Frame::GSM> ref = opd_field_at(pts[i], std::sin(0.28), std::cos(0.28), 18.0, 430.0, -55.0);
        if (r.value) {
            for (int c = 0; c < 3; ++c) EXPECT_NEAR(out[i].v[c], ref.v[c], 2e-2) << "device lane, point " << i;
        } else {
            EXPECT_EQ(out[i].v[0], ref.v[0]);
            EXPECT_EQ(out[i].v[1], ref.v[1]);
            EXPECT_EQ(out[i].v[2], ref.v[2]);
        }
    }
}

TEST(IrbemOpd, BatchRejectsMismatchedSpans) {
    const std::vector<Position<Frame::GSM>> pts = scatter(4);
    std::vector<ir::FieldVector<Frame::GSM>> shorter(3);
    EXPECT_EQ(opd_field_batch(pts, 0.1, 10.0, 400.0, 0.0, shorter).status, Status::DomainError);
    std::vector<ir::FieldVector<Frame::GSM>> right(4);
    EXPECT_EQ(opd_field_batch(pts, std::nan(""), 10.0, 400.0, 0.0, right).status, Status::DomainError);
    EXPECT_EQ(opd_field_batch(pts, 0.1, std::nan(""), 400.0, 0.0, right).status, Status::DomainError);
    EXPECT_EQ(opd_field_batch(pts, 0.1, 10.0, std::nan(""), 0.0, right).status, Status::DomainError);
    EXPECT_EQ(opd_field_batch(pts, 0.1, 10.0, 400.0, std::nan(""), right).status, Status::DomainError);
    EXPECT_EQ(opd_field_batch(pts, std::numbers::pi / 2.0, 10.0, 400.0, 0.0, right).status, Status::DomainError);
    EXPECT_EQ(opd_field_batch(pts, 0.1, 0.0, 400.0, 0.0, right).status, Status::DomainError);
    EXPECT_EQ(opd_field_batch(pts, 0.1, 10.0, -1.0, 0.0, right).status, Status::DomainError);
    const ir::Result<bool> empty = opd_field_batch({}, 0.1, 10.0, 400.0, 0.0, {});
    EXPECT_EQ(empty.status, Status::Ok);
    EXPECT_FALSE(empty.value);
    EXPECT_EQ(opd_field_batch(pts, 0.1, 10.0, 400.0, -300.0, right).status, Status::OutOfValidityRange);
}

TEST(IrbemOpd, NothingOnTheHeapInTheHotPath) {
    const std::vector<Position<Frame::GSM>> pts = scatter(256);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    (void)opd_field_batch(pts, 0.2, 10.0, 400.0, -30.0, out);
    (void)opd_field(at(5.0, 1.0, 1.0), 0.2, 10.0, 400.0, -30.0);

    const std::size_t before = cheatah_space_test::allocation_count();
    for (int i = 0; i < 64; ++i) {
        sink = sink + opd_field(at(5.0 + (0.01 * i), 1.0, 1.0), 0.2, 10.0 + i, 400.0, -30.0).value.v[2];
    }
    (void)opd_field_batch(pts, 0.2, 10.0, 400.0, -30.0, out);   // below the crossover: the host lane
    sink = sink + out[0].v[2];
    EXPECT_EQ(before, cheatah_space_test::allocation_count());
}

// ================================================================================================
// The total field
// ================================================================================================

TEST(IrbemOpd, TotalFieldSuperposesInternalAndExternal) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.0).value();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::TotalFieldOpd<10> total(igrf, rot, 12.0, 450.0, -40.0);
    static_assert(ir::GeoFieldModel<ir::TotalFieldOpd<10>>, "a tracer must be able to follow it");
    static_assert(ir::TotalFieldOpd<10>::degree == 10);

    const Position<Frame::GEO> p{fx::vec3d{6.0, 0.0, 0.0}};
    const double internal_only = igrf.evaluate(p).magnitude();
    const double with_external = total.evaluate(p).magnitude();
    EXPECT_EQ(igrf.g(1, 0), total.g(1, 0));
    EXPECT_EQ(igrf.h(2, 1), total.h(2, 1));
    EXPECT_EQ(&total.internal(), &igrf);
    EXPECT_EQ(&total.rotations(), &rot);
    EXPECT_GT(std::abs(with_external - internal_only) / internal_only, 1e-3);
    EXPECT_LT(std::abs(with_external - internal_only) / internal_only, 0.5);

    // The superposition is the rotated external field added to the internal one, exactly.
    const Position<Frame::GSM> p_gsm = ir::transform<Frame::GSM>(p, rot);
    const double tilt = rot.dipole_tilt_deg * (std::numbers::pi / 180.0);
    const ir::Result<ir::FieldVector<Frame::GSM>> ext = opd_field(p_gsm, tilt, 12.0, 450.0, -40.0);
    ASSERT_EQ(ext.status, Status::Ok);
    const fx::vec3d expect = igrf.evaluate(p).v + ir::transform<Frame::GEO>(ext.value, rot).v;
    const fx::vec3d got = total.evaluate(p).v;
    for (int c = 0; c < 3; ++c) EXPECT_EQ(got[c], expect[c]);

    // Activity propagates into the invariants: a deeper Dst and a stronger wind move B_min.
    const ir::TotalFieldOpd<10> quiet(igrf, rot, 5.0, 380.0, -8.0);
    const ir::TotalFieldOpd<10> storm(igrf, rot, 20.0, 500.0, -100.0);
    const ir::Result<ir::FieldLine> lq = ir::trace_invariant(quiet, p, 60.0);
    const ir::Result<ir::FieldLine> ls = ir::trace_invariant(storm, p, 60.0);
    ASSERT_EQ(lq.status, Status::Ok);
    ASSERT_EQ(ls.status, Status::Ok);
    std::printf("[ MEASURED ] B_min at GEO (6, 0, 0): quiet %.3f nT, storm %.3f nT\n", lq.value.b_min, ls.value.b_min);
    EXPECT_NE(lq.value.b_min, ls.value.b_min);
    EXPECT_GT(std::fabs(lq.value.b_min - ls.value.b_min), 1.0);
}

TEST(IrbemOpd, TotalFieldTracesAndReportsWhenTheExternalModelDeclines) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.0).value();
    const ir::Rotations rot = epoch_rotations(igrf);
    // Drivers outside the envelope: the external field is still summed in, and says so.
    const ir::TotalFieldOpd<10> extrapolated(igrf, rot, 60.0, 450.0, -200.0);
    const Position<Frame::GEO> p{fx::vec3d{5.0, 1.0, 0.5}};
    EXPECT_EQ(extrapolated.external_status(p), Status::OutOfValidityRange);
    EXPECT_NE(extrapolated.evaluate(p).magnitude(), igrf.evaluate(p).magnitude());
    // A point the external model has no answer for (inside the Earth): the internal field alone.
    const Position<Frame::GEO> deep{fx::vec3d{0.5, 0.1, 0.1}};
    EXPECT_EQ(extrapolated.external_status(deep), Status::DomainError);
    EXPECT_EQ(extrapolated.evaluate(deep).v[0], igrf.evaluate(deep).v[0]);
    EXPECT_EQ(extrapolated.evaluate(deep).v[2], igrf.evaluate(deep).v[2]);
    // And with the tilt pushed to a right angle the whole field declines everywhere: the trace
    // then follows IGRF alone and still closes.
    ir::Rotations flat = rot;
    flat.dipole_tilt_deg = 90.0;
    const ir::TotalFieldOpd<10> declined(igrf, flat, 10.0, 400.0, 0.0);
    EXPECT_EQ(declined.external_status(p), Status::DomainError);
    EXPECT_EQ(declined.evaluate(p).v[1], igrf.evaluate(p).v[1]);
    // In range, a normal trace through the total field closes and reports Ok.
    const ir::TotalFieldOpd<10> nominal(igrf, rot, 10.0, 400.0, -30.0);
    EXPECT_EQ(nominal.external_status(p), Status::Ok);
    const ir::Result<ir::FieldLine> line = ir::trace_invariant(nominal, p, 45.0);
    EXPECT_EQ(line.status, Status::Ok);
    EXPECT_GT(line.value.b_min, 0.0);
}

// ================================================================================================
// The differential against the IRBEM oracle
// ================================================================================================

TEST(IrbemOpd, DiffersFromTheIrbemOracleByTheMeasuredEnvelope) {
    const Oracle& o = oracle();
    if (!o.usable()) {
        GTEST_SKIP() << "IRBEM oracle not present (set CHEATAH_SPACE_IRBEM_ORACLE to its .so); "
                        "the oracle is a dev-only black box and is never linked";
    }
    // The radiation-belt region of tools/oracle/opd_diff.cpp section 3, exactly: three epochs
    // spanning the tilt range, the four regimes clipped to the envelope, the drivers swept
    // CONTINUOUSLY around each (nine triples), 60 scattered points at 3..10 Re per triple, inside
    // the Shue magnetopause, external field isolated as kext=6 minus kext=0.
    //
    // The caps are a MEASUREMENT, not an agreement target: the model's form and coefficients were
    // never published, this implements the documented structure with published constants, and the
    // harness measures the gap at RMS 15.2 / 26.4 / 49.7 / 67.7 nT (52 / 47 / 47 / 51 % of the
    // oracle's external field) across the four regimes. Caps ~1.4x above, so a regression in THIS
    // implementation fails while the documented model-family gap does not.
    struct Epoch {
        int doy;
        double ut;
    };
    const std::array<Epoch, 3> epochs{{{80, 39183.0}, {180, 43200.0}, {355, 7200.0}}};
    const std::array<Drivers, 4> regimes{{{5.0, 380.0, -8.0}, {8.0, 450.0, -42.0}, {20.0, 500.0, -100.0}, {45.0, 500.0, -100.0}}};
    const std::array<double, 4> rms_cap{{21.5, 37.0, 70.0, 95.0}};

    for (std::size_t ri = 0; ri < regimes.size(); ++ri) {
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
            for (int k = 0; k < 9; ++k) {
                const double dsw = std::clamp(regimes[ri].n * (0.8 + (0.05 * k)), 5.0, 50.0);
                const double vsw = std::clamp(regimes[ri].v * (0.92 + (0.02 * k)), 300.0, 500.0);
                const double dst = std::clamp(regimes[ri].dst + (-15.0 + (3.75 * k)), -100.0, 20.0);
                std::uint64_t s = 0xA5A5A5A5DEADBEEFULL + static_cast<std::uint64_t>(k);
                const auto next = [&s] {
                    s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
                    return static_cast<double>(s >> 11) / 9007199254740992.0;
                };
                for (int i = 0; i < 60; ++i) {
                    const double r = 3.0 + (7.0 * next());
                    const double th = std::acos(1.0 - (2.0 * next()));
                    const double ph = 6.283185307179586 * next();
                    const double x = r * std::sin(th) * std::cos(ph);
                    const double y = r * std::sin(th) * std::sin(ph);
                    const double z = r * std::cos(th);
                    if (!inside_magnetopause(x, y, z, dsw, vsw)) continue;
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
                    int k6 = 6;
                    std::vector<double> mag(25, 0.0);
                    mag[1] = dst;
                    mag[2] = dsw;
                    mag[3] = vsw;
                    std::array<double, 3> b0{};
                    std::array<double, 3> b6{};
                    double m0 = 0.0;
                    double m6 = 0.0;
                    double x1 = geo[0];
                    double x2 = geo[1];
                    double x3 = geo[2];
                    o.get_field(&k0, options.data(), &sysaxes, &iyear, &idoy, &ut, &x1, &x2, &x3, mag.data(), b0.data(), &m0);
                    o.get_field(&k6, options.data(), &sysaxes, &iyear, &idoy, &ut, &x1, &x2, &x3, mag.data(), b6.data(), &m6);
                    if (b6[0] < -1e30) continue;   // the oracle's refusal; never inside the envelope
                    std::array<double, 3> dgeo{b6[0] - b0[0], b6[1] - b0[1], b6[2] - b0[2]};
                    std::array<double, 3> ora{};
                    {
                        int si = 1;
                        int so = 2;
                        o.coord_trans(&one, &si, &so, &iyear, &idoy, &ut, dgeo.data(), ora.data());
                    }
                    const ir::FieldVector<Frame::GSM> mine = opd_field_at(at(x, y, z), std::sin(ps), std::cos(ps), dsw, vsw, dst);
                    for (int c = 0; c < 3; ++c) {
                        const std::size_t kc = static_cast<std::size_t>(c);
                        const double d = mine.v[c] - ora[kc];
                        sum2 += d * d;
                        sig2 += ora[kc] * ora[kc];
                    }
                    ++n;
                }
            }
        }
        ASSERT_GT(n, 1000U);
        const double rms = std::sqrt(sum2 / static_cast<double>(n));
        std::printf("[ MEASURED ] regime %zu: RMS |dB| vs IRBEM kext=6 = %7.3f nT (%.1f%% of the external "
                    "field, %zu points)\n",
                    ri, rms, 100.0 * std::sqrt(sum2 / sig2), n);
        EXPECT_LT(rms, rms_cap[ri]) << "regime " << ri << ": the gap against the oracle has GROWN beyond "
                                       "the documented structure-vs-distributed-code difference";
        // Same sign and order of magnitude everywhere: the weak claim no re-parameterization breaks.
        EXPECT_LT(std::sqrt(sum2 / sig2), 0.75) << "regime " << ri;
    }
}

#if CHEATAH_SPACE_IRBEM_OPD_GPU
TEST(IrbemOpd, BatchUsesTheDeviceWhenOneIsAvailable) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    if (!std::filesystem::exists(ir::gpu::shader_path("irbem_opd_f32"))) {
        GTEST_SKIP() << "irbem_opd_f32.spv not built under " << ir::gpu::shader_path("irbem_opd_f32");
    }
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_opd_f32");
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    const ir::Result<bool> r = opd_field_batch(pts, 0.28, 18.0, 430.0, -55.0, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_TRUE(r.value) << "the batch fell back to the host with a device present";

    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = opd_field_at(pts[i], std::sin(0.28), std::cos(0.28), 18.0, 430.0, -55.0);
        for (int c = 0; c < 3; ++c) worst = std::max(worst, std::fabs(out[i].v[c] - ref.v[c]));
    }
    std::printf("[ MEASURED ] device batch of %zu vs fp64 reference: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 2e-2);

    // Below the crossover the same call runs on the host, bit-identical to the reference.
    const std::vector<Position<Frame::GSM>> few = scatter(16);
    std::vector<ir::FieldVector<Frame::GSM>> few_out(few.size());
    const ir::Result<bool> host = opd_field_batch(few, 0.28, 18.0, 430.0, -55.0, few_out);
    EXPECT_FALSE(host.value);
    for (std::size_t i = 0; i < few.size(); ++i) {
        const ir::FieldVector<Frame::GSM> ref = opd_field_at(few[i], std::sin(0.28), std::cos(0.28), 18.0, 430.0, -55.0);
        EXPECT_EQ(few_out[i].v[2], ref.v[2]);
    }
}

TEST(IrbemOpd, TheDeviceLaneRefusesABadPointBeforeItDispatches) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    if (!std::filesystem::exists(ir::gpu::shader_path("irbem_opd_f32"))) GTEST_SKIP() << "no irbem_opd_f32.spv";
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_opd_f32");
    std::vector<Position<Frame::GSM>> pts = scatter(n);
    pts[n / 3] = at(0.25, 0.0, 0.0);
    std::vector<ir::FieldVector<Frame::GSM>> out(n, ir::FieldVector<Frame::GSM>{fx::vec3d{7.0, 7.0, 7.0}});
    const ir::Result<bool> r = opd_field_batch(pts, 0.28, 18.0, 430.0, -55.0, out);
    EXPECT_EQ(r.status, Status::DomainError);
    EXPECT_FALSE(r.value) << "a refused batch must not have reached the device";
    for (const ir::FieldVector<Frame::GSM>& b : out) {
        ASSERT_EQ(b.v[0], 0.0);
        ASSERT_EQ(b.v[1], 0.0);
        ASSERT_EQ(b.v[2], 0.0);
    }
    std::vector<Position<Frame::GSM>> far = scatter(n);
    far[n / 2] = at(-70.0, 0.0, 0.0);
    const ir::Result<bool> f = opd_field_batch(far, 0.28, 18.0, 430.0, -55.0, out);
    EXPECT_EQ(f.status, Status::OutOfValidityRange);
    EXPECT_TRUE(f.value) << "an out-of-validity batch must still be computed on the device";
}

TEST(IrbemOpd, DeviceKernelAgreesWithTheHostLane) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    if (!std::filesystem::exists(ir::gpu::shader_path("irbem_opd_f32"))) GTEST_SKIP() << "no irbem_opd_f32.spv";
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
    const OpdParameters<float> p = opd_parameters<float>(20.0, 450.0, -60.0);
    std::vector<float> host(3 * n);
    std::vector<float> device(3 * n);
    ASSERT_TRUE(opd_field_host(pos, host, sp, cp, p));
    const std::array<float, opd_param_count> block = opd_param_block(sp, cp, p);
    ir::gpu::dispatch_batch("irbem_opd_f32", pos, device, std::span<const float>(block));
    double worst = 0.0;
    for (std::size_t i = 0; i < 3 * n; ++i) worst = std::max(worst, std::fabs(static_cast<double>(device[i]) - host[i]));
    std::printf("[ MEASURED ] device vs host, %zu points: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-3) << "the device is not evaluating the same expressions";
}
#endif
