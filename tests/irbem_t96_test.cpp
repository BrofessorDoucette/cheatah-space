/// @file irbem_t96_test.cpp
/// @brief The suite for `space/irbem/ext_t96.hpp` — Tsyganenko (1995/1996), "T96".
///
/// A fitted model cannot be checked against an analytic truth, and THIS fitted model cannot even be
/// checked against its own published coefficients, because there are none (see the header's brief).
/// What can be answered, and is answered here, is every question about whether the implementation is
/// what it claims to be:
///
///  - **`div B = 0`, with no floor.** Every basis field is exactly divergence-free — the discs are
///    curls of potentials, the toroidal modules are `∇T × r`, the box harmonics are gradients of
///    harmonic functions — so a central-difference divergence of the FULL field must fall as `h^2`
///    over three decades and does. Unlike T89, whose eq. (20) table rounding sets a floor, nothing
///    here is rounded: every coefficient is used exactly as fitted. A sign error in any one of the
///    ~thirty analytic derivatives in the basis makes the residual stop scaling.
///  - **The driver structure the oracle was measured to have.** Exactly affine in Dst; exactly
///    homogeneous of degree one in the IMF; continuous everywhere, kinked only on the northward ray;
///    smooth in pressure. Those are theorems about eq. (1) of the header, asserted as such.
///  - **The exact symmetries** the disc and box blocks carry, bitwise.
///  - **The physics the model exists for.** Dayside compression, nightside stretching, a ring-current
///    depression that deepens linearly with Dst, a tail that loads with southward Bz.
///  - **The envelope from both sides**, driver by driver and bound by bound.
///  - **The differential against the IRBEM oracle** across the four corpus regimes and a dense
///    southward-Bz sweep, when the oracle is present. Those caps are MEASUREMENTS of the documented
///    published-form gap, ~1.3x above the measured RMS, so a regression fails while the known gap
///    does not.
///
/// The oracle tests `dlopen` IRBEM at runtime rather than linking it: this binary must build and
/// pass on a machine that has never heard of IRBEM.

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
#include "irbem_domain_corpus.hpp"
#include "space/irbem/api.hpp"
#include "space/irbem/ext_t96.hpp"
#include "space/irbem/lstar.hpp"

namespace {

namespace ir = cheatah::space::irbem;
namespace corpus = cheatah_space_test;

using ir::Frame;
using ir::Position;
using ir::Status;
using ir::T96Amplitudes;
using ir::T96Basis;
using ir::t96_amplitudes;
using ir::t96_basis;
using ir::t96_basis_count;
using ir::t96_clock_driver;
using ir::t96_components;
using ir::t96_field;
using ir::t96_field_at;
using ir::t96_field_batch;
using ir::t96_field_host;
using ir::t96_geometry;
using ir::t96_param_block;
using ir::t96_param_count;
using ir::t96_tables;

/// A sink the optimizer cannot see through, so the allocation test's calls actually happen.
volatile double sink = 0.0;

/// A GSM point, spelled so a test reads as coordinates rather than as a constructor call.
Position<Frame::GSM> at(double x, double y, double z) {
    return Position<Frame::GSM>{cheatah::fixarray::vec3d(x, y, z)};
}

/// The tilts the suite sweeps: zero, and both signs of a realistic seasonal-diurnal excursion.
constexpr std::array<double, 5> kTilts{-0.55, -0.21, 0.0, 0.21, 0.55};

/// A representative disturbed driver set — the moderate corpus regime — used wherever a test needs
/// "some real storm" rather than a specific one.
struct Drivers {
    double dst;
    double pdyn;
    double by;
    double bz;
};
constexpr Drivers kModerate{-42.0, 3.2, -4.0, -5.0};

/// `div B` at one point by central differences of the FULL field for a given amplitude set.
double divergence(const T96Amplitudes<double>& a, double x, double y, double z, double h) {
    const auto b = [&](double p, double q, double r) { return t96_components<double>(a, p, q, r); };
    return ((b(x + h, y, z)[0] - b(x - h, y, z)[0]) + (b(x, y + h, z)[1] - b(x, y - h, z)[1]) +
            (b(x, y, z + h)[2] - b(x, y, z - h)[2])) /
           (2.0 * h);
}

/// The largest `|div B|` over the sampled box at difference step @p h, over tilts and drivers.
double worst_divergence(double h) {
    double worst = 0.0;
    for (double ps : kTilts) {
        for (const Drivers& d : {Drivers{0.0, 2.0, 0.0, 0.0}, kModerate, Drivers{-100.0, 9.0, 8.0, -10.0}}) {
            const T96Amplitudes<double> a = t96_amplitudes<double>(ps, d.dst, d.pdyn, d.by, d.bz);
            // Integer induction, coordinates derived per iteration (cert-flp30).
            for (int ix = 0; ix <= 8; ++ix)
                for (int iy = 0; iy <= 4; ++iy)
                    for (int iz = 0; iz <= 4; ++iz) {
                        const double x = -25.0 + (4.1 * ix);
                        const double y = -9.0 + (4.5 * iy);
                        const double z = -7.0 + (3.5 * iz);
                        if (std::sqrt((x * x) + (y * y) + (z * z)) < 1.5) continue;
                        worst = std::max(worst, std::fabs(divergence(a, x, y, z, h)));
                    }
        }
    }
    return worst;
}

/// A deterministic scatter of GSM points over the inner magnetosphere and near tail, never a
/// lattice, so no component is systematically zero.
std::vector<Position<Frame::GSM>> scatter(std::size_t n) {
    std::vector<Position<Frame::GSM>> out;
    out.reserve(n);
    std::uint64_t s = 0x9E3779B97F4A7C15ULL;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) / 9007199254740992.0;
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double r = 2.5 + (12.5 * next());
        const double th = std::acos(1.0 - (2.0 * next()));
        const double ph = 6.283185307179586 * next();
        out.push_back(at(r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph),
                         r * std::cos(th)));
    }
    return out;
}

// ---- the oracle, opened at runtime and never linked --------------------------------------------

using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

struct Oracle {
    void* handle = nullptr;
    GetField1 get_field = nullptr;
    CoordTransVec1 coord_trans = nullptr;
    [[nodiscard]] bool usable() const { return get_field != nullptr && coord_trans != nullptr; }

    /// The external field at a GSM point in GSM nT, kext = 7 minus kext = 0; NaN when refused.
    [[nodiscard]] std::array<double, 3> ext(int doy, double ut, double x, double y, double z,
                                            double dst, double pdyn, double by, double bz) const {
        int one = 1;
        int iyear = 2015;
        int idoy = doy;
        double t = ut;
        std::array<double, 3> gsm{x, y, z};
        std::array<double, 3> geo{};
        {
            int si = 2;
            int so = 1;
            coord_trans(&one, &si, &so, &iyear, &idoy, &t, gsm.data(), geo.data());
        }
        std::array<int, 5> options{0, 0, 0, 0, 0};
        int sysaxes = 1;
        int k0 = 0;
        int k7 = 7;
        std::vector<double> mag(25, 0.0);
        mag[1] = dst;
        mag[4] = pdyn;
        mag[5] = by;
        mag[6] = bz;
        std::array<double, 3> b0{};
        std::array<double, 3> b7{};
        double m0 = 0.0;
        double m7 = 0.0;
        get_field(&k0, options.data(), &sysaxes, &iyear, &idoy, &t, &geo[0], &geo[1], &geo[2],
                  mag.data(), b0.data(), &m0);
        get_field(&k7, options.data(), &sysaxes, &iyear, &idoy, &t, &geo[0], &geo[1], &geo[2],
                  mag.data(), b7.data(), &m7);
        if (b7[0] < -1e30) return {std::nan(""), std::nan(""), std::nan("")};
        std::array<double, 3> dgeo{b7[0] - b0[0], b7[1] - b0[1], b7[2] - b0[2]};
        std::array<double, 3> out{};
        int si = 1;
        int so = 2;
        coord_trans(&one, &si, &so, &iyear, &idoy, &t, dgeo.data(), out.data());
        return out;
    }

    /// The oracle's own dipole tilt at the epoch.
    [[nodiscard]] double tilt(int doy, double ut) const {
        int one = 1;
        int si = 4;
        int so = 2;
        int iyear = 2015;
        int idoy = doy;
        double t = ut;
        std::array<double, 3> in{0.0, 0.0, 1.0};
        std::array<double, 3> out{};
        coord_trans(&one, &si, &so, &iyear, &idoy, &t, in.data(), out.data());
        return std::atan2(out[0], out[2]);
    }
};

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

/// The differential grid: the same magnetopause-shaped scatter the harness reports on, but its
/// own seed, so this test is not the fit's own training data.
std::vector<Position<Frame::GSM>> differential_points() {
    std::vector<Position<Frame::GSM>> out;
    std::uint64_t s = 0x7C4A7F9B3E1D5A2FULL;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) / 9007199254740992.0;
    };
    while (out.size() < 200) {
        const double r = 1.5 * std::pow(30.0 / 1.5, next());
        const double c = 1.0 - (2.0 * next());
        const double snt = std::sqrt(1.0 - (c * c));
        const double ph = 6.283185307179586 * next();
        const double x = r * snt * std::cos(ph);
        const double y = r * snt * std::sin(ph);
        const double z = r * c;
        if (x < -30.0 || x > 9.0) continue;
        if ((y * y) + (z * z) > 180.0 - (15.0 * x)) continue;
        out.push_back(at(x, y, z));
    }
    return out;
}

struct Epoch {
    int doy;
    double ut;
};
constexpr std::array<Epoch, 3> kEpochs{{{172, 61200.0}, {80, 39183.0}, {355, 20500.0}}};

ir::Rotations epoch_rotations(const ir::Igrf<10>& m) {
    const auto r = ir::api::rotations_at(2015, 180, 43200.0, m);
    EXPECT_EQ(Status::Ok, r.status);
    return r.value;
}

}  // namespace

// ================================================================================================
// The tables and the geometry
// ================================================================================================

TEST(IrbemT96, GeometryIsTheFrozenFitGeometry) {
    // The header's geometry is the output of the harness's Nelder-Mead pass, frozen. Pinning the
    // ladders' ORDER and sanity is what a transposed or shuffled regeneration would break; the
    // exact values are the harness's and are asserted bitwise against themselves so a hand edit
    // of one entry fails here rather than silently changing the model.
    EXPECT_GT(t96_geometry.d_0, 0.3);
    EXPECT_GT(t96_geometry.d_y, 3.0);
    for (double a : t96_geometry.a_t) EXPECT_GT(a, 2.0);
    for (double a : t96_geometry.a_rc) EXPECT_GT(a, 1.5);
    for (std::size_t i = 1; i < t96_geometry.fac_shell.size(); ++i) {
        EXPECT_LT(t96_geometry.fac_shell[i - 1], t96_geometry.fac_shell[i]) << "shell ladder order";
    }
    for (double w : t96_geometry.fac_width) {
        EXPECT_GT(w, 0.03);
        EXPECT_LT(w, 0.8);
    }
    for (double p : t96_geometry.box_p) EXPECT_GT(p, 3.0);
    for (double r : t96_geometry.box_r) EXPECT_GT(r, 2.0);
    EXPECT_EQ(t96_basis_count, 128U);
    EXPECT_EQ(ir::t96_column_count, 768U);
}

TEST(IrbemT96, TablesAreNotEmpty) {
    // Every family carries a fitted amplitude on every block of the basis; a regeneration that
    // dropped a block (or a placeholder that was never replaced) is caught here.
    for (std::size_t f = 0; f < ir::t96_family_count; ++f) {
        double discs = 0.0;
        double fac = 0.0;
        double box = 0.0;
        for (std::size_t j = 0; j < t96_basis_count; ++j) {
            for (std::size_t k = 0; k < ir::t96_factor_count; ++k) {
                const double v = std::fabs(t96_tables[f][(j * ir::t96_factor_count) + k]);
                if (j < ir::t96_fac_first) discs += v;
                else if (j < ir::t96_box_first) fac += v;
                else box += v;
            }
        }
        EXPECT_GT(discs, 0.0) << "family " << f;
        EXPECT_GT(fac, 0.0) << "family " << f;
        EXPECT_GT(box, 0.0) << "family " << f;
    }
}

TEST(IrbemT96, ClockDriverIsTheMeasuredThirdFunction) {
    // h = B_t sin(theta/2): |Bz| southward, 0 northward, |By|/sqrt 2 for pure By, homogeneous of
    // degree one, even in By. Exactly-representable inputs so the first three are exact.
    EXPECT_EQ(t96_clock_driver(0.0, -8.0), 8.0);
    EXPECT_EQ(t96_clock_driver(0.0, 8.0), 0.0);
    EXPECT_EQ(t96_clock_driver(0.0, 0.0), 0.0);
    EXPECT_NEAR(t96_clock_driver(8.0, 0.0), 8.0 / std::sqrt(2.0), 1e-14);
    EXPECT_EQ(t96_clock_driver(6.0, -2.5), t96_clock_driver(-6.0, -2.5));
    EXPECT_NEAR(t96_clock_driver(3.0, -4.5), 0.5 * t96_clock_driver(6.0, -9.0), 1e-14);
    // The general form: B_t sin(theta/2) with cos theta = Bz / B_t.
    const double bt = 5.0;
    const double th = 2.0;
    EXPECT_NEAR(t96_clock_driver(bt * std::sin(th), bt * std::cos(th)), bt * std::sin(th / 2.0), 1e-12);
}

TEST(IrbemT96, AmplitudesCollapseTheTablesLinearly) {
    // Eq. (1): the amplitudes are B_0 + Dst B_D + Bz B_Z + h B_A + By B_Y, so they must be exactly
    // affine in Dst and exactly homogeneous in the IMF, at any tilt and pressure.
    const double ps = 0.31;
    const double pd = 4.0;
    const T96Amplitudes<double> quiet = t96_amplitudes<double>(ps, 0.0, pd, 0.0, 0.0);
    const T96Amplitudes<double> d100 = t96_amplitudes<double>(ps, -100.0, pd, 0.0, 0.0);
    const T96Amplitudes<double> d25 = t96_amplitudes<double>(ps, -25.0, pd, 0.0, 0.0);
    const T96Amplitudes<double> imf = t96_amplitudes<double>(ps, 0.0, pd, 6.0, -8.0);
    const T96Amplitudes<double> imf_half = t96_amplitudes<double>(ps, 0.0, pd, 3.0, -4.0);
    const T96Amplitudes<double> both = t96_amplitudes<double>(ps, -25.0, pd, 6.0, -8.0);
    for (std::size_t j = 0; j < t96_basis_count; ++j) {
        const double scale = std::fabs(quiet.c[j]) + std::fabs(d100.c[j]) + std::fabs(imf.c[j]) + 1e-300;
        EXPECT_NEAR(d25.c[j], quiet.c[j] + (0.25 * (d100.c[j] - quiet.c[j])), 1e-12 * scale) << j;
        EXPECT_NEAR(imf_half.c[j] - quiet.c[j], 0.5 * (imf.c[j] - quiet.c[j]), 1e-12 * scale) << j;
        EXPECT_NEAR(both.c[j], d25.c[j] + imf.c[j] - quiet.c[j], 1e-12 * scale) << j;
    }
    EXPECT_EQ(quiet.sin_tilt, std::sin(ps));
    EXPECT_EQ(quiet.cos_tilt, std::cos(ps));
}

TEST(IrbemT96, AmplitudesRoundTripThroughFloat) {
    // The float amplitudes are the double ones rounded ONCE — not a float recomputation.
    const T96Amplitudes<double> d = t96_amplitudes<double>(0.4, -60.0, 2.5, 3.0, -7.0);
    const T96Amplitudes<float> f = t96_amplitudes<float>(0.4, -60.0, 2.5, 3.0, -7.0);
    for (std::size_t j = 0; j < t96_basis_count; ++j) EXPECT_EQ(f.c[j], static_cast<float>(d.c[j])) << j;
    EXPECT_EQ(f.sin_tilt, static_cast<float>(d.sin_tilt));
    EXPECT_EQ(f.cos_tilt, static_cast<float>(d.cos_tilt));
}

// ================================================================================================
// The mathematics: div B = 0, and the symmetries
// ================================================================================================

TEST(IrbemT96, DivergenceVanishesEverywhere) {
    // Three decades of h. Every coefficient is used exactly as fitted, so unlike T89 there is no
    // table-rounding floor: the residual must keep falling as h^2 until fp64 roundoff takes over.
    const double d2 = worst_divergence(1e-2);
    const double d3 = worst_divergence(1e-3);
    const double d4 = worst_divergence(1e-4);
    std::printf("[ MEASURED ] worst |div B|: h=1e-2 %.3e  h=1e-3 %.3e  h=1e-4 %.3e nT/Re\n", d2, d3, d4);
    EXPECT_GT(d2 / d3, 50.0) << "divergence is not falling as h^2 — an analytic derivative is wrong";
    EXPECT_GT(d3 / d4, 20.0) << "divergence is not falling as h^2 — an analytic derivative is wrong";
    EXPECT_LT(d4, 1e-5);
}

TEST(IrbemT96, EveryBasisFieldIsDivergenceFreeOnItsOwn) {
    // The full-field test above can be fooled by two errors that cancel in one driver set; this
    // one cannot. Each of the 128 basis fields is divergence-free by itself — discs as curls,
    // toroidal terms as grad T x r, box terms as gradients of harmonic functions — so each one's
    // central-difference divergence, RELATIVE to the size of the field's own gradient, must be
    // second-order small.
    const double h = 1e-3;
    double worst = 0.0;
    std::size_t worst_j = 0;
    for (double ps : {-0.4, 0.0, 0.35}) {
        const double s = std::sin(ps);
        const double c = std::cos(ps);
        for (const Position<Frame::GSM>& p : {at(-6.5, 2.25, 1.5), at(4.0, -3.0, 2.0), at(-15.0, 6.0, -3.0),
                                              at(2.5, 0.75, -1.25), at(-3.0, -4.0, 4.5)}) {
            T96Basis<double> xp{};
            T96Basis<double> xm{};
            T96Basis<double> yp{};
            T96Basis<double> ym{};
            T96Basis<double> zp{};
            T96Basis<double> zm{};
            t96_basis<double>(t96_geometry, s, c, p.v[0] + h, p.v[1], p.v[2], xp);
            t96_basis<double>(t96_geometry, s, c, p.v[0] - h, p.v[1], p.v[2], xm);
            t96_basis<double>(t96_geometry, s, c, p.v[0], p.v[1] + h, p.v[2], yp);
            t96_basis<double>(t96_geometry, s, c, p.v[0], p.v[1] - h, p.v[2], ym);
            t96_basis<double>(t96_geometry, s, c, p.v[0], p.v[1], p.v[2] + h, zp);
            t96_basis<double>(t96_geometry, s, c, p.v[0], p.v[1], p.v[2] - h, zm);
            for (std::size_t j = 0; j < t96_basis_count; ++j) {
                const double dxx = (xp[j][0] - xm[j][0]) / (2.0 * h);
                const double dyy = (yp[j][1] - ym[j][1]) / (2.0 * h);
                const double dzz = (zp[j][2] - zm[j][2]) / (2.0 * h);
                const double div = dxx + dyy + dzz;
                const double grad = std::fabs(dxx) + std::fabs(dyy) + std::fabs(dzz) + 1e-30;
                const double rel = std::fabs(div) / grad;
                if (rel > worst) {
                    worst = rel;
                    worst_j = j;
                }
            }
        }
    }
    std::printf("[ MEASURED ] worst per-basis |div B| / |grad terms| at h = 1e-3: %.3e (basis %zu)\n", worst, worst_j);
    EXPECT_LT(worst, 1e-4) << "basis field " << worst_j << " is not divergence-free";
}

TEST(IrbemT96, DawnDuskSymmetryOfTheDiscAndBoxBlocks) {
    // The disc block depends on y only through y^2 and y^4 (and linearly in B_y), so y -> -y gives
    // (B_x, -B_y, B_z) BITWISE at any tilt. The cos(y/p) box families do the same; the sin(y/p)
    // families are the opposite parity, (-B_x, B_y, -B_z). This is the test that catches a sign
    // slip in dW/dy, dz_s/dy or dD/dy that the divergence test happens not to see.
    for (double ps : kTilts) {
        const double s = std::sin(ps);
        const double c = std::cos(ps);
        for (double x : {-16.0, -4.5, 2.0, 8.75}) {
            for (double y : {0.75, 3.5, 9.25}) {
                for (double z : {-5.5, -1.0, 0.0, 4.25}) {
                    T96Basis<double> dusk{};
                    T96Basis<double> dawn{};
                    t96_basis<double>(t96_geometry, s, c, x, y, z, dusk);
                    t96_basis<double>(t96_geometry, s, c, x, -y, z, dawn);
                    for (std::size_t j = 0; j < ir::t96_disc_count; ++j) {
                        EXPECT_EQ(dusk[j][0], dawn[j][0]) << "disc " << j;
                        EXPECT_EQ(dusk[j][1], -dawn[j][1]) << "disc " << j;
                        EXPECT_EQ(dusk[j][2], dawn[j][2]) << "disc " << j;
                    }
                    for (std::size_t j = ir::t96_box_first; j < t96_basis_count; ++j) {
                        const std::size_t fam = (j - ir::t96_box_first) / 9;
                        const double sign = fam < 2 ? 1.0 : -1.0;  // cos(y/p) families even in y
                        EXPECT_EQ(dusk[j][0], sign * dawn[j][0]) << "box " << j;
                        EXPECT_EQ(dusk[j][1], -sign * dawn[j][1]) << "box " << j;
                        EXPECT_EQ(dusk[j][2], sign * dawn[j][2]) << "box " << j;
                    }
                }
            }
        }
    }
}

TEST(IrbemT96, ZeroTiltDiscsAreMirrorSymmetricAboutTheEquator) {
    // At psi = 0 the warped sheet is flat (tan psi = sin psi = 0 exactly), so the disc block is even
    // about z = 0 in B_z and odd in B_x, B_y — bitwise.
    for (double x : {-18.0, -6.5, -1.5, 3.25, 9.0}) {
        for (double y : {-7.5, -1.25, 0.0, 2.5, 8.0}) {
            for (double z : {0.5, 2.25, 6.0}) {
                T96Basis<double> up{};
                T96Basis<double> dn{};
                t96_basis<double>(t96_geometry, 0.0, 1.0, x, y, z, up);
                t96_basis<double>(t96_geometry, 0.0, 1.0, x, y, -z, dn);
                for (std::size_t j = 0; j < ir::t96_disc_count; ++j) {
                    EXPECT_EQ(up[j][0], -dn[j][0]) << "disc " << j;
                    EXPECT_EQ(up[j][1], -dn[j][1]) << "disc " << j;
                    EXPECT_EQ(up[j][2], dn[j][2]) << "disc " << j;
                }
            }
        }
    }
}

TEST(IrbemT96, ToroidalModulesHaveNoRadialComponentInSm) {
    // B = grad T x r is tangential to spheres about the DIPOLE: in SM, B . r = 0 identically, for
    // every one of the 72 field-aligned-current terms, at any tilt. Evaluated at a tilt, the GSM
    // vector must be rotated to SM before the dot product — which is also a check on the rotation.
    for (double ps : kTilts) {
        const double s = std::sin(ps);
        const double c = std::cos(ps);
        for (const Position<Frame::GSM>& p : {at(-5.0, 2.0, 3.0), at(3.0, -1.5, 2.5), at(1.0, 0.5, 1.5), at(-9.0, 4.0, -6.0)}) {
            T96Basis<double> b{};
            t96_basis<double>(t96_geometry, s, c, p.v[0], p.v[1], p.v[2], b);
            const double xs = (p.v[0] * c) - (p.v[2] * s);
            const double zs = (p.v[0] * s) + (p.v[2] * c);
            for (std::size_t j = ir::t96_fac_first; j < ir::t96_box_first; ++j) {
                const double bxs = (b[j][0] * c) - (b[j][2] * s);
                const double bzs = (b[j][0] * s) + (b[j][2] * c);
                const double radial = (bxs * xs) + (b[j][1] * p.v[1]) + (bzs * zs);
                const double mag = std::sqrt((bxs * bxs) + (b[j][1] * b[j][1]) + (bzs * bzs));
                EXPECT_LT(std::fabs(radial), 1e-12 * (mag * std::sqrt((xs * xs) + (p.v[1] * p.v[1]) + (zs * zs)) + 1e-30)) << "fac " << j;
            }
        }
    }
}

TEST(IrbemT96, OnTheNoonMidnightMeridianTheFieldStaysInThatPlaneWithoutBy) {
    // With By = 0 every remaining driver family is even in y, so at y = 0 exactly B_y is exactly
    // zero — the property a field-line trace in the meridian plane relies on.
    for (double ps : kTilts) {
        for (const Drivers& d : {Drivers{0.0, 2.0, 0.0, 0.0}, Drivers{-80.0, 5.0, 0.0, -9.0}, Drivers{-30.0, 1.0, 0.0, 6.0}}) {
            const T96Amplitudes<double> a = t96_amplitudes<double>(ps, d.dst, d.pdyn, d.by, d.bz);
            for (double x : {-20.0, -6.0, 2.5, 8.0}) {
                for (double z : {-6.0, 0.0, 3.5}) {
                    EXPECT_EQ(t96_components<double>(a, x, 0.0, z)[1], 0.0) << "x " << x << " z " << z;
                }
            }
        }
    }
}

// ================================================================================================
// The driver structure: continuous, affine in Dst, homogeneous in the IMF
// ================================================================================================

TEST(IrbemT96, TheFieldIsExactlyAffineInDst) {
    const Position<Frame::GSM> p = at(-6.5, 1.25, 0.75);
    for (double ps : kTilts) {
        const auto b0 = t96_field(p, ps, 0.0, 3.0, 2.0, -3.0).value.v;
        const auto b1 = t96_field(p, ps, -100.0, 3.0, 2.0, -3.0).value.v;
        // Dst swept continuously across and beyond the envelope; every value is the affine one.
        for (int k = -300; k <= 40; k += 17) {
            const double dst = static_cast<double>(k);
            const auto b = t96_field(p, ps, dst, 3.0, 2.0, -3.0).value.v;
            for (int c = 0; c < 3; ++c) {
                const double expect = b0[c] + ((dst / -100.0) * (b1[c] - b0[c]));
                EXPECT_NEAR(b[c], expect, 1e-9 * (std::fabs(expect) + 1.0)) << "Dst " << dst;
            }
        }
    }
}

TEST(IrbemT96, TheImfResponseIsHomogeneousAndContinuous) {
    // Along every ray from the origin of the (By, Bz) plane the response is linear; across the
    // plane it is continuous, including at the origin and on the northward ray where h is kinked.
    const Position<Frame::GSM> p = at(-6.6, 0.5, 0.25);
    const double ps = 0.3;
    const auto quiet = t96_field(p, ps, -20.0, 2.0, 0.0, 0.0).value.v;
    for (double by : {-8.0, -2.0, 0.0, 3.0, 9.0}) {
        for (double bz : {-10.0, -6.0, -1.0, 0.0, 2.0, 7.0}) {
            if (by == 0.0 && bz == 0.0) continue;
            const auto full = t96_field(p, ps, -20.0, 2.0, by, bz).value.v;
            for (double l : {0.25, 0.5, 0.75}) {
                const auto part = t96_field(p, ps, -20.0, 2.0, l * by, l * bz).value.v;
                for (int c = 0; c < 3; ++c) {
                    const double expect = quiet[c] + (l * (full[c] - quiet[c]));
                    EXPECT_NEAR(part[c], expect, 1e-9 * (std::fabs(expect) + 1.0)) << by << " " << bz << " " << l;
                }
            }
        }
    }
    // Continuity: a step of 1e-3 nT in any driver moves the field by no more than the slope allows.
    const auto near = [&](double by, double bz) { return t96_field(p, ps, -20.0, 2.0, by, bz).value.v; };
    for (const std::array<double, 2>& imf : {std::array<double, 2>{0.0, 0.0}, {0.0, 5.0}, {0.0, -5.0}, {4.0, 0.0}, {-3.0, 2.0}}) {
        const auto a = near(imf[0], imf[1]);
        for (const std::array<double, 2>& d : {std::array<double, 2>{1e-3, 0.0}, {-1e-3, 0.0}, {0.0, 1e-3}, {0.0, -1e-3}}) {
            const auto b = near(imf[0] + d[0], imf[1] + d[1]);
            for (int c = 0; c < 3; ++c) EXPECT_LT(std::fabs(a[c] - b[c]), 0.05) << "IMF step at " << imf[0] << "," << imf[1];
        }
    }
}

TEST(IrbemT96, SouthwardBzIsSweptDenselyAndTheFieldMovesSmoothly) {
    // The storm-time reconfiguration the model exists to capture: Bz from +5 down to -10 in
    // 0.25 nT steps, at the storm regime's other drivers. On the southward branch the field must
    // be EXACTLY linear in Bz (a theorem of eq. (1)), and the tail's B_x at midnight must load
    // monotonically as Bz turns south.
    const Position<Frame::GSM> p = at(-8.0, 0.0, 1.0);
    for (double ps : {0.0, 0.35}) {
        std::vector<double> bx;
        for (int k = 20; k >= -40; --k) {
            const double bz = 0.25 * k;
            const auto r = t96_field(p, ps, -100.0, 6.0, 5.0, bz);
            EXPECT_EQ(r.status, Status::Ok) << "Bz " << bz;
            bx.push_back(r.value.v[0]);
        }
        // Linear on the south branch: second differences vanish for Bz < 0.
        for (std::size_t i = 21; i + 1 < bx.size(); ++i) {
            EXPECT_NEAR((bx[i + 1] - bx[i]) - (bx[i] - bx[i - 1]), 0.0, 1e-9) << "second difference at index " << i;
        }
        // Continuous across Bz = 0 and never a jump anywhere.
        for (std::size_t i = 1; i < bx.size(); ++i) EXPECT_LT(std::fabs(bx[i] - bx[i - 1]), 2.0) << "step at index " << i;
    }
}

TEST(IrbemT96, PressureEntersContinuously) {
    const Position<Frame::GSM> p = at(6.0, 1.0, 1.5);
    for (double ps : {-0.3, 0.2}) {
        double prev = t96_field(p, ps, -30.0, 0.5, 2.0, -4.0).value.v[2];
        for (int k = 1; k <= 95; ++k) {
            const double pd = 0.5 + (0.1 * k);
            const double bz = t96_field(p, ps, -30.0, pd, 2.0, -4.0).value.v[2];
            EXPECT_LT(std::fabs(bz - prev), 1.5) << "Pdyn " << pd;
            prev = bz;
        }
    }
}

// ================================================================================================
// The physics the model exists for
// ================================================================================================

TEST(IrbemT96, TheDaysideFieldIsCompressedAndTheNightsideStretched) {
    // On the dayside the magnetopause currents ADD to the northward internal field; on the
    // nightside the tail current subtracts from it. At quiet drivers, zero tilt.
    const T96Amplitudes<double> a = t96_amplitudes<double>(0.0, 0.0, 2.0, 0.0, 0.0);
    const std::array<double, 3> day = t96_components<double>(a, 8.0, 0.0, 0.0);
    const std::array<double, 3> night = t96_components<double>(a, -8.0, 0.0, 0.0);
    std::printf("[ MEASURED ] quiet external Bz at x = +8: %+.3f nT; at x = -8: %+.3f nT\n", day[2], night[2]);
    EXPECT_GT(day[2], 0.0);
    EXPECT_LT(night[2], 0.0);
}

TEST(IrbemT96, TheRingCurrentDepressionDeepensLinearlyWithDst) {
    // The Dst family's signature: a southward external field at the inner-belt equator that grows
    // with |Dst|. Exactly linear, and of the right sign and order — a Dst of -100 nT depresses the
    // equatorial field at 3 Re by tens of nT in T96.
    const Position<Frame::GSM> p = at(-3.0, 0.0, 0.0);
    const double b0 = t96_field(p, 0.0, 0.0, 2.0, 0.0, 0.0).value.v[2];
    const double b100 = t96_field(p, 0.0, -100.0, 2.0, 0.0, 0.0).value.v[2];
    std::printf("[ MEASURED ] external Bz at (-3, 0, 0): Dst 0 %+.3f nT, Dst -100 %+.3f nT\n", b0, b100);
    EXPECT_LT(b100, b0);
    EXPECT_LT(b100 - b0, -15.0);
    EXPECT_GT(b100 - b0, -80.0);
}

TEST(IrbemT96, SouthwardBzLoadsTheTail) {
    // Southward IMF strengthens the tail current: the near-tail B_z depression at midnight deepens
    // from Bz = 0 to Bz = -10, and northward Bz does much less (the measured kink).
    const Position<Frame::GSM> p = at(-10.0, 0.0, 0.0);
    const double b0 = t96_field(p, 0.0, -20.0, 2.0, 0.0, 0.0).value.v[2];
    const double bs = t96_field(p, 0.0, -20.0, 2.0, 0.0, -10.0).value.v[2];
    const double bn = t96_field(p, 0.0, -20.0, 2.0, 0.0, 10.0).value.v[2];
    std::printf("[ MEASURED ] external Bz at (-10, 0, 0): Bz_IMF 0 %+.3f, -10 %+.3f, +10 %+.3f nT\n", b0, bs, bn);
    EXPECT_LT(bs, b0);
    EXPECT_LT(std::fabs(bn - b0), std::fabs(bs - b0)) << "the northward response must be the weaker branch";
}

TEST(IrbemT96, ByPenetratesAsADawnDuskField) {
    // The interconnection field: a By of +10 nT produces an external B_y of the same sign near
    // the Earth on the noon-midnight meridian, at a fraction of the IMF value.
    const Position<Frame::GSM> p = at(-4.0, 0.0, 0.0);
    const double by0 = t96_field(p, 0.0, -20.0, 2.0, 0.0, 0.0).value.v[1];
    const double by10 = t96_field(p, 0.0, -20.0, 2.0, 10.0, 0.0).value.v[1];
    std::printf("[ MEASURED ] external By at (-4, 0, 0) for IMF By = +10: %+.3f nT\n", by10 - by0);
    EXPECT_GT(by10 - by0, 0.5);
    EXPECT_LT(by10 - by0, 10.0);
}

// ================================================================================================
// The API surface and the envelope
// ================================================================================================

TEST(IrbemT96, ReferenceLaneMatchesTheComponentForm) {
    for (double ps : kTilts) {
        const T96Amplitudes<double> a = t96_amplitudes<double>(ps, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
        const std::array<double, 3> raw = t96_components<double>(a, 4.5, -2.25, 1.5);
        const cheatah::fixarray::vec3d wrapped = t96_field_at(at(4.5, -2.25, 1.5), a).v;
        EXPECT_EQ(raw[0], wrapped[0]);
        EXPECT_EQ(raw[1], wrapped[1]);
        EXPECT_EQ(raw[2], wrapped[2]);
        const auto viaField = t96_field(at(4.5, -2.25, 1.5), ps, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
        EXPECT_EQ(viaField.status, Status::Ok);
        EXPECT_EQ(viaField.value.v[2], raw[2]);
    }
}

TEST(IrbemT96, ValidityIsCheckedFromBothSidesOfEveryBound) {
    // The published envelope: -100 <= Dst <= 20, 0.5 <= Pdyn <= 10, |By| <= 10, |Bz| <= 10. Closed
    // intervals: ON a bound is Ok, one ulp-ish beyond is OutOfValidityRange, and the value is still
    // computed and non-zero on both sides.
    const Position<Frame::GSM> p = at(-6.6, 1.0, 0.5);
    const double tilt = 0.2;
    const auto check = [&](Drivers inside, Drivers outside, const char* what) {
        const auto in = t96_field(p, tilt, inside.dst, inside.pdyn, inside.by, inside.bz);
        const auto out = t96_field(p, tilt, outside.dst, outside.pdyn, outside.by, outside.bz);
        EXPECT_EQ(in.status, Status::Ok) << what << " inside";
        EXPECT_EQ(out.status, Status::OutOfValidityRange) << what << " outside";
        EXPECT_NE(in.value.v[2], 0.0) << what;
        EXPECT_NE(out.value.v[2], 0.0) << what << ": the value must still be computed";
    };
    check({-100.0, 2.0, 0.0, 0.0}, {-100.001, 2.0, 0.0, 0.0}, "Dst low");
    check({20.0, 2.0, 0.0, 0.0}, {20.001, 2.0, 0.0, 0.0}, "Dst high");
    check({-20.0, 0.5, 0.0, 0.0}, {-20.0, 0.499, 0.0, 0.0}, "Pdyn low");
    check({-20.0, 10.0, 0.0, 0.0}, {-20.0, 10.001, 0.0, 0.0}, "Pdyn high");
    check({-20.0, 2.0, -10.0, 0.0}, {-20.0, 2.0, -10.001, 0.0}, "By low");
    check({-20.0, 2.0, 10.0, 0.0}, {-20.0, 2.0, 10.001, 0.0}, "By high");
    check({-20.0, 2.0, 0.0, -10.0}, {-20.0, 2.0, 0.0, -10.001}, "Bz low");
    check({-20.0, 2.0, 0.0, 10.0}, {-20.0, 2.0, 0.0, 10.001}, "Bz high");
}

TEST(IrbemT96, OutOfRangeDriversAreReportedButStillEvaluated) {
    // The corpus's storm and extreme regimes are past the envelope on purpose: the model must say
    // so and still answer, and the answer must be the affine/homogeneous extrapolation.
    const Position<Frame::GSM> p = at(-6.6, 0.0, 0.0);
    for (const corpus::MagInput& m : {corpus::regime_drivers[2], corpus::regime_drivers[3]}) {
        const auto r = t96_field(p, 0.2, m.dst, m.pdyn, m.by_imf, m.bz_imf);
        EXPECT_EQ(r.status, Status::OutOfValidityRange);
        EXPECT_TRUE(std::isfinite(r.value.v[2]));
        EXPECT_NE(r.value.v[2], 0.0);
    }
    for (const corpus::MagInput& m : {corpus::regime_drivers[0], corpus::regime_drivers[1]}) {
        EXPECT_EQ(t96_field(p, 0.2, m.dst, m.pdyn, m.by_imf, m.bz_imf).status, Status::Ok);
    }
}

TEST(IrbemT96, PositionOutsideTheFittedRegionIsReported) {
    // T96's published spatial envelope is r_GEO <= 40 R_E — the radius the oracle itself refuses
    // beyond. Just inside is Ok, just outside is a caveat with a finite value, inside the Earth is
    // a domain error.
    EXPECT_EQ(t96_field(at(-39.9, 0.0, 0.0), 0.1, -20.0, 2.0, 0.0, 0.0).status, Status::Ok);
    const auto far = t96_field(at(-40.1, 0.0, 0.0), 0.1, -20.0, 2.0, 0.0, 0.0);
    EXPECT_EQ(far.status, Status::OutOfValidityRange);
    EXPECT_TRUE(std::isfinite(far.value.v[2]));
    EXPECT_EQ(t96_field(at(0.5, 0.0, 0.0), 0.1, -20.0, 2.0, 0.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(t96_field(at(0.0, 0.0, 0.0), 0.1, -20.0, 2.0, 0.0, 0.0).status, Status::DomainError);
}

TEST(IrbemT96, NonFiniteInputIsADomainError) {
    const double nan = std::nan("");
    const double inf = std::numeric_limits<double>::infinity();
    for (const Position<Frame::GSM>& p : {at(nan, 0.0, 0.0), at(0.0, inf, 0.0), at(4.0, 0.0, -nan)}) {
        const auto r = t96_field(p, 0.2, -20.0, 2.0, 0.0, 0.0);
        EXPECT_EQ(r.status, Status::DomainError);
        EXPECT_EQ(r.value.v[0], 0.0);
        EXPECT_EQ(r.value.v[1], 0.0);
        EXPECT_EQ(r.value.v[2], 0.0);
    }
    const Position<Frame::GSM> p = at(5.0, 0.0, 0.0);
    EXPECT_EQ(t96_field(p, nan, -20.0, 2.0, 0.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(t96_field(p, 0.2, nan, 2.0, 0.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(t96_field(p, 0.2, -20.0, inf, 0.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(t96_field(p, 0.2, -20.0, 2.0, nan, 0.0).status, Status::DomainError);
    EXPECT_EQ(t96_field(p, 0.2, -20.0, 2.0, 0.0, -inf).status, Status::DomainError);
}

TEST(IrbemT96, RightAngleTiltAndNegativePressureAreDomainErrors) {
    const double quarter_turn = std::numbers::pi / 2.0;
    const Position<Frame::GSM> p = at(5.0, 0.0, 0.0);
    EXPECT_EQ(t96_field(p, quarter_turn, -20.0, 2.0, 0.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(t96_field(p, -quarter_turn, -20.0, 2.0, 0.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(t96_field(p, std::nextafter(quarter_turn, 0.0), -20.0, 2.0, 0.0, 0.0).status, Status::Ok);
    // A negative pressure has no square root: refused, not extrapolated.
    EXPECT_EQ(t96_field(p, 0.2, -20.0, -1.0, 0.0, 0.0).status, Status::DomainError);
    EXPECT_EQ(t96_field(p, 0.2, -20.0, 0.0, 0.0, 0.0).status, Status::OutOfValidityRange);
}

TEST(IrbemT96, AnOverflowingExtrapolationIsADomainErrorNotANaN) {
    // The box harmonics carry exp(kappa x); far enough sunward they overflow and the assembly makes
    // a NaN. Refused at the boundary, as T89's eq. (20) overflow is.
    for (double x : {1.0e5, 1.0e6, 1.0e30}) {
        const auto r = t96_field(at(x, 0.0, 0.0), 0.2, -20.0, 2.0, 0.0, 0.0);
        EXPECT_EQ(r.status, Status::DomainError) << "x = " << x;
        EXPECT_EQ(r.value.v[0], 0.0);
        EXPECT_EQ(r.value.v[2], 0.0);
    }
    // Far out but representable stays what it was: an extrapolation, reported and returned.
    const auto big = t96_field(at(1.0e2, 0.0, 0.0), 0.2, -20.0, 2.0, 0.0, 0.0);
    EXPECT_EQ(big.status, Status::OutOfValidityRange);
    EXPECT_TRUE(std::isfinite(big.value.v[2]));
    EXPECT_EQ(t96_field(at(-1.0e5, 0.0, 0.0), 0.2, -20.0, 2.0, 0.0, 0.0).status, Status::OutOfValidityRange);
}

TEST(IrbemT96, ContextOverloadAgreesWithTheExplicitOne) {
    const ir::Epoch epoch{2015.5, 43200.0, 2015, 180};
    ir::RotationTable identity{};
    for (cheatah::fixarray::mat3d& m : identity) m = cheatah::fixarray::mat3d::identity();
    ir::DriverSet drivers{};
    drivers[static_cast<std::size_t>(ir::Driver::Dst)] = -47.0;
    drivers[static_cast<std::size_t>(ir::Driver::Pdyn)] = 3.5;
    drivers[static_cast<std::size_t>(ir::Driver::ByIMF)] = 2.5;
    drivers[static_cast<std::size_t>(ir::Driver::BzIMF)] = -6.5;
    const ir::ContextResult built = ir::make_field_context(epoch, -0.42, identity, drivers);
    ASSERT_TRUE(built.has_value()) << ir::describe(built.error());
    const Position<Frame::GSM> p = at(3.75, -1.5, 2.25);
    const auto viaContext = t96_field(p, built.value());
    const auto viaScalars = t96_field(p, -0.42, -47.0, 3.5, 2.5, -6.5);
    EXPECT_EQ(viaContext.status, viaScalars.status);
    EXPECT_EQ(viaContext.value.v[0], viaScalars.value.v[0]);
    EXPECT_EQ(viaContext.value.v[1], viaScalars.value.v[1]);
    EXPECT_EQ(viaContext.value.v[2], viaScalars.value.v[2]);
}

// ================================================================================================
// The batch lanes
// ================================================================================================

TEST(IrbemT96, HostFloatLaneTracksTheReferenceLane) {
    const std::vector<Position<Frame::GSM>> pts = scatter(4096);
    std::vector<float> pos(3 * pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    std::vector<float> out(3 * pts.size());
    const double ps = 0.35;
    const T96Amplitudes<float> af = t96_amplitudes<float>(ps, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
    const T96Amplitudes<double> ad = t96_amplitudes<double>(ps, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
    ASSERT_TRUE(t96_field_host(pos, out, af));
    double worst_abs = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const std::array<double, 3> ref = t96_components<double>(ad, pos[(3 * i) + 0], pos[(3 * i) + 1], pos[(3 * i) + 2]);
        for (std::size_t c = 0; c < 3; ++c) worst_abs = std::max(worst_abs, std::fabs(out[(3 * i) + c] - ref[c]));
    }
    std::printf("[ MEASURED ] fp32 host lane vs fp64 reference over 4096 points: max |dB| %.3e nT\n", worst_abs);
    EXPECT_LT(worst_abs, 5e-2);
}

TEST(IrbemT96, HostFloatLaneRejectsMismatchedSpans) {
    const T96Amplitudes<float> a = t96_amplitudes<float>(0.1, -10.0, 2.0, 0.0, 0.0);
    std::vector<float> pos(7, 0.0F);
    std::vector<float> out(7, 0.0F);
    EXPECT_FALSE(t96_field_host(pos, out, a));
    std::vector<float> pos6(6, 1.0F);
    std::vector<float> out3(3, 0.0F);
    EXPECT_FALSE(t96_field_host(pos6, out3, a));
    std::vector<float> out6(6, 0.0F);
    EXPECT_TRUE(t96_field_host(pos6, out6, a));
}

TEST(IrbemT96, ParameterBlockCarriesTheTiltThenTheAmplitudes) {
    const T96Amplitudes<float> a = t96_amplitudes<float>(0.25, -33.0, 2.5, 1.5, -4.0);
    const std::array<float, t96_param_count> block = t96_param_block(a);
    EXPECT_EQ(block.size(), 130U);
    EXPECT_EQ(block[0], a.sin_tilt);
    EXPECT_EQ(block[1], a.cos_tilt);
    for (std::size_t j = 0; j < t96_basis_count; ++j) EXPECT_EQ(block[2 + j], a.c[j]) << j;
}

TEST(IrbemT96, BatchAgreesWithTheReferenceLane) {
    const std::vector<Position<Frame::GSM>> pts = scatter(1000);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    const ir::Result<bool> r = t96_field_batch(pts, 0.28, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz, out);
    EXPECT_EQ(r.status, Status::Ok);
    const T96Amplitudes<double> a = t96_amplitudes<double>(0.28, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const ir::FieldVector<Frame::GSM> ref = t96_field_at(pts[i], a);
        if (r.value) {
            for (int c = 0; c < 3; ++c) EXPECT_NEAR(out[i].v[c], ref.v[c], 5e-2) << "device lane, point " << i;
        } else {
            EXPECT_EQ(out[i].v[0], ref.v[0]);
            EXPECT_EQ(out[i].v[1], ref.v[1]);
            EXPECT_EQ(out[i].v[2], ref.v[2]);
        }
    }
}

TEST(IrbemT96, BatchRejectsMismatchedSpans) {
    const std::vector<Position<Frame::GSM>> pts = scatter(4);
    std::vector<ir::FieldVector<Frame::GSM>> shorter(3);
    EXPECT_EQ(t96_field_batch(pts, 0.1, -20.0, 2.0, 0.0, 0.0, shorter).status, Status::DomainError);
    std::vector<ir::FieldVector<Frame::GSM>> right(4);
    EXPECT_EQ(t96_field_batch(pts, std::nan(""), -20.0, 2.0, 0.0, 0.0, right).status, Status::DomainError);
    EXPECT_EQ(t96_field_batch(pts, 0.1, std::nan(""), 2.0, 0.0, 0.0, right).status, Status::DomainError);
    EXPECT_EQ(t96_field_batch(pts, 0.1, -20.0, -2.0, 0.0, 0.0, right).status, Status::DomainError);
    EXPECT_EQ(t96_field_batch(pts, std::numbers::pi / 2.0, -20.0, 2.0, 0.0, 0.0, right).status, Status::DomainError);
    const ir::Result<bool> empty = t96_field_batch({}, 0.1, -20.0, 2.0, 0.0, 0.0, {});
    EXPECT_EQ(empty.status, Status::Ok);
    EXPECT_FALSE(empty.value);
    EXPECT_EQ(t96_field_batch(pts, 0.1, -300.0, 2.0, 0.0, 0.0, right).status, Status::OutOfValidityRange);
}

TEST(IrbemT96, BatchReportsTheSameEnvelopeTheScalarLaneDoes) {
    const std::vector<Position<Frame::GSM>> good{at(5.0, 1.0, 1.0), at(-8.0, 2.0, -1.0)};
    std::vector<ir::FieldVector<Frame::GSM>> out(good.size());
    EXPECT_EQ(t96_field_batch(good, 0.2, -20.0, 2.0, 0.0, 0.0, out).status, Status::Ok);

    const std::vector<Position<Frame::GSM>> far{at(5.0, 1.0, 1.0), at(-41.0, 0.0, 0.0)};
    EXPECT_EQ(t96_field_batch(far, 0.2, -20.0, 2.0, 0.0, 0.0, out).status, Status::OutOfValidityRange);
    EXPECT_EQ(t96_field(far[1], 0.2, -20.0, 2.0, 0.0, 0.0).status, Status::OutOfValidityRange);
    for (const ir::FieldVector<Frame::GSM>& b : out) EXPECT_NE(b.v[2], 0.0);

    for (const Position<Frame::GSM>& bad : {at(0.5, 0.0, 0.0), at(std::nan(""), 0.0, 0.0),
                                            at(std::numeric_limits<double>::infinity(), 0.0, 0.0)}) {
        const std::vector<Position<Frame::GSM>> mixed{at(5.0, 1.0, 1.0), bad};
        std::vector<ir::FieldVector<Frame::GSM>> mixed_out(mixed.size(), ir::FieldVector<Frame::GSM>{cheatah::fixarray::vec3d{1.0, 1.0, 1.0}});
        EXPECT_EQ(t96_field_batch(mixed, 0.2, -20.0, 2.0, 0.0, 0.0, mixed_out).status, Status::DomainError);
        EXPECT_EQ(t96_field(bad, 0.2, -20.0, 2.0, 0.0, 0.0).status, Status::DomainError);
        for (const ir::FieldVector<Frame::GSM>& b : mixed_out) {
            EXPECT_EQ(b.v[0], 0.0);
            EXPECT_EQ(b.v[1], 0.0);
            EXPECT_EQ(b.v[2], 0.0);
        }
    }
    const std::vector<Position<Frame::GSM>> two_bad{at(std::nan(""), 0.0, 0.0), at(std::numeric_limits<double>::infinity(), 1.0, 0.0), at(5.0, 1.0, 1.0)};
    std::vector<ir::FieldVector<Frame::GSM>> two_out(two_bad.size());
    EXPECT_EQ(t96_field_batch(two_bad, 0.2, -20.0, 2.0, 0.0, 0.0, two_out).status, Status::DomainError);
    // An out-of-range driver is still reported on an EMPTY batch.
    const std::vector<Position<Frame::GSM>> none;
    std::vector<ir::FieldVector<Frame::GSM>> none_out;
    EXPECT_EQ(t96_field_batch(none, 0.2, -20.0, 30.0, 0.0, 0.0, none_out).status, Status::OutOfValidityRange);
    EXPECT_EQ(t96_field_batch(none, 0.2, -20.0, 2.0, 0.0, 0.0, none_out).status, Status::Ok);
}

TEST(IrbemT96, NothingOnTheHeapInTheHotPath) {
    const std::vector<Position<Frame::GSM>> pts = scatter(256);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    (void)t96_field_batch(pts, 0.2, -30.0, 2.0, 1.0, -2.0, out);
    (void)t96_field(at(5.0, 1.0, 1.0), 0.2, -30.0, 2.0, 1.0, -2.0);
    const T96Amplitudes<double> a = t96_amplitudes<double>(0.2, -30.0, 2.0, 1.0, -2.0);

    const std::size_t before = cheatah_space_test::allocation_count();
    for (int i = 0; i < 64; ++i) {
        sink = sink + t96_field(at(5.0 + (0.01 * i), 1.0, 1.0), 0.2, -30.0, 2.0, 1.0, -2.0).value.v[2];
        sink = sink + t96_field_at(at(5.0 - (0.01 * i), 1.0, 1.0), a).v[2];
    }
    (void)t96_field_batch(pts, 0.2, -30.0, 2.0, 1.0, -2.0, out);  // below the crossover: the host lane
    sink = sink + out[0].v[2];
    EXPECT_EQ(before, cheatah_space_test::allocation_count());
}

// ================================================================================================
// The total field
// ================================================================================================

TEST(IrbemT96, TotalFieldSuperposesInternalAndExternal) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.0).value();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::TotalFieldT96<10> total(igrf, rot, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
    static_assert(ir::GeoFieldModel<ir::TotalFieldT96<10>>, "a tracer must be able to follow it");
    EXPECT_EQ(ir::TotalFieldT96<10>::degree, 10);
    EXPECT_EQ(igrf.g(1, 0), total.g(1, 0));
    EXPECT_EQ(igrf.h(2, 1), total.h(2, 1));
    EXPECT_EQ(igrf.g(1, 0), total.internal().g(1, 0));
    EXPECT_EQ(&total.rotations(), &rot);

    const ir::Position<Frame::GEO> p{cheatah::fixarray::vec3d{6.0, 0.0, 0.0}};
    const double internal_only = igrf.evaluate(p).magnitude();
    const double with_external = total.evaluate(p).magnitude();
    EXPECT_GT(std::fabs(with_external - internal_only) / internal_only, 1e-3);
    EXPECT_LT(std::fabs(with_external - internal_only) / internal_only, 0.5);
    EXPECT_EQ(total.external_status(p), Status::Ok);

    // The superposition is exactly the sum of the two parts in GEO: recompute it by hand.
    const ir::Position<Frame::GSM> p_gsm = ir::transform<Frame::GSM>(p, rot);
    const ir::FieldVector<Frame::GSM> ext = t96_field_at(p_gsm, total.amplitudes());
    const ir::FieldVector<Frame::GEO> ext_geo = ir::transform<Frame::GEO>(ext, rot);
    const ir::FieldVector<Frame::GEO> sum{igrf.evaluate(p).v + ext_geo.v};
    EXPECT_EQ(total.evaluate(p).v[0], sum.v[0]);
    EXPECT_EQ(total.evaluate(p).v[2], sum.v[2]);
}

TEST(IrbemT96, TotalFieldTracesAndReportsWhenTheExternalModelDeclines) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.0).value();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::Position<Frame::GEO> p{cheatah::fixarray::vec3d{6.0, 0.0, 0.0}};

    // A trace through the total field closes and the invariants move with activity.
    std::vector<double> bmin;
    for (const corpus::MagInput& m : corpus::regime_drivers) {
        const ir::TotalFieldT96<10> total(igrf, rot, m.dst, m.pdyn, m.by_imf, m.bz_imf);
        const auto t = ir::trace_invariant(total, p, 45.0);
        ASSERT_EQ(Status::Ok, t.status) << "Dst = " << m.dst;
        EXPECT_GT(t.value.invariant_i, 0.0);
        bmin.push_back(t.value.b_min);
    }
    EXPECT_NE(bmin.front(), bmin[2]);
    EXPECT_LT(bmin[2], bmin.front()) << "the ring current should reduce Bmin during a storm";

    // The extreme regime is outside the envelope: the field still evaluates, the status says so.
    const corpus::MagInput& x = corpus::regime_drivers[3];
    const ir::TotalFieldT96<10> extreme(igrf, rot, x.dst, x.pdyn, x.by_imf, x.bz_imf);
    EXPECT_EQ(extreme.external_status(p), Status::OutOfValidityRange);
    EXPECT_NE(extreme.evaluate(p).v[2], igrf.evaluate(p).v[2]);

    // Inside the Earth the external model declines and the internal field is returned alone.
    const ir::Position<Frame::GEO> deep{cheatah::fixarray::vec3d{0.5, 0.0, 0.0}};
    const ir::TotalFieldT96<10> total(igrf, rot, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
    EXPECT_EQ(total.external_status(deep), Status::DomainError);
    EXPECT_EQ(total.evaluate(deep).v[2], igrf.evaluate(deep).v[2]);
    // A non-finite point likewise.
    const ir::Position<Frame::GEO> bad{cheatah::fixarray::vec3d{std::nan(""), 0.0, 0.0}};
    EXPECT_EQ(total.external_status(bad), Status::DomainError);
    // A negative pressure at construction: the external model declines everywhere.
    const ir::TotalFieldT96<10> vacuum(igrf, rot, 0.0, -1.0, 0.0, 0.0);
    EXPECT_EQ(vacuum.external_status(p), Status::DomainError);
    EXPECT_EQ(vacuum.evaluate(p).v[2], igrf.evaluate(p).v[2]);
    // Far beyond the envelope the value is still an extrapolation, not a refusal.
    const ir::Position<Frame::GEO> far{cheatah::fixarray::vec3d{-35.0, 0.0, 0.0}};
    EXPECT_EQ(total.external_status(far), Status::Ok);
    const ir::Position<Frame::GEO> farther{cheatah::fixarray::vec3d{-45.0, 0.0, 0.0}};
    EXPECT_EQ(total.external_status(farther), Status::OutOfValidityRange);
}

// ================================================================================================
// The differential against the IRBEM oracle
// ================================================================================================

TEST(IrbemT96, DiffersFromTheIrbemOracleByTheMeasuredEnvelope) {
    const Oracle& o = oracle();
    if (!o.usable()) GTEST_SKIP() << "IRBEM oracle not present (set CHEATAH_SPACE_IRBEM_ORACLE to its .so)";
    // The caps below are MEASUREMENTS of the documented published-form gap (see the header's file
    // brief and tools/oracle/t96_diff.cpp `report`), ~1.3x above the measured relative RMS in the
    // belts (3-8 Re), so a regression in THIS implementation fails while the known gap does not.
    // The drivers are the corpus regimes, the two beyond the envelope clamped to it.
    struct Regime {
        const char* name;
        Drivers d;
        double cap_rel;
    };
    const std::array<Regime, 4> regimes{{
        {"quiet", {-8.0, 1.8, 1.0, 2.0}, T96_CAP_QUIET},
        {"moderate", {-42.0, 3.2, -4.0, -5.0}, T96_CAP_MODERATE},
        {"storm (clamped)", {-100.0, 9.0, 8.0, -10.0}, T96_CAP_STORM},
        {"extreme (clamped)", {-100.0, 10.0, -10.0, -10.0}, T96_CAP_EXTREME},
    }};
    const std::vector<Position<Frame::GSM>> pts = differential_points();
    for (const Regime& rg : regimes) {
        double sum2 = 0.0;
        double sig2 = 0.0;
        std::size_t n = 0;
        for (const Epoch& e : kEpochs) {
            const double ps = o.tilt(e.doy, e.ut);
            const T96Amplitudes<double> a = t96_amplitudes<double>(ps, rg.d.dst, rg.d.pdyn, rg.d.by, rg.d.bz);
            for (const Position<Frame::GSM>& p : pts) {
                const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
                if (r < 3.0 || r > 8.0) continue;
                const std::array<double, 3> ora = o.ext(e.doy, e.ut, p.v[0], p.v[1], p.v[2], rg.d.dst, rg.d.pdyn, rg.d.by, rg.d.bz);
                if (!std::isfinite(ora[0])) continue;
                const std::array<double, 3> mine = t96_components<double>(a, p.v[0], p.v[1], p.v[2]);
                for (std::size_t c = 0; c < 3; ++c) {
                    sum2 += (mine[c] - ora[c]) * (mine[c] - ora[c]);
                    sig2 += ora[c] * ora[c];
                }
                ++n;
            }
        }
        ASSERT_GT(n, 0U);
        const double rel = std::sqrt(sum2 / sig2);
        std::printf("[ MEASURED ] %-18s belts 3-8 Re: RMS |dB| = %6.3f nT, %.1f%% of the oracle's external field (%zu points)\n",
                    rg.name, std::sqrt(sum2 / static_cast<double>(n)), 100.0 * rel, n);
        EXPECT_LT(rel, rg.cap_rel) << rg.name << ": the gap against the oracle has GROWN beyond the documented published-form floor";
    }
}

TEST(IrbemT96, StormSweepAgainstTheOracleStaysInsideTheMeasuredEnvelope) {
    const Oracle& o = oracle();
    if (!o.usable()) GTEST_SKIP() << "IRBEM oracle not present";
    // Southward Bz densely, at Dst -100, Pdyn 6, By +5: the belts-region gap must stay inside the
    // measured cap at every step, and the two models must move the SAME WAY as Bz turns south —
    // the midnight tail B_x at (-8, 0, 1) must load with the same sign in both.
    const std::vector<Position<Frame::GSM>> pts = differential_points();
    double prev_mine = 0.0;
    double prev_ora = 0.0;
    for (int k = 8; k >= -40; k -= 4) {
        const double bz = 0.25 * k;
        double sum2 = 0.0;
        double sig2 = 0.0;
        for (const Epoch& e : kEpochs) {
            const double ps = o.tilt(e.doy, e.ut);
            const T96Amplitudes<double> a = t96_amplitudes<double>(ps, -100.0, 6.0, 5.0, bz);
            for (const Position<Frame::GSM>& p : pts) {
                const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
                if (r < 3.0 || r > 8.0) continue;
                const std::array<double, 3> ora = o.ext(e.doy, e.ut, p.v[0], p.v[1], p.v[2], -100.0, 6.0, 5.0, bz);
                if (!std::isfinite(ora[0])) continue;
                const std::array<double, 3> mine = t96_components<double>(a, p.v[0], p.v[1], p.v[2]);
                for (std::size_t c = 0; c < 3; ++c) {
                    sum2 += (mine[c] - ora[c]) * (mine[c] - ora[c]);
                    sig2 += ora[c] * ora[c];
                }
            }
        }
        const double rel = std::sqrt(sum2 / sig2);
        EXPECT_LT(rel, T96_CAP_STORM_SWEEP) << "Bz = " << bz;
        const double ps0 = o.tilt(kEpochs[1].doy, kEpochs[1].ut);
        const double mine_bx = t96_components<double>(t96_amplitudes<double>(ps0, -100.0, 6.0, 5.0, bz), -8.0, 0.0, 1.0)[0];
        const double ora_bx = o.ext(kEpochs[1].doy, kEpochs[1].ut, -8.0, 0.0, 1.0, -100.0, 6.0, 5.0, bz)[0];
        if (k < 8 && bz < 0.0) {
            EXPECT_EQ(mine_bx - prev_mine > 0.0, ora_bx - prev_ora > 0.0) << "tail loading sense differs at Bz = " << bz;
        }
        prev_mine = mine_bx;
        prev_ora = ora_bx;
    }
}

#if CHEATAH_SPACE_IRBEM_T96_GPU
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

}  // namespace

TEST(IrbemT96, BatchFallsBackToTheHostWhenTheShaderWasNeverBuilt) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_t96_f32");
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    {
        const SpvDirScope nowhere(std::filesystem::temp_directory_path().string() + "/cheatah-space-no-such-shader-dir");
        const ir::Result<bool> r = t96_field_batch(pts, 0.28, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz, out);
        EXPECT_EQ(r.status, Status::Ok);
        EXPECT_FALSE(r.value) << "with no compiled shader the batch must run on the host";
    }
    const T96Amplitudes<double> a = t96_amplitudes<double>(0.28, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = t96_field_at(pts[i], a);
        ASSERT_EQ(out[i].v[0], ref.v[0]) << "point " << i;
        ASSERT_EQ(out[i].v[2], ref.v[2]) << "point " << i;
    }
}

TEST(IrbemT96, BatchUsesTheDeviceWhenOneIsAvailable) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_t96_f32");
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    const ir::Result<bool> r = t96_field_batch(pts, 0.28, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_TRUE(r.value) << "the batch fell back to the host with a device present";
    const T96Amplitudes<double> a = t96_amplitudes<double>(0.28, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = t96_field_at(pts[i], a);
        for (int c = 0; c < 3; ++c) worst = std::max(worst, std::fabs(out[i].v[c] - ref.v[c]));
    }
    std::printf("[ MEASURED ] device batch of %zu vs fp64 reference: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 5e-2);
    const std::vector<Position<Frame::GSM>> few = scatter(16);
    std::vector<ir::FieldVector<Frame::GSM>> few_out(few.size());
    EXPECT_FALSE(t96_field_batch(few, 0.28, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz, few_out).value);
}

TEST(IrbemT96, TheDeviceLaneRefusesABadPointBeforeItDispatches) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_t96_f32");
    std::vector<Position<Frame::GSM>> pts = scatter(n);
    pts[n / 3] = at(0.25, 0.0, 0.0);
    std::vector<ir::FieldVector<Frame::GSM>> out(n, ir::FieldVector<Frame::GSM>{cheatah::fixarray::vec3d{7.0, 7.0, 7.0}});
    const ir::Result<bool> r = t96_field_batch(pts, 0.28, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz, out);
    EXPECT_EQ(r.status, Status::DomainError);
    EXPECT_FALSE(r.value);
    for (const ir::FieldVector<Frame::GSM>& b : out) {
        ASSERT_EQ(b.v[0], 0.0);
        ASSERT_EQ(b.v[2], 0.0);
    }
    std::vector<Position<Frame::GSM>> far = scatter(n);
    far[n / 2] = at(-45.0, 0.0, 0.0);
    const ir::Result<bool> f = t96_field_batch(far, 0.28, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz, out);
    EXPECT_EQ(f.status, Status::OutOfValidityRange);
    EXPECT_TRUE(f.value) << "an out-of-validity batch must still be computed on the device";
}

TEST(IrbemT96, DeviceKernelAgreesWithTheHostLane) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 1 << 16;
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<float> pos(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    const T96Amplitudes<float> a = t96_amplitudes<float>(0.31, kModerate.dst, kModerate.pdyn, kModerate.by, kModerate.bz);
    std::vector<float> host(3 * n);
    std::vector<float> device(3 * n);
    ASSERT_TRUE(t96_field_host(pos, host, a));
    const std::array<float, t96_param_count> block = t96_param_block(a);
    ir::gpu::dispatch_batch("irbem_t96_f32", pos, device, std::span<const float>(block));
    double worst = 0.0;
    for (std::size_t i = 0; i < 3 * n; ++i) worst = std::max(worst, std::fabs(static_cast<double>(device[i]) - host[i]));
    std::printf("[ MEASURED ] device vs host, %zu points: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 2e-2) << "the device is not evaluating the same expressions";
}
#endif
