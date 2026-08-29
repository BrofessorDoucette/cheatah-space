/// @file irbem_opq_test.cpp
/// @brief The suite for `space/irbem/ext_opq.hpp` — Olson & Pfitzer (1977), IRBEM's `kext = 5`.
///
/// This model is a least-squares power-series fit, so three different kinds of question are asked
/// of the implementation, and none of them is "does it look plausible":
///
///  - **Is it the published model?** The tables are checked for shape and for the mirror-symmetry
///    rule that fixes every tilt parity; the enumeration of monomials is checked by one-hot
///    evaluation against an independently derived exponent list; the fold is checked for its
///    even/odd structure bit for bit; the three radial rules are checked at their edges.
///  - **Is the arithmetic right?** The field is a polynomial times an exponential, so its
///    divergence exists in closed form; a second-order stencil of the implemented field must fall
///    as `h^2` onto it — three decades — and a deliberately broken coefficient must stop it doing
///    so. `div B` is NOT zero for this model, and the suite measures and states by how much rather
///    than asserting a property the published fit does not have.
///  - **Is it the oracle?** When IRBEM is present, `kext = 5` minus `kext = 0` is compared at three
///    tilts across the belts AND the taper region, with the corpus's four activity regimes swept to
///    show the oracle ignores every driver exactly as this model does, and the drift-shell chain is
///    run through the total field against IRBEM's own `make_lstar`. The provenance study in the
///    header's brief says the two are the same model coefficient for coefficient; these tests hold
///    the implementation to that at a tolerance three orders tighter than the 1e-6 budget.
///
/// The oracle tests `dlopen` IRBEM at runtime rather than linking it: this binary must build and
/// pass on a machine that has never heard of IRBEM. They carry the `HeavyDifferential` prefix the
/// Valgrind gate excludes by policy.

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
#include "space/irbem/driftshell.hpp"
#include "space/irbem/ext_opq.hpp"

namespace {

namespace ir = cheatah::space::irbem;
namespace corpus = cheatah_space_test;

using ir::Frame;
using ir::OpqParameters;
using ir::opq_components;
using ir::opq_field;
using ir::opq_field_at;
using ir::opq_field_batch;
using ir::opq_field_host;
using ir::opq_param_block;
using ir::opq_param_count;
using ir::opq_parameters;
using ir::opq_table;
using ir::opq_xz_count;
using ir::opq_y_count;
using ir::Position;
using ir::Status;

/// A sink the optimizer cannot see through, so the allocation test's calls actually happen.
volatile double sink = 0.0;

/// Degrees per radian, spelled once.
constexpr double kDegPerRad = 180.0 / std::numbers::pi;

/// A GSM point, spelled so a test reads as coordinates rather than as a constructor call.
Position<Frame::GSM> at(double x, double y, double z) {
    return Position<Frame::GSM>{cheatah::fixarray::vec3d(x, y, z)};
}

/// The tilts the suite sweeps, radians: zero, and both signs of a realistic seasonal-diurnal
/// excursion. Every odd-in-tilt pair is invisible at zero tilt, so a single tilt tests half the
/// tables.
constexpr std::array<double, 5> kTilts{-0.55, -0.21, 0.0, 0.21, 0.55};

/// One monomial's exponents: `x^p y^{2q} z^s`.
struct Mono {
    int p, q, s;
};

/// The published enumeration, re-derived from the report's loop rules INDEPENDENTLY of the header
/// (the header keeps the report's counters; this keeps the exponents), so the two can disagree.
void enumerate(std::vector<Mono>& xz, std::vector<Mono>& y) {
    for (int i = 1; i <= 5; ++i) {
        for (int j = 1; j <= 3; ++j) {
            if (i + (2 * j) > 8) break;
            int ijk = i + (2 * j) + 1;
            int k = 0;
            for (;;) {
                xz.push_back({i - 1, j - 1, k});
                if (ijk > 8) break;
                y.push_back({i - 1, j - 1, k});
                ++ijk;
                ++k;
                if (ijk > 9 || k > 4) break;
            }
        }
    }
}

/// `b^e` by repeated multiplication — exact for the small integer exponents here.
double ipow(double b, int e) {
    double r = 1.0;
    for (int i = 0; i < e; ++i) r *= b;
    return r;
}

/// The closed-form divergence of the model in SM, valid where the taper is 1 (2.5 < r < 15).
double analytic_divergence(const OpqParameters<double>& p, double x, double y, double z) {
    static const std::vector<std::vector<Mono>> tables = [] {
        std::vector<std::vector<Mono>> t(2);
        enumerate(t[0], t[1]);
        return t;
    }();
    const std::vector<Mono>& xz = tables[0];
    const std::vector<Mono>& yy = tables[1];
    const double r2 = (x * x) + (y * y) + (z * z);
    const double env = std::exp(-ir::opq_exp_rate * r2);
    const double denv = -2.0 * ir::opq_exp_rate * env;  // d(env)/dx = denv * x, etc.
    double d = 0.0;
    for (std::size_t i = 0; i < xz.size(); ++i) {
        const Mono& m = xz[i];
        const double mono = ipow(x, m.p) * ipow(y, 2 * m.q) * ipow(z, m.s);
        const double dx = m.p > 0 ? m.p * ipow(x, m.p - 1) * ipow(y, 2 * m.q) * ipow(z, m.s) : 0.0;
        const double dz = m.s > 0 ? m.s * ipow(x, m.p) * ipow(y, 2 * m.q) * ipow(z, m.s - 1) : 0.0;
        d += (p.a[i] * dx) + (p.b[i] * ((dx * env) + (mono * denv * x)));
        d += (p.e[i] * dz) + (p.f[i] * ((dz * env) + (mono * denv * z)));
    }
    for (std::size_t j = 0; j < yy.size(); ++j) {
        const Mono& m = yy[j];
        const double mono = ipow(x, m.p) * ipow(y, 2 * m.q) * ipow(z, m.s);
        const double dy =
            m.q > 0 ? 2 * m.q * ipow(x, m.p) * ipow(y, (2 * m.q) - 1) * ipow(z, m.s) : 0.0;
        const double c = p.c[j] + (p.d[j] * env);
        d += (c * mono) + (y * ((c * dy) + (p.d[j] * mono * denv * y)));
    }
    return d;
}

/// `div B` at one SM point by central differences of the implemented field.
double stencil_divergence(const OpqParameters<double>& p, double x, double y, double z, double h) {
    const auto b = [&](double a, double c, double e) { return opq_components<double>(p, a, c, e); };
    return ((b(x + h, y, z)[0] - b(x - h, y, z)[0]) + (b(x, y + h, z)[1] - b(x, y - h, z)[1]) +
            (b(x, y, z + h)[2] - b(x, y, z - h)[2])) /
           (2.0 * h);
}

/// The box the divergence tests sample: inside the published region, clear of the taper. Integer
/// induction with coordinates derived per iteration, never a floating accumulator.
template <class F>
void for_each_box_point(const F& f) {
    for (int ix = 0; ix <= 8; ++ix)
        for (int iy = 0; iy <= 6; ++iy)
            for (int iz = 0; iz <= 5; ++iz) {
                const double x = -13.0 + (3.1 * ix);
                const double y = -9.0 + (3.0 * iy);
                const double z = -7.0 + (2.8 * iz);
                const double r = std::sqrt((x * x) + (y * y) + (z * z));
                if (r < 2.8 || r > 14.5) continue;
                f(x, y, z, r);
            }
}

/// A deterministic scatter of GSM points over the model's region. A 64-bit LCG so a disagreement
/// is reproducible, and never a lattice, so no component is systematically zero.
std::vector<Position<Frame::GSM>> scatter(std::size_t n, double r_lo = 2.6, double r_hi = 14.8) {
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

// ---- the oracle, opened at runtime and never linked --------------------------------------------

/// `get_field1_`, `coord_trans_vec1_` and `make_lstar1_`, as the vendored `matlab/libirbem.h`
/// documents them.
using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);
using MakeLstar1 = void (*)(int*, int*, int*, int*, int*, int*, double*, double*, double*, double*,
                            double*, double*, double*, double*, double*, double*, double*);

/// The oracle handle, or nulls when IRBEM is not on this machine.
struct Oracle {
    void* handle = nullptr;
    GetField1 get_field = nullptr;
    CoordTransVec1 coord_trans = nullptr;
    MakeLstar1 make_lstar = nullptr;

    /// Whether the oracle can be called.
    [[nodiscard]] bool usable() const {
        return get_field != nullptr && coord_trans != nullptr && make_lstar != nullptr;
    }

    /// The oracle's own dipole tilt at an epoch, radians: the SM z axis in GSM is (sin, 0, cos).
    [[nodiscard]] double tilt(int year, int doy, double ut) const {
        int one = 1;
        int si = 4;
        int so = 2;
        std::array<double, 3> in{0.0, 0.0, 1.0};
        std::array<double, 3> out{};
        coord_trans(&one, &si, &so, &year, &doy, &ut, in.data(), out.data());
        return std::atan2(out[0], out[2]);
    }

    /// The oracle's external field (kext = 5 minus kext = 0) at a GSM point, in GSM, nT, under
    /// the given driver set. False when the oracle refused the point.
    [[nodiscard]] bool external(int year, int doy, double ut, const Position<Frame::GSM>& p,
                                const corpus::MagInput& m, std::array<double, 3>& out) const {
        int one = 1;
        std::array<double, 3> gsm{p.v[0], p.v[1], p.v[2]};
        std::array<double, 3> geo{};
        {
            int si = 2;
            int so = 1;
            coord_trans(&one, &si, &so, &year, &doy, &ut, gsm.data(), geo.data());
        }
        std::array<int, 5> options{0, 0, 0, 0, 0};
        int sysaxes = 1;
        int k0 = 0;
        int k5 = 5;
        std::vector<double> mag(25, 0.0);
        mag[0] = m.kp * 10.0;
        mag[1] = m.dst;
        mag[2] = m.density;
        mag[3] = m.velocity;
        mag[4] = m.pdyn;
        mag[5] = m.by_imf;
        mag[6] = m.bz_imf;
        mag[7] = m.g1;
        mag[8] = m.g2;
        mag[9] = m.g3;
        std::array<double, 3> b0{};
        std::array<double, 3> b5{};
        double m0 = 0.0;
        double m5 = 0.0;
        double x1 = geo[0];
        double x2 = geo[1];
        double x3 = geo[2];
        get_field(&k0, options.data(), &sysaxes, &year, &doy, &ut, &x1, &x2, &x3, mag.data(),
                  b0.data(), &m0);
        get_field(&k5, options.data(), &sysaxes, &year, &doy, &ut, &x1, &x2, &x3, mag.data(),
                  b5.data(), &m5);
        if (b5[0] < -1e30 || b0[0] < -1e30) return false;
        std::array<double, 3> dgeo{b5[0] - b0[0], b5[1] - b0[1], b5[2] - b0[2]};
        int si = 1;
        int so = 2;
        coord_trans(&one, &si, &so, &year, &doy, &ut, dgeo.data(), out.data());
        return true;
    }
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
        o->make_lstar = reinterpret_cast<MakeLstar1>(dlsym(o->handle, "make_lstar1_"));
        return o;
    }();
    return *opened;
}

/// The three epochs the oracle differentials sample: tilts of ~0, +25.6 and -30.4 degrees.
struct Epoch {
    int year;
    int doy;
    double ut;
};
constexpr std::array<Epoch, 3> kEpochs{{{2015, 80, 39183.0}, {2015, 180, 43200.0}, {2015, 355, 7200.0}}};

}  // namespace

// ================================================================================================
// The published tables
// ================================================================================================

TEST(IrbemOpq, TablesHaveThePublishedShape) {
    // Olson & Pfitzer (1977) p. 64: AA(64), BB(64), CC(44), DD(44), EE(64), FF(64), ITA(32),
    // ITB(22), ITC(32). The counts are the type; the corners are the printed values.
    static_assert(opq_xz_count == 32 && opq_y_count == 22);
    static_assert(opq_table.aa.size() == 64 && opq_table.cc.size() == 44);
    static_assert(opq_param_count == 174, "sin, cos, then 172 folded coefficients");
    EXPECT_EQ(opq_table.aa[0], -2.26836e-02);
    EXPECT_EQ(opq_table.aa[63], 1.25818e-07);
    EXPECT_EQ(opq_table.bb[0], 9.47753e-02);
    EXPECT_EQ(opq_table.cc[0], -1.88177e-02);
    EXPECT_EQ(opq_table.dd[2], 3.23821e+00);
    EXPECT_EQ(opq_table.ee[0], -2.77924e+01);
    EXPECT_EQ(opq_table.ff[0], -5.07092e+00);
    EXPECT_EQ(opq_table.ff[63], -6.10021e-08);
    // The report states the fit's leading term in so many words: the largest coefficient is the
    // constant B_z depression, -27.79 nT at zero tilt — a quiet-time ring current.
    double largest = 0.0;
    for (double v : opq_table.ee) largest = std::max(largest, std::fabs(v));
    EXPECT_EQ(largest, 2.77924e+01);
    // Every parity table entry is one of the report's two selector values, and the odd terms are
    // the majority in B_x/B_y (which vanish at the equator at zero tilt) and the minority in B_z.
    int odd_x = 0;
    int odd_y = 0;
    int odd_z = 0;
    for (bool b : opq_table.odd_x) odd_x += b ? 1 : 0;
    for (bool b : opq_table.odd_y) odd_y += b ? 1 : 0;
    for (bool b : opq_table.odd_z) odd_z += b ? 1 : 0;
    EXPECT_EQ(odd_x, 19);
    EXPECT_EQ(odd_y, 14);
    EXPECT_EQ(odd_z, 13);
}

TEST(IrbemOpq, ParityTablesFollowTheEquatorialMirror) {
    // At zero tilt the quiet magnetosphere is a mirror about the SM equator: B_x and B_y are odd
    // in z, B_z is even. A monomial z^s with the WRONG parity for its component can therefore only
    // enter multiplied by an odd power of the tilt. That rule fixes every one of the 86 selector
    // entries from the enumeration alone — and it is exactly what the printed tables say. A
    // transcription check that cannot pass by accident: one flipped entry fails it.
    std::vector<Mono> xz;
    std::vector<Mono> y;
    enumerate(xz, y);
    ASSERT_EQ(xz.size(), opq_xz_count);
    ASSERT_EQ(y.size(), opq_y_count);
    for (std::size_t i = 0; i < xz.size(); ++i) {
        const bool s_even = (xz[i].s % 2) == 0;
        EXPECT_EQ(opq_table.odd_x[i], s_even) << "B_x term " << i;   // B_x odd in z
        EXPECT_EQ(opq_table.odd_z[i], !s_even) << "B_z term " << i;  // B_z even in z
    }
    for (std::size_t j = 0; j < y.size(); ++j) {
        const bool s_even = (y[j].s % 2) == 0;
        EXPECT_EQ(opq_table.odd_y[j], s_even) << "B_y term " << j;   // B_y odd in z
    }
}

TEST(IrbemOpq, MonomialEnumerationHasThePublishedCounts) {
    // The evaluator keeps the report's counters; this suite keeps the exponents. One-hot
    // coefficient vectors make the evaluator print its k-th monomial, which must be the k-th
    // exponent triple here — order included, because that order is how the tables are laid out.
    std::vector<Mono> xz;
    std::vector<Mono> y;
    enumerate(xz, y);
    ASSERT_EQ(xz.size(), 32U);
    ASSERT_EQ(y.size(), 22U);
    const double x = 1.5;
    const double yy = -2.25;
    const double z = 0.75;
    for (std::size_t k = 0; k < xz.size(); ++k) {
        OpqParameters<double> one{};
        one.a[k] = 1.0;
        one.e[k] = 1.0;
        const std::array<double, 3> b = opq_components<double>(one, x, yy, z);
        const double expect = ipow(x, xz[k].p) * ipow(yy, 2 * xz[k].q) * ipow(z, xz[k].s);
        EXPECT_EQ(b[0], expect) << "B_x term " << k;
        EXPECT_EQ(b[2], expect) << "B_z term " << k;
        EXPECT_EQ(b[1], 0.0);
    }
    for (std::size_t k = 0; k < y.size(); ++k) {
        OpqParameters<double> one{};
        one.c[k] = 1.0;
        const std::array<double, 3> b = opq_components<double>(one, x, yy, z);
        const double expect = yy * ipow(x, y[k].p) * ipow(yy, 2 * y[k].q) * ipow(z, y[k].s);
        EXPECT_EQ(b[1], expect) << "B_y term " << k;
        EXPECT_EQ(b[0], 0.0);
    }
    // And the enveloped set is the same monomial under exp(-0.06 r^2).
    OpqParameters<double> env{};
    env.f[5] = 1.0;
    const double r2 = (x * x) + (yy * yy) + (z * z);
    EXPECT_EQ(opq_components<double>(env, x, yy, z)[2],
              std::exp(-0.06 * r2) * ipow(x, xz[5].p) * ipow(yy, 2 * xz[5].q) * ipow(z, xz[5].s));
}

// ================================================================================================
// The fold
// ================================================================================================

TEST(IrbemOpq, FoldedCoefficientsAreEvenOrOddInTheTilt) {
    // An even pair folds to a + b t^2, an odd one to a t + b t^3: under t -> -t the first is
    // unchanged and the second negated, BITWISE, because negating t negates t^3 exactly.
    for (double t : {7.5, 25.636, 34.0}) {
        const OpqParameters<double> plus = opq_parameters<double>(t);
        const OpqParameters<double> minus = opq_parameters<double>(-t);
        for (std::size_t i = 0; i < opq_xz_count; ++i) {
            const double sx = opq_table.odd_x[i] ? -1.0 : 1.0;
            const double sz = opq_table.odd_z[i] ? -1.0 : 1.0;
            EXPECT_EQ(minus.a[i], sx * plus.a[i]) << i;
            EXPECT_EQ(minus.b[i], sx * plus.b[i]) << i;
            EXPECT_EQ(minus.e[i], sz * plus.e[i]) << i;
            EXPECT_EQ(minus.f[i], sz * plus.f[i]) << i;
        }
        for (std::size_t j = 0; j < opq_y_count; ++j) {
            const double sy = opq_table.odd_y[j] ? -1.0 : 1.0;
            EXPECT_EQ(minus.c[j], sy * plus.c[j]) << j;
            EXPECT_EQ(minus.d[j], sy * plus.d[j]) << j;
        }
    }
    // At zero tilt every odd pair is exactly zero and every even pair is its first entry.
    const OpqParameters<double> flat = opq_parameters<double>(0.0);
    for (std::size_t i = 0; i < opq_xz_count; ++i) {
        EXPECT_EQ(flat.a[i], opq_table.odd_x[i] ? 0.0 : opq_table.aa[2 * i]);
        EXPECT_EQ(flat.e[i], opq_table.odd_z[i] ? 0.0 : opq_table.ee[2 * i]);
    }
    // And a worked value: the first B_x pair is odd, so at 10 degrees it is AA1*10 + AA2*1000.
    EXPECT_DOUBLE_EQ(opq_parameters<double>(10.0).a[0], (-2.26836e-02 * 10.0) + (-1.01863e-04 * 1000.0));
}

TEST(IrbemOpq, ParametersRoundTripThroughFloat) {
    // The float fold is the double fold rounded element by element — never a fold done in float.
    for (double t : kTilts) {
        const OpqParameters<double> d = opq_parameters<double>(t * kDegPerRad);
        const OpqParameters<float> f = opq_parameters<float>(t * kDegPerRad);
        for (std::size_t i = 0; i < opq_xz_count; ++i) {
            EXPECT_EQ(f.a[i], static_cast<float>(d.a[i]));
            EXPECT_EQ(f.b[i], static_cast<float>(d.b[i]));
            EXPECT_EQ(f.e[i], static_cast<float>(d.e[i]));
            EXPECT_EQ(f.f[i], static_cast<float>(d.f[i]));
        }
        for (std::size_t j = 0; j < opq_y_count; ++j) {
            EXPECT_EQ(f.c[j], static_cast<float>(d.c[j]));
            EXPECT_EQ(f.d[j], static_cast<float>(d.d[j]));
        }
    }
}

// ================================================================================================
// The mathematics: the symmetries, the divergence, the radial rules
// ================================================================================================

TEST(IrbemOpq, ZeroTiltIsMirrorSymmetricAboutTheEquator) {
    // z -> -z flips B_x and B_y and keeps B_z, bitwise: the odd pairs are exactly zero at zero
    // tilt and (-z)^s is exactly -+z^s.
    const OpqParameters<double> p = opq_parameters<double>(0.0);
    for (const Position<Frame::GSM>& q : scatter(64)) {
        const std::array<double, 3> up = opq_components<double>(p, q.v[0], q.v[1], q.v[2]);
        const std::array<double, 3> down = opq_components<double>(p, q.v[0], q.v[1], -q.v[2]);
        EXPECT_EQ(down[0], -up[0]);
        EXPECT_EQ(down[1], -up[1]);
        EXPECT_EQ(down[2], up[2]);
    }
    // And the symmetry is BROKEN by a tilt, which is the whole point of the 1977 model over the
    // 1974 one: at 25 degrees the equatorial mirror no longer holds.
    const OpqParameters<double> tilted = opq_parameters<double>(25.0);
    const std::array<double, 3> up = opq_components<double>(tilted, 6.0, 2.0, 3.0);
    const std::array<double, 3> down = opq_components<double>(tilted, 6.0, 2.0, -3.0);
    EXPECT_NE(down[2], up[2]);
}

TEST(IrbemOpq, DawnDuskSymmetryHoldsAtEveryTilt) {
    // y -> -y flips B_y and keeps B_x, B_z, bitwise and at EVERY tilt: y enters B_x and B_z only as
    // y^2 and B_y as y times a series in y^2.
    for (double t : kTilts) {
        const OpqParameters<double> p = opq_parameters<double>(t * kDegPerRad);
        for (const Position<Frame::GSM>& q : scatter(48)) {
            const std::array<double, 3> dusk = opq_components<double>(p, q.v[0], q.v[1], q.v[2]);
            const std::array<double, 3> dawn = opq_components<double>(p, q.v[0], -q.v[1], q.v[2]);
            EXPECT_EQ(dawn[0], dusk[0]);
            EXPECT_EQ(dawn[1], -dusk[1]);
            EXPECT_EQ(dawn[2], dusk[2]);
        }
        // On the noon-midnight meridian B_y is therefore exactly zero — in GSM too, since the
        // tilt rotation is about y.
        for (double x : {-12.0, -6.0, 3.0, 9.0}) {
            EXPECT_EQ(opq_field_at(at(x, 0.0, 2.0), t).v[1], 0.0);
        }
    }
}

TEST(IrbemOpq, StencilDivergenceConvergesToTheAnalyticOne) {
    // The field is a polynomial times exp(-0.06 r^2), so its divergence is a closed form the suite
    // computes from an independently derived exponent list. A central difference of the
    // implemented field must fall onto it as h^2 — three decades — which verifies every term of
    // the evaluator against its own derivative with no reference model in the loop.
    std::array<double, 3> worst_gap{};
    const std::array<double, 3> steps{1e-2, 1e-3, 1e-4};
    double worst_intrinsic_belts = 0.0;     // r <= 8: where every shell this library traces lives
    double worst_intrinsic_region = 0.0;    // the whole published region
    double worst_ratio_belts = 0.0;         // |div B| against the dipole's own gradient |B_dip|/r
    for (double t : kTilts) {
        const OpqParameters<double> p = opq_parameters<double>(t * kDegPerRad);
        for_each_box_point([&](double x, double y, double z, double r) {
            const double a = analytic_divergence(p, x, y, z);
            for (std::size_t k = 0; k < steps.size(); ++k) {
                worst_gap[k] = std::max(worst_gap[k], std::fabs(stencil_divergence(p, x, y, z, steps[k]) - a));
            }
            worst_intrinsic_region = std::max(worst_intrinsic_region, std::fabs(a));
            if (r <= 8.0) {
                worst_intrinsic_belts = std::max(worst_intrinsic_belts, std::fabs(a));
                worst_ratio_belts = std::max(worst_ratio_belts, std::fabs(a) / (30000.0 / (r * r * r * r)));
            }
        });
    }
    std::printf("[ MEASURED ] worst |stencil - analytic div B|: h=1e-2 %.3e  h=1e-3 %.3e  h=1e-4 %.3e"
                " nT/Re\n",
                worst_gap[0], worst_gap[1], worst_gap[2]);
    EXPECT_GT(worst_gap[0] / worst_gap[1], 50.0) << "not falling as h^2 — a term is not its derivative";
    EXPECT_GT(worst_gap[1] / worst_gap[2], 50.0) << "not falling as h^2 — a term is not its derivative";
    EXPECT_LT(worst_gap[2], 1e-6);

    // What the published fit's divergence actually IS, stated rather than hidden: this model is
    // not the curl of anything, and the figure matters to anyone who integrates it.
    std::printf("[ MEASURED ] intrinsic |div B| of the published fit: worst %.2f nT/Re inside 8 Re "
                "(%.3f of the dipole's |B|/r there), worst %.2f nT/Re over the whole region to 14.5 Re\n",
                worst_intrinsic_belts, worst_ratio_belts, worst_intrinsic_region);
    EXPECT_LT(worst_intrinsic_belts, 20.0);
    EXPECT_GT(worst_intrinsic_belts, 1.0) << "the fit is not divergence-free; a zero here means "
                                             "the analytic form was zeroed, not that the model improved";

    // The perturbation control: break one coefficient in the EVALUATOR's parameter set but not in
    // the closed form, and the gap must stop falling — it plateaus at the broken term's own
    // divergence. A test that cannot fail proves nothing.
    OpqParameters<double> broken = opq_parameters<double>(0.21 * kDegPerRad);
    broken.b[7] += 1.0e-3;
    const OpqParameters<double> intact = opq_parameters<double>(0.21 * kDegPerRad);
    double gap3 = 0.0;
    double gap4 = 0.0;
    for_each_box_point([&](double x, double y, double z, double) {
        const double a = analytic_divergence(intact, x, y, z);
        gap3 = std::max(gap3, std::fabs(stencil_divergence(broken, x, y, z, 1e-3) - a));
        gap4 = std::max(gap4, std::fabs(stencil_divergence(broken, x, y, z, 1e-4) - a));
    });
    EXPECT_LT(gap3 / gap4, 5.0) << "a broken coefficient must be detected as a floor";
    EXPECT_GT(gap4, 1e-4);
}

TEST(IrbemOpq, TheThreeRadialRulesArePublishedOnes) {
    const OpqParameters<double> p = opq_parameters<double>(12.0);
    // Inside 2 Re: identically zero (report p. 66, "set to zero below 2 earth radii").
    for (const Position<Frame::GSM>& q : scatter(32, 1.0, 1.999)) {
        const std::array<double, 3> b = opq_components<double>(p, q.v[0], q.v[1], q.v[2]);
        EXPECT_EQ(b[0], 0.0);
        EXPECT_EQ(b[1], 0.0);
        EXPECT_EQ(b[2], 0.0);
    }
    // The taper: linear in r^2 from 0 at r = 2 to 1 at r = 2.5. Along one ray the untapered series
    // varies slowly, so the ratio of the field at r = 2.25 to the field just past 2.5 is the taper
    // factor (2.25^2 - 4) / 2.25 = 0.4722 to within the series' own few-per-cent variation, and the
    // field is continuous through both ends.
    const double ux = 0.36;
    const double uy = 0.48;
    const double uz = 0.80;
    const auto along = [&](double r) { return opq_components<double>(p, r * ux, r * uy, r * uz); };
    const auto mag = [](const std::array<double, 3>& b) {
        return std::sqrt((b[0] * b[0]) + (b[1] * b[1]) + (b[2] * b[2]));
    };
    const double taper_25 = ((2.25 * 2.25) - 4.0) / 2.25;
    EXPECT_NEAR(mag(along(2.25)) / mag(along(2.5)), taper_25, 0.05);
    EXPECT_LT(mag(along(2.0 + 1e-6)), 1e-4);                    // continuous into the inner zero
    EXPECT_NEAR(mag(along(2.5 - 1e-9)), mag(along(2.5 + 1e-9)), 1e-6);  // continuous at 2.5
    EXPECT_GT(mag(along(2.5)), 10.0);                            // and not small: it is a field
    // The template: exactly zero beyond 15 Re, a field just inside.
    const std::array<double, 3> outside = along(15.0 + 1e-9);
    EXPECT_EQ(outside[0], 0.0);
    EXPECT_EQ(outside[1], 0.0);
    EXPECT_EQ(outside[2], 0.0);
    EXPECT_GT(mag(along(15.0 - 1e-6)), 1.0);
    // A NaN radius takes the same zero exit, so no NaN can propagate out of the series.
    const std::array<double, 3> nan = opq_components<double>(p, std::nan(""), 1.0, 1.0);
    EXPECT_EQ(nan[0], 0.0);
    EXPECT_EQ(nan[2], 0.0);
    // The float lane applies the same three rules.
    const OpqParameters<float> pf = opq_parameters<float>(12.0);
    EXPECT_EQ(opq_components<float>(pf, 1.5F, 0.0F, 0.0F)[2], 0.0F);
    EXPECT_EQ(opq_components<float>(pf, 16.0F, 0.0F, 0.0F)[2], 0.0F);
    EXPECT_NE(opq_components<float>(pf, 6.0F, 0.0F, 0.0F)[2], 0.0F);
}

TEST(IrbemOpq, TiltIsAContinuousParameter) {
    // Unlike T89, the tilt here is a real polynomial parameter, so nearby tilts give nearby fields:
    // the finite-difference slope in the tilt is bounded and its second difference is small. A
    // step of 0.5 degrees is a few minutes of UT.
    double worst_slope = 0.0;
    double worst_curvature = 0.0;
    for (const Position<Frame::GSM>& q : scatter(24)) {
        for (int k = -70; k < 70; ++k) {
            const double t0 = 0.5 * k;
            const OpqParameters<double> lo = opq_parameters<double>(t0 - 0.5);
            const OpqParameters<double> mid = opq_parameters<double>(t0);
            const OpqParameters<double> hi = opq_parameters<double>(t0 + 0.5);
            for (int c = 0; c < 3; ++c) {
                const double a = opq_components<double>(lo, q.v[0], q.v[1], q.v[2])[static_cast<std::size_t>(c)];
                const double b = opq_components<double>(mid, q.v[0], q.v[1], q.v[2])[static_cast<std::size_t>(c)];
                const double d = opq_components<double>(hi, q.v[0], q.v[1], q.v[2])[static_cast<std::size_t>(c)];
                worst_slope = std::max(worst_slope, std::fabs(d - a));
                worst_curvature = std::max(worst_curvature, std::fabs((d - (2.0 * b)) + a));
            }
        }
    }
    std::printf("[ MEASURED ] tilt sweep -35..35 deg: worst |dB| per degree %.3f nT, worst second "
                "difference over 0.5 deg %.4f nT\n",
                worst_slope, worst_curvature);
    EXPECT_LT(worst_slope, 5.0);
    EXPECT_LT(worst_curvature, 0.05);
    EXPECT_GT(worst_slope, 0.1) << "the tilt must actually move the field";
}

// ================================================================================================
// The scalar entry points and their verdicts
// ================================================================================================

TEST(IrbemOpq, ReferenceLaneMatchesTheComponentForm) {
    // The GSM entry point is the SM series inside one rotation about y: rotate by hand, evaluate,
    // rotate back, and the two must agree to roundoff; the folding overload agrees exactly.
    const double t = 0.37;
    const OpqParameters<double> p = opq_parameters<double>(t * kDegPerRad);
    for (const Position<Frame::GSM>& q : scatter(64)) {
        const double s = std::sin(t);
        const double c = std::cos(t);
        const std::array<double, 3> sm =
            opq_components<double>(p, (q.v[0] * c) - (q.v[2] * s), q.v[1], (q.v[0] * s) + (q.v[2] * c));
        const ir::FieldVector<Frame::GSM> b = opq_field_at(q, p, s, c);
        EXPECT_NEAR(b.v[0], (sm[0] * c) + (sm[2] * s), 1e-12);
        EXPECT_EQ(b.v[1], sm[1]);
        EXPECT_NEAR(b.v[2], (-sm[0] * s) + (sm[2] * c), 1e-12);
        const ir::FieldVector<Frame::GSM> folded = opq_field_at(q, t);
        EXPECT_EQ(folded.v[0], b.v[0]);
        EXPECT_EQ(folded.v[1], b.v[1]);
        EXPECT_EQ(folded.v[2], b.v[2]);
    }
    // At zero tilt GSM and SM coincide and the entry point IS the series.
    const std::array<double, 3> sm = opq_components<double>(opq_parameters<double>(0.0), -7.0, 2.0, 1.5);
    const ir::FieldVector<Frame::GSM> b = opq_field_at(at(-7.0, 2.0, 1.5), 0.0);
    EXPECT_EQ(b.v[0], sm[0]);
    EXPECT_EQ(b.v[2], sm[2]);
}

TEST(IrbemOpq, ValidityIsReportedFromBothSides) {
    // Just inside 15 Re: Ok with a real field. Just outside: OutOfValidityRange, and the value is
    // the published template's ZERO — returned, not suppressed, and stated as what it is.
    const ir::Result<ir::FieldVector<Frame::GSM>> in = opq_field(at(-14.99, 0.0, 0.0), 0.2);
    EXPECT_EQ(in.status, Status::Ok);
    EXPECT_GT(in.value.magnitude(), 0.5);
    const ir::Result<ir::FieldVector<Frame::GSM>> out = opq_field(at(-15.01, 0.0, 0.0), 0.2);
    EXPECT_EQ(out.status, Status::OutOfValidityRange);
    EXPECT_EQ(out.value.v[0], 0.0);
    EXPECT_EQ(out.value.v[2], 0.0);
    // Exactly on the bound is inside (the envelope's intervals are closed).
    EXPECT_EQ(opq_field(at(15.0, 0.0, 0.0), 0.2).status, Status::Ok);
    // Below 2 Re is inside the envelope and the published answer there is zero: Ok, zero.
    const ir::Result<ir::FieldVector<Frame::GSM>> low = opq_field(at(1.5, 0.5, 0.5), 0.2);
    EXPECT_EQ(low.status, Status::Ok);
    EXPECT_EQ(low.value.v[2], 0.0);
    // Inside the Earth is not a point at all.
    EXPECT_EQ(opq_field(at(0.5, 0.0, 0.5), 0.2).status, Status::DomainError);
    EXPECT_EQ(ir::opq_status(at(0.99, 0.0, 0.0), 0.0), Status::DomainError);
    EXPECT_EQ(ir::opq_status(at(1.0, 0.0, 0.0), 0.0), Status::Ok);
}

TEST(IrbemOpq, NonFiniteInputIsADomainError) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const Position<Frame::GSM>& p : {at(nan, 0.0, 0.0), at(0.0, inf, 0.0), at(0.0, 0.0, -inf)}) {
        const ir::Result<ir::FieldVector<Frame::GSM>> r = opq_field(p, 0.1);
        EXPECT_EQ(r.status, Status::DomainError);
        EXPECT_EQ(r.value.v[0], 0.0);
        EXPECT_EQ(r.value.v[1], 0.0);
        EXPECT_EQ(r.value.v[2], 0.0);
    }
    EXPECT_EQ(opq_field(at(5.0, 0.0, 0.0), nan).status, Status::DomainError);
    EXPECT_EQ(opq_field(at(5.0, 0.0, 0.0), inf).status, Status::DomainError);
}

TEST(IrbemOpq, ATiltBeyondARightAngleIsADomainError) {
    // The model has no tan(psi), so a right angle itself is a defined point; beyond it the "tilt"
    // is not a dipole tilt, which is the bound FieldContext already guarantees.
    const double half_pi = std::numbers::pi / 2.0;
    EXPECT_EQ(opq_field(at(5.0, 0.0, 0.0), half_pi).status, Status::Ok);
    EXPECT_EQ(opq_field(at(5.0, 0.0, 0.0), -half_pi).status, Status::Ok);
    EXPECT_EQ(opq_field(at(5.0, 0.0, 0.0), half_pi + 1e-6).status, Status::DomainError);
    EXPECT_EQ(opq_field(at(5.0, 0.0, 0.0), -2.0).status, Status::DomainError);
    std::array<ir::FieldVector<Frame::GSM>, 1> out{};
    const std::array<Position<Frame::GSM>, 1> pts{at(5.0, 0.0, 0.0)};
    EXPECT_EQ(opq_field_batch(pts, 2.0, out).status, Status::DomainError);
    EXPECT_EQ(opq_field_batch(pts, std::nan(""), out).status, Status::DomainError);
}

TEST(IrbemOpq, ContextOverloadAgreesWithTheExplicitOne) {
    const ir::Epoch epoch{2015.5, 43200.0, 2015, 180};
    ir::RotationTable identity{};
    for (cheatah::fixarray::mat3d& m : identity) m = cheatah::fixarray::mat3d::identity();
    ir::DriverSet drivers{};
    drivers[static_cast<std::size_t>(ir::Driver::Kp)] = 47.0;
    const ir::ContextResult built = ir::make_field_context(epoch, -0.42, identity, drivers);
    ASSERT_TRUE(built.has_value()) << ir::describe(built.error());

    const Position<Frame::GSM> p = at(3.75, -1.5, 2.25);
    const ir::Result<ir::FieldVector<Frame::GSM>> via_context = opq_field(p, built.value());
    const ir::Result<ir::FieldVector<Frame::GSM>> via_scalar = opq_field(p, -0.42);
    EXPECT_EQ(via_context.status, via_scalar.status);
    EXPECT_EQ(via_context.value.v[0], via_scalar.value.v[0]);
    EXPECT_EQ(via_context.value.v[1], via_scalar.value.v[1]);
    EXPECT_EQ(via_context.value.v[2], via_scalar.value.v[2]);
}

TEST(IrbemOpq, TheModelReadsNoDrivers) {
    // A quiet-time model: the four corpus regimes, quiet to extreme, must give BIT-identical
    // fields, because nothing in the driver vector is read. This is the regime sweep for a model
    // with no drivers — and the envelope table says so too.
    EXPECT_EQ(ir::envelope_of(ir::ExternalModel::OlsonPfitzerQuiet1977).bound_count, 0U);
    EXPECT_EQ(ir::envelope_of(ir::ExternalModel::OlsonPfitzerQuiet1977).max_r_geo, 15.0);
    const ir::Epoch epoch{2015.5, 43200.0, 2015, 180};
    ir::RotationTable identity{};
    for (cheatah::fixarray::mat3d& m : identity) m = cheatah::fixarray::mat3d::identity();
    std::vector<ir::Result<ir::FieldVector<Frame::GSM>>> answers;
    for (const corpus::MagInput& m : corpus::regime_drivers) {
        ir::DriverSet drivers{};
        drivers[static_cast<std::size_t>(ir::Driver::Kp)] = m.kp * 10.0;
        drivers[static_cast<std::size_t>(ir::Driver::Dst)] = m.dst;
        drivers[static_cast<std::size_t>(ir::Driver::Pdyn)] = m.pdyn;
        drivers[static_cast<std::size_t>(ir::Driver::ByIMF)] = m.by_imf;
        drivers[static_cast<std::size_t>(ir::Driver::BzIMF)] = m.bz_imf;
        const ir::ContextResult built = ir::make_field_context(epoch, 0.3, identity, drivers);
        ASSERT_TRUE(built.has_value()) << ir::describe(built.error());
        answers.push_back(opq_field(at(6.6, 1.0, -0.5), built.value()));
    }
    ASSERT_EQ(answers.size(), corpus::regime_drivers.size());
    for (const ir::Result<ir::FieldVector<Frame::GSM>>& r : answers) {
        EXPECT_EQ(r.status, Status::Ok);
        EXPECT_EQ(r.value.v[0], answers[0].value.v[0]);
        EXPECT_EQ(r.value.v[1], answers[0].value.v[1]);
        EXPECT_EQ(r.value.v[2], answers[0].value.v[2]);
    }
    EXPECT_NE(answers[0].value.v[2], 0.0);
}

// ================================================================================================
// The batch lanes
// ================================================================================================

TEST(IrbemOpq, HostFloatLaneTracksTheReferenceLane) {
    const std::size_t n = 4096;
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<float> pos(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    const double t = 0.31;
    const float sp = static_cast<float>(std::sin(t));
    const float cp = static_cast<float>(std::cos(t));
    std::vector<float> out(3 * n);
    ASSERT_TRUE(opq_field_host(pos, out, sp, cp, opq_parameters<float>(t * kDegPerRad)));
    const OpqParameters<double> p = opq_parameters<double>(t * kDegPerRad);
    double worst = 0.0;
    double biggest = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref =
            opq_field_at(at(pos[3 * i], pos[(3 * i) + 1], pos[(3 * i) + 2]), p, std::sin(t), std::cos(t));
        for (std::size_t c = 0; c < 3; ++c) {
            worst = std::max(worst, std::fabs(static_cast<double>(out[(3 * i) + c]) - ref.v[c]));
            biggest = std::max(biggest, std::fabs(ref.v[c]));
        }
    }
    std::printf("[ MEASURED ] fp32 host lane vs fp64 reference, %zu points: max |dB| = %.3e nT of "
                "%.1f nT\n",
                n, worst, biggest);
    // Fifth-order monomials at 15 Re reach 1e6 before their coefficients tame them, so fp32
    // carries ~1e-3 nT of cancellation error at the outer edge — three orders inside the field.
    EXPECT_LT(worst, 5e-3);
}

TEST(IrbemOpq, HostFloatLaneRejectsMismatchedSpans) {
    const OpqParameters<float> p = opq_parameters<float>(5.0);
    std::vector<float> pos(7, 1.0F);
    std::vector<float> out(7, 0.0F);
    EXPECT_FALSE(opq_field_host(pos, out, 0.1F, 0.99F, p));
    std::vector<float> pos6(6, 3.0F);
    std::vector<float> out3(3, 0.0F);
    EXPECT_FALSE(opq_field_host(pos6, out3, 0.1F, 0.99F, p));
    std::vector<float> out6(6, 0.0F);
    EXPECT_TRUE(opq_field_host(pos6, out6, 0.1F, 0.99F, p));
}

TEST(IrbemOpq, ParameterBlockCarriesTheTiltThenTheCoefficients) {
    const std::array<float, opq_param_count> block = opq_param_block(0.25F, 0.968F, 14.5);
    const OpqParameters<float> p = opq_parameters<float>(14.5);
    EXPECT_EQ(block[0], 0.25F);
    EXPECT_EQ(block[1], 0.968F);
    for (std::size_t k = 0; k < opq_xz_count; ++k) {
        EXPECT_EQ(block[2 + k], p.a[k]);
        EXPECT_EQ(block[34 + k], p.b[k]);
        EXPECT_EQ(block[110 + k], p.e[k]);
        EXPECT_EQ(block[142 + k], p.f[k]);
    }
    for (std::size_t k = 0; k < opq_y_count; ++k) {
        EXPECT_EQ(block[66 + k], p.c[k]);
        EXPECT_EQ(block[88 + k], p.d[k]);
    }
#if CHEATAH_SPACE_IRBEM_OPQ_GPU
    EXPECT_EQ(ir::gpu::kernel_info("irbem_opq_f32").params, opq_param_count);
    EXPECT_EQ(ir::gpu::kernel_info("irbem_opq_f32").bindings, 4U);
#endif
}

TEST(IrbemOpq, BatchAgreesWithTheReferenceLane) {
    // Below the crossover the batch is the host lane and must be BIT-identical to the reference.
    const std::vector<Position<Frame::GSM>> pts = scatter(256);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    const ir::Result<bool> r = opq_field_batch(pts, 0.28, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_FALSE(r.value);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const ir::FieldVector<Frame::GSM> ref = opq_field_at(pts[i], 0.28);
        EXPECT_EQ(out[i].v[0], ref.v[0]) << i;
        EXPECT_EQ(out[i].v[1], ref.v[1]) << i;
        EXPECT_EQ(out[i].v[2], ref.v[2]) << i;
    }
    // An empty batch is Ok and touches nothing.
    const std::vector<Position<Frame::GSM>> none;
    std::vector<ir::FieldVector<Frame::GSM>> nothing;
    EXPECT_EQ(opq_field_batch(none, 0.28, nothing).status, Status::Ok);
}

TEST(IrbemOpq, BatchRejectsMismatchedSpans) {
    const std::vector<Position<Frame::GSM>> pts = scatter(8);
    std::vector<ir::FieldVector<Frame::GSM>> out(7);
    const ir::Result<bool> r = opq_field_batch(pts, 0.1, out);
    EXPECT_EQ(r.status, Status::DomainError);
    EXPECT_FALSE(r.value);
}

TEST(IrbemOpq, BatchReportsTheSameEnvelopeTheScalarLaneDoes) {
    // One point beyond 15 Re makes the batch a caveat, still computed in full; one point inside
    // the Earth makes it a refusal, zeroed in full; a NaN anywhere is a refusal too.
    std::vector<Position<Frame::GSM>> pts = scatter(64);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    pts[10] = at(-16.0, 0.0, 0.0);
    const ir::Result<bool> far = opq_field_batch(pts, 0.1, out);
    EXPECT_EQ(far.status, Status::OutOfValidityRange);
    EXPECT_EQ(out[10].v[2], 0.0) << "the template's zero, computed and returned";
    EXPECT_NE(out[11].v[2], 0.0) << "every other point is still computed";
    for (std::size_t i = 0; i < pts.size(); ++i) {
        EXPECT_EQ(opq_field(pts[i], 0.1).status, i == 10 ? Status::OutOfValidityRange : Status::Ok);
    }

    pts[10] = at(0.3, 0.0, 0.0);
    const ir::Result<bool> inside = opq_field_batch(pts, 0.1, out);
    EXPECT_EQ(inside.status, Status::DomainError);
    for (const ir::FieldVector<Frame::GSM>& b : out) {
        EXPECT_EQ(b.v[0], 0.0);
        EXPECT_EQ(b.v[1], 0.0);
        EXPECT_EQ(b.v[2], 0.0);
    }

    pts[10] = at(std::nan(""), 0.0, 0.0);
    EXPECT_EQ(opq_field_batch(pts, 0.1, out).status, Status::DomainError);

    // The fold itself, in isolation: it folds squares and reports through check_position.
    ir::OpqPositionFold fold;
    EXPECT_EQ(fold.verdict(), Status::DomainError) << "no points folded: r2_lo is +inf";
    fold.add(at(3.0, 4.0, 0.0));
    EXPECT_EQ(fold.r2_lo, 25.0);
    EXPECT_EQ(fold.r2_hi, 25.0);
    EXPECT_EQ(fold.verdict(), Status::Ok);
    fold.add(at(0.0, 0.0, 15.5));
    EXPECT_EQ(fold.verdict(), Status::OutOfValidityRange);
}

TEST(IrbemOpq, NothingOnTheHeapInTheHotPath) {
    // The claim is "no heap in a hot path", so it gets a counter. The SECOND call is the one that
    // matters: a routine that allocates a workspace per invocation passes a single-call check.
    const std::vector<Position<Frame::GSM>> pts = scatter(256);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    (void)opq_field_batch(pts, 0.2, out);
    (void)opq_field(at(5.0, 1.0, 1.0), 0.2);
    const OpqParameters<double> p = opq_parameters<double>(11.5);

    const std::size_t before = cheatah_space_test::allocation_count();
    for (int i = 0; i < 64; ++i) {
        sink = sink + opq_field(at(5.0 + (0.01 * i), 1.0, 1.0), 0.2).value.v[2];
        sink = sink + opq_field_at(at(4.0, 1.0 + (0.01 * i), 1.0), p, 0.2, 0.98).v[0];
        sink = sink + opq_components<double>(p, 3.0, 2.0, 1.0 + (0.01 * i))[1];
    }
    (void)opq_field_batch(pts, 0.2, out);   // below the crossover: the host lane
    sink = sink + out[0].v[2];
    const std::array<float, opq_param_count> block = opq_param_block(0.2F, 0.98F, 11.5);
    sink = sink + block[3];
    EXPECT_EQ(before, cheatah_space_test::allocation_count());
}

// ================================================================================================
// The total field
// ================================================================================================

namespace {
ir::Rotations epoch_rotations(const ir::Igrf<10>& m) {
    const ir::Result<ir::Rotations> r = ir::api::rotations_at(2015, 180, 43200.0, m);
    EXPECT_EQ(r.status, Status::Ok);
    return r.value;
}
}  // namespace

TEST(IrbemOpq, TotalFieldSuperposesInternalAndExternal) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.5).value();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::TotalFieldOpq<10> total(igrf, rot);
    static_assert(ir::GeoFieldModel<ir::TotalFieldOpq<10>>,
                  "a tracer must be able to follow it without knowing it is a sum");
    static_assert(ir::TotalFieldOpq<10>::degree == 10);

    const ir::Position<Frame::GEO> p{cheatah::fixarray::vec3d{6.0, 0.0, 0.0}};
    const double internal_only = igrf.evaluate(p).magnitude();
    const double with_external = total.evaluate(p).magnitude();
    EXPECT_EQ(igrf.g(1, 0), total.g(1, 0));
    EXPECT_EQ(igrf.h(2, 1), total.h(2, 1));
    EXPECT_EQ(&total.internal(), &igrf);
    EXPECT_EQ(&total.rotations(), &rot);
    EXPECT_EQ(total.parameters().e[0], opq_parameters<double>(rot.dipole_tilt_deg).e[0]);
    // A real contribution at L = 6 — a few per cent — and a perturbation, not a replacement.
    EXPECT_GT(std::fabs(with_external - internal_only) / internal_only, 1e-3);
    EXPECT_LT(std::fabs(with_external - internal_only) / internal_only, 0.5);

    // The superposition, by hand: rotate the external field into GEO and add it.
    const ir::Position<Frame::GSM> p_gsm = ir::transform<Frame::GSM>(p, rot);
    const ir::FieldVector<Frame::GSM> ext =
        opq_field_at(p_gsm, rot.dipole_tilt_deg / kDegPerRad);
    const ir::FieldVector<Frame::GEO> ext_geo = ir::transform<Frame::GEO>(ext, rot);
    const ir::FieldVector<Frame::GEO> sum = total.evaluate(p);
    for (std::size_t c = 0; c < 3; ++c) {
        EXPECT_NEAR(sum.v[c], igrf.evaluate(p).v[c] + ext_geo.v[c], 1e-9);
    }
    // Inside the Earth the external model declines and the internal field is returned alone.
    const ir::Position<Frame::GEO> deep{cheatah::fixarray::vec3d{0.5, 0.0, 0.0}};
    EXPECT_EQ(total.evaluate(deep).v[0], igrf.evaluate(deep).v[0]);
}

TEST(IrbemOpq, TotalFieldReportsWhenTheExternalModelDeclines) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.5).value();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::TotalFieldOpq<10> total(igrf, rot);
    EXPECT_EQ(total.external_status(ir::Position<Frame::GEO>{cheatah::fixarray::vec3d{6.0, 0.0, 0.0}}),
              Status::Ok);
    EXPECT_EQ(total.external_status(ir::Position<Frame::GEO>{cheatah::fixarray::vec3d{16.0, 0.0, 0.0}}),
              Status::OutOfValidityRange);
    EXPECT_EQ(total.external_status(ir::Position<Frame::GEO>{cheatah::fixarray::vec3d{0.5, 0.0, 0.0}}),
              Status::DomainError);
    // Beyond 15 Re the total IS the internal field: the published external part is zero there.
    const ir::Position<Frame::GEO> far{cheatah::fixarray::vec3d{-16.0, 0.0, 0.0}};
    EXPECT_EQ(total.evaluate(far).v[2], igrf.evaluate(far).v[2]);
}

TEST(IrbemOpq, LstarRunsThroughTheTotalField) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.5).value();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::TotalFieldOpq<10> total(igrf, rot);
    const ir::Position<Frame::GEO> p{cheatah::fixarray::vec3d{6.0, 0.0, 0.0}};
    const ir::Result<ir::FieldLine> line = ir::trace_invariant(total, p, 45.0);
    ASSERT_EQ(line.status, Status::Ok);
    EXPECT_GT(line.value.invariant_i, 0.0);
    // The physics the model exists for, at 12 UT: GEO +x is local noon, where the magnetopause
    // current COMPRESSES the field, and GEO -x is midnight, where the tail current STRETCHES it.
    // B_min through the total field is therefore above the internal one on the dayside and below
    // it on the nightside — the signature no dipole-only trace can show.
    const ir::Result<ir::FieldLine> internal = ir::trace_invariant(igrf, p, 45.0);
    ASSERT_EQ(internal.status, Status::Ok);
    EXPECT_GT(line.value.b_min, internal.value.b_min) << "noon: compressed";
    const ir::Position<Frame::GEO> night{cheatah::fixarray::vec3d{-6.0, 0.0, 0.0}};
    const ir::Result<ir::FieldLine> night_total = ir::trace_invariant(total, night, 45.0);
    const ir::Result<ir::FieldLine> night_internal = ir::trace_invariant(igrf, night, 45.0);
    ASSERT_EQ(night_total.status, Status::Ok);
    ASSERT_EQ(night_internal.status, Status::Ok);
    EXPECT_LT(night_total.value.b_min, night_internal.value.b_min) << "midnight: stretched";

    const ir::Result<ir::DriftShell> shell = ir::make_lstar(total, rot, p, 90.0);
    ASSERT_EQ(shell.status, Status::Ok) << ir::describe(shell.status);
    // The oracle says 5.4482 here at matched default resolution (tools/oracle/opq_diff.cpp);
    // a loose band, since this test's job is that the chain runs, not the differential below.
    EXPECT_GT(shell.value.lstar, 5.3);
    EXPECT_LT(shell.value.lstar, 5.6);
}

TEST(IrbemOpq, TotalFieldBatchTraceIsTheHostLaneAndSaysSo) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.5).value();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::TotalFieldOpq<10> total(igrf, rot);
    const std::array<ir::Position<Frame::GEO>, 3> starts{
        ir::Position<Frame::GEO>{cheatah::fixarray::vec3d{4.0, 0.0, 0.0}},
        ir::Position<Frame::GEO>{cheatah::fixarray::vec3d{0.0, 5.5, 0.5}},
        ir::Position<Frame::GEO>{cheatah::fixarray::vec3d{-6.0, 1.0, 0.0}}};
    const std::array<double, 3> pitch{30.0, 60.0, 90.0};
    std::array<ir::FieldLine, 3> lines{};
    std::array<Status, 3> st{};
    const ir::Result<bool> r = ir::trace_invariant_batch(total, starts, pitch, lines, st);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_FALSE(r.value) << "no device tracer composes OP-77; the answer must say host";
    for (std::size_t i = 0; i < starts.size(); ++i) {
        const ir::Result<ir::FieldLine> one = ir::trace_invariant(total, starts[i], pitch[i]);
        EXPECT_EQ(st[i], one.status);
        EXPECT_EQ(lines[i].invariant_i, one.value.invariant_i);
        EXPECT_EQ(lines[i].b_min, one.value.b_min);
    }
    std::array<double, 2> short_pitch{30.0, 60.0};
    EXPECT_EQ(ir::trace_invariant_batch(total, starts, short_pitch, lines, st).status,
              Status::DomainError);
    // A step cap no line can close under makes the batch report OpenFieldLine, per line and in
    // the aggregate — the same verdict the scalar tracer gives, never a silently truncated I.
    ir::TraceOptions starved;
    starved.max_steps = 2;
    const ir::Result<bool> open = ir::trace_invariant_batch(total, starts, pitch, lines, st, starved);
    EXPECT_EQ(open.status, Status::OpenFieldLine);
    EXPECT_FALSE(open.value);
    for (std::size_t i = 0; i < starts.size(); ++i) {
        EXPECT_EQ(st[i], ir::trace_invariant(total, starts[i], pitch[i], starved).status);
    }
}

// ================================================================================================
// The differentials against the IRBEM oracle
// ================================================================================================

TEST(IrbemOpq, HeavyDifferentialAgreesWithTheIrbemOracle) {
    const Oracle& o = oracle();
    if (!o.usable()) {
        GTEST_SKIP() << "IRBEM oracle not present (set CHEATAH_SPACE_IRBEM_ORACLE to its .so); "
                        "the oracle is a dev-only black box and is never linked";
    }
    // Three tilts, 300 scattered points each from the inner zero through the taper to the edge of
    // the published region. The budget is 1e-6 relative on B; the measurement is fp64 roundoff
    // (tools/oracle/opq_diff.cpp: RMS 3e-14 nT), so the cap here sits three orders inside the
    // budget and a single transcription error in any table would fail it.
    const corpus::MagInput quiet = corpus::regime_drivers[0];
    double worst_abs = 0.0;
    double worst_rel = 0.0;
    std::size_t n = 0;
    for (const Epoch& e : kEpochs) {
        const double ps = o.tilt(e.year, e.doy, e.ut);
        const OpqParameters<double> p = opq_parameters<double>(ps * kDegPerRad);
        for (const Position<Frame::GSM>& q : scatter(300, 1.2, 14.9)) {
            std::array<double, 3> ora{};
            ASSERT_TRUE(o.external(e.year, e.doy, e.ut, q, quiet, ora));
            const ir::FieldVector<Frame::GSM> mine = opq_field_at(q, p, std::sin(ps), std::cos(ps));
            double d2 = 0.0;
            double o2 = 0.0;
            for (std::size_t c = 0; c < 3; ++c) {
                d2 += (mine.v[c] - ora[c]) * (mine.v[c] - ora[c]);
                o2 += ora[c] * ora[c];
            }
            worst_abs = std::max(worst_abs, std::sqrt(d2));
            if (o2 > 1e-6) worst_rel = std::max(worst_rel, std::sqrt(d2 / o2));
            ++n;
        }
    }
    std::printf("[ MEASURED ] vs IRBEM kext=5, %zu points, three tilts: max |dB| = %.3e nT, max rel = %.3e\n",
                n, worst_abs, worst_rel);
    EXPECT_LT(worst_rel, 1e-9);
    EXPECT_LT(worst_abs, 1e-9);
}

TEST(IrbemOpq, HeavyDifferentialOracleIgnoresEveryDriverRegime) {
    const Oracle& o = oracle();
    if (!o.usable()) GTEST_SKIP() << "IRBEM oracle not present";
    // The regime sweep for a model with no drivers: quiet, moderate, storm and extreme must give
    // the ORACLE the same field to the last bit — and ours too — or "no drivers" would be a claim
    // about this implementation rather than about the model.
    const Epoch e = kEpochs[1];
    const double ps = o.tilt(e.year, e.doy, e.ut);
    for (const Position<Frame::GSM>& q : scatter(24, 2.6, 14.5)) {
        std::array<double, 3> base{};
        ASSERT_TRUE(o.external(e.year, e.doy, e.ut, q, corpus::regime_drivers[0], base));
        for (std::size_t g = 1; g < corpus::regime_drivers.size(); ++g) {
            std::array<double, 3> again{};
            ASSERT_TRUE(o.external(e.year, e.doy, e.ut, q, corpus::regime_drivers[g], again));
            EXPECT_EQ(again[0], base[0]) << "regime " << g;
            EXPECT_EQ(again[1], base[1]) << "regime " << g;
            EXPECT_EQ(again[2], base[2]) << "regime " << g;
        }
        for (const corpus::StormEvent& s : corpus::storm_events) {
            std::array<double, 3> storm{};
            ASSERT_TRUE(o.external(e.year, e.doy, e.ut, q, s.mag, storm));
            EXPECT_EQ(storm[2], base[2]) << s.name;
        }
        const ir::FieldVector<Frame::GSM> mine = opq_field_at(q, ps);
        EXPECT_NEAR(mine.v[2], base[2], 1e-9);
    }
    // And beyond 15 Re the oracle refuses outright where the published template says zero — the
    // one place the two disagree, by design on both sides, and stated.
    std::array<double, 3> far{};
    EXPECT_FALSE(o.external(e.year, e.doy, e.ut, at(-16.0, 0.0, 0.0), corpus::regime_drivers[0], far));
    EXPECT_EQ(opq_field(at(-16.0, 0.0, 0.0), ps).status, Status::OutOfValidityRange);
}

TEST(IrbemOpq, HeavyDifferentialLstarMatchesTheOracleThroughTheTotalField) {
    const Oracle& o = oracle();
    if (!o.usable()) GTEST_SKIP() << "IRBEM oracle not present";
    // The drift-shell chain through IGRF + OP-77 against IRBEM's make_lstar at kext = 5 — the model
    // its LANDI2LSTAR path hard-wires — at matched default resolution (options(3,4) = 0). IRBEM's
    // own default-resolution discretization error is 0.010-0.017 at L ~ 6 (docs/ERROR_BUDGET.md),
    // and at matched resolution 9 the harness measures 0.0062, the same figure the internal field
    // alone gives: the two agree to the numerics.
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.5).value();
    const ir::Rotations rot = epoch_rotations(igrf);
    const ir::TotalFieldOpq<10> total(igrf, rot);
    double worst = 0.0;
    for (double x : {3.0, 4.5, 6.0}) {
        int ntime = 1;
        int kext = 5;
        int sysaxes = 1;
        std::array<int, 5> options{1, 0, 0, 0, 0};
        std::array<int, 1> iyear{2015};
        std::array<int, 1> idoy{180};
        std::array<double, 1> ut{43200.0};
        std::array<double, 1> x1{x};
        std::array<double, 1> x2{0.0};
        std::array<double, 1> x3{0.0};
        std::vector<double> mag(25, 0.0);
        std::array<double, 1> lm{};
        std::array<double, 1> lstar{};
        std::array<double, 1> blocal{};
        std::array<double, 1> bmin{};
        std::array<double, 1> xj{};
        std::array<double, 1> mlt{};
        o.make_lstar(&ntime, &kext, options.data(), &sysaxes, iyear.data(), idoy.data(), ut.data(),
                     x1.data(), x2.data(), x3.data(), mag.data(), lm.data(), lstar.data(),
                     blocal.data(), bmin.data(), xj.data(), mlt.data());
        ASSERT_GT(lstar[0], 0.0);
        const ir::Result<ir::DriftShell> ours =
            ir::make_lstar(total, rot, ir::Position<Frame::GEO>{cheatah::fixarray::vec3d{x, 0.0, 0.0}}, 90.0);
        ASSERT_EQ(ours.status, Status::Ok);
        worst = std::max(worst, std::fabs(ours.value.lstar - lstar[0]));
        std::printf("[ MEASURED ] L* at (%.1f, 0, 0): oracle kext=5 %.5f, ours %.5f\n", x, lstar[0],
                    ours.value.lstar);
    }
    EXPECT_LT(worst, 0.02);
}

// ================================================================================================
// The device lane
// ================================================================================================

#if CHEATAH_SPACE_IRBEM_OPQ_GPU
namespace {

/// Point `CHEATAH_SPACE_IRBEM_SPV_DIR` somewhere for the life of the object, and put it back —
/// value or absence — afterwards; the seam resolves the shader directory on every launch.
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

TEST(IrbemOpq, BatchFallsBackToTheHostWhenTheShaderWasNeverBuilt) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_opq_f32");
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    {
        const SpvDirScope nowhere(std::filesystem::temp_directory_path().string() +
                                  "/cheatah-space-no-such-shader-dir");
        const ir::Result<bool> r = opq_field_batch(pts, 0.28, out);
        EXPECT_EQ(r.status, Status::Ok);
        EXPECT_FALSE(r.value) << "with no compiled shader the batch must run on the host";
    }
    const OpqParameters<double> p = opq_parameters<double>(0.28 * kDegPerRad);
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = opq_field_at(pts[i], p, std::sin(0.28), std::cos(0.28));
        ASSERT_EQ(out[i].v[0], ref.v[0]) << "point " << i;
        ASSERT_EQ(out[i].v[1], ref.v[1]) << "point " << i;
        ASSERT_EQ(out[i].v[2], ref.v[2]) << "point " << i;
    }
}

TEST(IrbemOpq, BatchUsesTheDeviceWhenOneIsAvailable) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_opq_f32");
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    const ir::Result<bool> r = opq_field_batch(pts, 0.28, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_TRUE(r.value) << "the batch fell back to the host with a device present";
    const OpqParameters<double> p = opq_parameters<double>(0.28 * kDegPerRad);
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = opq_field_at(pts[i], p, std::sin(0.28), std::cos(0.28));
        for (std::size_t c = 0; c < 3; ++c) worst = std::max(worst, std::fabs(out[i].v[c] - ref.v[c]));
    }
    std::printf("[ MEASURED ] device batch of %zu vs fp64 reference: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-2);
    // Below the crossover the same call runs on the host, bit-identical to the reference.
    const std::vector<Position<Frame::GSM>> few = scatter(16);
    std::vector<ir::FieldVector<Frame::GSM>> few_out(few.size());
    EXPECT_FALSE(opq_field_batch(few, 0.28, few_out).value);
}

TEST(IrbemOpq, TheDeviceLaneRefusesABadPointBeforeItDispatches) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_opq_f32");
    std::vector<Position<Frame::GSM>> pts = scatter(n);
    pts[n / 3] = at(0.25, 0.0, 0.0);
    std::vector<ir::FieldVector<Frame::GSM>> out(
        n, ir::FieldVector<Frame::GSM>{cheatah::fixarray::vec3d{7.0, 7.0, 7.0}});
    const ir::Result<bool> r = opq_field_batch(pts, 0.28, out);
    EXPECT_EQ(r.status, Status::DomainError);
    EXPECT_FALSE(r.value) << "a refused batch must not have reached the device";
    for (const ir::FieldVector<Frame::GSM>& b : out) {
        ASSERT_EQ(b.v[0], 0.0);
        ASSERT_EQ(b.v[1], 0.0);
        ASSERT_EQ(b.v[2], 0.0);
    }
    // One point past 15 Re is a caveat, not a refusal: the device still runs, and that point's
    // answer is the published zero.
    std::vector<Position<Frame::GSM>> far = scatter(n);
    far[n / 2] = at(-16.0, 0.0, 0.0);
    const ir::Result<bool> f = opq_field_batch(far, 0.28, out);
    EXPECT_EQ(f.status, Status::OutOfValidityRange);
    EXPECT_TRUE(f.value) << "an out-of-validity batch must still be computed on the device";
    EXPECT_EQ(out[n / 2].v[2], 0.0);
}

TEST(IrbemOpq, DeviceKernelAgreesWithTheHostLane) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 1 << 16;
    const std::vector<Position<Frame::GSM>> pts = scatter(n, 1.5, 15.5);   // taper and template too
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
    ASSERT_TRUE(opq_field_host(pos, host, sp, cp, opq_parameters<float>(ps * kDegPerRad)));
    const std::array<float, opq_param_count> block = opq_param_block(sp, cp, ps * kDegPerRad);
    ir::gpu::dispatch_batch("irbem_opq_f32", pos, device, std::span<const float>(block));
    double worst = 0.0;
    for (std::size_t i = 0; i < 3 * n; ++i) {
        worst = std::max(worst, std::fabs(static_cast<double>(device[i]) - host[i]));
    }
    std::printf("[ MEASURED ] device vs host, %zu points: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-3) << "the device is not evaluating the same expressions";
}
#endif
