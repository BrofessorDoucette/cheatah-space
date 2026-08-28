// Tests for space/irbem/field.hpp — GET_FIELD_MULTI, GET_BDERIVS, COMPUTE_GRAD_CURV_CURL and
// GET_HEMI_MULTI.
//
// Three kinds of assertion, and they are deliberately different in kind:
//
//   1. ANALYTIC. `Igrf<1>` is a tilted centred dipole and has a closed-form Jacobian, so the
//      finite-difference machinery is checked against exact mathematics rather than against
//      another approximation. It also gives `div B = 0` and `curl B = 0` for free, which is the
//      residual `COMPUTE_GRAD_CURV_CURL` carries precisely so a differencing error cannot hide.
//   2. GOLDEN. Numbers printed by `tools/oracle/field_differential.cpp`, which runs the compiled
//      IRBEM `-O2` build as a black box. Every golden below is quoted with the exact settings that
//      produced it — kext=0, options={0,0,0,0,0}, 2015-180 12:00 UT (so epoch 2015.5 by the
//      published options(2)=0 rule), sysaxes=1, dX=1e-3, IGRF truncated at degree 10 — because a
//      comparison at unmatched settings measures resolutions rather than algorithms.
//   3. SELF-CONSISTENCY. The batch lanes against the pointwise lane, and the device lane against
//      the host lane, at a matched step. A silent fallback to the host is the failure mode that
//      makes a speedup claim worthless, so the tests that mean to exercise the device assert that
//      it actually ran.
//
// Build standalone (the GPU lane is opt-in by include path, so this file compiles and passes both
// with and without cheatah-gpu-linalg on it):
//   g++ -std=c++20 -O3 -march=native -ffp-contract=off tests/irbem_field_test.cpp -I.
//       -I$CHEATAH_DIR/stdlib/{ndarray,builtins,fixarray} -Ibuild/debug/_deps/.../include
//       -Lbuild/debug/lib -lgtest_main -lgtest -pthread
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "space/irbem/field.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;
using ib::Frame;

namespace {

/// The epoch every oracle comparison uses. IRBEM's `options(2) = 0` is documented as "initialize
/// IGRF field once per year (year.5)", so a call stamped 2015 day 180 at 12:00 UT evaluates the
/// coefficients at exactly 2015.5 — not at the day-of-year fraction 2015.4918. Getting that wrong
/// costs ~3e-6 in |B|, which is three times the Bgeo budget and looks exactly like a bad kernel.
constexpr double kEpoch = 2015.5;

/// IRBEM truncates its internal IGRF at degree 10. Established black-box: sweeping the degree over
/// 8..13 against `GET_FIELD_MULTI` gives 1.5e-4, 4.2e-5, **1.7e-15**, 1.2e-5, 1.5e-5. Ten orders
/// of magnitude in one column is not a coincidence of tolerance.
using OracleIgrf = ib::Igrf<10>;

/// The full IGRF-14 truncation, for the tests that are about this file's arithmetic rather than
/// about agreeing with IRBEM.
using FullIgrf = ib::Igrf<13>;

/// A tilted centred dipole: IGRF truncated at degree 1. Its potential is `V = (m·r)/r³` with the
/// moment vector `m = (g₁¹, h₁¹, g₁⁰)` in GEO — read straight off the degree-1 terms of the IAGA
/// series, since `P₁⁰ = cos θ = z/r` and `P₁¹ cos φ = sin θ cos φ = x/r`. Everything about this
/// field is closed form, which is what makes it a real check on a finite difference.
struct AnalyticDipole {
    fx::vec3d m{};

    /// The moment of a model's degree-1 part.
    template <int NMAX>
    static AnalyticDipole of(const ib::Igrf<NMAX>& model) {
        return AnalyticDipole{fx::vec3d{model.g(1, 1), model.h(1, 1), model.g(1, 0)}};
    }

    /// `B = 3(m·r)r/r⁵ − m/r³`.
    [[nodiscard]] fx::vec3d b(const fx::vec3d& r) const {
        const double rn = fx::norm(r);
        const double r3 = rn * rn * rn;
        const double r5 = r3 * rn * rn;
        const double mr = fx::dot(m, r);
        return (r * (3.0 * mr / r5)) - (m / r3);
    }

    /// `∂Bᵢ/∂xⱼ`, by differentiating the line above: three terms and a `−15(m·r)rᵢrⱼ/r⁷`.
    [[nodiscard]] double jacobian(const fx::vec3d& r, std::size_t i, std::size_t j) const {
        const double rn = fx::norm(r);
        const double r5 = std::pow(rn, 5.0);
        const double r7 = std::pow(rn, 7.0);
        const double mr = fx::dot(m, r);
        const double delta = (i == j) ? 1.0 : 0.0;
        return ((3.0 * ((m[j] * r[i]) + (mr * delta))) / r5) - ((15.0 * mr * r[i] * r[j]) / r7) +
               ((3.0 * m[i] * r[j]) / r5);
    }

    /// `∂|B|/∂xⱼ = Σᵢ Bᵢ ∂Bᵢ/∂xⱼ / |B|`.
    [[nodiscard]] fx::vec3d grad_magnitude(const fx::vec3d& r) const {
        const fx::vec3d bv = b(r);
        const double bm = fx::norm(bv);
        fx::vec3d g{};
        for (std::size_t j = 0; j < 3; ++j) {
            double s = 0.0;
            for (std::size_t i = 0; i < 3; ++i) { s += bv[i] * jacobian(r, i, j); }
            g[j] = s / bm;
        }
        return g;
    }

    /// The @ref ib::GeoFieldModel surface, so this can be handed to `bderivs` directly.
    [[nodiscard]] ib::FieldVector<Frame::GEO> evaluate(const ib::Position<Frame::GEO>& p) const {
        return ib::FieldVector<Frame::GEO>{b(p.v)};
    }
};
/// A perfectly uniform field. Physically it is the one field with no magnetic equator anywhere —
/// `|B|` is the same in every direction, so `d|B|/ds` is exactly zero and neither hemisphere is the
/// answer. It exists here because that exact-zero case is a real branch of @ref ib::hemisphere and
/// the only way to reach it deliberately.
struct UniformField {
    fx::vec3d b{0.0, 0.0, 100.0};

    [[nodiscard]] ib::FieldVector<Frame::GEO> evaluate(const ib::Position<Frame::GEO>& /*p*/) const {
        return ib::FieldVector<Frame::GEO>{b};
    }
};

static_assert(ib::GeoFieldModel<AnalyticDipole>,
              "the analytic dipole must satisfy the same concept the real models do");
static_assert(ib::GeoFieldModel<ib::Igrf<13>>, "Igrf is the model the concept was written for");
static_assert(!ib::GeoFieldModel<double>, "the concept must actually exclude something");

/// A reproducible spread of points from LEO to beyond geosynchronous, off the equator.
std::vector<ib::Position<Frame::GEO>> sample_points(std::size_t n, unsigned seed = 7) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<ib::Position<Frame::GEO>> pts(n);
    for (auto& p : pts) {
        const double a = u(rng) * 3.14159265358979;
        const double b = u(rng) * 1.45;
        const double r = 1.5 + (std::abs(u(rng)) * 7.0);
        p = ib::Position<Frame::GEO>{
            fx::vec3d{r * std::cos(b) * std::cos(a), r * std::cos(b) * std::sin(a), r * std::sin(b)}};
    }
    return pts;
}

/// Relative difference, guarded so a zero reference reads as an absolute one.
double rel(double a, double b) { return std::abs(a - b) / std::max(1e-12, std::abs(b)); }

/// One oracle golden: everything `tools/oracle/field_differential.cpp` prints for one point.
struct Golden {
    fx::vec3d p;
    fx::vec3d b;
    double b_mag;
    fx::vec3d grad_b;
    std::array<double, 9> diff_b;  ///< flat `3j + i` = `∂Bᵢ/∂xⱼ`, IRBEM's own column-major order
    double grad_par;
    double r_curv;
    double div_b;
    int hemi;
    fx::vec3d curl_b;
    fx::vec3d curvature;
    /// The three outputs whose MAGNITUDE nothing else in this file pins. Without them a scale
    /// factor on any of the three — a missing `1/Bmag` on `grad_drift`, say — passes every other
    /// assertion here, because the rest are orthogonality relations and a scale factor leaves
    /// those exactly intact. Verified: doubling `grad_perp`, `grad_drift` or `curv_drift` passed
    /// 24 of 24 tests before these were added.
    fx::vec3d grad_perp;
    fx::vec3d grad_drift;
    fx::vec3d curv_drift;
};

/// Printed by `tools/oracle/field_differential.cpp` against `/tmp/irbem-builds/libirbem-O2.so`
/// (IRBEM built `-O2 -ffp-contract=off -fno-fast-math`) at kext=0, options={0,0,0,0,0},
/// sysaxes=1, iyear=2015, idoy=180, UT=43200, dX=1e-3.
const std::array<Golden, 4> kGoldens{{
    {fx::vec3d{2.0, 0.0, 0.0},
     fx::vec3d{402.430556, -483.318073, 3392.41413},
     3450.220329,
     fx::vec3d{-5129.02045, 160.002969, 158.158525},
     {-1152.5443, 646.042321, -4987.70051, 645.991425, 381.3266, 140.343255, -4993.82164,
      141.689953, 769.678032},
     -465.149328,
     0.674063942,
     -1.53967192,
     -1,
     fx::vec3d{-1.34669861, -6.121123, 0.0508960765},
     fx::vec3d{-1.47248554, 0.0277037521, 0.178622915},
     fx::vec3d{-5074.76586, 94.8433397, 615.514577},
     fx::vec3d{-0.0520191657, -1.46701776, -0.202835524},
     fx::vec3d{-0.052261672, -1.46864942, -0.203039219}},
    {fx::vec3d{3.0, 1.0, 0.5},
     fx::vec3d{-274.607367, -208.669366, 837.355493},
     905.6027446,
     fx::vec3d{-808.08802, -259.840188, -18.5595417},
     {298.927174, 214.732412, -722.409087, 214.865141, 9.86461851, -208.10877, -722.589311,
      -207.954384, -309.259294},
     287.749433,
     1.13332509,
     -0.467501519,
     1,
     fx::vec3d{-0.154385684, -0.180224269, -0.132729456},
     fx::vec3d{-0.796302935, -0.213606822, -0.314375262},
     fx::vec3d{-720.8333, -193.536842, -284.623851},
     fx::vec3d{0.270024395, -0.83128881, -0.118604131},
     fx::vec3d{0.26994765, -0.831621154, -0.118712119}},
    {fx::vec3d{-4.0, 2.0, 1.5},
     fx::vec3d{169.459647, -118.026119, 239.435133},
     316.1896888,
     fx::vec3d{182.018301, -93.2461982, -32.0132132},
     {121.553579, -96.0266279, 106.997983, -95.9655765, 19.0728138, -45.8228182, 106.912997,
      -45.7835158, -140.578657},
     108.115973,
     1.79186761,
     0.0477357051,
     1,
     fx::vec3d{-0.0393024063, -0.084985231, -0.0610514323},
     fx::vec3d{0.392183297, -0.167340157, -0.360054731},
     fx::vec3d{124.074303, -52.8890672, -113.884202},
     fx::vec3d{0.261110989, 0.490183278, 0.0568281433},
     fx::vec3d{0.261118494, 0.489950219, 0.0567079486}},
    {fx::vec3d{5.5, 0.0, 2.0},
     fx::vec3d{-144.308017, -21.3782364, 91.759924},
     172.3418591,
     fx::vec3d{-89.5738108, -2.64332972, -9.05548778},
     {87.3332661, 9.67438723, -28.6384718, 9.68737299, -26.2878631, 4.14149803, -28.6191911,
      4.14430746, -61.0752745},
     70.5098456,
     3.07548428,
     -0.0298714284,
     1,
     fx::vec3d{-0.002809428, 0.01928069, -0.0129857529},
     fx::vec3d{-0.177152226, 0.0354672736, -0.270338701},
     fx::vec3d{-30.5333894, 6.10310111, -46.5970235},
     fx::vec3d{0.0146840703, -0.320724413, -0.0516292195},
     fx::vec3d{0.0146504764, -0.320685392, -0.0516729605}},
}};

/// Whether a device is present AND its kernel is compiled — the tests that assert the device ran
/// skip themselves rather than fail on a machine that has none.
bool device_ready() {
#ifdef CHEATAH_SPACE_IRBEM_FIELD_GPU
    return ib::gpu::available() && std::filesystem::exists(ib::gpu::shader_path("irbem_igrf_f32"));
#else
    return false;
#endif
}

/// Set `CHEATAH_SPACE_IRBEM_GPU_CROSSOVER` for a scope and put it back — the seam the batch
/// routines consult, and the only way to drive the device lane at a batch small enough to compare
/// element by element in a unit test.
class CrossoverOverride {
  public:
    explicit CrossoverOverride(const char* value) {
        if (const char* prev = std::getenv(kVar)) {
            had_ = true;
            prev_ = prev;
        }
        ::setenv(kVar, value, 1);
    }
    CrossoverOverride(const CrossoverOverride&) = delete;
    CrossoverOverride& operator=(const CrossoverOverride&) = delete;
    CrossoverOverride(CrossoverOverride&&) = delete;
    CrossoverOverride& operator=(CrossoverOverride&&) = delete;
    ~CrossoverOverride() {
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

/// Point `CHEATAH_SPACE_IRBEM_SPV_DIR` somewhere for a scope and put it back. `gpu::shader_dir`
/// consults that variable first, so this is how a "the shaders were never built" deployment is
/// reproduced on a machine where they were.
class ShaderDirOverride {
  public:
    explicit ShaderDirOverride(const char* value) {
        if (const char* prev = std::getenv(kVar)) {
            had_ = true;
            prev_ = prev;
        }
        ::setenv(kVar, value, 1);
    }
    ShaderDirOverride(const ShaderDirOverride&) = delete;
    ShaderDirOverride& operator=(const ShaderDirOverride&) = delete;
    ShaderDirOverride(ShaderDirOverride&&) = delete;
    ShaderDirOverride& operator=(ShaderDirOverride&&) = delete;
    ~ShaderDirOverride() {
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

// ---------------------------------------------------------------------------------------------
// The step
// ---------------------------------------------------------------------------------------------

// The whole design of GET_BDERIVS is in these two constants, so they are pinned rather than left to
// drift: the ratio scales with radius (both error terms do), the device's is four orders larger
// than the host's (that ratio is sqrt of the ratio of the two lanes' effective epsilons), and a
// radius at or inside 1 Re floors rather than collapsing the step to zero.
TEST(IrbemField, AutoStepTracksTheRadiusAndTheLane) {
    EXPECT_EQ(ib::auto_step(1.0, ib::DifferenceLane::Fp64Host), ib::host_step_ratio);
    EXPECT_EQ(ib::auto_step(4.0, ib::DifferenceLane::Fp64Host), 4.0 * ib::host_step_ratio);
    EXPECT_EQ(ib::auto_step(4.0, ib::DifferenceLane::Fp32Device), 4.0 * ib::device_step_ratio);
    // Inside the Earth, at the origin, and NaN all floor at r = 1 rather than producing a zero or
    // negative step that would divide by zero three lines later.
    EXPECT_EQ(ib::auto_step(0.0, ib::DifferenceLane::Fp64Host), ib::host_step_ratio);
    EXPECT_EQ(ib::auto_step(-3.0, ib::DifferenceLane::Fp64Host), ib::host_step_ratio);
    EXPECT_EQ(ib::auto_step(std::nan(""), ib::DifferenceLane::Fp32Device), ib::device_step_ratio);
    // The device needs a far larger step, and by a factor that is a measurement, not a taste.
    EXPECT_GT(ib::device_step_ratio / ib::host_step_ratio, 1000.0);
}

#ifdef CHEATAH_SPACE_IRBEM_FIELD_GPU
// The device divides by the step it ACTUALLY took, not the one it was asked for -- a correction
// the header devotes a paragraph to and which, measured over 4096 points at the automatic device
// step, is worth 2.3e-4 relative on every derivative. From inside the device lane the correction
// is smaller than the fp32 field noise it sits under, so deleting it there is invisible: it was.
// Pinning it here, on exactly representable values, is what makes the claim falsifiable.
TEST(IrbemField, TheRealisedStepIsTheOneTheFloatLaneCanTake) {
    // 4.0 is exact in float and its float neighbour is 4 + 2^-21, so the step lands on one of
    // exactly three answers and each is a different branch of the routine.
    constexpr double kUlp = 1.0 / 2097152.0;   // 2^-21, the float spacing at 4.0
    EXPECT_EQ(ib::detail::realised_step(4.0, kUlp), kUlp);          // representable: taken exactly
    EXPECT_EQ(ib::detail::realised_step(4.0, 0.75 * kUlp), kUlp);   // rounds UP to one ulp, and
                                                                    // dividing by 0.75*kUlp would
                                                                    // be a 33% error on the
                                                                    // derivative
    EXPECT_EQ(ib::detail::realised_step(4.0, 0.25 * kUlp), 0.25 * kUlp);  // collapses onto x, so
                                                                          // the intended step is
                                                                          // returned rather than
                                                                          // a division by zero
    // And a realistic case: the automatic device step at r = 8 is not representable, and the
    // realised step differs from it by the ~1e-4 relative the header claims.
    const double h = ib::auto_step(8.0, ib::DifferenceLane::Fp32Device);
    const double taken = ib::detail::realised_step(8.0, h);
    EXPECT_NE(taken, h);
    EXPECT_LT(std::abs(taken - h) / h, 1.0e-3);
    EXPECT_GT(std::abs(taken - h) / h, 1.0e-6) << "if this is negligible the correction is dead code";
    // It never returns zero, whatever it is handed -- the caller divides by it.
    EXPECT_NE(ib::detail::realised_step(1.0e30, 1.0e-30), 0.0);
}
#endif

// ---------------------------------------------------------------------------------------------
// GET_BDERIVS, against exact mathematics
// ---------------------------------------------------------------------------------------------

// A tilted centred dipole has a closed-form Jacobian, so this compares the finite difference to the
// answer rather than to another approximation. The tolerance is the truncation error the file
// brief predicts — about 2h/r for a one-sided difference — and it is asserted as a BOUND on the
// error, so a step-selection regression that made the derivative worse would fail here.
TEST(IrbemField, BderivsMatchesTheAnalyticDipoleJacobian) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const AnalyticDipole dip = AnalyticDipole::of(model);
    for (const auto& p : sample_points(40)) {
        const ib::Result<ib::BDerivatives> d = ib::bderivs(dip, p);
        ASSERT_EQ(d.status, ib::Status::Ok);
        const double r = fx::norm(p.v);
        const double scale = fx::norm(dip.grad_magnitude(p.v));
        EXPECT_LT(fx::norm(d.value.grad_b_mag - dip.grad_magnitude(p.v)) / scale, 1.0e-6)
            << "at r = " << r;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                EXPECT_LT(std::abs(d.value.diff_b(i, j) - dip.jacobian(p.v, i, j)) / scale, 1.0e-6)
                    << "component (" << i << "," << j << ") at r = " << r;
            }
        }
    }
}

// The reference differences FORWARD, and that is not a stylistic detail: at a given step a central
// difference is two orders more accurate, so an implementation that quietly used one would agree
// with IRBEM far LESS well while looking better against truth. This pins the scheme by showing the
// error is first order in the step -- halving the step halves the error, which a central difference
// would quarter.
TEST(IrbemField, BderivsIsAForwardDifferenceNotACentralOne) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const AnalyticDipole dip = AnalyticDipole::of(model);
    const ib::Position<Frame::GEO> p{fx::vec3d{3.0, 1.0, 0.5}};
    const fx::vec3d truth = dip.grad_magnitude(p.v);
    const double coarse = fx::norm(ib::bderivs(dip, p, 1.0e-3).value.grad_b_mag - truth);
    const double fine = fx::norm(ib::bderivs(dip, p, 5.0e-4).value.grad_b_mag - truth);
    // First order: halving the step halves the error, to within a few percent. A central difference
    // would land near 0.25 and fail both bounds.
    EXPECT_GT(coarse / fine, 1.9);
    EXPECT_LT(coarse / fine, 2.1);
    // And the difference really is one-sided: a forward difference of a function with positive
    // curvature lands on the same side of the truth at every step, which a central one does not.
    const fx::vec3d e1 = ib::bderivs(dip, p, 1.0e-3).value.grad_b_mag - truth;
    const fx::vec3d e2 = ib::bderivs(dip, p, 5.0e-4).value.grad_b_mag - truth;
    EXPECT_GT(e1[0] * e2[0], 0.0);
}

// Every failure path of the pointwise routine, because a derivative routine that silently returns
// zeros at a bad input is worse than one that refuses.
TEST(IrbemField, BderivsRefusesInputsItCannotAnswer) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const ib::Position<Frame::GEO> origin{fx::vec3d{0.0, 0.0, 0.0}};
    EXPECT_EQ(ib::bderivs(model, origin).status, ib::Status::DomainError);
    const ib::Position<Frame::GEO> nan_point{fx::vec3d{std::nan(""), 0.0, 0.0}};
    EXPECT_EQ(ib::bderivs(model, nan_point).status, ib::Status::DomainError);
    const ib::Position<Frame::GEO> good{fx::vec3d{3.0, 0.0, 0.0}};
    EXPECT_EQ(ib::bderivs(model, good, -1.0).status, ib::Status::DomainError);
    EXPECT_EQ(ib::bderivs(model, good, std::nan("")).status, ib::Status::DomainError);
    // The value is zero-filled on refusal rather than left indeterminate.
    EXPECT_EQ(ib::bderivs(model, origin).value.b_mag, 0.0);
}

// ---------------------------------------------------------------------------------------------
// COMPUTE_GRAD_CURV_CURL
// ---------------------------------------------------------------------------------------------

// The eight outputs, reproduced from the oracle's own inputs. `tools/oracle/field_differential.cpp`
// runs this comparison over 46 storm configurations (Kp 0-9 under T89; Dst 0..-400 nT under T96 and
// T01-storm; Pdyn 0.5..40 nPa; Bz 0..-30 nT southward) x 300 points and reports grad_par,
// grad_perp, grad_drift, curl_b and div_b at exactly 0.0, r_curv at 5.9e-14, curvature at 7.7e-12
// and curv_drift at 9.1e-12. Here the same algebra is checked against four goldens captured from
// that run, so the property survives in the gate without linking IRBEM.
TEST(IrbemField, GradCurvCurlReproducesTheOracleAlgebra) {
    for (const Golden& g : kGoldens) {
        ib::BDerivatives d{};
        d.b = ib::FieldVector<Frame::GEO>{g.b};
        d.b_mag = g.b_mag;
        d.grad_b_mag = g.grad_b;
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t i = 0; i < 3; ++i) { d.diff_b(i, j) = g.diff_b[(3 * j) + i]; }
        }
        const ib::Result<ib::GradCurvCurl> r = ib::grad_curv_curl(d);
        ASSERT_EQ(r.status, ib::Status::Ok);
        // The goldens carry nine significant figures, so 1e-8 is the print precision and not a
        // slack tolerance: the underlying agreement measured against unrounded oracle output is
        // exactly 0.0 for five of the eight outputs and 1e-12 for the rest.
        EXPECT_LT(rel(r.value.grad_par, g.grad_par), 1.0e-8);
        EXPECT_LT(rel(r.value.r_curv, g.r_curv), 1.0e-8);
        for (std::size_t k = 0; k < 3; ++k) {
            EXPECT_LT(rel(r.value.curvature[k], g.curvature[k]), 1.0e-7) << "curvature " << k;
        }
        // div B and curl B must be judged ABSOLUTELY, against the size of the Jacobian they are
        // built from. Both are near-total cancellations: the diffB entries run to ~5000 nT/Re and
        // their signed sum is ~1, so the goldens' 1e-5 absolute print truncation is already a 1e-5
        // RELATIVE error in div B before either implementation does anything. That is not a
        // weakness of the test — it is the same cancellation that makes div B a useful residual,
        // measured from the other side.
        double jacobian_scale = 0.0;
        for (const double e : g.diff_b) jacobian_scale = std::max(jacobian_scale, std::abs(e));
        EXPECT_LT(std::abs(r.value.div_b - g.div_b), 1.0e-8 * jacobian_scale);
        for (std::size_t k = 0; k < 3; ++k) {
            EXPECT_LT(std::abs(r.value.curl_b[k] - g.curl_b[k]), 1.0e-8 * jacobian_scale)
                << "curl component " << k;
        }
        // grad_perp, grad_drift and curv_drift are pinned by VALUE, not only by their
        // orthogonality relations. That distinction is the whole point: an orthogonality
        // assertion is invariant under a scale factor, so before these three goldens existed,
        // doubling any of the three passed all 24 tests. The oracle's own numbers are the only
        // thing that pins the 1/Bmag in `grad_drift` and the |B| that is absent from `curv_drift`.
        for (std::size_t k = 0; k < 3; ++k) {
            const double gp_scale = fx::norm(g.grad_perp);
            EXPECT_LT(std::abs(r.value.grad_perp[k] - g.grad_perp[k]) / gp_scale, 1.0e-8)
                << "grad_perp " << k;
            EXPECT_LT(std::abs(r.value.grad_drift[k] - g.grad_drift[k]) / fx::norm(g.grad_drift),
                      1.0e-8)
                << "grad_drift " << k;
            EXPECT_LT(std::abs(r.value.curv_drift[k] - g.curv_drift[k]) / fx::norm(g.curv_drift),
                      1.0e-7)
                << "curv_drift " << k;
        }
        // And then the definitions on top of the values: the perpendicular part is perpendicular,
        // and the drift is perpendicular to both B and the gradient that drives it.
        // The tolerance is set by the goldens, not by the arithmetic: nine printed significant
        // figures means the reconstructed B̂ is itself only good to ~1e-9, so the perpendicularity
        // residual cannot be smaller than that. Against unrounded inputs it is 1e-16 — which is
        // what CurvatureIsPerpendicularToTheField below actually measures.
        const fx::vec3d bhat = g.b / g.b_mag;
        EXPECT_LT(std::abs(fx::dot(r.value.grad_perp, bhat)) / fx::norm(g.grad_b), 1.0e-8);
        EXPECT_LT(std::abs(fx::dot(r.value.grad_drift, bhat)), 1.0e-12);
        EXPECT_LT(std::abs(fx::dot(r.value.grad_drift, r.value.grad_perp)), 1.0e-12);
    }
}

// `(B̂·∇)B̂` is perpendicular to `B̂` identically, because `B̂·B̂ = 1`. The implementation gets that
// by construction — it projects out the parallel part of `Â` rather than subtracting
// `grad_par·B̂/|B|`, which is only equal when the gradient and the Jacobian are exactly consistent.
// This is the test that would fail if that choice were reverted, and it is also the reason the
// oracle probe could tell the two candidates apart.
TEST(IrbemField, CurvatureIsPerpendicularToTheField) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    for (const auto& p : sample_points(30, 11)) {
        const ib::Result<ib::BDerivatives> d = ib::bderivs(model, p);
        const ib::Result<ib::GradCurvCurl> g = ib::grad_curv_curl(d.value);
        ASSERT_EQ(g.status, ib::Status::Ok);
        const fx::vec3d bhat = d.value.b.v / d.value.b_mag;
        const double kappa = fx::norm(g.value.curvature);
        EXPECT_LT(std::abs(fx::dot(g.value.curvature, bhat)) / kappa, 1.0e-14);
        EXPECT_LT(std::abs(fx::dot(g.value.curv_drift, bhat)), 1.0e-14 * kappa);
        EXPECT_NEAR(g.value.r_curv * kappa, 1.0, 1.0e-12);
    }
}

// div B = 0 is Maxwell, so whatever this routine reports is the differencing residual and nothing
// else — which is what makes it the free diagnostic the struct documents. On a dipole the field is
// exactly divergence- and curl-free, so the residual must sit at the truncation error of the step
// and not above it. A test that only checked "small" would pass on a broken Jacobian; this one
// checks the residual is small RELATIVE to the derivatives it is built from.
TEST(IrbemField, DivergenceAndCurlAreTheDifferencingResidual) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const AnalyticDipole dip = AnalyticDipole::of(model);
    for (const auto& p : sample_points(25, 13)) {
        const ib::Result<ib::BDerivatives> d = ib::bderivs(dip, p, 1.0e-5);
        const ib::Result<ib::GradCurvCurl> g = ib::grad_curv_curl(d.value);
        ASSERT_EQ(g.status, ib::Status::Ok);
        const double scale = fx::norm(d.value.grad_b_mag);
        EXPECT_LT(std::abs(g.value.div_b) / scale, 1.0e-4);
        EXPECT_LT(fx::norm(g.value.curl_b) / scale, 1.0e-4);
    }
}

// A field of exactly zero has no direction, so every output is meaningless; and a straight field
// line has an infinite radius of curvature, which is the physics rather than an overflow.
TEST(IrbemField, GradCurvCurlHandlesTheDegenerateCases) {
    ib::BDerivatives zero{};
    EXPECT_EQ(ib::grad_curv_curl(zero).status, ib::Status::DomainError);
    ib::BDerivatives nan_mag{};
    nan_mag.b_mag = std::nan("");
    EXPECT_EQ(ib::grad_curv_curl(nan_mag).status, ib::Status::DomainError);

    // A uniform field along z: no gradient, no Jacobian, so no curvature at all.
    ib::BDerivatives uniform{};
    uniform.b = ib::FieldVector<Frame::GEO>{fx::vec3d{0.0, 0.0, 100.0}};
    uniform.b_mag = 100.0;
    const ib::Result<ib::GradCurvCurl> g = ib::grad_curv_curl(uniform);
    ASSERT_EQ(g.status, ib::Status::Ok);
    EXPECT_EQ(g.value.grad_par, 0.0);
    EXPECT_EQ(fx::norm(g.value.curvature), 0.0);
    EXPECT_TRUE(std::isinf(g.value.r_curv));
    EXPECT_EQ(g.value.div_b, 0.0);
}

// ---------------------------------------------------------------------------------------------
// GET_HEMI_MULTI
// ---------------------------------------------------------------------------------------------

// The four goldens' hemispheres, and the reason the first one is southern is worth stating: (2,0,0)
// is on the GEOGRAPHIC equator, and it comes back SOUTH because the dipole is tilted and its
// magnetic equator there lies to the north of it. A geographic test would get this backwards, which
// is exactly why the routine exists.
TEST(IrbemField, HemisphereAgreesWithTheOracleGoldens) {
    const OracleIgrf model = *OracleIgrf::at(kEpoch);
    for (const Golden& g : kGoldens) {
        const ib::Result<ib::Hemisphere> h =
            ib::hemisphere(model, ib::Position<Frame::GEO>{g.p});
        EXPECT_EQ(h.status, ib::Status::Ok);
        EXPECT_EQ(static_cast<int>(h.value), g.hemi)
            << "at (" << g.p[0] << ", " << g.p[1] << ", " << g.p[2] << ")";
        // The hemisphere is the sign of grad_par -- the same scalar COMPUTE_GRAD_CURV_CURL
        // reports -- which is why the two live in one file. The goldens pin both, so this also
        // pins that they agree.
        EXPECT_EQ(g.grad_par > 0.0 ? 1 : -1, g.hemi);
    }
}

// Walking a meridian of the tilted dipole from deep south to deep north must cross exactly once,
// and the crossing must be where the parallel gradient changes sign. A routine that returned a
// constant, or that keyed off geographic latitude, passes neither half.
TEST(IrbemField, HemisphereFlipsAcrossTheDipoleEquator) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const AnalyticDipole dip = AnalyticDipole::of(model);
    int flips = 0;
    ib::Hemisphere previous = ib::Hemisphere::Invalid;
    for (int k = -60; k <= 60; ++k) {
        const double z = 0.05 * k;
        const ib::Position<Frame::GEO> p{fx::vec3d{4.0, 0.0, z}};
        const ib::Result<ib::Hemisphere> h = ib::hemisphere(dip, p);
        ASSERT_EQ(h.status, ib::Status::Ok);
        ASSERT_NE(h.value, ib::Hemisphere::Invalid);
        if (previous != ib::Hemisphere::Invalid && h.value != previous) ++flips;
        previous = h.value;
    }
    EXPECT_EQ(flips, 1);
    // Deep in each hemisphere the answer is unambiguous.
    EXPECT_EQ(ib::hemisphere(dip, ib::Position<Frame::GEO>{fx::vec3d{4.0, 0.0, 3.0}}).value,
              ib::Hemisphere::North);
    EXPECT_EQ(ib::hemisphere(dip, ib::Position<Frame::GEO>{fx::vec3d{4.0, 0.0, -3.0}}).value,
              ib::Hemisphere::South);
}

// A field that vanishes has no hemisphere, and a point with no radius has no field. IRBEM spells
// the first `0`, and so does this.
TEST(IrbemField, HemisphereRefusesInputsItCannotAnswer) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    EXPECT_EQ(ib::hemisphere(model, ib::Position<Frame::GEO>{fx::vec3d{0.0, 0.0, 0.0}}).status,
              ib::Status::DomainError);
    EXPECT_EQ(ib::hemisphere(model, ib::Position<Frame::GEO>{fx::vec3d{3.0, 0.0, 0.0}}, -1.0).status,
              ib::Status::DomainError);
    const AnalyticDipole null_field{fx::vec3d{0.0, 0.0, 0.0}};
    const ib::Result<ib::Hemisphere> h =
        ib::hemisphere(null_field, ib::Position<Frame::GEO>{fx::vec3d{3.0, 0.0, 0.0}});
    EXPECT_EQ(h.status, ib::Status::Ok);
    EXPECT_EQ(h.value, ib::Hemisphere::Invalid);
    EXPECT_EQ(static_cast<int>(ib::Hemisphere::Invalid), 0);   // IRBEM's own xHEMI value
    EXPECT_EQ(static_cast<int>(ib::Hemisphere::North), 1);
    EXPECT_EQ(static_cast<int>(ib::Hemisphere::South), -1);

    // A uniform field has no magnetic equator at all, so the parallel derivative is EXACTLY zero
    // and neither hemisphere is the answer. This is the one input that separates "the field
    // vanishes" from "the field has no structure", and both come back as IRBEM's 0.
    const UniformField flat{};
    const ib::Result<ib::Hemisphere> u =
        ib::hemisphere(flat, ib::Position<Frame::GEO>{fx::vec3d{4.0, 1.0, -2.0}});
    EXPECT_EQ(u.status, ib::Status::Ok);
    EXPECT_EQ(u.value, ib::Hemisphere::Invalid);
}

// ---------------------------------------------------------------------------------------------
// GET_FIELD_MULTI
// ---------------------------------------------------------------------------------------------

// The four goldens' field vectors and magnitudes, at IRBEM's own truncation and epoch. Nine
// significant figures, and the underlying agreement measured over 300 points is 1.7e-15 relative —
// nine orders inside the 1e-6 Bgeo budget of docs/ERROR_BUDGET.md §4.
TEST(IrbemField, FieldBatchMatchesTheOracleGoldens) {
    const OracleIgrf model = *OracleIgrf::at(kEpoch);
    std::vector<ib::Position<Frame::GEO>> pts;
    for (const Golden& g : kGoldens) pts.emplace_back(ib::Position<Frame::GEO>{g.p});
    std::vector<ib::FieldVector<Frame::GEO>> b(pts.size());
    std::vector<double> bm(pts.size());
    const ib::Result<bool> r = ib::field_batch(model, pts, b, bm);
    ASSERT_EQ(r.status, ib::Status::Ok);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        EXPECT_LT(rel(bm[i], kGoldens[i].b_mag), 1.0e-9);
        for (std::size_t k = 0; k < 3; ++k) {
            EXPECT_LT(rel(b[i].v[k], kGoldens[i].b[k]), 1.0e-8) << "point " << i << " comp " << k;
        }
    }
    // And the truncation matters: at degree 13 the same comparison is four orders looser, because
    // the two libraries are then evaluating DIFFERENT models. Asserting that keeps a future
    // "upgrade" to degree 13 from silently loosening this test.
    const FullIgrf full = *FullIgrf::at(kEpoch);
    std::vector<ib::FieldVector<Frame::GEO>> b13(pts.size());
    std::vector<double> bm13(pts.size());
    ASSERT_EQ(ib::field_batch(full, pts, b13, bm13).status, ib::Status::Ok);
    double worst = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) worst = std::max(worst, rel(bm13[i], kGoldens[i].b_mag));
    EXPECT_GT(worst, 1.0e-7) << "degree 13 should NOT match IRBEM's degree-10 field to noise";
}

// The batch entry point and the model's own pointwise evaluation must be the same computation on
// the host lane -- not close, identical. A batch routine that reordered or refactored the
// arithmetic would show up here as a last-bit difference.
TEST(IrbemField, FieldBatchAgreesWithTheReferenceLane) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto pts = sample_points(64);
    std::vector<ib::FieldVector<Frame::GEO>> b(pts.size());
    std::vector<double> bm(pts.size());
    const ib::Result<bool> r = ib::field_batch(model, pts, b, bm);
    ASSERT_EQ(r.status, ib::Status::Ok);
    EXPECT_FALSE(r.value) << "64 points is below every crossover; this must be the host lane";
    for (std::size_t i = 0; i < pts.size(); ++i) {
        EXPECT_EQ(b[i], model.evaluate(pts[i]));
        EXPECT_EQ(bm[i], model.evaluate(pts[i]).magnitude());
    }
}

// Malformed batches are refused, and an empty one is a no-op rather than an error -- an ephemeris
// with no points is a legitimate thing to hand a batch routine.
TEST(IrbemField, BatchRoutinesRefuseMismatchedSpans) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto pts = sample_points(4);
    std::vector<ib::FieldVector<Frame::GEO>> b(3);
    std::vector<double> bm(4);
    EXPECT_EQ(ib::field_batch(model, pts, b, bm).status, ib::Status::DomainError);
    std::vector<ib::FieldVector<Frame::GEO>> b4(4);
    std::vector<double> bm3(3);
    EXPECT_EQ(ib::field_batch(model, pts, b4, bm3).status, ib::Status::DomainError);

    std::vector<ib::BDerivatives> d(3);
    EXPECT_EQ(ib::bderivs_batch(model, pts, d).status, ib::Status::DomainError);
    std::vector<ib::BDerivatives> d4(4);
    EXPECT_EQ(ib::bderivs_batch(model, pts, d4, -1.0).status, ib::Status::DomainError);
    EXPECT_EQ(ib::bderivs_batch(model, pts, d4, std::nan("")).status, ib::Status::DomainError);

    std::vector<ib::Hemisphere> h(3);
    EXPECT_EQ(ib::hemisphere_batch(model, pts, h).status, ib::Status::DomainError);
    std::vector<ib::Hemisphere> h4(4);
    EXPECT_EQ(ib::hemisphere_batch(model, pts, h4, -1.0).status, ib::Status::DomainError);

    std::vector<ib::GradCurvCurl> g(3);
    std::vector<ib::Status> st(4);
    EXPECT_EQ(ib::grad_curv_curl_batch(d4, g, st), ib::Status::DomainError);

    const std::span<const ib::Position<Frame::GEO>> none;
    std::vector<ib::FieldVector<Frame::GEO>> nb;
    std::vector<double> nm;
    std::vector<ib::BDerivatives> nd;
    std::vector<ib::Hemisphere> nh;
    EXPECT_EQ(ib::field_batch(model, none, nb, nm).status, ib::Status::Ok);
    EXPECT_EQ(ib::bderivs_batch(model, none, nd).status, ib::Status::Ok);
    EXPECT_EQ(ib::hemisphere_batch(model, none, nh).status, ib::Status::Ok);
}

// ---------------------------------------------------------------------------------------------
// The batch lanes against the pointwise lanes
// ---------------------------------------------------------------------------------------------

// The four goldens' gradBmag and diffB, at IRBEM's dX = 1e-3. The full comparison in
// tools/oracle/field_differential.cpp is 4.5e-12 relative over 300 points; the tolerance here is
// the goldens' own print precision.
TEST(IrbemField, BderivsBatchMatchesTheOracleGoldens) {
    const OracleIgrf model = *OracleIgrf::at(kEpoch);
    std::vector<ib::Position<Frame::GEO>> pts;
    for (const Golden& g : kGoldens) pts.emplace_back(ib::Position<Frame::GEO>{g.p});
    std::vector<ib::BDerivatives> d(pts.size());
    ASSERT_EQ(ib::bderivs_batch(model, pts, d, 1.0e-3).status, ib::Status::Ok);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const double scale = fx::norm(kGoldens[i].grad_b);
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_LT(std::abs(d[i].grad_b_mag[j] - kGoldens[i].grad_b[j]) / scale, 1.0e-8)
                << "gradBmag " << i << "," << j;
            for (std::size_t c = 0; c < 3; ++c) {
                EXPECT_LT(std::abs(d[i].diff_b(c, j) - kGoldens[i].diff_b[(3 * j) + c]) / scale,
                          1.0e-8)
                    << "diffB (" << c << "," << j << ") at point " << i;
            }
        }
    }
}

// The batch lane must be the pointwise lane, point for point, on the host. Exact equality, because
// it IS the same code path -- and if it ever stops being, that is the thing to find out.
TEST(IrbemField, BderivsBatchAgreesWithTheReferenceLane) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto pts = sample_points(48, 21);
    std::vector<ib::BDerivatives> d(pts.size());
    const ib::Result<bool> r = ib::bderivs_batch(model, pts, d, 1.0e-4);
    ASSERT_EQ(r.status, ib::Status::Ok);
    EXPECT_FALSE(r.value);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const ib::BDerivatives one = ib::bderivs(model, pts[i], 1.0e-4).value;
        EXPECT_EQ(d[i].b, one.b);
        EXPECT_EQ(d[i].b_mag, one.b_mag);
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_EQ(d[i].grad_b_mag[j], one.grad_b_mag[j]);
            for (std::size_t c = 0; c < 3; ++c) { EXPECT_EQ(d[i].diff_b(c, j), one.diff_b(c, j)); }
        }
    }
    // And the auto step really is used when none is given: the two must differ, because the auto
    // step is 5e-8*r and not 1e-4.
    std::vector<ib::BDerivatives> automatic(pts.size());
    ASSERT_EQ(ib::bderivs_batch(model, pts, automatic).status, ib::Status::Ok);
    EXPECT_NE(automatic[0].grad_b_mag[0], d[0].grad_b_mag[0]);
}

// The batch form of the pure-algebra routine, and the per-point statuses that let one degenerate
// point report itself instead of spoiling the batch.
TEST(IrbemField, GradCurvCurlBatchMatchesThePointwiseRoutine) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto pts = sample_points(32, 31);
    std::vector<ib::BDerivatives> d(pts.size());
    ASSERT_EQ(ib::bderivs_batch(model, pts, d).status, ib::Status::Ok);
    std::vector<ib::GradCurvCurl> g(pts.size());
    std::vector<ib::Status> st(pts.size());
    EXPECT_EQ(ib::grad_curv_curl_batch(d, g, st), ib::Status::Ok);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        EXPECT_EQ(st[i], ib::Status::Ok);
        const ib::GradCurvCurl one = ib::grad_curv_curl(d[i]).value;
        EXPECT_EQ(g[i].grad_par, one.grad_par);
        EXPECT_EQ(g[i].r_curv, one.r_curv);
        EXPECT_EQ(g[i].div_b, one.div_b);
        for (std::size_t k = 0; k < 3; ++k) {
            EXPECT_EQ(g[i].curvature[k], one.curvature[k]);
            EXPECT_EQ(g[i].curl_b[k], one.curl_b[k]);
            EXPECT_EQ(g[i].curv_drift[k], one.curv_drift[k]);
            EXPECT_EQ(g[i].grad_perp[k], one.grad_perp[k]);
            EXPECT_EQ(g[i].grad_drift[k], one.grad_drift[k]);
        }
    }
    // One dead point must be reported, not hidden: the batch's status carries it and the point's
    // own slot names it.
    d[3] = ib::BDerivatives{};
    EXPECT_EQ(ib::grad_curv_curl_batch(d, g, st), ib::Status::DomainError);
    EXPECT_EQ(st[3], ib::Status::DomainError);
    EXPECT_EQ(st[4], ib::Status::Ok);
}

// The batch hemisphere is the pointwise hemisphere, point for point, on the host lane.
TEST(IrbemField, HemisphereBatchAgreesWithTheReferenceLane) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto pts = sample_points(64, 41);
    std::vector<ib::Hemisphere> h(pts.size());
    const ib::Result<bool> r = ib::hemisphere_batch(model, pts, h);
    ASSERT_EQ(r.status, ib::Status::Ok);
    EXPECT_FALSE(r.value);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        EXPECT_EQ(h[i], ib::hemisphere(model, pts[i]).value);
        EXPECT_NE(h[i], ib::Hemisphere::Invalid);
    }
}

// ---------------------------------------------------------------------------------------------
// Lane selection and the device
// ---------------------------------------------------------------------------------------------

// The registry's crossover is the KERNEL's, derived from an assumed ~30 us submit floor; this
// launcher's measured floor is ~115 us, so each routine carries its own. The ordering is the
// physics: four evaluations per point reach the device at a quarter of the batch size, and the
// routine that needs two dispatches sits between.
TEST(IrbemField, CrossoverIsTheRoutinesNotTheKernels) {
    EXPECT_EQ(ib::field_batch_crossover, 4 * ib::bderivs_batch_crossover);
    EXPECT_LT(ib::bderivs_batch_crossover, ib::hemisphere_batch_crossover);
    EXPECT_LT(ib::hemisphere_batch_crossover, ib::field_batch_crossover);
#ifdef CHEATAH_SPACE_IRBEM_FIELD_GPU
    // Without a device nothing prefers one, whatever the batch size -- so this half of the
    // assertion is meaningful on every machine.
    if (!device_ready()) {
        EXPECT_FALSE(ib::detail::prefer_device(1u << 20, ib::field_batch_crossover));
        return;
    }
    EXPECT_FALSE(ib::detail::prefer_device(8, ib::field_batch_crossover));
    EXPECT_TRUE(ib::detail::prefer_device(ib::field_batch_crossover, ib::field_batch_crossover));
    // THE BAND THAT MATTERS. Below the registry's own crossover every batch is rejected by
    // gpu::prefer_gpu before this function's threshold is ever consulted, so a test that only
    // probes 8 points proves nothing about the routine's own number -- deleting `points >=
    // crossover` entirely passed the whole suite before these three lines existed. 256 points is
    // above the registry's 128 and below field_batch's measured 512, and it is exactly where the
    // header says the device LOSES (0.63x measured): only the routine's threshold can refuse it.
    ASSERT_LT(ib::gpu::gpu_crossover("irbem_igrf_f32"), ib::field_batch_crossover)
        << "the band this asserts over must be non-empty";
    EXPECT_FALSE(ib::detail::prefer_device(256, ib::field_batch_crossover));
    EXPECT_FALSE(ib::detail::prefer_device(ib::field_batch_crossover - 1,
                                           ib::field_batch_crossover));
    EXPECT_TRUE(ib::detail::prefer_device(256, ib::hemisphere_batch_crossover));
    // ... and end to end, because a threshold nothing consults is a threshold that does not exist.
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto mid = sample_points(256, 97);
    std::vector<ib::FieldVector<Frame::GEO>> b(mid.size());
    std::vector<double> bm(mid.size());
    EXPECT_FALSE(ib::field_batch(model, mid, b, bm).value)
        << "256 points is inside the band where the device is measured at 0.63x";
    std::vector<ib::BDerivatives> d(mid.size());
    EXPECT_TRUE(ib::bderivs_batch(model, mid, d, 1.0e-3).value)
        << "but four evaluations per point put the SAME batch over bderivs' crossover";
    {   // An explicit operator override replaces BOTH thresholds, in both directions.
        const CrossoverOverride force("1");
        EXPECT_TRUE(ib::detail::prefer_device(1, ib::field_batch_crossover));
    }
    {
        const CrossoverOverride never("999999999");
        EXPECT_FALSE(ib::detail::prefer_device(1u << 20, ib::field_batch_crossover));
    }
#endif
}

// A test that means to exercise the device must assert the device ran. Everything else in this file
// would pass unchanged on a machine with no GPU, which is exactly why these three exist.
TEST(IrbemField, FieldBatchUsesTheDeviceWhenOneIsAvailable) {
    if (!device_ready()) GTEST_SKIP() << "no device or no compiled irbem_igrf_f32.spv";
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto pts = sample_points(4096, 51);
    std::vector<ib::FieldVector<Frame::GEO>> dev(pts.size());
    std::vector<double> dev_mag(pts.size());
    const ib::Result<bool> r = ib::field_batch(model, pts, dev, dev_mag);
    ASSERT_EQ(r.status, ib::Status::Ok);
    ASSERT_TRUE(r.value) << "4096 points is above every crossover; the device must have run";

    double worst = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        worst = std::max(worst, rel(dev_mag[i], model.evaluate(pts[i]).magnitude()));
    }
    // The fp32 lane is held to the Blocal budget of docs/ERROR_BUDGET.md §4, with the honest
    // margin the header records: 7.3e-7 over 2000 points, 1.1e-6 over 2^20. At 4096 points a 2e-6
    // cap is a real constraint and not a rubber stamp -- the measured value here is ~7e-7.
    EXPECT_LT(worst, 2.0e-6) << "fp32 device lane vs fp64 host lane";
    // Run to run, the same inputs must give the same bits. A device that did not would make every
    // tolerance above meaningless.
    std::vector<ib::FieldVector<Frame::GEO>> again(pts.size());
    std::vector<double> again_mag(pts.size());
    ASSERT_TRUE(ib::field_batch(model, pts, again, again_mag).value);
    for (std::size_t i = 0; i < pts.size(); ++i) EXPECT_EQ(again_mag[i], dev_mag[i]);

    // The truncation degree is a RUNTIME argument to the kernel (dims[1]), not a compile-time one,
    // so a kernel that ignored it would return the degree-13 field for a degree-10 model and every
    // test above would still pass. Running IRBEM's own truncation through the device and requiring
    // it to differ — by the ~1e-5 the two models differ by on the host — is what pins that.
    const OracleIgrf truncated = *OracleIgrf::at(kEpoch);
    std::vector<ib::FieldVector<Frame::GEO>> dev10(pts.size());
    std::vector<double> dev10_mag(pts.size());
    ASSERT_TRUE(ib::field_batch(truncated, pts, dev10, dev10_mag).value);
    double model_gap = 0.0;
    double against_host = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        model_gap = std::max(model_gap, rel(dev10_mag[i], dev_mag[i]));
        against_host = std::max(against_host, rel(dev10_mag[i], truncated.evaluate(pts[i]).magnitude()));
    }
    EXPECT_GT(model_gap, 1.0e-6) << "the device ignored the truncation degree";
    EXPECT_LT(against_host, 2.0e-6) << "and it must still match ITS OWN host lane";
}

// The derivative device lane, at a MATCHED step. The two lanes' automatic steps differ by four
// orders of magnitude on purpose, so comparing them without pinning the step would be comparing
// resolutions -- the same discipline docs/ERROR_BUDGET.md §2(a) imposes on L*.
TEST(IrbemField, BderivsBatchUsesTheDeviceWhenOneIsAvailable) {
    if (!device_ready()) GTEST_SKIP() << "no device or no compiled irbem_igrf_f32.spv";
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto pts = sample_points(4096, 61);
    const double step = 2.0e-3;   // in the truncation-dominated regime, where both lanes agree
    std::vector<ib::BDerivatives> dev(pts.size());
    const ib::Result<bool> r = ib::bderivs_batch(model, pts, dev, step);
    ASSERT_EQ(r.status, ib::Status::Ok);
    ASSERT_TRUE(r.value);
    std::vector<ib::BDerivatives> host(pts.size());
    {
        const CrossoverOverride never("999999999");
        ASSERT_FALSE(ib::bderivs_batch(model, pts, host, step).value);
    }
    double worst = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const double scale = fx::norm(host[i].grad_b_mag);
        worst = std::max(worst, fx::norm(dev[i].grad_b_mag - host[i].grad_b_mag) / scale);
    }
    // At dX = 2e-3 truncation dominates on both lanes, so they agree to ~1e-3 -- the file brief's
    // table says 3.5e-3 at 1e-3 and better at larger steps. This is the number that would blow up
    // if the device stopped recovering the step it actually took.
    EXPECT_LT(worst, 5.0e-3) << "fp32 vs fp64 derivative lane at a matched step";
}

// Only the SIGN survives here, so fp32 costs almost nothing -- and what it does cost is points
// within a step of the magnetic equator, where the two hemispheres are genuinely adjacent.
TEST(IrbemField, HemisphereBatchUsesTheDeviceWhenOneIsAvailable) {
    if (!device_ready()) GTEST_SKIP() << "no device or no compiled irbem_igrf_f32.spv";
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto pts = sample_points(4096, 71);
    std::vector<ib::Hemisphere> dev(pts.size());
    const ib::Result<bool> r = ib::hemisphere_batch(model, pts, dev);
    ASSERT_EQ(r.status, ib::Status::Ok);
    ASSERT_TRUE(r.value);
    std::size_t disagree = 0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (dev[i] != ib::hemisphere(model, pts[i]).value) ++disagree;
    }
    // Measured over 2^18 points: 11 disagreements, 4e-5. A cap of one in a thousand is loose enough
    // to survive a different driver and tight enough that a broken sign convention fails it.
    EXPECT_LE(disagree, pts.size() / 1000) << disagree << " of " << pts.size();
}

// A registered kernel whose SPIR-V was never built is a DEPLOYMENT problem, not a reason to refuse
// to compute — the header says so and this is what makes that true. Pointing the shader directory
// at nothing must silently produce the host lane's answer, bit for bit, not an exception and not a
// buffer of zeros.
TEST(IrbemField, AMissingShaderFallsBackToTheHostLane) {
#ifndef CHEATAH_SPACE_IRBEM_FIELD_GPU
    GTEST_SKIP() << "built without cheatah-gpu-linalg";
#else
    if (!device_ready()) GTEST_SKIP() << "no device";
    const FullIgrf model = *FullIgrf::at(kEpoch);
    const auto pts = sample_points(1024, 81);
    std::vector<ib::FieldVector<Frame::GEO>> expected(pts.size());
    std::vector<double> expected_mag(pts.size());
    {
        const CrossoverOverride never("999999999");
        ASSERT_FALSE(ib::field_batch(model, pts, expected, expected_mag).value);
    }

    const CrossoverOverride force("1");
    const ShaderDirOverride nowhere("/nonexistent/space-irbem-shaders");
    ASSERT_FALSE(device_ready()) << "the override must make the kernel unresolvable";
    std::vector<ib::FieldVector<Frame::GEO>> got(pts.size());
    std::vector<double> got_mag(pts.size());
    const ib::Result<bool> r = ib::field_batch(model, pts, got, got_mag);
    EXPECT_EQ(r.status, ib::Status::Ok);
    EXPECT_FALSE(r.value) << "with no shader the device lane must report that it did not run";
    for (std::size_t i = 0; i < pts.size(); ++i) EXPECT_EQ(got_mag[i], expected_mag[i]);

    std::vector<ib::BDerivatives> d(pts.size());
    EXPECT_FALSE(ib::bderivs_batch(model, pts, d, 1.0e-3).value);
    EXPECT_GT(d[0].b_mag, 0.0);
    std::vector<ib::Hemisphere> h(pts.size());
    EXPECT_FALSE(ib::hemisphere_batch(model, pts, h).value);
    EXPECT_NE(h[0], ib::Hemisphere::Invalid);
#endif
}

// The two derivative lanes must refuse the same point. They did not: before the guard in
// bderivs_batch's device lane, a batch containing the origin came back with b_mag = 0 and
// grad_b_mag = 0 on the host and b_mag = 0 with grad_b_mag = NaN on the device -- the device
// differencing a field that is infinite one step inside the Earth. A NaN that only appears above
// the crossover is the worst shape a bug can have, because the small batches a test writes are
// exactly the ones that never see it.
TEST(IrbemField, BothDerivativeLanesRefuseTheSameDegeneratePoint) {
    const FullIgrf model = *FullIgrf::at(kEpoch);
    auto pts = sample_points(1024, 93);
    pts[17] = ib::Position<Frame::GEO>{fx::vec3d{0.0, 0.0, 0.0}};
    std::vector<ib::BDerivatives> host(pts.size());
    {
        const CrossoverOverride never("999999999");
        ASSERT_FALSE(ib::bderivs_batch(model, pts, host, 1.0e-3).value);
    }
    EXPECT_EQ(host[17].b_mag, 0.0);
    EXPECT_EQ(host[17].grad_b_mag[0], 0.0);
    EXPECT_GT(host[18].b_mag, 0.0);

    if (!device_ready()) GTEST_SKIP() << "no device or no compiled irbem_igrf_f32.spv";
    std::vector<ib::BDerivatives> dev(pts.size());
    ASSERT_TRUE(ib::bderivs_batch(model, pts, dev, 1.0e-3).value);
    EXPECT_EQ(dev[17].b_mag, 0.0);
    for (std::size_t j = 0; j < 3; ++j) {
        EXPECT_EQ(dev[17].grad_b_mag[j], 0.0) << "component " << j;
        for (std::size_t k = 0; k < 3; ++k) { EXPECT_EQ(dev[17].diff_b(k, j), 0.0); }
    }
    // ... and the guard must not have swallowed the rest of the batch with it.
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (i == 17) continue;
        ASSERT_TRUE(std::isfinite(dev[i].b_mag)) << "point " << i;
        EXPECT_GT(dev[i].b_mag, 0.0) << "point " << i;
        EXPECT_TRUE(std::isfinite(dev[i].grad_b_mag[0])) << "point " << i;
    }
}

// The origin is the one point where the series has no answer, and the device lane must report it
// the way the host lane does — as IRBEM's 0 — rather than as a NaN that propagates into whatever
// the caller does next. Buried in a batch of good points, because that is how it would actually
// arrive.
TEST(IrbemField, TheDeviceLaneReportsAPointWithNoField) {
    if (!device_ready()) GTEST_SKIP() << "no device or no compiled irbem_igrf_f32.spv";
    const FullIgrf model = *FullIgrf::at(kEpoch);
    auto pts = sample_points(1024, 91);
    pts[17] = ib::Position<Frame::GEO>{fx::vec3d{0.0, 0.0, 0.0}};
    std::vector<ib::Hemisphere> h(pts.size());
    const ib::Result<bool> r = ib::hemisphere_batch(model, pts, h);
    ASSERT_EQ(r.status, ib::Status::Ok);
    ASSERT_TRUE(r.value);
    EXPECT_EQ(h[17], ib::Hemisphere::Invalid);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (i != 17) { EXPECT_NE(h[i], ib::Hemisphere::Invalid) << "point " << i; }
    }
}

// A model defined in a test file is still separately-compiled code: llvm-cov instantiates
// `bderivs<AnalyticDipole>` and `hemisphere<AnalyticDipole>` independently of the `Igrf` ones, so
// guards driven only through `Igrf` above leave this instantiation's copies untested. The analytic
// model is not a throwaway either — it is the closed-form truth the Jacobian test measures against,
// so a guard misbehaving here would corrupt the reference rather than the subject.
TEST(IrbemField, TheAnalyticModelDrivesTheSameGuardsAndTheSameBody) {
    const ib::Igrf<13> model = ib::Igrf<13>::at(2015.0).value();
    const AnalyticDipole dip = AnalyticDipole::of(model);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const ib::Position<ib::Frame::GEO> good{fx::vec3d{4.0, 1.0, 0.5}};
    const ib::Position<ib::Frame::GEO> origin{fx::vec3d{0.0, 0.0, 0.0}};

    EXPECT_EQ(ib::Status::DomainError, ib::bderivs(dip, origin).status);
    EXPECT_EQ(ib::Status::DomainError, ib::bderivs(dip, good, nan).status);
    EXPECT_EQ(ib::Status::DomainError, ib::bderivs(dip, good, -1.0).status);
    EXPECT_EQ(ib::Status::Ok, ib::bderivs(dip, good, 1.0e-3).status);

    EXPECT_EQ(ib::Status::DomainError, ib::hemisphere(dip, origin).status);
    EXPECT_EQ(ib::Status::DomainError, ib::hemisphere(dip, good, -1.0).status);

    // The BODY, not just the guards: an explicit step and the automatic one, north and south of the
    // magnetic equator. For a centred dipole |B| grows away from the equator along a field line, so
    // the two hemispheres must come back opposite — which also shows the slope sign is not
    // arbitrarily fixed.
    const ib::Position<ib::Frame::GEO> north{fx::vec3d{2.0, 0.0, 2.0}};
    const ib::Position<ib::Frame::GEO> south{fx::vec3d{2.0, 0.0, -2.0}};
    const auto hn = ib::hemisphere(dip, north);
    const auto hs = ib::hemisphere(dip, south);
    EXPECT_EQ(ib::Status::Ok, hn.status);
    EXPECT_EQ(ib::Status::Ok, hs.status);
    EXPECT_NE(hn.value, hs.value) << "opposite sides of the equator must not agree";
    EXPECT_EQ(ib::Status::Ok, ib::hemisphere(dip, north, 1.0e-3).status);
}
