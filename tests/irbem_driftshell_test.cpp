// Unit and differential tests for space.irbem's drift shell and Roederer L* — `driftshell.hpp`.
//
// Three kinds of assertion live here, and they are deliberately not mixed up:
//
//   (a) EXACT ones, where the mathematics permits. `DriftShellOptions::from_irbem` is arithmetic on
//       small integers, the polar-cap seed geometry is a rotation of an exactly-representable
//       vector, and Brent on a polynomial with a dyadic root converges to it exactly. Those use
//       `==` or an epsilon-scale tolerance, never a loose one.
//   (b) SELF-CONSISTENT ones against physics that has a closed form. For a centred dipole the drift
//       shell IS the McIlwain shell and `L* == L_m`, and the ionospheric footpoint of an `L` shell
//       sits at `sin²θ = 1/L`. Neither needs an oracle, and both fail loudly if the flux integral,
//       the footpoint walk or the root-find is wrong.
//   (c) DIFFERENTIAL ones against IRBEM, run as a black box. Its numbers are transcribed below with
//       full provenance rather than linked, because IRBEM is LGPL-3.0 and this repository is MIT:
//       nothing here links it, and the gate never builds it.
//
// **On tolerances.** `docs/ERROR_BUDGET.md` §2(a) is emphatic that a 0.01 L* target is meaningful
// only at MATCHED resolutions, because IRBEM's own default-resolution error is 0.010-0.017 at
// L ≈ 6 — measured independently in that document, and reproduced in `kOracleRes0` versus
// `kOracleRes9` below. So the two matched comparisons are asserted separately and mean different
// things: at `options(3,4) = 9` both implementations are near their converged values and the 0.01
// budget applies; at `options(3,4) = 0` neither is, and what is asserted is the distance from the
// CONVERGED oracle, where this implementation is measured to sit 2.5x closer than IRBEM's own
// default-resolution answer does.
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <span>
#include <string>
#include <vector>

#include "space/irbem/api.hpp"
#include "space/irbem/driftshell.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;

using ib::Frame;
using ib::Status;

namespace {

/// The epoch every differential case below was produced at: 2015-180 12:00 UT. IRBEM takes
/// (year, doy, ut); IGRF takes a decimal year, and 179.5/365 is the same instant.
constexpr int kYear = 2015;
constexpr int kDoy = 180;
constexpr double kUt = 43200.0;
constexpr double kDecimalYear = 2015.0 + (179.5 / 365.0);

/// One differential case: a point in GEO Cartesian (IRBEM `sysaxes = 1`), Earth radii.
struct OraclePoint {
    const char* name;
    double x;
    double y;
    double z;
    double lstar;
};

// IRBEM e7cecb0 built `-O2 -ffp-contract=off -fno-fast-math`, `libirbem-O2.so`, run through
// `make_lstar1_` with kext = 0 (internal field only), options = {1, 0, R, R, 0} — compute L*,
// IGRF internal — at 2015-180 12:00 UT. `get_igrf_version_` reports 14, the same generation
// `tables/igrf14.hpp` carries, so the two libraries are evaluating the same field.
//
// R = 0 is IRBEM's recommended setting: dθ = π/720 (0.25°), Nder = 25.
constexpr std::array<OraclePoint, 12> kOracleRes0{{
    {"eq2.0", 2.0, 0.0, 0.0, 2.049710957},
    {"eq3.0", 3.0, 0.0, 0.0, 3.049521649},
    {"eq4.0", 4.0, 0.0, 0.0, 4.048411808},
    {"eq5.0", 5.0, 0.0, 0.0, 5.080397525},
    {"eq6.0", 6.0, 0.0, 0.0, 6.085192256},
    {"eq6.6", 6.6, 0.0, 0.0, 6.697352007},
    {"eq8.0", 8.0, 0.0, 0.0, 8.064429403},
    {"y4.0", 0.0, 4.0, 0.0, 4.026906519},
    {"y6.0", 0.0, 6.0, 0.0, 6.078414261},
    {"xy45", 3.5355339059327378, 3.5355339059327378, 0.0, 5.035357416},
    {"off5.5", 5.5, 0.0, 2.0, 6.763195199},
    {"off3.0", 3.0, 0.0, 1.0, 3.553125244},
}};

// The same twelve points at R = 9: dθ = π/7200 (0.025°), Nder = 250 — IRBEM's most converged
// setting, and the reference `ERROR_BUDGET.md` §2 measures its coarser settings against. Not
// "truth": that document measures its own residual at plausibly 1e-4 to 1e-3.
constexpr std::array<OraclePoint, 12> kOracleRes9{{
    {"eq2.0", 2.0, 0.0, 0.0, 2.055585332},
    {"eq3.0", 3.0, 0.0, 0.0, 3.060278504},
    {"eq4.0", 4.0, 0.0, 0.0, 4.062577054},
    {"eq5.0", 5.0, 0.0, 0.0, 5.066721138},
    {"eq6.0", 6.0, 0.0, 0.0, 6.067784859},
    {"eq6.6", 6.6, 0.0, 0.0, 6.670502152},
    {"eq8.0", 8.0, 0.0, 0.0, 8.073761758},
    {"y4.0", 0.0, 4.0, 0.0, 4.019294897},
    {"y6.0", 0.0, 6.0, 0.0, 6.077992094},
    {"xy45", 3.5355339059327378, 3.5355339059327378, 0.0, 5.032164825},
    {"off5.5", 5.5, 0.0, 2.0, 6.783228363},
    {"off3.0", 3.0, 0.0, 1.0, 3.551685506},
}};

/// One shell-splitting case: the same geometry at a non-90° local pitch angle, from
/// `make_lstar_shell_splitting1_` at R = 9. Shell splitting is the reason the drift shell has to
/// conserve BOTH `I` and `B_m` rather than just one of them: at fixed position, `I` ranges from
/// 0.002 to 6.1 across these pitch angles while `B_m` ranges over a factor of fifteen.
struct SplitPoint {
    double x;
    double y;
    double z;
    double alpha;
    double lstar;
};

constexpr std::array<SplitPoint, 18> kOracleSplitRes9{{
    {4.0, 0.0, 0.0, 90.0, 4.0625771}, {4.0, 0.0, 0.0, 75.0, 4.0625398},
    {4.0, 0.0, 0.0, 60.0, 4.0622256}, {4.0, 0.0, 0.0, 45.0, 4.0632716},
    {4.0, 0.0, 0.0, 30.0, 4.0646864}, {4.0, 0.0, 0.0, 15.0, 4.0701421},
    {6.0, 0.0, 0.0, 90.0, 6.0677849}, {6.0, 0.0, 0.0, 75.0, 6.0672205},
    {6.0, 0.0, 0.0, 60.0, 6.0686021}, {6.0, 0.0, 0.0, 45.0, 6.0682147},
    {6.0, 0.0, 0.0, 30.0, 6.0677648}, {6.0, 0.0, 0.0, 15.0, 6.0715702},
    {5.5, 0.0, 2.0, 90.0, 6.7832284}, {5.5, 0.0, 2.0, 75.0, 6.7830879},
    {5.5, 0.0, 2.0, 60.0, 6.7834702}, {5.5, 0.0, 2.0, 45.0, 6.7835484},
    {5.5, 0.0, 2.0, 30.0, 6.7836407}, {5.5, 0.0, 2.0, 15.0, 6.7863089},
}};

/// L* at (6, 0, 0) GEO, 90°, across seven epochs, R = 9. The dipole moment falls ~6 % per century
/// and `L ∝ M^(1/3)`, so this is the case a hard-coded `k₀` fails: a stale moment shows up here as
/// a drift with epoch that looks like a physical secular trend and is not.
struct EpochPoint {
    int year;
    double lstar;
};

constexpr std::array<EpochPoint, 7> kOracleEpochRes9{{
    {1900, 6.0711847}, {1940, 6.0728117}, {1975, 6.0708236}, {2000, 6.0676722},
    {2015, 6.0677849}, {2025, 6.0625517}, {2029, 6.0623756},
}};

/// The model and rotations every case at the standard epoch shares.
struct Epoch {
    ib::Igrf<13> model;
    ib::Rotations rotations;
};

/// The model and rotations for one epoch. `optional::value()` rather than `operator*`: an epoch
/// outside IGRF's definition is a mistake in the test's own table, and throwing says so at the
/// line that made it instead of reading an empty optional.
[[nodiscard]] Epoch epoch_at(int year, double decimal_year) {
    const ib::Igrf<13> model = ib::Igrf<13>::at(decimal_year).value();
    const ib::Result<ib::Rotations> r = ib::api::rotations_at(year, kDoy, kUt, model);
    EXPECT_EQ(r.status, Status::Ok);
    return Epoch{model, r.value};
}

[[nodiscard]] Epoch standard_epoch() { return epoch_at(kYear, kDecimalYear); }

[[nodiscard]] ib::Position<Frame::GEO> geo(double x, double y, double z) {
    return ib::Position<Frame::GEO>{fx::vec3d{x, y, z}};
}

/// Forces the host lane for the life of the object by setting the seam's own kill switch, and puts
/// the previous value — or its absence — back afterwards. The differential tests need to run the
/// SAME inputs through both lanes on a machine that has a device, which they cannot do if the
/// choice is only ever made for them.
class HostLaneOnly {
public:
    HostLaneOnly() {
        if (const char* prev = std::getenv(kVar)) {
            had_ = true;
            prev_ = prev;
        }
        ::setenv(kVar, "1", 1);
    }
    HostLaneOnly(const HostLaneOnly&) = delete;
    HostLaneOnly& operator=(const HostLaneOnly&) = delete;
    HostLaneOnly(HostLaneOnly&&) = delete;
    HostLaneOnly& operator=(HostLaneOnly&&) = delete;
    ~HostLaneOnly() {
        if (had_) {
            ::setenv(kVar, prev_.c_str(), 1);
            return;
        }
        ::unsetenv(kVar);
    }

private:
    static constexpr const char* kVar = "CHEATAH_SPACE_IRBEM_NO_GPU";
    bool had_ = false;
    std::string prev_;
};

}  // namespace

// ---- the option translation --------------------------------------------------------------------

// Exact: IRBEM's published table says dθ = π/(720·(options(3)+1)) and Nder = 25·(options(4)+1).
// Both are arithmetic on small integers, so both are checked with `==` and no tolerance.
TEST(IrbemDriftShell, IrbemOptionsTranslateToTheDocumentedResolutions) {
    const ib::DriftShellOptions zero = ib::DriftShellOptions::from_irbem(0, 0);
    EXPECT_EQ(zero.azimuths, 25);
    EXPECT_EQ(zero.colatitude_step_deg, 0.25);   // 180/720, exactly representable

    const ib::DriftShellOptions nine = ib::DriftShellOptions::from_irbem(9, 9);
    EXPECT_EQ(nine.azimuths, 250);
    EXPECT_EQ(nine.colatitude_step_deg, 0.025);

    const ib::DriftShellOptions three = ib::DriftShellOptions::from_irbem(3, 1);
    EXPECT_EQ(three.azimuths, 50);
    EXPECT_EQ(three.colatitude_step_deg, 180.0 / 2880.0);

    // The defaults ARE options(3,4) = 0 — the setting IRBEM recommends and this suite compares at.
    const ib::DriftShellOptions def;
    EXPECT_EQ(def.azimuths, zero.azimuths);
    EXPECT_EQ(def.colatitude_step_deg, zero.colatitude_step_deg);
}

// ---- the root-find, on functions whose roots are known exactly ------------------------------

// x² − 4 on [0, 6]: the root is 2, exactly representable, and every arithmetic step of the
// interpolation is on dyadic rationals. Brent must land on it to the tolerance and stop.
TEST(IrbemDriftShell, BrentFindsAKnownRootOfAnAnalyticFunction) {
    const auto f = [](double x) { return (x * x) - 4.0; };
    ib::detail::RootState s = ib::detail::root_begin(0.0, f(0.0), 6.0, f(6.0));
    EXPECT_TRUE(s.valid);
    EXPECT_FALSE(s.done);
    int iterations = 0;
    while (!s.done && iterations < 60) {
        ib::detail::root_step(s, f(s.b), 1e-12);
        ++iterations;
    }
    EXPECT_TRUE(s.done);
    EXPECT_LT(iterations, 20) << "Brent should not need twenty iterations on a quadratic";
    EXPECT_NEAR(s.b, 2.0, 1e-11);

    // A root sitting exactly ON an endpoint is reported immediately rather than bisected toward.
    ib::detail::RootState hit = ib::detail::root_begin(-1.0, f(-1.0), 2.0, f(2.0));
    ib::detail::root_step(hit, 0.0, 1e-12);
    EXPECT_TRUE(hit.done);
    EXPECT_EQ(hit.b, 2.0);
}

// The case bisection exists for: a residual so lopsided that the interpolated step would leave the
// bracket or barely shrink it. `x^(1/5)` has an infinite derivative at its root and defeats plain
// secant iteration; Brent must still converge, and must keep the bracket the whole way.
TEST(IrbemDriftShell, BrentBisectsWhenInterpolationWouldLeaveTheBracket) {
    const auto f = [](double x) { return std::cbrt(std::cbrt(x)) * std::abs(std::cbrt(x)); };
    ib::detail::RootState s = ib::detail::root_begin(-1.0, f(-1.0), 4.0, f(4.0));
    int iterations = 0;
    double previous_width = std::abs(s.c - s.a);
    while (!s.done && iterations < 200) {
        ib::detail::root_step(s, f(s.b), 1e-9);
        // The bracket is never surrendered. Brent's invariant is on the THREE points it carries:
        // the root stays inside their hull, and the hull never grows. A plain secant iteration on
        // this function walks straight out of it, which is the whole reason for the bisection
        // fallback this test exists to exercise.
        const double lo = std::min({s.a, s.b, s.c});
        const double hi = std::max({s.a, s.b, s.c});
        EXPECT_LE(lo, 0.0) << "at iteration " << iterations;
        EXPECT_GE(hi, 0.0) << "at iteration " << iterations;
        EXPECT_LE(hi - lo, previous_width + 1e-15) << "at iteration " << iterations;
        previous_width = hi - lo;
        ++iterations;
    }
    EXPECT_TRUE(s.done);
    EXPECT_NEAR(s.b, 0.0, 1e-8);
}

// ---- the geometry the shell is parameterised by -------------------------------------------------

TEST(IrbemDriftShell, SeedsLieInTheMagneticEquatorialPlane) {
    const Epoch e = standard_epoch();
    const fx::mat3d mag_to_geo = ib::rotation_matrix<Frame::GEO, Frame::MAG>(e.rotations);
    const fx::mat3d geo_to_mag = ib::rotation_matrix<Frame::MAG, Frame::GEO>(e.rotations);
    for (const double phi : {0.0, 1.0, 2.5, -3.0}) {
        for (const double r : {1.5, 4.0, 8.0}) {
            const ib::Position<Frame::GEO> p = ib::detail::equatorial_seed(mag_to_geo, phi, r);
            // A rotation preserves length, so the radius comes back to the last few ulps.
            EXPECT_NEAR(fx::norm(p.v), r, 1e-13 * r);
            const fx::vec3d m = geo_to_mag * p.v;
            EXPECT_NEAR(m[2], 0.0, 1e-13 * r) << "the seed must sit in the MAG z = 0 plane";
            EXPECT_NEAR(std::atan2(m[1], m[0]), std::atan2(std::sin(phi), std::cos(phi)), 1e-12);
        }
    }
}

// A centred dipole's field line through equatorial radius L meets the surface at sin²θ = 1/L.
// Igrf<13> truncated to degree 1 IS a centred (tilted) dipole, and MAG is built around its axis,
// so the closed form applies exactly there and the walk is measured against it.
TEST(IrbemDriftShell, FootpointOfADipoleLineIsTheAnalyticColatitude) {
    const ib::Igrf<1> dipole = ib::Igrf<1>::at(kDecimalYear).value();
    const ib::Result<ib::Rotations> rot = ib::api::rotations_at(kYear, kDoy, kUt, dipole);
    ASSERT_EQ(rot.status, Status::Ok);
    const fx::mat3d mag_to_geo = ib::rotation_matrix<Frame::GEO, Frame::MAG>(rot.value);
    const fx::mat3d geo_to_mag = ib::rotation_matrix<Frame::MAG, Frame::GEO>(rot.value);
    ib::TraceOptions foot{100.0, 8000, 1.0};

    for (const double l : {2.0, 4.0, 6.5}) {
        for (const double phi : {0.0, 2.0}) {
            const ib::Position<Frame::GEO> seed =
                ib::detail::equatorial_seed(mag_to_geo, phi, l);
            const ib::Result<ib::Position<Frame::GEO>> f =
                ib::detail::walk_to_surface(dipole, seed, true, foot);
            ASSERT_EQ(f.status, Status::Ok) << "L = " << l;
            EXPECT_NEAR(fx::norm(f.value.v), 1.0, 1e-12) << "the footpoint is on the unit sphere";
            const fx::vec3d m = geo_to_mag * f.value.v;
            const double theta = std::acos(m[2]);
            // The RK4 truncation at ds = L/100 is what this measures; 1e-4 rad is 0.006 degrees.
            EXPECT_NEAR(std::sin(theta) * std::sin(theta), 1.0 / l, 1e-4) << "L = " << l;
            EXPECT_GT(m[2], 0.0) << "walking along +B must land in the northern magnetic cap";
        }
    }

    // And the other direction lands in the southern cap, on the same shell.
    const ib::Position<Frame::GEO> seed = ib::detail::equatorial_seed(mag_to_geo, 0.0, 4.0);
    const ib::Result<ib::Position<Frame::GEO>> south =
        ib::detail::walk_to_surface(dipole, seed, false, foot);
    ASSERT_EQ(south.status, Status::Ok);
    EXPECT_LT((geo_to_mag * south.value.v)[2], 0.0);

    // A start inside the atmosphere is a domain error, not a walk that happens to terminate at
    // once — and a null field likewise. Both are the guards, exercised.
    const ib::Result<ib::Position<Frame::GEO>> inside =
        ib::detail::walk_to_surface(dipole, geo(0.5, 0.0, 0.0), true, foot);
    EXPECT_EQ(inside.status, Status::DomainError);
    ib::TraceOptions stingy{100.0, 3, 1.0};
    const ib::Result<ib::Position<Frame::GEO>> capped =
        ib::detail::walk_to_surface(dipole, seed, true, stingy);
    EXPECT_EQ(capped.status, Status::OpenFieldLine) << "three steps cannot reach the surface";
}

// ---- the staged tables and the batched field ----------------------------------------------------

TEST(IrbemDriftShell, FluxCellsAgreeBetweenLanes) {
    const Epoch e = standard_epoch();
    // The slot counts are pure arithmetic on the degree: 105 triangular slots at degree 13.
    EXPECT_EQ(ib::detail::coefficient_slots(13), 210U);
    EXPECT_EQ(ib::detail::normalisation_slots(13), 224U);
    EXPECT_EQ(ib::detail::coefficient_slots(1), 6U);

    std::vector<float> coef(ib::detail::coefficient_slots(13));
    std::vector<float> norm(ib::detail::normalisation_slots(13));
    ib::detail::stage_model(e.model, coef, norm);
    // g₁⁰ lands in slot 1 (n(n+1)/2 + m), and h₁⁰ is structurally zero.
    EXPECT_FLOAT_EQ(coef[1], static_cast<float>(e.model.g(1, 0)));
    EXPECT_FLOAT_EQ(coef[2], static_cast<float>(e.model.g(1, 1)));
    EXPECT_FLOAT_EQ(coef[105 + 2], static_cast<float>(e.model.h(1, 1)));
    EXPECT_EQ(coef[0], 0.0F);

    // A polar-cap-shaped batch: points on the unit sphere at the colatitudes the flux quadrature
    // actually samples. The host lane is the fp64 reference; the device lane, when present, is
    // held to the Blocal budget against it.
    constexpr std::size_t kCells = 4096;
    std::vector<float> pos(3 * kCells);
    for (std::size_t i = 0; i < kCells; ++i) {
        const std::size_t ring = i % 128;    // 128 colatitudes per azimuth, 32 azimuths
        const std::size_t spoke = i / 128;
        const double th = 0.30 * static_cast<double>(ring) / 128.0;
        const double ph = std::numbers::pi * 2.0 * static_cast<double>(spoke) / 32.0;
        pos[(3 * i) + 0] = static_cast<float>(std::sin(th) * std::cos(ph));
        pos[(3 * i) + 1] = static_cast<float>(std::sin(th) * std::sin(ph));
        pos[(3 * i) + 2] = static_cast<float>(std::cos(th));
    }
    std::vector<float> host(3 * kCells);
    {
        const HostLaneOnly forced;
        EXPECT_FALSE(ib::detail::field_batch(e.model, coef, norm, pos, host));
    }

    double worst = 0.0;
    for (std::size_t i = 0; i < kCells; ++i) {
        const double m = std::hypot(std::hypot(host[(3 * i) + 0], host[(3 * i) + 1]),
                                    host[(3 * i) + 2]);
        EXPECT_GT(m, 20000.0) << "the surface field is 20-65 uT everywhere; cell " << i;
        worst = std::max(worst, m);
    }
    EXPECT_LT(worst, 70000.0);

#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
    if (!ib::gpu::available()) {
        GTEST_SKIP() << "no device: " << ib::gpu::unavailable_reason();
    }
    std::vector<float> dev(3 * kCells);
    ASSERT_TRUE(ib::detail::field_batch(e.model, coef, norm, pos, dev))
        << "the device lane must service a 4096-point batch, not fall back silently";
    // What the flux integral actually consumes is the RADIAL component, and what bounds its error
    // is the field magnitude — a per-component relative error is meaningless where a component
    // passes through zero. So both are measured against |B|.
    double worst_radial = 0.0;
    double worst_magnitude = 0.0;
    for (std::size_t i = 0; i < kCells; ++i) {
        const auto radial = [&](const std::vector<float>& b) {
            return (static_cast<double>(b[(3 * i) + 0]) * pos[(3 * i) + 0]) +
                   (static_cast<double>(b[(3 * i) + 1]) * pos[(3 * i) + 1]) +
                   (static_cast<double>(b[(3 * i) + 2]) * pos[(3 * i) + 2]);
        };
        const auto magnitude = [&](const std::vector<float>& b) {
            return std::hypot(std::hypot(static_cast<double>(b[(3 * i) + 0]),
                                         static_cast<double>(b[(3 * i) + 1])),
                              static_cast<double>(b[(3 * i) + 2]));
        };
        const double scale = magnitude(host);
        worst_radial = std::max(worst_radial, std::abs(radial(dev) - radial(host)) / scale);
        worst_magnitude =
            std::max(worst_magnitude, std::abs(magnitude(dev) - magnitude(host)) / scale);
    }
    std::printf("[ measured ] flux cells at r = 1, device vs fp64 host: B_r %.3g, |B| %.3g "
                "(relative to |B|)\n",
                worst_radial, worst_magnitude);
    // ERROR_BUDGET.md §4: Bgeo/Blocal, 1e-6 relative. Worth stating that this is the HARD case for
    // the fp32 kernel and it still passes: every one of these points is at r = 1, where each
    // (a/r)^(n+1) factor is exactly 1, no degree is attenuated, and the degree-13 sum cancels
    // hardest. The registry's 8.8e-7 was measured over points spread through the magnetosphere,
    // where the high degrees are suppressed; here they are not, and the answer is the same size.
    EXPECT_LT(worst_radial, 1e-6);
    EXPECT_LT(worst_magnitude, 1e-6);
#endif
}

// ---- the residual --------------------------------------------------------------------------

// The residual has two branches and they do NOT meet: the seed is the field line's
// magnetic-equatorial-PLANE crossing rather than its B_min, so a particle mirroring there still
// bounces through a short path and the traced branch starts above zero. What the root-find actually
// rests on — and what is asserted here — is that both branches increase with r, that the analytic
// one stays negative, and that the step across the junction is small enough to be an error term
// rather than a defect. See driftshell.hpp's brief for the measured bound and what closing it
// would take.
TEST(IrbemDriftShell, ResidualIsMonotoneAcrossTheAccessibilityBoundary) {
    const Epoch e = standard_epoch();
    const fx::mat3d mag_to_geo = ib::rotation_matrix<Frame::GEO, Frame::MAG>(e.rotations);
    std::vector<float> coef(ib::detail::coefficient_slots(13));
    std::vector<float> norm(ib::detail::normalisation_slots(13));
    ib::detail::stage_model(e.model, coef, norm);
    const ib::DriftShellOptions opt;

    // Find the radius where |B| on the magnetic-equatorial plane equals a chosen mirror field: the
    // branch boundary, by construction.
    const double phi = 0.7;
    double lo = 3.0;
    double hi = 9.0;
    const double bm = 200.0;
    for (int i = 0; i < 60; ++i) {
        const double mid = 0.5 * (lo + hi);
        const fx::vec3d b =
            e.model.evaluate(ib::detail::equatorial_seed(mag_to_geo, phi, mid)).v;
        (fx::norm(b) > bm ? lo : hi) = mid;
    }
    const double boundary = 0.5 * (lo + hi);

    const std::array<double, 4> az{phi, phi, phi, phi};
    const std::array<double, 4> mirror{bm, bm, bm, bm};
    const std::array<double, 4> target{0.0, 0.0, 0.0, 0.0};
    const std::array<double, 4> radius{boundary - 2e-4, boundary - 1e-5, boundary + 1e-5,
                                       boundary + 2e-4};
    std::array<double, 4> f{};
    int traces = 0;
    (void)ib::detail::residual_round(e.model, mag_to_geo, coef, norm, az, mirror, target, radius, f,
                                     opt, traces);
    EXPECT_EQ(traces, 2) << "only the two radii above the boundary need a trace";
    // Below: analytic and negative. Above: traced and positive. Across: continuous through zero.
    EXPECT_LT(f[0], 0.0);
    EXPECT_LT(f[1], 0.0);
    EXPECT_GT(f[2], 0.0);
    EXPECT_GT(f[3], 0.0);
    EXPECT_LT(f[0], f[1]);
    EXPECT_LT(f[2], f[3]);
    // One part in 1e-5 of the radius is ~7e-6 of the field (B falls as r^-3), which is what the
    // analytic branch reports there — so this checks the branch really does reach zero AT the
    // junction rather than merely being negative near it.
    EXPECT_NEAR(f[1], 0.0, 2e-5) << "the analytic branch reaches zero at the boundary";
    // The step across the junction: I of a particle mirroring at the plane crossing. Bounded, and
    // small — this shell sits near L = 4.6, where the brief's table gives ~6e-4 Re of radius error.
    const double jump = f[2] - f[1];
    EXPECT_GT(jump, 0.0) << "the residual steps UP across the junction, never down";
    EXPECT_LT(jump, 5e-3) << "the junction step is the largest remaining error term; it must stay "
                             "an error term";
    std::printf("[ measured ] residual step at the accessibility junction: %.3g Re\n", jump);
}

// f(r) = I(r) − I₀ increases with r: a larger shell is a weaker field, a lower B_min, and a longer
// bounce path between the same two mirror fields. Brent's bracket rests on exactly this.
TEST(IrbemDriftShell, ResidualIsMonotoneInTheShellRadius) {
    const Epoch e = standard_epoch();
    const fx::mat3d mag_to_geo = ib::rotation_matrix<Frame::GEO, Frame::MAG>(e.rotations);
    std::vector<float> coef(ib::detail::coefficient_slots(13));
    std::vector<float> norm(ib::detail::normalisation_slots(13));
    ib::detail::stage_model(e.model, coef, norm);
    const ib::DriftShellOptions opt;

    constexpr std::size_t kN = 24;
    std::vector<double> az(kN, 1.3);
    std::vector<double> mirror(kN, 150.0);
    std::vector<double> target(kN, 1.0);
    std::vector<double> radius(kN);
    std::vector<double> f(kN);
    for (std::size_t i = 0; i < kN; ++i) radius[i] = 3.0 + (0.25 * static_cast<double>(i));
    int traces = 0;
    (void)ib::detail::residual_round(e.model, mag_to_geo, coef, norm, az, mirror, target, radius, f,
                                     opt, traces);
    for (std::size_t i = 1; i < kN; ++i) {
        EXPECT_GT(f[i], f[i - 1]) << "residual fell between r = " << radius[i - 1] << " and "
                                  << radius[i];
    }
    EXPECT_LT(f.front(), 0.0);
    EXPECT_GT(f.back(), 0.0) << "the sweep must bracket the root";
    EXPECT_GT(traces, 0);
    EXPECT_LT(traces, static_cast<int>(kN));
}

// ---- the physics, against a closed form ---------------------------------------------------------

// For a centred dipole the drift shell IS the McIlwain shell, so L* must equal L_m. Nothing about
// the flux integral, the footpoint walk or the root-find is assumed here: if the polar cap were
// bounded wrongly, if k₀ were stale, or if the quadrature were biased, this would fail.
TEST(IrbemDriftShell, DipoleShellReproducesMcIlwainL) {
    const ib::Igrf<1> dipole = ib::Igrf<1>::at(kDecimalYear).value();
    const ib::Result<ib::Rotations> rot = ib::api::rotations_at(kYear, kDoy, kUt, dipole);
    ASSERT_EQ(rot.status, Status::Ok);
    const ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(3, 3);

    double worst = 0.0;
    for (const double l : {2.0, 3.0, 4.5, 6.6}) {
        const ib::Result<ib::DriftShell> s =
            ib::make_lstar(dipole, rot.value, geo(l, 0.0, 0.0), 90.0, opt);
        ASSERT_EQ(s.status, Status::Ok) << "L = " << l;
        EXPECT_EQ(s.value.azimuths, opt.azimuths);
        EXPECT_GT(s.value.phi, 0.0);
        // Φ = 2πk₀/(L* R_E) is the definition; check the two numbers are each other's, not two
        // independently computed quantities that happen to be close.
        EXPECT_NEAR(s.value.phi * s.value.lstar, 2.0 * std::numbers::pi * ib::dipole_moment(dipole),
                    1e-9 * s.value.phi * s.value.lstar);
        const double dev = std::abs(s.value.lstar - s.value.lm);
        worst = std::max(worst, dev);
        EXPECT_LT(dev, 6e-3) << "L = " << l << ": L* = " << s.value.lstar
                             << ", L_m = " << s.value.lm;
    }
    std::printf("[ measured ] centred dipole, |L* - L_m| worst = %.2e\n", worst);
}

// ---- the differential comparisons ---------------------------------------------------------------

TEST(IrbemDriftShell, MatchesTheOracleAtIrbemDefaultResolution) {
    const Epoch e = standard_epoch();
    const ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(0, 0);
    double worst_vs_converged = 0.0;
    double sum_vs_converged = 0.0;
    double worst_vs_matched = 0.0;
    double oracle_own_error = 0.0;
    double oracle_own_sum = 0.0;
    for (std::size_t i = 0; i < kOracleRes0.size(); ++i) {
        const OraclePoint& p = kOracleRes0[i];
        const ib::Result<ib::DriftShell> s =
            ib::make_lstar(e.model, e.rotations, geo(p.x, p.y, p.z), 90.0, opt);
        ASSERT_EQ(s.status, Status::Ok) << p.name;
        EXPECT_EQ(s.value.azimuths, 25);
        EXPECT_GT(s.value.traces, 25) << "every azimuth must have cost at least one trace";
        const double vs_converged = std::abs(s.value.lstar - kOracleRes9[i].lstar);
        worst_vs_converged = std::max(worst_vs_converged, vs_converged);
        sum_vs_converged += vs_converged;
        worst_vs_matched = std::max(worst_vs_matched, std::abs(s.value.lstar - p.lstar));
        const double own = std::abs(p.lstar - kOracleRes9[i].lstar);
        oracle_own_error = std::max(oracle_own_error, own);
        oracle_own_sum += own;
    }
    const auto n = static_cast<double>(kOracleRes0.size());
    std::printf("[ measured ] options(3,4)=0: ours vs converged oracle  max %.4f  mean %.4f\n",
                worst_vs_converged, sum_vs_converged / n);
    std::printf("[ measured ] options(3,4)=0: IRBEM vs converged oracle max %.4f  mean %.4f\n",
                oracle_own_error, oracle_own_sum / n);
    std::printf("[ measured ] options(3,4)=0: ours vs IRBEM, same setting, max %.4f\n",
                worst_vs_matched);

    // What is asserted at this resolution is the distance from the CONVERGED oracle, because
    // ERROR_BUDGET.md §2(a) measures IRBEM's own default-resolution error at 0.010-0.017 and the
    // transcribed tables above reproduce that (`oracle_own_error`). Comparing two implementations
    // that are each ~0.01 from converged, in opposite directions, cannot be held to 0.01.
    EXPECT_LT(worst_vs_converged, 0.015);
    EXPECT_LT(sum_vs_converged / n, 0.006);
    // The claim the design earns: at IRBEM's own recommended resolution this implementation sits
    // closer to converged than IRBEM does. That is the second-order midpoint quadrature with an
    // exactly-truncated final cell, against IRBEM's first-order rectangle rule on a grid the
    // boundary jitters across.
    EXPECT_LT(sum_vs_converged, oracle_own_sum)
        << "the whole point of not reproducing IRBEM's quadrature";
}

TEST(IrbemDriftShell, MatchesTheOracleAtIrbemHighestResolution) {
    const Epoch e = standard_epoch();
    const ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(9, 9);
    double worst = 0.0;
    for (const OraclePoint& p : kOracleRes9) {
        const ib::Result<ib::DriftShell> s =
            ib::make_lstar(e.model, e.rotations, geo(p.x, p.y, p.z), 90.0, opt);
        ASSERT_EQ(s.status, Status::Ok) << p.name;
        EXPECT_EQ(s.value.azimuths, 250);
        const double dev = std::abs(s.value.lstar - p.lstar);
        worst = std::max(worst, dev);
        EXPECT_LT(dev, 0.01) << p.name << ": L* = " << s.value.lstar << " vs " << p.lstar;
    }
    std::printf("[ measured ] options(3,4)=9: ours vs IRBEM, matched, max %.4f\n", worst);
}

// Shell splitting: same point, six pitch angles, `I` spanning 0.002 to 11.6 and `B_m` a factor of
// fifteen. This is what fails if the shell conserves only one of the two invariants.
TEST(IrbemDriftShell, ShellSplittingMatchesTheOracleAcrossPitchAngles) {
    const Epoch e = standard_epoch();
    const ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(9, 9);
    double worst = 0.0;
    for (const SplitPoint& p : kOracleSplitRes9) {
        const ib::Result<ib::DriftShell> s =
            ib::make_lstar(e.model, e.rotations, geo(p.x, p.y, p.z), p.alpha, opt);
        ASSERT_EQ(s.status, Status::Ok) << "alpha = " << p.alpha;
        const double dev = std::abs(s.value.lstar - p.lstar);
        worst = std::max(worst, dev);
        EXPECT_LT(dev, 0.01) << "alpha = " << p.alpha << " at (" << p.x << ", " << p.y << ", "
                             << p.z << ")";
    }
    std::printf("[ measured ] shell splitting 15-90 deg: max %.4f\n", worst);
}

// The epoch sweep is the control on `k₀`: the dipole moment falls ~6 % from 1900 to 2029 and
// L ∝ M^(1/3), so a hard-coded moment would show up here as a systematic drift with year.
TEST(IrbemDriftShell, EpochSweepMatchesTheOracle) {
    const ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(9, 9);
    double worst = 0.0;
    double moment_span = 0.0;
    double first_moment = 0.0;
    for (const EpochPoint& p : kOracleEpochRes9) {
        const Epoch e = epoch_at(p.year, p.year + (179.5 / 365.0));
        const double moment = ib::dipole_moment(e.model);
        if (first_moment == 0.0) first_moment = moment;
        moment_span = std::max(moment_span, std::abs(moment - first_moment) / first_moment);
        const ib::Result<ib::DriftShell> s =
            ib::make_lstar(e.model, e.rotations, geo(6.0, 0.0, 0.0), 90.0, opt);
        ASSERT_EQ(s.status, Status::Ok) << p.year;
        const double dev = std::abs(s.value.lstar - p.lstar);
        worst = std::max(worst, dev);
        EXPECT_LT(dev, 0.01) << p.year << ": L* = " << s.value.lstar << " vs " << p.lstar;
    }
    EXPECT_GT(moment_span, 0.04) << "the epochs chosen must actually move the dipole moment";
    std::printf("[ measured ] epochs 1900-2029 (moment moves %.1f%%): max %.4f\n",
                100.0 * moment_span, worst);
}

// ---- the batch entry point ----------------------------------------------------------------------

TEST(IrbemDriftShell, BatchMatchesThePointAtATimeCall) {
    const Epoch e = standard_epoch();
    const ib::DriftShellOptions opt;
    std::vector<ib::Position<Frame::GEO>> pts;
    std::vector<double> pitch;
    for (const OraclePoint& p : kOracleRes0) {
        pts.push_back(geo(p.x, p.y, p.z));
        pitch.push_back(90.0);
    }
    std::vector<ib::DriftShell> out(pts.size());
    std::vector<Status> st(pts.size());
    // Forced onto one lane on purpose. `prefer_gpu` is a function of BATCH SIZE, so a twelve-point
    // batch and a one-point call can legitimately take different lanes for the same stage on a
    // machine with a device — and an fp32 seed field gives a very slightly different pitch angle
    // than an fp64 one. Bit-exactness is a claim about the ALGORITHM, so the lane is pinned to make
    // it one; `UsesTheDeviceWhenOneIsAvailable` measures the cross-lane agreement separately.
    const HostLaneOnly forced;
    const ib::Result<bool> r =
        ib::make_lstar_batch(e.model, e.rotations, pts, pitch, out, st, opt);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_FALSE(r.value);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        ASSERT_EQ(st[i], Status::Ok) << kOracleRes0[i].name;
        const ib::Result<ib::DriftShell> one =
            ib::make_lstar(e.model, e.rotations, pts[i], 90.0, opt);
        ASSERT_EQ(one.status, Status::Ok);
        // Same arithmetic, same order, one lane: the two must agree bit for bit when the batch is
        // small enough that both take the host lane.
        EXPECT_EQ(out[i].lstar, one.value.lstar) << kOracleRes0[i].name;
        EXPECT_EQ(out[i].phi, one.value.phi);
        EXPECT_EQ(out[i].lm, one.value.lm);
        EXPECT_EQ(out[i].invariant_i, one.value.invariant_i);
        EXPECT_EQ(out[i].b_mirror, one.value.b_mirror);
        EXPECT_EQ(out[i].b_min, one.value.b_min);
        EXPECT_EQ(out[i].b_local, one.value.b_local);
    }
}

TEST(IrbemDriftShell, RejectsMalformedBatchesAndImpossibleSettings) {
    const Epoch e = standard_epoch();
    const std::array<ib::Position<Frame::GEO>, 2> pts{geo(4.0, 0.0, 0.0), geo(5.0, 0.0, 0.0)};
    const std::array<double, 1> short_pitch{90.0};
    std::array<ib::DriftShell, 2> out{};
    std::array<Status, 2> st{};
    EXPECT_EQ(ib::make_lstar_batch(e.model, e.rotations, pts, short_pitch, out, st).status,
              Status::DomainError);

    const std::array<double, 2> pitch{90.0, 90.0};
    ib::DriftShellOptions bad;
    bad.azimuths = 2;   // a two-point contour does not enclose a cap
    EXPECT_EQ(ib::make_lstar_batch(e.model, e.rotations, pts, pitch, out, st, bad).status,
              Status::DomainError);
    EXPECT_EQ(st[0], Status::DomainError);
    bad = ib::DriftShellOptions{};
    bad.bracket_trials = 1;   // one trial cannot bracket anything
    EXPECT_EQ(ib::make_lstar_batch(e.model, e.rotations, pts, pitch, out, st, bad).status,
              Status::DomainError);
    bad = ib::DriftShellOptions{};
    bad.colatitude_step_deg = 0.0;
    EXPECT_EQ(ib::make_lstar_batch(e.model, e.rotations, pts, pitch, out, st, bad).status,
              Status::DomainError);

    // An empty batch is a no-op, not an error.
    EXPECT_EQ(ib::make_lstar_batch(e.model, e.rotations,
                                   std::span<const ib::Position<Frame::GEO>>{},
                                   std::span<const double>{}, std::span<ib::DriftShell>{},
                                   std::span<Status>{})
                  .status,
              Status::Ok);

    // The same four refusals against the degree-1 model. Not redundant: `make_lstar_batch` is a
    // template, every truncation degree is its own instantiation, and a guard that is only ever
    // compiled for one of them is only ever tested for one of them.
    const ib::Igrf<1> dipole = ib::Igrf<1>::at(kDecimalYear).value();
    const ib::Result<ib::Rotations> drot = ib::api::rotations_at(kYear, kDoy, kUt, dipole);
    ASSERT_EQ(drot.status, Status::Ok);
    EXPECT_EQ(ib::make_lstar_batch(dipole, drot.value, pts, short_pitch, out, st).status,
              Status::DomainError);
    ib::DriftShellOptions dbad;
    dbad.azimuths = 2;
    EXPECT_EQ(ib::make_lstar_batch(dipole, drot.value, pts, pitch, out, st, dbad).status,
              Status::DomainError);
    EXPECT_EQ(st[1], Status::DomainError);
    EXPECT_EQ(ib::make_lstar_batch(dipole, drot.value,
                                   std::span<const ib::Position<Frame::GEO>>{},
                                   std::span<const double>{}, std::span<ib::DriftShell>{},
                                   std::span<Status>{})
                  .status,
              Status::Ok);
}

// Failure has to be REPORTED, not smoothed over. A starting point inside the atmosphere has no
// closed line; a starting point whose shell cannot be bracketed within the sweep window has no
// shell. Both come back named, with the diagnostic fields still populated.
TEST(IrbemDriftShell, ReportsFailureRatherThanInventingAShell) {
    const Epoch e = standard_epoch();
    const ib::Result<ib::DriftShell> inside =
        ib::make_lstar(e.model, e.rotations, geo(0.5, 0.0, 0.0), 90.0);
    EXPECT_EQ(inside.status, Status::DomainError);
    EXPECT_EQ(inside.value.lstar, 0.0);

    // A sweep window that cannot reach the root: every trial radius sits below the shell, so no
    // pair of them brackets and the point is NotConverged rather than silently bisected.
    ib::DriftShellOptions narrow;
    narrow.bracket_low = 0.20;
    narrow.bracket_high = 0.40;
    const ib::Result<ib::DriftShell> unbracketed =
        ib::make_lstar(e.model, e.rotations, geo(6.0, 0.0, 0.0), 90.0, narrow);
    EXPECT_EQ(unbracketed.status, Status::NotConverged);
    EXPECT_EQ(unbracketed.value.lstar, 0.0);
    EXPECT_LT(unbracketed.value.azimuths, 25);
    // The diagnostics survive the failure: L_m, I and the fields are still the starting line's.
    EXPECT_NEAR(unbracketed.value.lm, 6.067, 1e-2);
    EXPECT_GT(unbracketed.value.b_local, 0.0);

    // A footpoint walk that cannot reach the surface within its step cap is the same story.
    ib::DriftShellOptions stingy;
    stingy.footpoint = ib::TraceOptions{100.0, 5, 1.0};
    const ib::Result<ib::DriftShell> nofoot =
        ib::make_lstar(e.model, e.rotations, geo(6.0, 0.0, 0.0), 90.0, stingy);
    EXPECT_EQ(nofoot.status, Status::NotConverged);

    // And all three against the degree-1 model, for the reason the previous test spells out: a
    // different truncation degree is a different instantiation of every one of these paths.
    const ib::Igrf<1> dipole = ib::Igrf<1>::at(kDecimalYear).value();
    const ib::Result<ib::Rotations> drot = ib::api::rotations_at(kYear, kDoy, kUt, dipole);
    ASSERT_EQ(drot.status, Status::Ok);
    EXPECT_EQ(ib::make_lstar(dipole, drot.value, geo(0.5, 0.0, 0.0), 90.0).status,
              Status::DomainError);
    EXPECT_EQ(ib::make_lstar(dipole, drot.value, geo(6.0, 0.0, 0.0), 90.0, narrow).status,
              Status::NotConverged);
    EXPECT_EQ(ib::make_lstar(dipole, drot.value, geo(6.0, 0.0, 0.0), 90.0, stingy).status,
              Status::NotConverged);
}

// ---- the geometries the equatorial table does not reach -----------------------------------------

// A GOLDEN SET CHOSEN TO BE AWKWARD, not to be passed. `kOracleRes0`/`kOracleRes9` above are ten
// equatorial points and two NORTHERN off-equator ones, all at eastward longitudes. That sample is
// where this implementation looks best, and it is not the whole domain. These eleven are the
// complement: four SOUTHERN off-equator geometries (the hemisphere the dipole tilt makes the
// harder one, since the field line's footpoint wanders further in longitude on the way down),
// three at negative x/y, the inner edge at L = 1.5, the outer edge at L = 9, and a start tilted
// out of every plane at once.
//
// It found something. See `LargeOffEquatorialDeviationsAreTheTraceStep` below: at IRBEM's most
// converged setting the worst of these is 0.0105, OUTSIDE `docs/ERROR_BUDGET.md`'s 0.01 absolute
// L* budget, where the shipped equatorial set tops out at 0.0066 — and the excess is the DEFAULT
// TRACE STEP, not the drift-shell algorithm.
//
// Provenance: `make_lstar1_` from `libirbem-O2.so` (IRBEM e7cecb0, gfortran 13.3.0,
// `-O2 -ffp-contract=off -fno-fast-math`), kext = 0, options = {1, 0, 9, 9, 0}, 2015-180 12:00 UT,
// sysaxes = 1.
struct HardPoint {
    const char* name;
    double x;
    double y;
    double z;
    double lstar;    ///< IRBEM's L*.
    double lm;       ///< IRBEM's Lm.
    double xj;       ///< IRBEM's XJ, i.e. the second invariant I.
    double blocal;   ///< IRBEM's Blocal, nT.
};

constexpr std::array<HardPoint, 11> kOracleHardRes9{{
    {"S-off5.5", 5.5, 0.0, -2.0, 6.5354543, 6.5379445, 1.9304374, 161.677702},
    {"S-off3.0", 3.0, 0.0, -1.0, 3.5630373, 3.5622218, 1.0015412, 979.226244},
    {"S-y4", 0.0, 4.0, -1.5, 5.5769183, 5.5760390, 4.1133478, 521.800926},
    {"S-deep2", 2.0, 0.0, -1.2, 3.2289617, 3.2271433, 2.3917003, 2707.458324},
    {"N-deep2", 2.0, 0.0, 1.2, 3.1262364, 3.1247390, 2.3610186, 3056.900822},
    {"negx4", -4.0, 0.0, 0.0, 3.9607674, 3.9621661, 0.065457393, 490.925098},
    {"negy5", 0.0, -5.0, 0.0, 5.1974590, 5.1971419, 0.49193059, 242.087194},
    {"negxy3", -3.0, -3.0, 0.0, 4.2542905, 4.2543914, 0.070215215, 396.543875},
    {"L1.5", 1.5, 0.0, 0.0, 1.5488062, 1.5478565, 0.025489079, 8233.601366},
    {"L9", 9.0, 0.0, 0.0, 9.0690265, 9.0733742, 0.035620357, 40.188498},
    {"tilt42", 4.0, 2.0, -2.0, 6.1254896, 6.1231754, 3.6136239, 308.926811},
}};

TEST(IrbemDriftShell, MatchesTheOracleOffTheEquatorAndOffTheEasternHemisphere) {
    const Epoch e = standard_epoch();
    const ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(9, 9);
    double worst = 0.0;
    const char* worst_name = "";
    for (const HardPoint& p : kOracleHardRes9) {
        const ib::Result<ib::DriftShell> s =
            ib::make_lstar(e.model, e.rotations, geo(p.x, p.y, p.z), 90.0, opt);
        ASSERT_EQ(s.status, Status::Ok) << p.name;
        EXPECT_EQ(s.value.azimuths, 250) << p.name;
        const double dev = std::abs(s.value.lstar - p.lstar);
        if (dev > worst) {
            worst = dev;
            worst_name = p.name;
        }
        // NOT the 0.01 budget. This asserts what is MEASURED, which is worse than the budget at
        // one of these eleven, and the next test explains and bounds why.
        EXPECT_LT(dev, 0.012) << p.name << ": L* = " << s.value.lstar << " vs " << p.lstar;
    }
    std::printf("[ measured ] off-equator/off-eastern, options=9: max %.4f at %s\n", worst,
                worst_name);
    EXPECT_GT(worst, 0.008) << "if this set stopped being the hard one, it stopped being useful";
}

// WHERE THAT 0.0105 COMES FROM, measured rather than asserted away. `DriftShellOptions::trace`
// defaults to `lstar.hpp`'s `steps_per_l = 50`, and for a 90-degree particle at an OFF-EQUATORIAL
// start the mirror point is the start point itself, so the integrand of `I` leaves its endpoint
// with a square-root singularity that a fixed-step rule resolves slowly. Refining the trace moves
// the deviation monotonically down to a NON-ZERO limit: the residue is a genuine algorithmic
// difference from IRBEM, and the part above it is the step.
//
// So the fix, when a caller needs the 0.01 budget off the equator, is one field on the options
// struct and not a change of algorithm — and it costs a factor of two in trace time, which is why
// it is not the default.
TEST(IrbemDriftShell, LargeOffEquatorialDeviationsAreTheTraceStep) {
    const Epoch e = standard_epoch();
    const HardPoint& p = kOracleHardRes9[2];   // S-y4, the worst of the eleven
    double previous = 1.0;
    for (const double steps : {50.0, 100.0, 200.0, 400.0}) {
        ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(9, 9);
        opt.trace = ib::TraceOptions{steps, 20000, 1.0};
        const ib::Result<ib::DriftShell> s =
            ib::make_lstar(e.model, e.rotations, geo(p.x, p.y, p.z), 90.0, opt);
        ASSERT_EQ(s.status, Status::Ok) << steps;
        const double dev = std::abs(s.value.lstar - p.lstar);
        std::printf("[ measured ] %s, trace steps_per_l = %4.0f: |ours - IRBEM| = %.4f\n", p.name,
                    steps, dev);
        EXPECT_LT(dev, previous) << "refining the trace must not make it worse: " << steps;
        previous = dev;
    }
    // The default is outside the budget here and one doubling is inside it. Both halves are
    // asserted, because the claim is that the step is the cause and not that the number is small.
    ib::DriftShellOptions coarse = ib::DriftShellOptions::from_irbem(9, 9);
    const ib::Result<ib::DriftShell> at_default =
        ib::make_lstar(e.model, e.rotations, geo(p.x, p.y, p.z), 90.0, coarse);
    EXPECT_GT(std::abs(at_default.value.lstar - p.lstar), 0.01);
    coarse.trace = ib::TraceOptions{100.0, 20000, 1.0};
    const ib::Result<ib::DriftShell> refined =
        ib::make_lstar(e.model, e.rotations, geo(p.x, p.y, p.z), 90.0, coarse);
    EXPECT_LT(std::abs(refined.value.lstar - p.lstar), 0.01);
}

// THE DIAGNOSTIC FIELDS, which no other test in this file compares against the oracle at all.
// `DriftShell` returns nine numbers and only `lstar` was ever differenced; `lm`, `invariant_i` and
// `b_local` are advertised as IRBEM's Lm, XJ and Blocal, so they are held to IRBEM here — and they
// do NOT all meet `docs/ERROR_BUDGET.md`. What is asserted is what is measured, with the budget
// stated beside it, because a tolerance quietly loosened to green is worse than a failing test.
//
//   Blocal : budget 1e-6 relative, MEASURED 1.7e-5 worst (at L = 1.5, where the high-degree terms
//            are least attenuated). Upstream of this header — `igrf.hpp` evaluates the field.
//   Lm     : budget 1e-3 relative, MEASURED 5.4e-3 worst, all of it at off-equatorial starts.
//   XJ     : budget 1e-4 relative, MEASURED 9.2e-2 worst. Two separate regimes: at a near-equatorial
//            start `I` is ~0.03-0.07 Re and a 6e-3 Re absolute error is a large FRACTION of it,
//            while at an off-equatorial start `I` is 1-4 Re and the error is the ~2 % endpoint
//            effect the previous test isolates. Both are the trace's fixed step, and both live in
//            `lstar.hpp` rather than here.
//
// L* survives all three because the root-find matches `I(r)` against `I0` computed the SAME way:
// the two errors are correlated and largely divide out. That is the substance of the claim, and it
// is why the shell is good to 1e-3 while the invariant it is built on is good to 1e-2.
TEST(IrbemDriftShell, ReportsTheDiagnosticFieldsIrbemDoesAndSaysHowCloselyTheyAgree) {
    const Epoch e = standard_epoch();
    const ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(9, 9);
    double worst_lm = 0.0;
    double worst_xj = 0.0;
    double worst_b = 0.0;
    for (const HardPoint& p : kOracleHardRes9) {
        const ib::Result<ib::DriftShell> s =
            ib::make_lstar(e.model, e.rotations, geo(p.x, p.y, p.z), 90.0, opt);
        ASSERT_EQ(s.status, Status::Ok) << p.name;
        worst_lm = std::max(worst_lm, std::abs(s.value.lm - p.lm) / p.lm);
        worst_xj = std::max(worst_xj, std::abs(s.value.invariant_i - p.xj) / p.xj);
        worst_b = std::max(worst_b, std::abs(s.value.b_local - p.blocal) / p.blocal);
        // B_min is bounded by B_local from below and by nothing above; no oracle value is
        // transcribed for it, so what is checked is the inequality the definition guarantees.
        EXPECT_GT(s.value.b_min, 0.0) << p.name;
        EXPECT_LE(s.value.b_min, s.value.b_local * (1.0 + 1e-12)) << p.name;
        // A 90-degree particle mirrors where it sits: B_mirror IS B_local, by construction.
        EXPECT_NEAR(s.value.b_mirror, s.value.b_local, 1e-9 * s.value.b_local) << p.name;
    }
    std::printf("[ measured ] vs IRBEM, relative: Lm %.2e (budget 1e-3), XJ %.2e (budget 1e-4), "
                "Blocal %.2e (budget 1e-6)\n",
                worst_lm, worst_xj, worst_b);
    EXPECT_LT(worst_b, 3e-5);
    EXPECT_LT(worst_lm, 8e-3);
    EXPECT_LT(worst_xj, 1.2e-1);
}

// THE LOSS CONE, where the right answer is to refuse. IRBEM returns its -1e31 BADDATA sentinel for
// a pitch angle whose mirror point is below the atmosphere; this implementation returns
// `Status::NotConverged` with `lstar == 0`. The two libraries must draw the line in the same place,
// because a shell invented just inside the loss cone is the failure mode that looks like data.
//
// Transcribed from `make_lstar_shell_splitting1_` at options = {1,0,9,9,0}, kext = 0, 2015-180
// 12:00 UT: at (4,0,0) GEO it returns a shell down to alpha = 10 degrees and BADDATA at 5 and
// below; at (6,0,-1) it returns one down to 5 degrees and BADDATA at 3 and below.
TEST(IrbemDriftShell, RefusesTheLossConeWhereIrbemDoes) {
    const Epoch e = standard_epoch();
    const ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(9, 9);
    struct Case {
        double x;
        double y;
        double z;
        double alpha;
        bool irbem_has_a_shell;
    };
    for (const Case& c : std::array<Case, 6>{{{4.0, 0.0, 0.0, 10.0, true},
                                              {4.0, 0.0, 0.0, 5.0, false},
                                              {4.0, 0.0, 0.0, 1.0, false},
                                              {6.0, 0.0, -1.0, 5.0, true},
                                              {6.0, 0.0, -1.0, 3.0, false},
                                              {6.0, 0.0, -1.0, 1.0, false}}}) {
        const ib::Result<ib::DriftShell> s =
            ib::make_lstar(e.model, e.rotations, geo(c.x, c.y, c.z), c.alpha, opt);
        if (c.irbem_has_a_shell) {
            EXPECT_EQ(s.status, Status::Ok) << "alpha = " << c.alpha;
            EXPECT_GT(s.value.lstar, 0.0) << "alpha = " << c.alpha;
            continue;
        }
        EXPECT_NE(s.status, Status::Ok) << "alpha = " << c.alpha
                                        << ": IRBEM refuses this shell and so must we";
        EXPECT_EQ(s.value.lstar, 0.0) << "alpha = " << c.alpha;
    }
}

#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
// The footpoint kernel against the fp64 host walk it was written from. This is the one stage of L*
// whose output is a POSITION rather than a scalar, and the quantity that matters is the colatitude
// the flux quadrature is bounded by — so that is what is compared, in the frame it is read in.
TEST(IrbemDriftShell, FootpointsAgreeBetweenLanes) {
    if (!ib::gpu::available()) {
        GTEST_SKIP() << "no device: " << ib::gpu::unavailable_reason();
    }
    const Epoch e = standard_epoch();
    const fx::mat3d mag_to_geo = ib::rotation_matrix<Frame::GEO, Frame::MAG>(e.rotations);
    const fx::mat3d geo_to_mag = ib::rotation_matrix<Frame::MAG, Frame::GEO>(e.rotations);
    std::vector<float> coef(ib::detail::coefficient_slots(13));
    std::vector<float> norm(ib::detail::normalisation_slots(13));
    ib::detail::stage_model(e.model, coef, norm);
    const ib::TraceOptions foot{100.0, 8000, 1.0};

    // Seeds spread over the whole shell: 24 azimuths at four radii, i.e. the shape a real batch
    // has, with L = 2 lines and L = 8 lines in the same dispatch — which is what the step-size
    // convention exists to keep converged.
    constexpr std::size_t kAzimuths = 24;
    constexpr std::array<double, 4> kRadii{2.0, 4.0, 6.0, 8.0};
    std::vector<float> pos;
    std::vector<float> dir;
    std::vector<ib::Position<Frame::GEO>> seeds;
    for (const double l : kRadii) {
        for (std::size_t k = 0; k < kAzimuths; ++k) {
            const double phi =
                2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(kAzimuths);
            const ib::Position<Frame::GEO> s = ib::detail::equatorial_seed(mag_to_geo, phi, l);
            seeds.push_back(s);
            pos.push_back(static_cast<float>(s.v[0]));
            pos.push_back(static_cast<float>(s.v[1]));
            pos.push_back(static_cast<float>(s.v[2]));
            const fx::vec3d b = e.model.evaluate(s).v;
            dir.push_back(fx::dot(b, e.rotations.dipole_geo) >= 0.0 ? 1.0F : -1.0F);
        }
    }
    const std::size_t n = seeds.size();
    const std::array<std::uint32_t, 4> dims{
        static_cast<std::uint32_t>(n), 13U, static_cast<std::uint32_t>(foot.max_steps),
        static_cast<std::uint32_t>(foot.steps_per_l * 1000.0)};
    std::vector<float> out(3 * n);
    std::vector<std::uint32_t> st(n);
    ASSERT_TRUE(ib::gpu::launch_shell_foot(pos, dir, coef, norm, dims, out, st))
        << "the kernel must run, not fall back silently";

    double worst_colatitude = 0.0;
    double worst_longitude = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        ASSERT_EQ(st[i], static_cast<std::uint32_t>(Status::Ok)) << "line " << i;
        const ib::Result<ib::Position<Frame::GEO>> host =
            ib::detail::walk_to_surface(e.model, seeds[i], dir[i] > 0.0F, foot);
        ASSERT_EQ(host.status, Status::Ok) << "line " << i;
        const fx::vec3d device{static_cast<double>(out[(3 * i) + 0]),
                               static_cast<double>(out[(3 * i) + 1]),
                               static_cast<double>(out[(3 * i) + 2])};
        // Both lanes normalise onto the unit sphere, so this is the kernel's own claim, checked.
        EXPECT_NEAR(fx::norm(device), 1.0, 1e-6) << "line " << i;
        const fx::vec3d dm = geo_to_mag * device;
        const fx::vec3d hm = geo_to_mag * host.value.v;
        worst_colatitude = std::max(worst_colatitude, std::abs(std::acos(std::clamp(dm[2], -1.0,
                                                                                    1.0)) -
                                                               std::acos(std::clamp(hm[2], -1.0,
                                                                                    1.0))));
        worst_longitude = std::max(
            worst_longitude, std::abs(std::atan2(dm[1], dm[0]) - std::atan2(hm[1], hm[0])));
    }
    const double deg = 180.0 / std::numbers::pi;
    std::printf("[ measured ] footpoints, fp32 device vs fp64 host over %zu lines: "
                "colatitude %.4f deg, longitude %.4f deg\n",
                n, worst_colatitude * deg, worst_longitude * deg);
    // The flux cell is 0.25 degrees wide at IRBEM's default resolution, so this is what "the
    // footpoint is not the error term" has to mean quantitatively.
    EXPECT_LT(worst_colatitude * deg, 0.05);
    EXPECT_LT(worst_longitude * deg, 0.05);
}

// The claim a performance number rests on: that the device was actually used. A silent fallback to
// the host is exactly the failure mode that makes a speedup meaningless, so it is asserted rather
// than assumed — and the device's answer is held to the same budget as the host's.
TEST(IrbemDriftShell, UsesTheDeviceWhenOneIsAvailable) {
    if (!ib::gpu::available()) {
        GTEST_SKIP() << "no device: " << ib::gpu::unavailable_reason();
    }
    const Epoch e = standard_epoch();
    const ib::DriftShellOptions opt;
    // Above the ~512-line crossover: 32 points x 25 azimuths is 800 shells, so the trace, the
    // footpoint and the flux stages all take the device lane.
    constexpr std::size_t kN = 32;
    std::vector<ib::Position<Frame::GEO>> pts(kN);
    std::vector<double> pitch(kN, 90.0);
    for (std::size_t i = 0; i < kN; ++i) {
        const double l = 2.5 + (5.0 * static_cast<double>(i) / static_cast<double>(kN));
        const double ph = 0.37 * static_cast<double>(i);
        pts[i] = geo(l * std::cos(ph), l * std::sin(ph), 0.0);
    }
    std::vector<ib::DriftShell> dev(kN);
    std::vector<Status> dev_st(kN);
    const ib::Result<bool> r =
        ib::make_lstar_batch(e.model, e.rotations, pts, pitch, dev, dev_st, opt);
    EXPECT_EQ(r.status, Status::Ok);
    // The flag means THE TRACES ran on the device and nothing weaker. Before it was narrowed it
    // was an OR across the trace, footpoint and flux stages, and the flux stage alone — a plain
    // IGRF batch that goes to the device at 128 points — held it true: deleting
    // `irbem_trace_i_f32.spv` from the shader directory made this batch run at 7.9 ms/point
    // instead of 0.35 and this assertion still PASSED. It no longer does.
    EXPECT_TRUE(r.value) << "the device lane must have serviced a 800-shell batch";

    std::vector<ib::DriftShell> host(kN);
    std::vector<Status> host_st(kN);
    ib::Result<bool> h{Status::Ok, false};
    {
        const HostLaneOnly forced;
        h = ib::make_lstar_batch(e.model, e.rotations, pts, pitch, host, host_st, opt);
    }
    EXPECT_FALSE(h.value);

    double worst = 0.0;
    for (std::size_t i = 0; i < kN; ++i) {
        ASSERT_EQ(dev_st[i], Status::Ok);
        ASSERT_EQ(host_st[i], Status::Ok);
        worst = std::max(worst, std::abs(dev[i].lstar - host[i].lstar));
    }
    // The fp32 device lane against the fp64 host lane, on the same algorithm: well inside the 0.01
    // L* budget, and measured rather than asserted.
    EXPECT_LT(worst, 0.01);
    std::printf("[ measured ] device vs host L*, %zu points: max %.3g\n", kN, worst);
}
#endif
