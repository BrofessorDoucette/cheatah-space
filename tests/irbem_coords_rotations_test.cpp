// Unit tests for space.irbem's Cartesian frame rotations — sidereal time, the solar ephemeris, the
// geodipole, the Hapgood 1992 matrix chain, and the frame-tagged transforms built on it.
//
// The suite is deliberately split in two, because two very different things are being checked and
// mixing them costs the exactness of one and the meaning of the other:
//
//   * The MACHINERY (which matrix is selected for which directed pair, that the reverse direction
//     is the transpose, that a round trip is a round trip) is tested against a synthetic epoch
//     whose eight matrices are hand-written signed permutations. Those have entries 0 and ±1, so
//     every product is exact in binary and every assertion is `==` with no tolerance at all. A
//     tolerance here would hide exactly the wire-crossing bug this layer exists to prevent.
//   * The PHYSICS (that GSM's Z lies in the dipole–Sun plane, that SM's Z is the dipole, that the
//     Sun is GSE's +X) is tested against a real epoch, where transcendentals make exactness
//     impossible. Those assertions are PROPERTIES, not stored numbers: a pair of matrices that are
//     mutually inverse but both wrong passes every round trip and fails every one of them.
//
// The reference epoch throughout is JD 2457203.0 = 2015-06-29 12:00 UT, with the DGRF 2015.0
// degree-1 coefficients, matching the setup already used by docs/ERROR_BUDGET.md.
#include <gtest/gtest.h>

#include "alloc_counter.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <new>

#include "space/irbem/coords_rotations.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;

using ib::Frame;

namespace {

// ---- the two epochs the suite runs against ------------------------------------------------------

/// 2015-06-29 12:00 UT. Chosen to match the oracle setup in docs/ERROR_BUDGET.md.
constexpr double kJdReference = 2457203.0;

/// 2000-01-01 12:00 UT — the J2000.0 epoch itself, where the sidereal-time series has a published
/// value that can be checked without an oracle.
constexpr double kJdJ2000 = 2451545.0;

/// DGRF 2015.0 degree-1 Gauss coefficients, nT (IGRF-13 definitive set for 2015.0).
constexpr ib::DipoleCoefficients kDipole2015{-29441.46, -1501.77, 4795.99};

/// IGRF-13 2020.0 degree-1 Gauss coefficients, nT.
constexpr ib::DipoleCoefficients kDipole2020{-29404.8, -1450.9, 4652.5};

ib::Rotations reference_rotations() { return ib::Rotations::at(kJdReference, kDipole2015); }

// ---- assertion helpers --------------------------------------------------------------------------

/// Componentwise closeness of two vectors. `tolerance == 0` means bit-exact equality.
void expect_vec_near(const fx::vec3d& actual, const fx::vec3d& expected, double tolerance,
                     const char* what) {
    for (std::size_t i = 0; i < 3; ++i) {
        if (tolerance == 0.0) {
            EXPECT_EQ(expected[i], actual[i]) << what << " component " << i;
        } else {
            EXPECT_NEAR(expected[i], actual[i], tolerance) << what << " component " << i;
        }
    }
}

/// A rotation must be orthogonal with unit determinant: no scale, no shear, no reflection. This is
/// the property that makes the reverse transform a transpose rather than an inversion, so it is
/// checked on every stored matrix rather than assumed.
void expect_proper_rotation(const fx::mat3d& m, const char* what) {
    const fx::mat3d product = m * fx::transpose(m);
    const fx::mat3d identity = fx::mat3d::identity();
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            // 1e-15 is ~4.5 machine epsilons: a 3-term dot product of unit-magnitude entries
            // carries at most a couple of ulps, and the composites are products of at most four
            // such matrices.
            EXPECT_NEAR(identity(r, c), product(r, c), 1e-15) << what << " (" << r << "," << c << ")";
        }
    }
    EXPECT_NEAR(1.0, fx::determinant(m), 1e-15) << what << " determinant";
}

// ---- the synthetic epoch: exact arithmetic for the dispatch machinery ---------------------------

/// Eight distinct proper signed permutation matrices — quarter and half turns about the axes. Their
/// entries are 0 and ±1, so a matrix-vector product is a permutation with sign flips and every
/// assertion downstream is exact rather than approximate.
constexpr fx::mat3d kQuarterZ{0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
constexpr fx::mat3d kQuarterX{1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0};
constexpr fx::mat3d kQuarterY{0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0};
constexpr fx::mat3d kHalfZ{-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0};
constexpr fx::mat3d kHalfX{1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0};
constexpr fx::mat3d kHalfY{-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0};
constexpr fx::mat3d kMinusQuarterZ{0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
constexpr fx::mat3d kCyclic{0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};

/// A Rotations whose every matrix is one of the permutations above. The scalars are placeholders —
/// nothing in the transform path reads them — and the eight matrices are deliberately mutually
/// inconsistent, because this fixture tests WHICH matrix a directed pair selects, not the physics
/// that relates them.
ib::Rotations permutation_rotations() {
    return ib::Rotations{kJdReference,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         fx::vec3d{0.0, 0.0, 1.0},
                         fx::vec3d{0.0, 0.0, 1.0},
                         kQuarterZ,       // gei_to_geo
                         kQuarterX,       // gei_to_gse
                         kQuarterY,       // gse_to_gsm
                         kHalfZ,          // gsm_to_sm
                         kHalfX,          // geo_to_mag
                         kHalfY,          // geo_to_gse
                         kMinusQuarterZ,  // geo_to_gsm
                         kCyclic};        // geo_to_sm
}

/// A probe with three distinct, exactly-representable, differently-signed components, so a
/// component swap or a sign flip anywhere in the dispatch cannot cancel out.
constexpr fx::vec3d kProbe{1.5, -2.25, 4.0};

// ---- the pair exerciser -------------------------------------------------------------------------

/// Exercises BOTH directions of one frame pair, for both tagged vector templates.
///
/// Every check here holds for any epoch, which is why the same function runs against the exact
/// permutation fixture (with `tolerance == 0`) and the real one.
///
/// @tparam A one frame of the pair. @tparam B the other.
template <Frame A, Frame B>
void exercise_pair(const ib::Rotations& rotations, double tolerance, const char* label) {
    const fx::mat3d forward = ib::rotation_matrix<B, A>(rotations);
    const fx::mat3d reverse = ib::rotation_matrix<A, B>(rotations);

    expect_proper_rotation(forward, label);
    expect_proper_rotation(reverse, label);

    // The reverse direction is not stored: it is the transpose, and that must be exact — a
    // reshuffle of the same nine doubles, never a re-derivation that could round differently.
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            EXPECT_EQ(forward(c, r), reverse(r, c)) << label << " transpose (" << r << "," << c << ")";
        }
    }

    const ib::Position<A> position{kProbe};
    const ib::Position<B> there = ib::transform<B>(position, rotations);
    const ib::Position<A> back = ib::transform<A>(there, rotations);
    expect_vec_near(back.v, kProbe, tolerance, label);
    // A rotation preserves length; a transposed-but-scaled matrix would not.
    EXPECT_NEAR(fx::norm(kProbe), fx::norm(there.v), 1e-15) << label << " position norm";

    const ib::FieldVector<A> field{kProbe};
    const ib::FieldVector<B> field_there = ib::transform<B>(field, rotations);
    const ib::FieldVector<A> field_back = ib::transform<A>(field_there, rotations);
    expect_vec_near(field_back.v, kProbe, tolerance, label);
    // A position and a field vector rotate identically — a rotation carries no translation, so the
    // two tagged templates must agree component for component.
    expect_vec_near(field_there.v, there.v, 0.0, label);

    // The matrix and the transform are the same operation, spelled two ways.
    expect_vec_near(there.v, forward * kProbe, 0.0, label);
}

/// Exercises the degenerate same-frame transform, which generic code parameterized on a target
/// frame needs to compile.
/// @tparam F the frame.
template <Frame F>
void exercise_identity(const ib::Rotations& rotations) {
    EXPECT_EQ(fx::mat3d::identity(), (ib::rotation_matrix<F, F>(rotations)));
    const ib::Position<F> position{kProbe};
    expect_vec_near(ib::transform<F>(position, rotations).v, kProbe, 0.0, "identity position");
    const ib::FieldVector<F> field{kProbe};
    expect_vec_near(ib::transform<F>(field, rotations).v, kProbe, 0.0, "identity field");
}

/// Runs @ref exercise_pair over every pair the header stores, in both directions.
void exercise_every_pair(const ib::Rotations& rotations, double tolerance) {
    exercise_pair<Frame::GEI, Frame::GEO>(rotations, tolerance, "GEI/GEO");
    exercise_pair<Frame::GEI, Frame::GSE>(rotations, tolerance, "GEI/GSE");
    exercise_pair<Frame::GSE, Frame::GSM>(rotations, tolerance, "GSE/GSM");
    exercise_pair<Frame::GSM, Frame::SM>(rotations, tolerance, "GSM/SM");
    exercise_pair<Frame::GEO, Frame::MAG>(rotations, tolerance, "GEO/MAG");
    exercise_pair<Frame::GEO, Frame::GSE>(rotations, tolerance, "GEO/GSE");
    exercise_pair<Frame::GEO, Frame::GSM>(rotations, tolerance, "GEO/GSM");
    exercise_pair<Frame::GEO, Frame::SM>(rotations, tolerance, "GEO/SM");
    exercise_identity<Frame::GEO>(rotations);
    exercise_identity<Frame::GSM>(rotations);
}

}  // namespace

// ---- compile-time surface -----------------------------------------------------------------------

/// Whether a directed frame pair has a rotation at all.
template <Frame To, Frame From>
concept HasRotation = requires(const ib::Rotations& r) { ib::rotation_matrix<To, From>(r); };

/// Whether a tagged value can be transformed into frame @p To.
template <Frame To, class V>
concept Transformable = requires(V v, const ib::Rotations& r) { ib::transform<To>(v, r); };

// The eight stored directions and their eight inverses exist...
static_assert(HasRotation<Frame::GEO, Frame::GEI> && HasRotation<Frame::GEI, Frame::GEO>);
static_assert(HasRotation<Frame::GSE, Frame::GEI> && HasRotation<Frame::GEI, Frame::GSE>);
static_assert(HasRotation<Frame::GSM, Frame::GSE> && HasRotation<Frame::GSE, Frame::GSM>);
static_assert(HasRotation<Frame::SM, Frame::GSM> && HasRotation<Frame::GSM, Frame::SM>);
static_assert(HasRotation<Frame::MAG, Frame::GEO> && HasRotation<Frame::GEO, Frame::MAG>);
static_assert(HasRotation<Frame::GSE, Frame::GEO> && HasRotation<Frame::GEO, Frame::GSE>);
static_assert(HasRotation<Frame::GSM, Frame::GEO> && HasRotation<Frame::GEO, Frame::GSM>);
static_assert(HasRotation<Frame::SM, Frame::GEO> && HasRotation<Frame::GEO, Frame::SM>);
static_assert(HasRotation<Frame::GSM, Frame::GSM>);  // the identity

// ...and nothing else does. An unstored pair is a compile error the caller must resolve by
// composing through GEO deliberately, not a silent chain picked for them.
static_assert(!HasRotation<Frame::SM, Frame::GEI>);
static_assert(!HasRotation<Frame::MAG, Frame::GSM>);
static_assert(!HasRotation<Frame::MAG, Frame::SM>);
static_assert(!HasRotation<Frame::SM, Frame::GSE>);

// The frame is in the type, so a transform that has no matrix for the frame it was handed cannot
// compile — the source frame is deduced from the argument and is never assumed...
static_assert(Transformable<Frame::GSM, ib::Position<Frame::GEO>>);
static_assert(Transformable<Frame::GSM, ib::FieldVector<Frame::GEO>>);
static_assert(!Transformable<Frame::SM, ib::Position<Frame::GEI>>);
static_assert(!Transformable<Frame::MAG, ib::FieldVector<Frame::GSM>>);
// ...and neither can a frame whose components are not (x, y, z) at all. Rotating a
// (radius, latitude, longitude) triple as though it were Cartesian is a classic silent disaster.
static_assert(!Transformable<Frame::GSM, ib::Position<Frame::GDZ>>);
static_assert(!Transformable<Frame::SPH, ib::Position<Frame::GEO>>);

// Rotations is a plain aggregate of doubles and fixarray values — no indirection, nothing to free,
// and cheap enough to hold one per epoch in an ephemeris-sized array.
static_assert(std::is_trivially_copyable_v<ib::Rotations>);
static_assert(std::is_aggregate_v<ib::Rotations>);

// ---- sidereal time ------------------------------------------------------------------------------

TEST(IrbemGmst, MatchesThePublishedValueAtJ2000) {
    // The IAU 1982 series evaluated at T = 0 is its own leading term, 67310.54841 s = 18h 41m
    // 50.54841s of sidereal time — the standard quoted GMST at the J2000.0 epoch. In degrees that
    // is the equally standard 280.46061837.
    EXPECT_NEAR(280.46061837, ib::gmst_iau1982_degrees(kJdJ2000), 1e-8);
    EXPECT_NEAR(67310.54841, ib::gmst_iau1982_degrees(kJdJ2000) * 240.0, 1e-6);
}

TEST(IrbemGmst, AgreesWithTheSeriesInItsOtherPublishedArrangement) {
    // The same IAU 1982 series is usually tabulated in degrees per day rather than seconds per
    // century: 280.46061837 + 360.98564736629 d + 3.87933e-4 T^2 - T^3/38710000. The header
    // evaluates the seconds-per-century form, but folds the huge linear coefficient into an exact
    // fractional-day term to avoid spending ten significant digits on whole turns about to be
    // discarded. Agreement between the two arrangements is what proves that trick sound.
    for (const double days : {-18262.0, -1000.0, 0.0, 0.375, 5478.5, 5658.0, 18262.0}) {
        const double jd = kJdJ2000 + days;
        const double t = days / 36525.0;
        const double reference =
            280.46061837 + (360.98564736629 * days) + (3.87933e-4 * t * t) - ((t * t * t) / 38710000.0);
        double wrapped = std::fmod(reference, 360.0);
        if (wrapped < 0.0) wrapped += 360.0;
        // 1e-7 deg = 0.36 milli-arcsecond, and the residual is entirely the REFERENCE
        // arrangement's own error: its constants are quoted to 11 decimals (360.98564736629
        // against the exact 360.985647366286...), which alone costs 7e-8 deg at +-50 years. That
        // is what limits this comparison, so it says nothing about the header's own arithmetic --
        // the test below does that instead.
        EXPECT_NEAR(wrapped, ib::gmst_iau1982_degrees(jd), 1e-7) << "day " << days;
    }
}

TEST(IrbemGmst, TheFractionalDayFoldingIsWorthAFactorOfHundreds) {
    // The header drops the whole days before multiplying by 86400, because only the fractional part
    // survives the reduction modulo a turn. This test is what makes that claim falsifiable: it
    // evaluates the identical series in extended precision and holds the header to 1e-10 deg. The
    // naive `86400 * days` form lands at 4.6e-9 and fails, so the folding is not decoration.
    //
    // long double is 80-bit on x86-64 and 128-bit on aarch64; where it is merely double (MSVC) the
    // reference carries no extra precision and the comparison would be vacuous, so it is skipped
    // rather than silently weakened.
    // `numeric_limits` describes the TYPE, which is not the same as what the hardware will
    // actually deliver at runtime. Under Valgrind the x87 unit is emulated in 64-bit, so
    // `long double` still reports 64 bits of mantissa and silently rounds like a double — the
    // reference then carries no extra precision, the comparison becomes vacuous, and this test
    // fails for a reason that has nothing to do with the header it is testing.
    //
    // So probe the property the test actually depends on: a value that survives only with more
    // than 53 bits of mantissa. This skips correctly under Valgrind, on MSVC where long double IS
    // double, and on any future target that quietly narrows it — without naming any of them.
    {
        volatile long double probe = 1.0L;
        probe += static_cast<long double>(std::numeric_limits<double>::epsilon()) / 2.0L;
        if (std::numeric_limits<long double>::digits < 64 || probe == 1.0L) {
            GTEST_SKIP() << "long double delivers no precision beyond double here, so the "
                            "extended-precision reference would be vacuous";
        }
    }
    double worst = 0.0;
    for (int step = -400; step <= 400; ++step) {
        const double jd = kJdJ2000 + (step * 91.7);  // ~+-100 years, off both midnight and noon
        const long double days = static_cast<long double>(jd) - 2451545.0L;
        const long double t = days / 36525.0L;
        const long double seconds = 67310.54841L + (86400.0L * days) + (8640184.812866L * t) +
                                    (0.093104L * t * t) - (6.2e-6L * t * t * t);
        long double reference = std::fmod(seconds / 240.0L, 360.0L);
        if (reference < 0.0L) reference += 360.0L;
        double difference = std::fabs(ib::gmst_iau1982_degrees(jd) - static_cast<double>(reference));
        if (difference > 180.0) difference = 360.0 - difference;
        worst = std::max(worst, difference);
    }
    EXPECT_LT(worst, 1e-10);
}

TEST(IrbemGmst, AdvancesByOneSiderealDayPerSolarDay) {
    // A sidereal day is shorter than a solar one by the Earth's daily orbital motion, so GMST gains
    // 360.98564736629 deg per day, not 360.
    const double before = ib::gmst_iau1982_degrees(kJdReference);
    const double after = ib::gmst_iau1982_degrees(kJdReference + 1.0);
    double gain = after - before;
    if (gain < 0.0) gain += 360.0;
    EXPECT_NEAR(0.98564736629, gain, 1e-8);
}

TEST(IrbemGmst, HapgoodsTruncationTracksTheIauSeriesToAFewArcseconds) {
    // Hapgood's eq. (2) drops the quadratic term of the same series. Over 1950-2050 the two must
    // therefore stay within a couple of arcseconds of each other; anything larger would mean one of
    // them is wrong rather than merely truncated. Measured worst case over that span: 1.38".
    double worst = 0.0;
    for (int day = 0; day <= 36525; day += 7) {
        const double jd = 2433282.5 + day;  // 1950-01-01 00:00 UT
        double difference = std::fabs(ib::gmst_iau1982_degrees(jd) - ib::gmst_hapgood_degrees(jd));
        if (difference > 180.0) difference = 360.0 - difference;
        worst = std::max(worst, difference);
    }
    EXPECT_LT(worst, 2.0 / 3600.0);
    EXPECT_GT(worst, 0.5 / 3600.0) << "the two series became identical — one of them is not the "
                                      "series it claims to be";
}

TEST(IrbemGmst, TheModelSelectorPicksTheNamedSeries) {
    EXPECT_EQ(ib::gmst_iau1982_degrees(kJdReference),
              ib::gmst_degrees(kJdReference, ib::GmstModel::Iau1982));
    EXPECT_EQ(ib::gmst_hapgood_degrees(kJdReference),
              ib::gmst_degrees(kJdReference, ib::GmstModel::Hapgood1992));
}

TEST(IrbemGmst, IsAlwaysReducedToOneTurn) {
    for (int day = 0; day < 400; ++day) {
        const double jd = kJdReference + (day * 0.917);
        for (const ib::GmstModel model : {ib::GmstModel::Iau1982, ib::GmstModel::Hapgood1992}) {
            const double theta = ib::gmst_degrees(jd, model);
            EXPECT_GE(theta, 0.0);
            EXPECT_LT(theta, 360.0);
        }
    }
    // Both signs of the wrap, exercised directly: an epoch before J2000 drives the raw series
    // negative, which is the branch that would otherwise never run.
    EXPECT_EQ(350.0, ib::detail::wrap_degrees(-10.0));
    EXPECT_EQ(10.0, ib::detail::wrap_degrees(370.0));
    EXPECT_EQ(0.0, ib::detail::wrap_degrees(720.0));
    EXPECT_GE(ib::gmst_iau1982_degrees(2400000.5), 0.0);
    EXPECT_GE(ib::gmst_hapgood_degrees(2400000.5), 0.0);
}

TEST(IrbemGmst, TheHapgoodEpochSplitsAtThePrecedingMidnight) {
    // JD 2457203.0 is 2015-06-29 12:00 UT, i.e. MJD 57202.5 — half a day past midnight.
    const ib::detail::HapgoodEpoch epoch = ib::detail::hapgood_epoch(kJdReference);
    EXPECT_EQ(12.0, epoch.ut_hours);
    EXPECT_EQ((57202.0 - 51544.5) / 36525.0, epoch.centuries);
    // Midnight itself: zero hours, and the same century count as the noon above.
    const ib::detail::HapgoodEpoch midnight = ib::detail::hapgood_epoch(kJdReference - 0.5);
    EXPECT_EQ(0.0, midnight.ut_hours);
    EXPECT_EQ(epoch.centuries, midnight.centuries);
}

// ---- the solar ephemeris ------------------------------------------------------------------------

TEST(IrbemSolarEphemeris, PutsTheSunWhereTheSeasonsSay) {
    // The Sun's ecliptic longitude IS the season: 0 deg at the March equinox, 90 at the June
    // solstice, 180 at the September equinox, 270 at the December one. The 2015 instants below are
    // the published ones; the series is a 0.01 deg approximation, so 0.05 deg is a tight check of
    // the right formula rather than of the almanac.
    struct Instant {
        double jd;
        double longitude_deg;
    };
    const std::array<Instant, 4> seasons{{
        {2457102.4479, 0.0},    // 2015-03-20 22:45 UT, March equinox
        {2457195.1931, 90.0},   // 2015-06-21 16:38 UT, June solstice
        {2457288.8479, 180.0},  // 2015-09-23 08:21 UT, September equinox
        {2457378.7000, 270.0},  // 2015-12-22 04:48 UT, December solstice
    }};
    for (const Instant& season : seasons) {
        const ib::SolarEphemeris sun = ib::solar_ephemeris(season.jd);
        double error = sun.ecliptic_longitude_deg - season.longitude_deg;
        if (error > 180.0) error -= 360.0;
        if (error < -180.0) error += 360.0;
        EXPECT_NEAR(0.0, error, 0.05) << "season at longitude " << season.longitude_deg;
    }
}

TEST(IrbemSolarEphemeris, ObliquityAndAnomalyAreInRange) {
    const ib::SolarEphemeris sun = ib::solar_ephemeris(kJdReference);
    // The obliquity is 23.44 deg and shrinking by 0.013 deg/century.
    EXPECT_NEAR(23.4369863792, sun.obliquity_deg, 1e-9);
    EXPECT_GT(ib::solar_ephemeris(kJdJ2000).obliquity_deg, sun.obliquity_deg);
    // The equation of centre never exceeds ~1.94 deg for the Earth's eccentricity.
    double centre = sun.ecliptic_longitude_deg - sun.mean_longitude_deg;
    if (centre > 180.0) centre -= 360.0;
    if (centre < -180.0) centre += 360.0;
    EXPECT_LT(std::fabs(centre), 1.95);
    EXPECT_GE(sun.mean_anomaly_deg, 0.0);
    EXPECT_LT(sun.mean_anomaly_deg, 360.0);
    EXPECT_GE(sun.mean_longitude_deg, 0.0);
    EXPECT_LT(sun.mean_longitude_deg, 360.0);
}

TEST(IrbemSolarEphemeris, TheSunDirectionIsAUnitVectorOnTheEcliptic) {
    const ib::SolarEphemeris sun = ib::solar_ephemeris(kJdReference);
    const fx::vec3d direction = ib::sun_direction_gei(sun);
    EXPECT_NEAR(1.0, fx::norm(direction), 1e-15);
    // Its ecliptic latitude is zero by construction: rotating by -eps about X must null the Z
    // component. This is the check that the obliquity is applied to the right axis with the right
    // sign — a sign error here would put the Sun 47 degrees from where it belongs at the solstices.
    const double eps = sun.obliquity_deg * ib::detail::kDegToRad;
    const double ecliptic_z =
        (-std::sin(eps) * direction[1]) + (std::cos(eps) * direction[2]);
    EXPECT_NEAR(0.0, ecliptic_z, 1e-16);
}

// ---- the geodipole ------------------------------------------------------------------------------

TEST(IrbemDipole, AnAxialDipoleGivesTheGeographicPoleExactly) {
    // g11 = h11 = 0 leaves a purely axial dipole, so every derived quantity is exact in binary and
    // can be asserted with ==. A sign error in axis_geo would put the pole at the SOUTH pole here,
    // which is precisely the mistake the derivation in the header exists to prevent.
    const ib::DipoleCoefficients axial{-30000.0, 0.0, 0.0};
    EXPECT_EQ(30000.0, axial.moment_nt());
    EXPECT_EQ((fx::vec3d{0.0, 0.0, 1.0}), axial.axis_geo());
    EXPECT_EQ(90.0, axial.north_pole_latitude_deg());
    // The longitude of a pole sitting exactly on the rotation axis is undefined, and IEEE signed
    // zero makes atan2(-0.0, -0.0) pick -180 rather than 0. Asserted as it is rather than papered
    // over: no real IGRF epoch has g11 and h11 both exactly zero, and inventing a branch for a case
    // that cannot occur would be untestable code in the hottest header in the module.
    EXPECT_EQ(-180.0, axial.north_pole_longitude_deg());
}

TEST(IrbemDipole, AThreeFourFiveDipoleIsExact) {
    // 3-4-5: |(g11, h11, g10)| = |(-3000, 0, -4000)| = 5000 exactly, and the axis is (0.6, 0, 0.8)
    // to the last bit the division allows.
    const ib::DipoleCoefficients tilted{-4000.0, -3000.0, 0.0};
    EXPECT_EQ(5000.0, tilted.moment_nt());
    EXPECT_EQ((fx::vec3d{3000.0 / 5000.0, 0.0, 4000.0 / 5000.0}), tilted.axis_geo());
    EXPECT_EQ(0.0, tilted.north_pole_longitude_deg());
    EXPECT_NEAR(53.13010235415598, tilted.north_pole_latitude_deg(), 1e-13);

    // An equatorial dipole: the axis lies in the equatorial plane, 90 deg east.
    const ib::DipoleCoefficients equatorial{0.0, 0.0, -5000.0};
    EXPECT_EQ((fx::vec3d{0.0, 1.0, 0.0}), equatorial.axis_geo());
    EXPECT_EQ(0.0, equatorial.north_pole_latitude_deg());
    EXPECT_EQ(90.0, equatorial.north_pole_longitude_deg());
}

TEST(IrbemDipole, ReproducesThePublishedIgrfPolePositions) {
    // The longitudes match the published geomagnetic pole directly. The latitudes are GEOCENTRIC
    // and the published figures are GEODETIC, which at 80 deg differ by 0.066 deg — so the check
    // is against the geocentric value, and the geodetic conversion is verified separately below so
    // the gap is demonstrated rather than merely asserted.
    EXPECT_NEAR(29867.313239, kDipole2015.moment_nt(), 1e-6);
    EXPECT_NEAR(80.313053, kDipole2015.north_pole_latitude_deg(), 1e-6);
    EXPECT_NEAR(-72.613078, kDipole2015.north_pole_longitude_deg(), 1e-6);
    EXPECT_NEAR(29805.924413, kDipole2020.moment_nt(), 1e-6);
    EXPECT_NEAR(80.589469, kDipole2020.north_pole_latitude_deg(), 1e-6);
    EXPECT_NEAR(-72.679710, kDipole2020.north_pole_longitude_deg(), 1e-6);

    // tan(geodetic) = tan(geocentric)/(1-f)^2, WGS-84 1/f = 298.257223563. Applying it to the
    // geocentric latitudes above lands on the published 80.37 (2015) and 80.65 (2020) N.
    constexpr double kOneMinusF = 1.0 - (1.0 / 298.257223563);
    for (const auto& [dipole, published] :
         std::array<std::pair<ib::DipoleCoefficients, double>, 2>{
             {{kDipole2015, 80.37}, {kDipole2020, 80.65}}}) {
        const double geocentric = dipole.north_pole_latitude_deg() * ib::detail::kDegToRad;
        const double geodetic =
            std::atan(std::tan(geocentric) / (kOneMinusF * kOneMinusF)) * ib::detail::kRadToDeg;
        // The published figures are quoted to 0.01 deg, which is the whole tolerance here.
        EXPECT_NEAR(published, geodetic, 0.01);
    }
    EXPECT_NEAR(1.0, fx::norm(kDipole2015.axis_geo()), 1e-16);
}

// ---- the machinery, exactly ---------------------------------------------------------------------

TEST(IrbemRotations, EveryDirectedPairSelectsTheRightMatrixExactly) {
    const ib::Rotations rotations = permutation_rotations();
    // tolerance 0: with signed permutation matrices every product is exact, so a round trip must
    // return the ORIGINAL BITS. Any wrong-matrix or transposed-matrix bug shows up as a permuted or
    // sign-flipped component, which no tolerance could absorb.
    exercise_every_pair(rotations, 0.0);

    // And each stored direction really is the matrix it was given, not a neighbour.
    EXPECT_EQ(kQuarterZ, (ib::rotation_matrix<Frame::GEO, Frame::GEI>(rotations)));
    EXPECT_EQ(kQuarterX, (ib::rotation_matrix<Frame::GSE, Frame::GEI>(rotations)));
    EXPECT_EQ(kQuarterY, (ib::rotation_matrix<Frame::GSM, Frame::GSE>(rotations)));
    EXPECT_EQ(kHalfZ, (ib::rotation_matrix<Frame::SM, Frame::GSM>(rotations)));
    EXPECT_EQ(kHalfX, (ib::rotation_matrix<Frame::MAG, Frame::GEO>(rotations)));
    EXPECT_EQ(kHalfY, (ib::rotation_matrix<Frame::GSE, Frame::GEO>(rotations)));
    EXPECT_EQ(kMinusQuarterZ, (ib::rotation_matrix<Frame::GSM, Frame::GEO>(rotations)));
    EXPECT_EQ(kCyclic, (ib::rotation_matrix<Frame::SM, Frame::GEO>(rotations)));

    // A worked example, so the sense of the rotation is pinned and not just its consistency: a
    // quarter turn of the frame about Z sends the GEI x axis to the GEO -y axis.
    const ib::Position<Frame::GEI> along_x{fx::vec3d{1.0, 0.0, 0.0}};
    EXPECT_EQ((fx::vec3d{0.0, -1.0, 0.0}), ib::transform<Frame::GEO>(along_x, rotations).v);
}

TEST(IrbemRotations, TransformIsLinear) {
    const ib::Rotations rotations = permutation_rotations();
    const ib::FieldVector<Frame::GEO> internal{fx::vec3d{100.0, -250.0, 3000.0}};
    const ib::FieldVector<Frame::GEO> external{fx::vec3d{-8.0, 16.5, -32.0}};
    // Superposition must survive the frame change: transforming the sum and summing the transforms
    // are the same thing, which is what lets a total field be assembled in whichever frame is
    // convenient. Exact, because the fixture is a permutation.
    const ib::FieldVector<Frame::GSM> sum_then_rotate =
        ib::transform<Frame::GSM>(internal + external, rotations);
    const ib::FieldVector<Frame::GSM> rotate_then_sum =
        ib::transform<Frame::GSM>(internal, rotations) + ib::transform<Frame::GSM>(external, rotations);
    EXPECT_EQ(sum_then_rotate, rotate_then_sum);
}

// ---- the physics --------------------------------------------------------------------------------

TEST(IrbemRotations, EveryPairRoundTripsAtARealEpoch) {
    // 1e-14 is two orders above the ~1e-16 actually observed; it is the ERROR_BUDGET line for
    // coordinate transforms with room for the four-matrix composites.
    exercise_every_pair(reference_rotations(), 1e-14);
}

TEST(IrbemRotations, TheCompositesAreTheProductsTheyClaimToBe) {
    const ib::Rotations r = reference_rotations();
    const fx::mat3d geo_to_gei = fx::transpose(r.gei_to_geo);
    const std::array<std::pair<fx::mat3d, fx::mat3d>, 3> claims{{
        {r.geo_to_gse, r.gei_to_gse * geo_to_gei},
        {r.geo_to_gsm, r.gse_to_gsm * (r.gei_to_gse * geo_to_gei)},
        {r.geo_to_sm, r.gsm_to_sm * (r.gse_to_gsm * (r.gei_to_gse * geo_to_gei))},
    }};
    for (const auto& [stored, product] : claims) {
        for (std::size_t row = 0; row < 3; ++row) {
            for (std::size_t col = 0; col < 3; ++col) {
                EXPECT_NEAR(product(row, col), stored(row, col), 1e-16);
            }
        }
    }
}

TEST(IrbemRotations, TheSmZAxisIsTheDipoleAxis) {
    const ib::Rotations r = reference_rotations();
    // The defining property of SM. A mutually-inverse-but-wrong matrix pair passes every round trip
    // and fails this.
    const ib::Position<Frame::GEO> dipole{r.dipole_geo};
    expect_vec_near(ib::transform<Frame::SM>(dipole, r).v, fx::vec3d{0.0, 0.0, 1.0}, 1e-15, "dipole in SM");
    // ...and MAG's, which shares the axis but keeps the Earth's rotation axis rather than the Sun
    // as its second reference direction.
    expect_vec_near(ib::transform<Frame::MAG>(dipole, r).v, fx::vec3d{0.0, 0.0, 1.0}, 1e-15,
                    "dipole in MAG");
}

TEST(IrbemRotations, TheGsmZAxisLiesInThePlaneOfTheDipoleAndTheSun) {
    const ib::Rotations r = reference_rotations();
    const fx::vec3d sun_geo =
        ib::transform<Frame::GEO>(ib::Position<Frame::GEI>{ib::sun_direction_gei(ib::solar_ephemeris(kJdReference))}, r).v;

    // Stated as GSM sees it: the dipole has no Y component there, so it lies in the GSM X-Z plane
    // together with the Sun, which is +X.
    const fx::vec3d dipole_gsm = ib::transform<Frame::GSM>(ib::Position<Frame::GEO>{r.dipole_geo}, r).v;
    EXPECT_NEAR(0.0, dipole_gsm[1], 1e-15);
    EXPECT_GT(dipole_gsm[2], 0.0) << "the dipole must point along +Z, not -Z";

    // Stated frame-independently, which is the stronger form: the GSM Z axis, the dipole and the
    // Sun direction are coplanar, so their scalar triple product vanishes. Expressed in GEO, so
    // none of the three vectors is trivially an axis of the frame it is written in.
    const fx::vec3d gsm_z_in_geo =
        ib::transform<Frame::GEO>(ib::Position<Frame::GSM>{fx::vec3d{0.0, 0.0, 1.0}}, r).v;
    EXPECT_NEAR(0.0, fx::dot(gsm_z_in_geo, fx::cross(r.dipole_geo, sun_geo)), 1e-15);

    // The other half of the GSM definition: +X is the Sun.
    expect_vec_near(ib::transform<Frame::GSM>(ib::Position<Frame::GEO>{sun_geo}, r).v,
                    fx::vec3d{1.0, 0.0, 0.0}, 1e-15, "sun in GSM");
    expect_vec_near(ib::transform<Frame::GSE>(ib::Position<Frame::GEO>{sun_geo}, r).v,
                    fx::vec3d{1.0, 0.0, 0.0}, 1e-15, "sun in GSE");
    // GSE and GSM share that axis and differ only by a rotation about it, which is what makes T3 a
    // pure X rotation.
    EXPECT_EQ(1.0, r.gse_to_gsm(0, 0));
    EXPECT_EQ(0.0, r.gse_to_gsm(0, 1));
    EXPECT_EQ(0.0, r.gse_to_gsm(1, 0));
}

TEST(IrbemRotations, TheMagYAxisIsPerpendicularToTheEarthsRotationAxis) {
    const ib::Rotations r = reference_rotations();
    // MAG's Y is (z_GEO x dipole) normalized, so written in GEO it has an exactly zero Z component
    // — exactly, because the construction never puts anything there.
    const fx::vec3d mag_y_in_geo =
        ib::transform<Frame::GEO>(ib::Position<Frame::MAG>{fx::vec3d{0.0, 1.0, 0.0}}, r).v;
    EXPECT_EQ(0.0, mag_y_in_geo[2]);
    const fx::vec3d expected = fx::normalize(fx::cross(fx::vec3d{0.0, 0.0, 1.0}, r.dipole_geo));
    expect_vec_near(mag_y_in_geo, expected, 1e-15, "MAG Y in GEO");
}

TEST(IrbemRotations, TheGeiToGeoRotationIsTheSiderealTime) {
    const ib::Rotations r = reference_rotations();
    // The Greenwich meridian sits at right ascension GMST, so an inertial direction at that right
    // ascension has zero geographic longitude. This is the only thing T1 asserts, and it fixes both
    // the magnitude and the sense of the rotation.
    const double theta = r.gmst_deg * ib::detail::kDegToRad;
    const ib::Position<Frame::GEI> greenwich{fx::vec3d{std::cos(theta), std::sin(theta), 0.0}};
    expect_vec_near(ib::transform<Frame::GEO>(greenwich, r).v, fx::vec3d{1.0, 0.0, 0.0}, 1e-15,
                    "Greenwich");
    // The rotation is about Z, so the polar axis is shared between GEI and GEO.
    expect_vec_near(ib::transform<Frame::GEO>(ib::Position<Frame::GEI>{fx::vec3d{0.0, 0.0, 1.0}}, r).v,
                    fx::vec3d{0.0, 0.0, 1.0}, 0.0, "polar axis");
}

TEST(IrbemRotations, TheDipoleTiltIsTheComplementOfTheDipoleSunAngle) {
    const ib::Rotations r = reference_rotations();
    const fx::vec3d sun_geo =
        ib::transform<Frame::GEO>(ib::Position<Frame::GEI>{ib::sun_direction_gei(ib::solar_ephemeris(kJdReference))}, r).v;
    // An independent derivation of mu that never touches T3 or T4: the tilt is by definition 90 deg
    // minus the angle between the dipole axis and the Earth-Sun line.
    const double angle = std::acos(fx::dot(r.dipole_geo, sun_geo)) * ib::detail::kRadToDeg;
    EXPECT_NEAR(90.0 - angle, r.dipole_tilt_deg, 1e-12);
    // Regression golden for the whole chain at the reference epoch, so an orchestrator can compare
    // it against the oracle.
    EXPECT_NEAR(25.6428938132, r.dipole_tilt_deg, 1e-9);
    EXPECT_NEAR(-13.5235408757, r.gsm_dipole_angle_deg, 1e-9);
    EXPECT_NEAR(97.2534261320, r.gmst_deg, 1e-9);
    EXPECT_NEAR(97.4472821808, r.sun_ecliptic_longitude_deg, 1e-9);
    EXPECT_EQ(kJdReference, r.jd_ut1);
    expect_vec_near(r.dipole_geo, kDipole2015.axis_geo(), 0.0, "stored dipole");
    // dipole_gse is the vector psi and mu are read off; it must be the dipole carried through.
    expect_vec_near(r.dipole_gse, r.geo_to_gse * r.dipole_geo, 0.0, "dipole in GSE");
}

TEST(IrbemRotations, TheTiltSwingsOverAYearByTheObliquityPlusThePoleColatitude) {
    // The dipole tilt is the sum of two independent leans: the pole's own colatitude, which swings
    // through +-itself once a day as the Earth turns, and the obliquity, which swings through
    // +-itself once a year. Their extremes coincide, so the annual maximum is exactly the sum. This
    // pins the tilt's amplitude against pure geometry, with no reference value involved.
    const double expected_extreme =
        (90.0 - kDipole2015.north_pole_latitude_deg()) + ib::solar_ephemeris(kJdReference).obliquity_deg;
    double lowest = 360.0;
    double highest = -360.0;
    for (int hour = 0; hour < 366 * 24; ++hour) {
        const double tilt = ib::Rotations::at(2457023.5 + (hour / 24.0), kDipole2015).dipole_tilt_deg;
        lowest = std::min(lowest, tilt);
        highest = std::max(highest, tilt);
    }
    // 0.03 deg covers the hourly sampling: the tilt moves ~0.8 deg/hour near the extreme, and the
    // extremum of a smooth maximum sampled at that step is missed by ~0.02 deg.
    EXPECT_NEAR(expected_extreme, highest, 0.03);
    EXPECT_NEAR(-expected_extreme, lowest, 0.03);
}

TEST(IrbemRotations, TheGmstModelChoiceMovesTheFramesButOnlyByArcseconds) {
    const ib::Rotations iau = ib::Rotations::at(kJdReference, kDipole2015, ib::GmstModel::Iau1982);
    const ib::Rotations hapgood =
        ib::Rotations::at(kJdReference, kDipole2015, ib::GmstModel::Hapgood1992);
    EXPECT_NE(iau.gmst_deg, hapgood.gmst_deg);
    // The two differ only through T1, so a GEO position lands within an arcsecond either way.
    const ib::Position<Frame::GEO> point{fx::vec3d{4.0, -3.0, 2.0}};
    const fx::vec3d a = ib::transform<Frame::GSM>(point, iau).v;
    const fx::vec3d b = ib::transform<Frame::GSM>(point, hapgood).v;
    const double separation_deg =
        std::asin(fx::norm(a - b) / fx::norm(point.v)) * ib::detail::kRadToDeg;
    EXPECT_LT(separation_deg, 2.0 / 3600.0);
    EXPECT_GT(separation_deg, 0.0);
}

// ---- allocation ---------------------------------------------------------------------------------

// The counter lives in tests/alloc_counter.cpp — one definition for the whole binary. The
// transforms are the hottest path in the library, and a heap touch there would dominate
// everything else they do.

TEST(IrbemRotations, BuildingAndApplyingRotationsTouchesNoHeap) {
    volatile double sink = 0.0;
    const std::size_t before = cheatah_space_test::allocation_count();
    for (int i = 0; i < 64; ++i) {
        const ib::Rotations r = ib::Rotations::at(kJdReference + (i * 0.25), kDipole2015);
        const ib::Position<Frame::GEO> p{fx::vec3d{1.0 + i, 2.0, 3.0}};
        const ib::FieldVector<Frame::GSM> f{
            ib::transform<Frame::GSM>(ib::FieldVector<Frame::GEO>{p.v}, r)};
        sink = sink + ib::transform<Frame::SM>(p, r).v[0] + f.v[2] + r.dipole_tilt_deg;
    }
    EXPECT_EQ(before, cheatah_space_test::allocation_count());
    EXPECT_NE(0.0, sink);
}
