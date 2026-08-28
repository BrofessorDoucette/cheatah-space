// Unit tests for space.irbem's non-Cartesian coordinate conversions — GDZ (geodetic), SPH
// (geographic spherical), RLL (geocentric radius with geodetic latitude), and the Cartesian pair.
//
// Three things this suite has to prove, because none of them is obvious from reading the code:
//
//   1. That the ellipsoid constants are the ones IRBEM uses. WGS84 is named outright in IRBEM's
//      documentation; the Re unit (6371.2 km) is not printed anywhere, and is pinned here by the
//      two closed-form identities it participates in.
//   2. That Bowring's iteration has actually converged — not "looks close", but that running it
//      four times longer changes no bit of the answer anywhere on the sweep.
//   3. That the singular directions behave. The poles, the equator, the antimeridian and the
//      geocentre are where every implementation of this conversion goes wrong, including IRBEM's
//      (see PoleWhereIrbemIsWrong below), so they are tested by name rather than left to a sweep.
//
// Assertions are exact `==` wherever the arithmetic permits it — the cardinal angles convert
// exactly (cos 0° = 1, atan2(1,0)·180/π = 90 with no rounding), so the 3-4-5-style cases and the
// equatorial/polar closed forms carry no tolerance at all. Where a transcendental makes that
// impossible the tolerance is stated in ulps of the quantity involved and justified on the line.
//
// The sweep assertions are deliberately just above the measured worst case, so any regression in
// the iteration count, the height formula, or the constants trips them.

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>

#include "alloc_counter.hpp"
#include "space/irbem/coords_geodetic.hpp"
#include "space/irbem/frames.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;
namespace det = cheatah::space::irbem::detail;

using ib::Frame;

namespace {

/// A GDZ position from its three named components — the aggregate's component order is
/// (altitude, latitude, longitude), which is easy to transpose by accident.
ib::Position<Frame::GDZ> gdz(double altitude_km, double latitude_deg, double longitude_deg) {
    return ib::Position<Frame::GDZ>{fx::vec3d{altitude_km, latitude_deg, longitude_deg}};
}

/// A GEO position from its three Cartesian components, in Earth radii.
ib::Position<Frame::GEO> geo(double x, double y, double z) {
    return ib::Position<Frame::GEO>{fx::vec3d{x, y, z}};
}

/// An SPH position — (radius Re, geocentric latitude deg, east longitude deg).
ib::Position<Frame::SPH> sph(double radius_re, double latitude_deg, double longitude_deg) {
    return ib::Position<Frame::SPH>{fx::vec3d{radius_re, latitude_deg, longitude_deg}};
}

/// An RLL position — (geocentric radius Re, *geodetic* latitude deg, east longitude deg).
ib::Position<Frame::RLL> rll(double radius_re, double latitude_deg, double longitude_deg) {
    return ib::Position<Frame::RLL>{fx::vec3d{radius_re, latitude_deg, longitude_deg}};
}

/// The altitudes the sweeps run over: below the ellipsoid, on it, and out past geosynchronous.
constexpr std::array<double, 10> sweep_altitudes_km{-100.0, 0.0,      1.0,     10.0,     100.0,
                                                    1000.0, 6371.2,  20000.0, 36000.0, 60000.0};

}  // namespace

// ---- the constants ------------------------------------------------------------------------------

// WGS84's two defining constants, and the three that follow from them. Written as literals so that
// a typo in the derivation cannot be self-consistent.
static_assert(det::wgs84_semi_major_km == 6378.137);
static_assert(det::wgs84_inverse_flattening == 298.257223563);
static_assert(det::wgs84_semi_minor_km == 6356.7523142451792);
static_assert(det::wgs84_eccentricity_squared == 0.0066943799901413165);
static_assert(det::wgs84_second_eccentricity_squared == 0.0067394967422764341);

// The Re unit is 6371.2 km, and it is neither axis of the ellipsoid. The two differences below are
// the ones the oracle pins (see the header): a point at 1 Re on the equator is 6.937 km *below* the
// ellipsoid, and one at 1 Re over the pole is 14.448 km above it.
static_assert(det::earth_radius_km == 6371.2);
static_assert(det::earth_radius_km - det::wgs84_semi_major_km == -6.9369999999998981);
static_assert(det::earth_radius_km - det::wgs84_semi_minor_km == 14.447685754820668);

static_assert(det::bowring_iterations == 4);
// Even, and that is load-bearing — see BowringIteration.OddPassCountsSitOnALimitCycle.
static_assert(det::bowring_iterations % 2 == 0);

// ---- exactly-representable cases: no tolerance --------------------------------------------------

TEST(GdzToGeo, EquatorIsTheSemiMajorAxisExactly) {
    // lat = 0 makes sin = 0 and cos = 1 with no rounding, so W = 1, N = a, and the whole transform
    // collapses to a/Re. Written as a LITERAL, not as `a/Re`: an expectation phrased in terms of
    // the constants under test cannot detect a wrong constant, and Re is a convention that nothing
    // internal can derive. (A perturbation run that changed Re to the semi-major axis passed every
    // self-referential assertion in this file — hence these three literals.)
    const auto p = ib::gdz_to_geo(gdz(0.0, 0.0, 0.0));
    EXPECT_EQ(p.v[0], 1.001088805876444);  // 6378.137 / 6371.2
    EXPECT_EQ(p.v[0], det::wgs84_semi_major_km / det::earth_radius_km);
    EXPECT_EQ(p.v[1], 0.0);
    EXPECT_EQ(p.v[2], 0.0);
}

TEST(GdzToGeo, AltitudeAddsAlongTheEquatorialNormalExactly) {
    // On the equator the ellipsoid normal is radial, so altitude adds straight onto a.
    const auto p = ib::gdz_to_geo(gdz(1000.0, 0.0, 0.0));
    EXPECT_EQ(p.v[0], (det::wgs84_semi_major_km + 1000.0) / det::earth_radius_km);
    EXPECT_EQ(p.v[2], 0.0);
}

TEST(GdzToGeo, LongitudeQuarterTurnsArePermutationsExactly) {
    // cos/sin of 0° and 90° are exact, so a quarter turn is a component swap and nothing else.
    const auto east = ib::gdz_to_geo(gdz(0.0, 0.0, 90.0));
    const auto reference = ib::gdz_to_geo(gdz(0.0, 0.0, 0.0));
    EXPECT_EQ(east.v[1], reference.v[0]);
    // cos(90° in radians) is 6.12e-17, not zero: the residual x is that times the equatorial radius.
    EXPECT_LT(std::abs(east.v[0]), 1e-16);
    EXPECT_EQ(east.v[2], 0.0);
}

TEST(GdzToGeo, PoleIsTheSemiMinorAxis) {
    // sin(90°) is exactly 1 in double, but W = sqrt(1-e²) is not bit-identical to b/a, so this one
    // carries a tolerance: 1e-16 Re is 0.6 nm, and ~2 ulp of the quantity being compared.
    const auto north = ib::gdz_to_geo(gdz(0.0, 90.0, 0.0));
    EXPECT_NEAR(north.v[2], det::wgs84_semi_minor_km / det::earth_radius_km, 1e-16);
    EXPECT_LT(std::abs(north.v[0]), 1e-16);
    EXPECT_EQ(north.v[1], 0.0);

    const auto south = ib::gdz_to_geo(gdz(0.0, -90.0, 0.0));
    EXPECT_EQ(south.v[2], -north.v[2]);
}

TEST(GeoToGdz, UnitEquatorialRadiusIsBelowTheEllipsoidExactly) {
    // The identity that pins the Re unit: 1 Re on the equator sits at altitude Re - a. Literal, for
    // the reason given in GdzToGeo.EquatorIsTheSemiMajorAxisExactly — and this is the exact number
    // IRBEM's geo_gdz(1,0,0) returns, which is where Re = 6371.2 was read off in the first place.
    const auto g = ib::geo_to_gdz(geo(1.0, 0.0, 0.0));
    EXPECT_EQ(g.radius(), -6.9369999999998981);
    EXPECT_EQ(g.radius(), det::earth_radius_km - det::wgs84_semi_major_km);
    EXPECT_EQ(g.latitude(), 0.0);
    EXPECT_EQ(g.longitude(), 0.0);
}

TEST(GeoToGdz, AntimeridianIsExactlyOneEightyNotMinusOneEighty) {
    // atan2(+0, -x) is +π, so a point on the -x axis reports +180. Worth pinning: the sign of the
    // zero decides it, and a "simplification" to atan(y/x) would lose it.
    const auto g = ib::geo_to_gdz(geo(-1.0, 0.0, 0.0));
    EXPECT_EQ(g.longitude(), 180.0);
    EXPECT_EQ(g.latitude(), 0.0);
    EXPECT_EQ(g.radius(), det::earth_radius_km - det::wgs84_semi_major_km);
}

TEST(GeoToGdz, SouthernAntimeridianReportsMinusOneEighty) {
    // ...and the mirrored point, with a negative zero in y, reports -180. Both are correct; the
    // pair is what makes the branch cut explicit.
    const auto g = ib::geo_to_gdz(geo(-1.0, -0.0, 0.0));
    EXPECT_EQ(g.longitude(), -180.0);
}

TEST(CarToSph, ThreeFourTwelveIsExact) {
    // 9 + 16 + 144 = 169 and sqrt(169) = 13, all integer-exact in double.
    const auto s = ib::car_to_sph(geo(3.0, 4.0, 12.0));
    EXPECT_EQ(s.radius(), 13.0);
    // The angles are transcendental; 1e-13 deg is 11 µm at the Earth's surface and ~500 ulp, which
    // is loose enough not to be an accident of libm and tight enough to catch a swapped component.
    EXPECT_NEAR(s.latitude(), std::asin(12.0 / 13.0) * det::deg_per_rad, 1e-13);
    EXPECT_NEAR(s.longitude(), std::atan2(4.0, 3.0) * det::deg_per_rad, 1e-13);
}

TEST(CarToSph, CardinalAxesAreExact) {
    const auto x_axis = ib::car_to_sph(geo(1.0, 0.0, 0.0));
    EXPECT_EQ(x_axis.radius(), 1.0);
    EXPECT_EQ(x_axis.latitude(), 0.0);   // 90 - atan2(1,0)·180/π, and atan2(1,0)·180/π is exactly 90
    EXPECT_EQ(x_axis.longitude(), 0.0);

    const auto north = ib::car_to_sph(geo(0.0, 0.0, 2.0));
    EXPECT_EQ(north.radius(), 2.0);
    EXPECT_EQ(north.latitude(), 90.0);
    EXPECT_EQ(north.longitude(), 0.0);

    const auto south = ib::car_to_sph(geo(0.0, 0.0, -2.0));
    EXPECT_EQ(south.radius(), 2.0);
    EXPECT_EQ(south.latitude(), -90.0);

    const auto west = ib::car_to_sph(geo(-1.0, 0.0, 0.0));
    EXPECT_EQ(west.longitude(), 180.0);
    EXPECT_EQ(west.latitude(), 0.0);
}

TEST(CarToSph, OriginTakesTheLimitAlongPlusZ) {
    // atan2(0, 0) is 0, so the colatitude is 0 and the latitude +90. A point with no direction has
    // no honest latitude; this pins the documented choice, which is also the one IRBEM makes.
    const auto s = ib::car_to_sph(geo(0.0, 0.0, 0.0));
    EXPECT_EQ(s.radius(), 0.0);
    EXPECT_EQ(s.latitude(), 90.0);
    EXPECT_EQ(s.longitude(), 0.0);
}

TEST(CarToSph, LongitudeOnThePolarAxisFollowsTheSignsOfTheZeros) {
    // The one place this file's convention differs from IRBEM's, pinned so the difference is a
    // decision rather than a surprise. Longitude is undefined on the axis; atan2 answers from the
    // signs of the zeros, so all four sign combinations are distinct and all four are legitimate.
    // IRBEM reports 0 for every one of them.
    EXPECT_EQ(ib::car_to_sph(geo(0.0, 0.0, 1.5)).longitude(), 0.0);
    EXPECT_EQ(ib::car_to_sph(geo(-0.0, -0.0, 1.5)).longitude(), -180.0);
    EXPECT_EQ(ib::car_to_sph(geo(-0.0, 0.0, 1.5)).longitude(), 180.0);
    EXPECT_EQ(ib::car_to_sph(geo(0.0, -0.0, 1.5)).longitude(), -0.0);
    // Latitude is unaffected — it is the one polar quantity that IS defined.
    EXPECT_EQ(ib::car_to_sph(geo(-0.0, -0.0, 1.5)).latitude(), 90.0);
}

TEST(SphToCar, PolarLongitudeSurvivesTheRoundTrip) {
    // ...and the reason atan2 is left alone: this library's own polar points carry x, y ~ 1e-17
    // rather than zero, so the longitude does come back. A hard-coded 0° would destroy it.
    for (const double longitude : {-180.0, -90.0, 0.0, 90.0, 179.5}) {
        const auto back = ib::car_to_sph(ib::sph_to_car(sph(1.5, 90.0, longitude)));
        EXPECT_EQ(back.longitude(), longitude);
        EXPECT_EQ(back.latitude(), 90.0);
    }
}

TEST(SphToCar, CardinalDirectionsAreExact) {
    EXPECT_EQ(ib::sph_to_car(sph(0.5, 0.0, 0.0)).v, (fx::vec3d{0.5, 0.0, 0.0}));
    EXPECT_EQ(ib::sph_to_car(sph(2.0, 0.0, 0.0)).v, (fx::vec3d{2.0, 0.0, 0.0}));

    // The pole: sin(90°) is exactly 1, cos is 6.12e-17, so z is exact and x/y are that residual.
    const auto north = ib::sph_to_car(sph(4.0, 90.0, 0.0));
    EXPECT_EQ(north.v[2], 4.0);
    EXPECT_LT(std::abs(north.v[0]), 1e-15);
    EXPECT_EQ(north.v[1], 0.0);
}

TEST(SphToCar, IsTheExactInverseOfCarToSphOnTheThreeFourTwelvePoint) {
    const auto s = ib::car_to_sph(geo(3.0, 4.0, 12.0));
    const auto back = ib::sph_to_car(s);
    // Four transcendentals round-trip; 1e-14 Re is 64 nm and ~50 ulp of 13.
    EXPECT_NEAR(back.v[0], 3.0, 1e-14);
    EXPECT_NEAR(back.v[1], 4.0, 1e-14);
    EXPECT_NEAR(back.v[2], 12.0, 1e-14);
}

TEST(RllToGdz, EquatorAndPoleAreTheClosedFormExactly) {
    // sinφ = 0 kills every ellipsoid term, so the equatorial case is r·Re - a with no rounding.
    const auto equator = ib::rll_to_gdz(rll(1.0, 0.0, 0.0));
    EXPECT_EQ(equator.radius(), det::earth_radius_km - det::wgs84_semi_major_km);
    EXPECT_EQ(equator.latitude(), 0.0);
    EXPECT_EQ(equator.longitude(), 0.0);

    // At the pole cosφ = 6.12e-17 so the offset term is ~1e-13 km and the subtraction under the
    // square root is exact to the last bit; 1e-11 km is 10 nm.
    const auto pole = ib::rll_to_gdz(rll(1.0, 90.0, 0.0));
    EXPECT_NEAR(pole.radius(), det::earth_radius_km - det::wgs84_semi_minor_km, 1e-11);
}

TEST(RllToGdz, PassesLatitudeAndLongitudeThrough) {
    // RLL and GDZ share both angles by definition — the conversion touches only the radial
    // component. Values chosen to be exactly representable so this is an `==`.
    const auto g = ib::rll_to_gdz(rll(1.5, -37.5, 175.25));
    EXPECT_EQ(g.latitude(), -37.5);
    EXPECT_EQ(g.longitude(), 175.25);
}

TEST(RllToGdz, DegenerateInteriorRadiusClampsToTheFootOfTheNormal) {
    // The discriminant r² - N²e⁴sin²φcos²φ goes negative inside ~0.0034 Re near 45°: no ellipsoid
    // normal at that latitude passes that close to the centre. Clamped, the answer is -aW, the
    // altitude of the foot of the normal — and it is the same for every radius inside the clamp,
    // which is exactly what makes it identifiable.
    const double w = det::ellipsoid_w(std::sin(45.0 * det::rad_per_deg));
    const auto tiny = ib::rll_to_gdz(rll(0.001, 45.0, 0.0));
    const auto centre = ib::rll_to_gdz(rll(0.0, 45.0, 0.0));
    EXPECT_EQ(tiny.radius(), -det::wgs84_semi_major_km * w);
    EXPECT_EQ(centre.radius(), tiny.radius());

    // ...and just outside the clamp the answer moves again, so the branch is not swallowing a
    // whole neighbourhood of real inputs.
    const auto outside = ib::rll_to_gdz(rll(0.01, 45.0, 0.0));
    EXPECT_GT(outside.radius(), tiny.radius());
}

// ---- the singular directions --------------------------------------------------------------------

TEST(GeoToGdz, PoleWhereIrbemIsWrong) {
    // 1 Re straight up the polar axis. p is exactly zero, which is where the `p/cosφ - N` height
    // form divides zero by zero — IRBEM returns 0 km here, and returns 3178.376 km for (0,0,1.5)
    // where the truth is 3200.048. The projection form used here is exact instead, and agrees with
    // IRBEM's own gdz_geo(90,0,0) inverse and with its p→0 limit (geo_gdz(1e-6,0,1) = 14.447690).
    const auto north = ib::geo_to_gdz(geo(0.0, 0.0, 1.0));
    EXPECT_EQ(north.latitude(), 90.0);
    EXPECT_EQ(north.longitude(), 0.0);
    EXPECT_NEAR(north.radius(), 14.447685754820668, 1e-11);  // Re - b, literal
    EXPECT_NEAR(north.radius(), det::earth_radius_km - det::wgs84_semi_minor_km, 1e-11);

    const auto south = ib::geo_to_gdz(geo(0.0, 0.0, -1.0));
    EXPECT_EQ(south.latitude(), -90.0);
    EXPECT_NEAR(south.radius(), det::earth_radius_km - det::wgs84_semi_minor_km, 1e-11);

    const auto high = ib::geo_to_gdz(geo(0.0, 0.0, 1.5));
    EXPECT_NEAR(high.radius(), (1.5 * det::earth_radius_km) - det::wgs84_semi_minor_km, 1e-11);
}

TEST(GeoToGdz, ApproachingThePoleIsContinuous) {
    // The pole is only interesting if the limit approaching it agrees. Step p down six decades and
    // watch the altitude converge on the polar value rather than wander.
    const double polar_altitude = det::earth_radius_km - det::wgs84_semi_minor_km;
    double previous = std::abs(ib::geo_to_gdz(geo(1e-3, 0.0, 1.0)).radius() - polar_altitude);
    for (const double p : {1e-4, 1e-6, 1e-9, 1e-12}) {
        const double error = std::abs(ib::geo_to_gdz(geo(p, 0.0, 1.0)).radius() - polar_altitude);
        EXPECT_LE(error, previous) << "at p = " << p;
        previous = error;
    }
    EXPECT_LT(previous, 1e-11);
}

TEST(GeoToGdz, GeocentreIsDegenerateAndSaysSo) {
    // The geocentre has no geodetic latitude at all, and the iteration says so out loud: with p and
    // z both zero, φ alternates between π and ~0 with the parity of the pass count. IRBEM, running
    // an odd count, reports 180°; this runs an even one and reports ~1e-46°. Neither is meaningful.
    //
    // The *altitude* is meaningful and is identical on both branches — h = p cosφ + z sinφ - aW
    // collapses to -a whichever φ comes out — so that is the assertion with teeth. The latitude is
    // only bounded, not pinned to a value, because pinning a meaningless number would make a change
    // to bowring_iterations look like a regression when it is not.
    const auto g = ib::geo_to_gdz(geo(0.0, 0.0, 0.0));
    EXPECT_EQ(g.radius(), -det::wgs84_semi_major_km);
    EXPECT_EQ(g.longitude(), 0.0);
    EXPECT_TRUE(std::isfinite(g.latitude()));
    // Whichever branch the parity lands on, the latitude is one of the two degenerate values.
    EXPECT_TRUE(std::abs(g.latitude()) < 1e-40 || std::abs(g.latitude()) == 180.0)
        << "latitude was " << g.latitude();
}

TEST(GdzToGeo, ExtremeLongitudesAreMirrorImages) {
    // ±180 must not collapse onto each other: the y components differ in sign even though both are
    // ~1e-16 of the equatorial radius.
    const auto east = ib::gdz_to_geo(gdz(0.0, 45.0, 180.0));
    const auto west = ib::gdz_to_geo(gdz(0.0, 45.0, -180.0));
    EXPECT_EQ(east.v[0], west.v[0]);
    EXPECT_EQ(east.v[2], west.v[2]);
    EXPECT_EQ(east.v[1], -west.v[1]);
    EXPECT_GT(east.v[1], 0.0);
}

TEST(GeoToGdz, PointsOnTheEllipsoidReturnZeroAltitude) {
    // The defining property: map the ellipsoid surface out and back, and the altitude must return
    // to zero. Swept over every half-degree of latitude and a spread of longitudes including the
    // antimeridian.
    double worst = 0.0;
    double worst_latitude = 0.0;
    for (int i = 0; i <= 360; ++i) {
        const double latitude = -90.0 + (0.5 * i);
        for (const double longitude : {-180.0, -73.25, 0.0, 37.5, 180.0}) {
            const double altitude =
                ib::geo_to_gdz(ib::gdz_to_geo(gdz(0.0, latitude, longitude))).radius();
            if (std::abs(altitude) > worst) {
                worst = std::abs(altitude);
                worst_latitude = latitude;
            }
        }
    }
    std::printf("[ on-ellipsoid ] worst |altitude| = %.3g km (%.3g mm) at latitude %g\n", worst,
                worst * 1e6, worst_latitude);
    // Measured worst is 2.7e-12 km (2.7 nm); the cap is one decade above it.
    EXPECT_LT(worst, 1e-11);
}

// ---- convergence ---------------------------------------------------------------------------------

TEST(BowringIteration, HasReachedItsFixedPointEverywhereOnTheSweep) {
    // The convergence criterion, made testable: the production pass count must give bit-identical
    // results to running the same iteration four times longer. Not "close" — identical. If a single
    // point on the sweep still moved, this fails.
    //
    // It also shows the criterion is not vacuous: one pass does NOT reach the fixed point (the
    // check below requires it to differ somewhere), so the loop is doing real work.
    std::size_t moved_after_one_pass = 0;
    double worst_one_pass_deg = 0.0;
    for (int i = 0; i <= 360; ++i) {
        const double latitude = -90.0 + (0.5 * i);
        for (const double altitude : sweep_altitudes_km) {
            const auto v = ib::gdz_to_geo(gdz(altitude, latitude, 37.5)).v;
            const double p_km = std::hypot(v[0], v[1]) * det::earth_radius_km;
            const double z_km = v[2] * det::earth_radius_km;

            const double converged = det::geodetic_latitude(p_km, z_km, 12);
            EXPECT_EQ(det::geodetic_latitude(p_km, z_km, det::bowring_iterations), converged)
                << "not converged at latitude " << latitude << ", altitude " << altitude;

            const double one_pass = det::geodetic_latitude(p_km, z_km, 1);
            if (one_pass != converged) {
                ++moved_after_one_pass;
                worst_one_pass_deg =
                    std::max(worst_one_pass_deg, std::abs(one_pass - converged) * det::deg_per_rad);
            }
        }
    }
    std::printf("[ convergence  ] single Bowring pass differs at %zu of %d points, worst %.3g deg\n",
                moved_after_one_pass, 361 * 10, worst_one_pass_deg);
    EXPECT_GT(moved_after_one_pass, 0u) << "a single pass already converged — the test proves "
                                           "nothing about the iteration";
}

TEST(BowringIteration, OddPassCountsSitOnALimitCycle) {
    // Why bowring_iterations is 4 and not 3. At a handful of points the iteration does not settle
    // on one double but alternates between two adjacent ones, so an odd count lands 1 ulp away from
    // where every even count lands. This measures that directly: three passes must differ somewhere
    // (otherwise the even/odd distinction is imaginary and the extra pass is waste), and the
    // difference must be at the 1-ulp level (otherwise it is not a limit cycle but a real
    // non-convergence, and four passes would not be enough either).
    std::size_t odd_differs = 0;
    double worst_odd = 0.0;
    for (int i = 0; i <= 360; ++i) {
        const double latitude = -90.0 + (0.5 * i);
        for (const double altitude : sweep_altitudes_km) {
            const auto v = ib::gdz_to_geo(gdz(altitude, latitude, 37.5)).v;
            const double p_km = std::hypot(v[0], v[1]) * det::earth_radius_km;
            const double z_km = v[2] * det::earth_radius_km;
            const double converged = det::geodetic_latitude(p_km, z_km, 16);
            const double odd = det::geodetic_latitude(p_km, z_km, 3);
            if (odd != converged) {
                ++odd_differs;
                worst_odd = std::max(worst_odd, std::abs(odd - converged));
            }
        }
    }
    std::printf("[ limit cycle  ] 3 passes differ from 16 at %zu of %d points, worst %.3g rad\n",
                odd_differs, 361 * 10, worst_odd);
    EXPECT_GT(odd_differs, 0u) << "no odd/even split — bowring_iterations could be 3";
    EXPECT_LT(worst_odd, 1e-15) << "the gap is larger than a rounding artefact";
}

TEST(BowringIteration, AlwaysRunsAtLeastOnePass) {
    // The loop is written `for (pass = 1; pass < iterations; ++pass)` after an unconditional first
    // pass, so zero and one behave identically. Pinned because the alternative — returning an
    // uninitialised or seed-only latitude for 0 — would be a silent wrong answer.
    const auto v = ib::gdz_to_geo(gdz(500.0, 51.5, -0.125)).v;
    const double p_km = std::hypot(v[0], v[1]) * det::earth_radius_km;
    const double z_km = v[2] * det::earth_radius_km;
    EXPECT_EQ(det::geodetic_latitude(p_km, z_km, 0), det::geodetic_latitude(p_km, z_km, 1));
}

TEST(EllipsoidW, IsBoundedByTheAxisRatio) {
    // W runs from b/a at the pole to 1 at the equator, and never reaches zero — which is why
    // dividing by it needs no guard.
    EXPECT_EQ(det::ellipsoid_w(0.0), 1.0);
    EXPECT_NEAR(det::ellipsoid_w(1.0), det::wgs84_semi_minor_km / det::wgs84_semi_major_km, 1e-16);
    EXPECT_EQ(det::ellipsoid_w(-1.0), det::ellipsoid_w(1.0));
    for (int i = 0; i <= 90; ++i) {
        const double w = det::ellipsoid_w(std::sin(i * det::rad_per_deg));
        EXPECT_GT(w, 0.99);
        EXPECT_LE(w, 1.0);
    }
}

TEST(BowringLatitude, SeededFromTheTrueParametricLatitudeIsExactInOnePass) {
    // Bowring's relation is an identity when β is the true parametric latitude, so feeding it that
    // β — obtained from the converged geodetic latitude via tan β = (1-f) tan φ — must reproduce
    // the geodetic latitude in a single evaluation. That is the claim the iteration rests on.
    for (const double latitude : {-89.0, -45.0, -1.0, 0.0, 12.5, 45.0, 78.25, 89.5}) {
        const auto v = ib::gdz_to_geo(gdz(250.0, latitude, 20.0)).v;
        const double p_km = std::hypot(v[0], v[1]) * det::earth_radius_km;
        const double z_km = v[2] * det::earth_radius_km;
        const double phi = det::geodetic_latitude(p_km, z_km, 12);
        const double beta =
            std::atan2((1.0 - det::wgs84_flattening) * std::sin(phi), std::cos(phi));
        EXPECT_EQ(det::bowring_latitude(beta, p_km, z_km), phi) << "at latitude " << latitude;
    }
}

// ---- the dense sweeps ---------------------------------------------------------------------------

TEST(RoundTrip, GeodeticToCartesianAndBackOverTheWholeRange) {
    // 361 latitudes × 145 longitudes × 10 altitudes = 523 450 points, including latitude exactly
    // ±90 and 0, longitude exactly ±180 and 0, altitude below the ellipsoid and out past
    // geosynchronous.
    double worst_altitude = 0.0;
    double worst_latitude = 0.0;
    double worst_longitude = 0.0;
    double at_latitude = 0.0;
    double at_altitude = 0.0;

    for (int i = 0; i <= 360; ++i) {
        const double latitude = -90.0 + (0.5 * i);
        for (int j = 0; j <= 144; ++j) {
            const double longitude = -180.0 + (2.5 * j);
            for (const double altitude : sweep_altitudes_km) {
                const auto back = ib::geo_to_gdz(ib::gdz_to_geo(gdz(altitude, latitude, longitude)));
                const double d_alt = std::abs(back.radius() - altitude);
                const double d_lat = std::abs(back.latitude() - latitude);
                // Longitude is undefined on the axis, so the pole rows are excluded from that one
                // statistic only — and they are the only rows excluded from anything.
                const double d_lon =
                    std::abs(latitude) == 90.0 ? 0.0 : std::abs(back.longitude() - longitude);
                if (d_alt > worst_altitude) {
                    worst_altitude = d_alt;
                    at_latitude = latitude;
                    at_altitude = altitude;
                }
                worst_latitude = std::max(worst_latitude, d_lat);
                worst_longitude = std::max(worst_longitude, d_lon);
            }
        }
    }
    std::printf(
        "[ GDZ round-trip ] worst |dalt| = %.3g km (at lat %g, alt %g km), |dlat| = %.3g deg, "
        "|dlon| = %.3g deg\n",
        worst_altitude, at_latitude, at_altitude, worst_latitude, worst_longitude);

    // Measured: 2.91e-11 km (29 nm) of altitude at 60 000 km, and 2.13e-14 deg in both angles
    // (2.4 nm on the ground). The caps are one decade above each. For scale, the footpoint budget
    // in docs/ERROR_BUDGET.md is 1e-3 deg, eleven orders of magnitude looser.
    EXPECT_LT(worst_altitude, 1e-10);
    EXPECT_LT(worst_latitude, 1e-13);
    EXPECT_LT(worst_longitude, 1e-12);
}

TEST(RoundTrip, CartesianToGeodeticAndBackOverTheWholeRange) {
    // The other direction, seeded in Cartesian so the sweep is not confined to the image of
    // gdz_to_geo. Radii span inside the Earth to well past geosynchronous.
    double worst = 0.0;
    double at_radius = 0.0;
    for (const double radius : {0.5, 0.98, 1.0, 1.5, 3.0, 6.6, 10.0}) {
        for (int i = 0; i <= 180; ++i) {
            const double latitude = -90.0 + i;
            for (int j = 0; j <= 72; ++j) {
                const double longitude = -180.0 + (5.0 * j);
                const auto start = ib::sph_to_car(sph(radius, latitude, longitude));
                const auto back = ib::gdz_to_geo(ib::geo_to_gdz(start));
                const double error = fx::norm(back.v - start.v);
                if (error > worst) {
                    worst = error;
                    at_radius = radius;
                }
            }
        }
    }
    std::printf("[ GEO round-trip ] worst |dx| = %.3g Re (%.3g mm) at r = %g Re\n", worst,
                worst * det::earth_radius_km * 1e6, at_radius);
    // Measured 3.97e-15 Re (25 nm) at the largest radius swept; cap one decade above.
    EXPECT_LT(worst, 1e-14);
}

TEST(RoundTrip, RllInvertsTheGeodeticRadiusExactly) {
    // rll_to_gdz solves for the altitude that reaches a given geocentric radius at a given geodetic
    // latitude. Compose it with the forward transform — take a GDZ point, measure its radius, ask
    // rll_to_gdz for the altitude — and the altitude must come back.
    double worst = 0.0;
    double at_latitude = 0.0;
    double at_altitude = 0.0;
    for (int i = 0; i <= 360; ++i) {
        const double latitude = -90.0 + (0.5 * i);
        for (const double altitude : sweep_altitudes_km) {
            const auto cartesian = ib::gdz_to_geo(gdz(altitude, latitude, 25.0));
            const double radius = fx::norm(cartesian.v);
            const double recovered = ib::rll_to_gdz(rll(radius, latitude, 25.0)).radius();
            if (std::abs(recovered - altitude) > worst) {
                worst = std::abs(recovered - altitude);
                at_latitude = latitude;
                at_altitude = altitude;
            }
        }
    }
    std::printf("[ RLL round-trip ] worst |dalt| = %.3g km at lat %g, alt %g km\n", worst,
                at_latitude, at_altitude);
    // Measured 2.2e-11 km (22 nm); cap one decade above.
    EXPECT_LT(worst, 1e-10);
}

TEST(RoundTrip, SphericalAndCartesianAreMutualInverses) {
    double worst_position = 0.0;
    double worst_radius = 0.0;
    for (const double radius : {0.25, 1.0, 4.5, 60.0}) {
        for (int i = 0; i <= 180; ++i) {
            const double latitude = -90.0 + i;
            for (int j = 0; j <= 72; ++j) {
                const double longitude = -180.0 + (5.0 * j);
                const auto cartesian = ib::sph_to_car(sph(radius, latitude, longitude));
                const auto back = ib::car_to_sph(cartesian);
                worst_radius = std::max(worst_radius, std::abs(back.radius() - radius));
                worst_position =
                    std::max(worst_position, fx::norm(ib::sph_to_car(back).v - cartesian.v));
            }
        }
    }
    std::printf("[ SPH round-trip ] worst |dr| = %.3g Re, worst |dx| = %.3g Re\n", worst_radius,
                worst_position);
    // Pure trigonometry, no ellipsoid: measured 1.42e-14 Re of radius and 3.5e-14 Re of position,
    // both at r = 60 where an ulp is largest. Caps one and a half decades above.
    EXPECT_LT(worst_radius, 1e-12);
    EXPECT_LT(worst_position, 1e-12);
}

TEST(GeodeticLatitude, DiffersFromGeocentricByUpToElevenArcminutes) {
    // The reason this file exists. If the two latitudes were interchangeable there would be no
    // ellipsoid code at all — so measure the gap and pin its size and its shape: zero at the
    // equator and the poles, maximal near 45°.
    double worst = 0.0;
    double at_latitude = 0.0;
    for (int i = 0; i <= 180; ++i) {
        const double geodetic = -90.0 + i;
        const auto cartesian = ib::gdz_to_geo(gdz(0.0, geodetic, 0.0));
        const double geocentric = ib::car_to_sph(cartesian).latitude();
        if (std::abs(geodetic - geocentric) > worst) {
            worst = std::abs(geodetic - geocentric);
            at_latitude = geodetic;
        }
    }
    std::printf("[ geodetic gap ] max |geodetic - geocentric| = %.6g deg (%.4g arcmin) at %g deg\n",
                worst, worst * 60.0, at_latitude);
    EXPECT_NEAR(worst, 0.19242, 1e-4);         // ~11.5 arcminutes
    EXPECT_EQ(std::abs(at_latitude), 45.0);    // and it peaks at 45°, as the theory says
    EXPECT_EQ(ib::car_to_sph(ib::gdz_to_geo(gdz(0.0, 0.0, 0.0))).latitude(), 0.0);
}

// ---- the heap tripwire --------------------------------------------------------------------------

TEST(Allocation, NoConversionTouchesTheHeap) {
    // Every routine in this header is a fixed-size computation over a vec3d, so the process-wide
    // counter must not move across a run of them. Per alloc_counter.hpp the assertion with teeth is
    // the one around the SECOND run: a routine that allocated a workspace once and reused it would
    // slip past a check that only wrapped the first.
    volatile double sink = 0.0;
    const auto exercise = [&sink]() {
        for (int i = 0; i < 1000; ++i) {
            const double t = 0.125 * i;
            const auto a = ib::gdz_to_geo(gdz(t, 0.05 * t, -180.0 + (0.36 * i)));
            const auto b = ib::geo_to_gdz(a);
            const auto c = ib::sph_to_car(sph(1.0 + t, 0.05 * t, 0.36 * i));
            const auto d = ib::car_to_sph(c);
            const auto e = ib::rll_to_gdz(rll(d.radius(), b.latitude(), b.longitude()));
            sink = sink + b.radius() + d.radius() + e.radius();
        }
    };

    exercise();
    const std::size_t before = cheatah_space_test::allocation_count();
    exercise();
    EXPECT_EQ(cheatah_space_test::allocation_count(), before);
    EXPECT_TRUE(std::isfinite(static_cast<double>(sink)));
}
