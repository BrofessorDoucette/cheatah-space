#pragma once

/**
 * @file coords_geodetic.hpp
 * @brief space.irbem — the conversions where the Earth stops being a sphere.
 *
 * Every other frame transform in this library is a rotation: an orthogonal 3×3 matrix, exact,
 * invertible, and cheap. This file holds the ones that are not. Geodetic coordinates (GDZ) are
 * referred to the WGS84 *ellipsoid*, so "latitude" means the angle of the local surface normal
 * rather than the angle of the position vector, and "altitude" means distance along that normal
 * rather than distance from the centre.
 *
 * The difference is not a rounding detail. The geodetic and geocentric latitudes of the same point
 * differ by up to 11.5 arcminutes near 45°, and the ellipsoid's own radius varies by 21.4 km from
 * equator to pole. A field-line footpoint reported with the two confused is wrong by roughly a
 * tenth of a degree — a hundred times the footpoint budget in `docs/ERROR_BUDGET.md`. So GDZ, SPH
 * and RLL are three *different* angular frames in @ref Frame, and this file is the only place that
 * knows how to move between them and Cartesian GEO.
 *
 * Three conversions here are exact closed forms. The fourth, Cartesian → geodetic, has no closed
 * form in elementary functions at all: it is the classic root of a quartic. This file uses
 * **Bowring's 1976 method**, iterated on the parametric latitude, which is a contraction whose
 * fixed point is the *exact* solution — see @ref detail::geodetic_latitude for why iterating it
 * (rather than taking Bowring's celebrated single step) is what makes it valid out to
 * geosynchronous altitude and beyond.
 *
 * @note Units are the frames' own, and they are not the same. GEO and SPH carry Earth radii; GDZ
 *       carries kilometres of altitude with degrees of latitude/longitude; RLL carries Earth radii
 *       of *geocentric* radius with *geodetic* latitude. That last combination is the one that
 *       surprises people — RLL is not SPH with a different name, and the two are not
 *       interchangeable (IRBEM's own `sysaxes` table says so explicitly).
 *
 * @note Where these routines diverge from IRBEM, the divergence is documented on the function that
 *       has it. Measured against IRBEM run as a black box over 132 830 points (latitude −90…90 in
 *       0.5° steps × longitude −180…180 in 10° steps × ten altitudes from −100 km to 60000 km):
 *       @ref gdz_to_geo agrees to **0.28 mm**, @ref sph_to_car to 6.7e-16 Re, and @ref rll_to_gdz
 *       to 0.25 mm — but @ref geo_to_gdz differs by up to **166 mm** of altitude and 1.8e-6° of
 *       latitude, because IRBEM's inverse does not iterate Bowring to convergence (its own
 *       GDZ→GEO→GDZ round-trip misses by the same 166 mm, off the polar axis; on it, by 2282 km).
 *       The only true convention difference is the longitude of the polar axis, which is undefined
 *       — see @ref car_to_sph.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

#include "fixarray.hpp"
#include "frames.hpp"

namespace cheatah::space::irbem {

namespace detail {

/**
 * The WGS84 semi-major axis, in kilometres.
 *
 * IRBEM's documentation names the ellipsoid outright — `docs/source/api/general_information.rst`,
 * the `sysaxes` table, entry 0: "Both the altitude and latitude depend on the ellipsoid used.
 * IRBEM uses the WGS84 reference ellipsoid." The defining constants are NIMA TR8350.2, 3rd ed.,
 * §3.2 (equivalently WGS84 as adopted by the IERS Conventions).
 */
inline constexpr double wgs84_semi_major_km = 6378.137;

/// The WGS84 inverse flattening, `1/f` — the second defining constant, NIMA TR8350.2 §3.2.
inline constexpr double wgs84_inverse_flattening = 298.257223563;

/// The WGS84 flattening `f`.
inline constexpr double wgs84_flattening = 1.0 / wgs84_inverse_flattening;

/// The WGS84 semi-minor axis `b = a(1 - f)`, in kilometres — 6356.7523142451792 km.
inline constexpr double wgs84_semi_minor_km = wgs84_semi_major_km * (1.0 - wgs84_flattening);

/// The first eccentricity squared, `e² = f(2 - f)` — 6.6943799901413165e-3.
inline constexpr double wgs84_eccentricity_squared =
    wgs84_flattening * (2.0 - wgs84_flattening);

/// The second eccentricity squared, `e'² = e²/(1 - e²)` — Bowring's method needs this one.
inline constexpr double wgs84_second_eccentricity_squared =
    wgs84_eccentricity_squared / (1.0 - wgs84_eccentricity_squared);

/**
 * The Earth radius that defines the "Re" unit of the Cartesian frames, in kilometres.
 *
 * This is the reference radius of the IGRF spherical-harmonic expansion (IAGA Working Group V-MOD;
 * every geopack-derived library inherits it), and it is *not* any axis of the WGS84 ellipsoid —
 * which is exactly why it needs stating: 1 Re on the equator is 6.937 km *below* the ellipsoid, and
 * 1 Re over the pole is 14.448 km above it.
 *
 * IRBEM's documentation never prints the number, so it was confirmed against the shipped library
 * run as a black box: `geo_gdz(1, 0, 0)` returns an altitude of -6.9369999999985916 km, i.e.
 * `Re - a` with `a = 6378.137`, giving `Re = 6371.2` exactly.
 */
inline constexpr double earth_radius_km = 6371.2;

/// Radians per degree — the inputs of GDZ/SPH/RLL are degrees, the trigonometry is radians.
inline constexpr double rad_per_deg = std::numbers::pi / 180.0;

/// Degrees per radian.
inline constexpr double deg_per_rad = 180.0 / std::numbers::pi;

/**
 * How many Bowring passes @ref geodetic_latitude takes.
 *
 * Measured, not guessed. Over a 36 010-point meridian sweep (latitude −90…90 in 0.05° steps ×
 * altitude −100 km … 60000 km), counting points whose latitude differs by even one bit from a
 * 16-pass run:
 *
 * | passes | points differing | worst |
 * |---|---|---|
 * | 2 | 6959 | 2.2e-16 rad |
 * | 3 | 40 | 2.2e-16 rad |
 * | **4** | **0** | **0** |
 * | 5 | 39 | 2.2e-16 rad |
 *
 * The odd/even pattern is real and is the reason the answer is four rather than three: at ~0.1% of
 * points the iteration does not settle on a single `double` but enters a **period-2 limit cycle
 * between two adjacent representable values**, 1 ulp (1.3e-14 deg, 1.4 nm on the ground) apart.
 * Neither element of that cycle is more correct than the other — it is a rounding artefact well
 * below every budget in `docs/ERROR_BUDGET.md` — but an even pass count lands consistently on one
 * of them, which makes the routine's output bit-reproducible and its convergence exactly testable
 * rather than testable to within a tolerance. That property is worth one extra sine and cosine.
 *
 * The count is fixed rather than checked at runtime so the routine stays branch-free — the batch
 * and GPU lanes want it that way — and the convergence claim is asserted by the test suite instead.
 */
inline constexpr std::size_t bowring_iterations = 4;

/**
 * The ellipsoid's auxiliary quantity `W = sqrt(1 - e² sin²φ)`, from which both radii of curvature
 * follow: the prime-vertical radius is `a/W` and the meridian radius is `a(1-e²)/W³`.
 * @param sin_lat the sine of the geodetic latitude.
 * @return `W`, in `[b/a, 1]` — never zero, so dividing by it is always safe.
 * @complexity O(1) — one multiply-chain and one square root.
 * @alloc none.
 */
inline double ellipsoid_w(double sin_lat) {
    return std::sqrt(1.0 - (wgs84_eccentricity_squared * sin_lat * sin_lat));
}

/**
 * The numerator and denominator of Bowring's latitude relation — a direction in the meridian plane
 * whose two-argument arctangent is the geodetic latitude.
 *
 * It exists so the iteration in @ref geodetic_latitude can carry the direction instead of the
 * angle. Nothing is normalised: only the ratio matters, and every use either divides the pair by
 * its own length or hands it to `atan2`, both of which are indifferent to a positive common factor.
 */
struct MeridianDirection {
    /// `z + e'² b sin³β` — the component along the polar axis.
    double numerator;
    /// `p - e² a cos³β` — the component along the equatorial plane. Negative only for points inside
    /// the ellipsoid near the axis, which is why the relation is an `atan2` and not an `atan`.
    double denominator;
};

/**
 * The meridian-plane direction whose two-argument arctangent is Bowring's geodetic latitude.
 *
 * Bowring, B. R., "Transformation from spatial to geographical coordinates", *Survey Review*
 * **23**(181), 323–327, 1976, eq. (1) is `tan φ = (z + e'² b sin³β) / (p - e² a cos³β)`. Written as
 * an `atan2` of numerator and denominator rather than an `atan` of their ratio it costs nothing
 * extra and removes both degeneracies: on the polar axis the denominator vanishes, and at the
 * geocentre both do.
 *
 * The relation is returned as that *pair* rather than as the angle because @ref geodetic_latitude
 * only ever needs the angle's sine and cosine, and forming the angle in between would mean an
 * `atan2` immediately undone by a `sin` and a `cos` — three libm calls per pass to recover a
 * direction the pair already holds exactly.
 *
 * @param sin_beta the sine of the parametric (reduced) latitude.
 * @param cos_beta its cosine.
 * @param p_km the distance from the polar axis, kilometres.
 * @param z_km the distance above the equatorial plane, kilometres.
 * @return the numerator and denominator of Bowring's ratio, unnormalised.
 * @complexity O(1) — nine multiplies and two adds, no transcendental and no division.
 * @alloc none.
 * @test IrbemPerf.BowringDirectionIsTheRelationBowringLatitudeTakesTheArctangentOf
 */
inline MeridianDirection bowring_direction(double sin_beta, double cos_beta, double p_km,
                                           double z_km) {
    return MeridianDirection{
        z_km + (wgs84_second_eccentricity_squared * wgs84_semi_minor_km * sin_beta * sin_beta *
                sin_beta),
        p_km - (wgs84_eccentricity_squared * wgs84_semi_major_km * cos_beta * cos_beta * cos_beta)};
}

/**
 * The sine and cosine of a direction, or `(0, 0)` when it has none.
 *
 * The one guard in the iteration, and it is for exactly one input: the geocentre, where the
 * direction is `(0, 0)` and has no angle. Every other point reaches this with a strictly positive
 * length, so the branch predicts perfectly and the compiler emits it as a select.
 *
 * @param y the component whose ratio to @p x is the tangent of the angle.
 * @param x the other component.
 * @return `(y, x)` divided by its length; `(0, 0)` when that length is zero.
 * @complexity O(1) — one square root and one division.
 * @alloc none.
 * @test IrbemPerf.NormalisingTheZeroDirectionYieldsZeroRatherThanNaN
 */
inline MeridianDirection normalised(double y, double x) {
    const double square = (y * y) + (x * x);
    const double inverse_length = square > 0.0 ? 1.0 / std::sqrt(square) : 0.0;
    return MeridianDirection{y * inverse_length, x * inverse_length};
}

/**
 * Bowring's relation evaluated at the parametric latitude a direction points along.
 *
 * The composition @ref normalised then @ref bowring_direction, which is the body of one pass of
 * @ref geodetic_direction — named so the pass reads as the single step it is.
 *
 * @param beta_y the numerator of `tan β`. @param beta_x its denominator.
 * @param p_km the distance from the polar axis, kilometres.
 * @param z_km the distance above the equatorial plane, kilometres.
 * @return the resulting geodetic direction, unnormalised.
 * @complexity O(1) — one square root, one division, and @ref bowring_direction's multiplies.
 * @alloc none.
 * @test IrbemPerf.TheDirectionCarryingBowringIterationMatchesTheAngleCarryingOne
 */
inline MeridianDirection bowring_direction_of(double beta_y, double beta_x, double p_km,
                                              double z_km) {
    const MeridianDirection beta = normalised(beta_y, beta_x);
    return bowring_direction(beta.numerator, beta.denominator, p_km, z_km);
}

/**
 * Bowring's exact relation between the parametric (reduced) latitude and the geodetic latitude,
 * as an angle.
 *
 * The identity holds *exactly* when @p beta is the true parametric latitude of the point; Bowring's
 * famous result is that seeding β with the geocentric direction already gives sub-millimetre φ near
 * the surface. See @ref bowring_direction for the relation itself; this is its arctangent, and is
 * the spelling the definition is stated in.
 *
 * @param beta the parametric latitude, radians.
 * @param p_km the distance from the polar axis, kilometres.
 * @param z_km the distance above the equatorial plane, kilometres.
 * @return the geodetic latitude, radians, in `[-π, π]`.
 * @complexity O(1) — one sine, one cosine, one `atan2`.
 * @alloc none.
 * @test BowringLatitude.SeededFromTheTrueParametricLatitudeIsExactInOnePass
 */
inline double bowring_latitude(double beta, double p_km, double z_km) {
    const MeridianDirection d = bowring_direction(std::sin(beta), std::cos(beta), p_km, z_km);
    return std::atan2(d.numerator, d.denominator);
}

/**
 * The converged Bowring direction of a meridian-plane point — the whole iteration, as a direction.
 *
 * The seed is Bowring's own: `tan β₀ = (z/p)(a/b)`, the parametric latitude of the geocentric
 * direction. Each pass then evaluates the exact relation @ref bowring_direction and maps the result
 * back through `tan β = (1-f) tan φ`, which is the *definition* of the parametric latitude and is
 * therefore also exact. Both halves being exact identities is the whole point: the composition is a
 * contraction whose fixed point is the true geodetic latitude, so iterating converges on the answer
 * instead of on Bowring's first-order approximation to it. That is what keeps the method valid at
 * 60000 km, where the single-step form is off by ~4.5e-7 deg.
 *
 * **The iteration carries a direction, not an angle**, which is what makes it cheap. β and φ enter
 * the recurrence only through their sines and cosines, so the classic spelling — `atan2` to get φ,
 * `sin`/`cos` to feed the parametric map, `atan2` again to get β, `sin`/`cos` again to feed
 * Bowring's relation — spends six libm calls a pass reconstructing a direction it already had.
 * Carrying `(y, x)` with `tan = y/x` instead, and normalising it once, replaces all six with a
 * single square root. Two further facts collapse it further: the parametric map `tan β = (1-f) tan φ`
 * applied to the direction `(num, den)` is just `((1-f)·num, den)`, and scaling BOTH components of
 * a direction by the same positive number cannot change it — so φ's own normalisation is redundant
 * and only β's survives. One `sqrt` and one divide per pass, and no `atan2` at all until a caller
 * asks for the angle.
 *
 * Measured over 65536 points at `-O3 -march=native`, against the angle-carrying form this replaced:
 * **@ref geodetic_latitude 275.0 ns → 74.1 ns per call, 3.71×**; @ref geo_to_gdz, which also spends
 * the direction rather than re-deriving `sin φ` and `cos φ`, **282.7 ns → 85.9 ns, 3.29×**. The two
 * forms agree to within 1 ulp — over a 21606-point meridian sweep they differ at 44 points, worst
 * 2.2e-16 rad, worst implied altitude difference 1.5e-11 km — and the four-pass convergence the
 * pass count rests on is unchanged: over that sweep four passes and sixteen give bit-identical
 * answers in both forms.
 *
 * @param p_km the distance from the polar axis, kilometres; `p ≥ 0`.
 * @param z_km the distance above the equatorial plane, kilometres.
 * @param iterations how many Bowring passes to run. One pass is always performed, so `0` and `1`
 *        behave identically; @ref bowring_iterations is the production value.
 * @return the converged direction, unnormalised. At the geocentre it is `(0, 0)`, which is the
 *         honest answer for a point that has no geodetic latitude at all.
 * @complexity O(@p iterations), each pass one square root and one division; no transcendental.
 * @alloc none.
 * @test IrbemPerf.TheDirectionCarryingBowringIterationMatchesTheAngleCarryingOne
 */
inline MeridianDirection geodetic_direction(double p_km, double z_km, std::size_t iterations) {
    // The seed, as a direction: tan β₀ = (z a)/(p b).
    double beta_y = z_km * wgs84_semi_major_km;
    double beta_x = p_km * wgs84_semi_minor_km;
    MeridianDirection phi{0.0, 0.0};
    for (std::size_t pass = 0;; ++pass) {
        phi = bowring_direction_of(beta_y, beta_x, p_km, z_km);
        if (pass + 1 >= iterations) break;
        // tan β = (1-f) tan φ, on the unnormalised direction: the length φ's direction would have
        // been divided by is a positive common factor of both components, so it cancels.
        beta_y = (1.0 - wgs84_flattening) * phi.numerator;
        beta_x = phi.denominator;
    }
    return phi;
}

/**
 * The geodetic latitude of a meridian-plane point, by Bowring's method iterated to its fixed point.
 *
 * The angle spelling of @ref geodetic_direction — see it for the method, the convergence evidence
 * and the measurements. The single `atan2` here is the only transcendental the whole iteration
 * costs.
 *
 * @param p_km the distance from the polar axis, kilometres; `p ≥ 0`.
 * @param z_km the distance above the equatorial plane, kilometres.
 * @param iterations how many Bowring passes to run; @ref bowring_iterations is the production value.
 * @return the geodetic latitude, radians, in `[-π, π]`. At the geocentre the direction is `(0, 0)`
 *         and its `atan2` is zero; the angle-carrying form this replaced landed on 1.8e-48 there,
 *         which is the same answer to every digit that means anything.
 * @complexity O(@p iterations) square roots plus one `atan2`.
 * @alloc none.
 * @test BowringIteration.HasReachedItsFixedPointEverywhereOnTheSweep
 */
inline double geodetic_latitude(double p_km, double z_km, std::size_t iterations) {
    const MeridianDirection phi = geodetic_direction(p_km, z_km, iterations);
    return std::atan2(phi.numerator, phi.denominator);
}

}  // namespace detail

/**
 * Geodetic → geocentric Cartesian: IRBEM's `GDZ2GEO` (`gdz_geo`).
 *
 * The exact closed form, e.g. Torge, *Geodesy* (3rd ed.) §4.1, or Hofmann-Wellenhof et al.,
 * *GPS: Theory and Practice*, eq. (10.1):
 * `x = (N+h)cosφ cosλ`, `y = (N+h)cosφ sinλ`, `z = (N(1-e²)+h) sinφ`, with the prime-vertical
 * radius `N = a/W`. The `(1-e²)` on `z` is the ellipsoid: the surface normal at latitude φ does not
 * pass through the centre, so the vertical drops short of the equatorial plane by exactly that
 * factor.
 *
 * @param gdz the geodetic position — altitude in km above the WGS84 ellipsoid, geodetic latitude in
 *        degrees (north positive), east longitude in degrees.
 * @return the geographic Cartesian position, in Earth radii.
 * @complexity O(1) — four transcendentals and a square root.
 * @alloc none.
 */
inline Position<Frame::GEO> gdz_to_geo(const Position<Frame::GDZ>& gdz) {
    const double lat = gdz.latitude() * detail::rad_per_deg;
    const double lon = gdz.longitude() * detail::rad_per_deg;
    const double sin_lat = std::sin(lat);
    const double cos_lat = std::cos(lat);
    const double prime_vertical = detail::wgs84_semi_major_km / detail::ellipsoid_w(sin_lat);

    // Work in km, divide once at the end: GEO is in Earth radii, the ellipsoid is in km.
    const double equatorial = (prime_vertical + gdz.radius()) * cos_lat;
    const double polar =
        ((prime_vertical * (1.0 - detail::wgs84_eccentricity_squared)) + gdz.radius()) * sin_lat;
    return Position<Frame::GEO>{fixarray::vec3d{equatorial * std::cos(lon) / detail::earth_radius_km,
                                                equatorial * std::sin(lon) / detail::earth_radius_km,
                                                polar / detail::earth_radius_km}};
}

/**
 * Geocentric Cartesian → geodetic: IRBEM's `GEO2GDZ` (`geo_gdz`).
 *
 * Latitude comes from @ref detail::geodetic_latitude (Bowring 1976, iterated). Altitude then uses
 * the projection form `h = p cosφ + z sinφ - aW` rather than Bowring's `h = p/cosφ - N`. The two
 * are algebraically identical — substituting the forward transform gives
 * `p cosφ + z sinφ = N(1 - e² sin²φ) + h = aW + h` — but this one has no cosine in a denominator,
 * so it is exact on the polar axis and loses no digits approaching it.
 *
 * @param geo the geographic Cartesian position, in Earth radii.
 * @return the geodetic position — altitude km, geodetic latitude deg, east longitude deg in
 *         `(-180, 180]`.
 * @complexity O(1) — @ref detail::bowring_iterations Bowring passes plus a handful of
 *             transcendentals.
 * @alloc none.
 *
 * @note **Diverges from IRBEM exactly on the polar axis, where IRBEM is wrong.** For `(0, 0, 1)`
 *       IRBEM returns altitude 0; this returns `Re - b = 14.447685754820668` km, which is what its
 *       own `gdz_geo(90, 0, 0)` inverts to and what IRBEM's own answer tends to as `p → 0`
 *       (`geo_gdz(1e-6, 0, 1)` gives 14.447690). A GDZ→GEO→GDZ sweep through IRBEM is off by up to
 *       2282 km at latitude −90; through this routine the same sweep is exact to 2.9e-11 km.
 * @note At the geocentre the answer is degenerate, because a point with no direction has no
 *       geodetic latitude: the iteration alternates between 0 and π with the parity of the pass
 *       count, so with @ref detail::bowring_iterations even the reported latitude is ~0 where IRBEM
 *       reports 180°. The **altitude is `-a` exactly either way**, which is the only part of the
 *       answer that means anything. Callers that care must reject `r = 0` themselves.
 */
inline Position<Frame::GDZ> geo_to_gdz(const Position<Frame::GEO>& geo) {
    const double p_km = std::hypot(geo.v[0], geo.v[1]) * detail::earth_radius_km;
    const double z_km = geo.v[2] * detail::earth_radius_km;

    // The converged direction already IS (sin phi, cos phi) up to its own length, so normalising it
    // costs one square root where std::sin and std::cos of the angle cost two libm calls — and the
    // angle itself is wanted only for the answer, not for the arithmetic. Measured 97.7 ns -> 85.9
    // ns per call; the two spellings differ by at most 1.4e-14 deg in latitude and 3.6e-11 km in
    // altitude over a 21606-point meridian sweep, four orders below `docs/ERROR_BUDGET.md`'s
    // coordinate-transform line.
    const detail::MeridianDirection phi =
        detail::geodetic_direction(p_km, z_km, detail::bowring_iterations);
    const detail::MeridianDirection unit = detail::normalised(phi.numerator, phi.denominator);
    const double sin_lat = unit.numerator;
    const double cos_lat = unit.denominator;
    const double altitude_km = (p_km * cos_lat) + (z_km * sin_lat) -
                               (detail::wgs84_semi_major_km * detail::ellipsoid_w(sin_lat));

    return Position<Frame::GDZ>{fixarray::vec3d{
        altitude_km, std::atan2(phi.numerator, phi.denominator) * detail::deg_per_rad,
        std::atan2(geo.v[1], geo.v[0]) * detail::deg_per_rad}};
}

/**
 * Geographic spherical → geocentric Cartesian: IRBEM's `SPH2CAR` (`sph_car`).
 *
 * The textbook spherical-polar relation with latitude in place of colatitude:
 * `x = r cosφ cosλ`, `y = r cosφ sinλ`, `z = r sinφ`. No ellipsoid is involved — SPH is GEO
 * re-expressed, nothing more — so this is exact and its inverse @ref car_to_sph is a true inverse.
 *
 * @param sph the spherical position — radius in Earth radii, geocentric latitude in degrees, east
 *        longitude in degrees.
 * @return the geographic Cartesian position, in Earth radii.
 * @complexity O(1) — four transcendentals.
 * @alloc none.
 */
inline Position<Frame::GEO> sph_to_car(const Position<Frame::SPH>& sph) {
    const double lat = sph.latitude() * detail::rad_per_deg;
    const double lon = sph.longitude() * detail::rad_per_deg;
    const double equatorial = sph.radius() * std::cos(lat);
    return Position<Frame::GEO>{fixarray::vec3d{equatorial * std::cos(lon),
                                                equatorial * std::sin(lon),
                                                sph.radius() * std::sin(lat)}};
}

/**
 * Geocentric Cartesian → geographic spherical: IRBEM's `CAR2SPH` (`car_sph`).
 *
 * Latitude is taken as `90° - θ` with the colatitude `θ = atan2(sqrt(x²+y²), z)`, and longitude as
 * `atan2(y, x)`. Both are `atan2` rather than `acos`/`atan` for the usual reason — full range, no
 * cancellation as the argument approaches ±1, and no division by a vanishing radius.
 *
 * @param car the geographic Cartesian position, in Earth radii.
 * @return the spherical position — radius in Earth radii, geocentric latitude in `[-90, 90]`, east
 *         longitude in `(-180, 180]`.
 * @complexity O(1) — a hypot and two `atan2`.
 * @alloc none.
 *
 * @note At the origin the colatitude `atan2(0, 0)` is 0, so the result is `(0, 90°, 0°)`. That is
 *       the limit along the `+z` axis and is what IRBEM reports there; a point with no direction
 *       has no honest latitude, and callers that care must reject `r = 0` themselves.
 * @note **Longitude on the polar axis is undefined, and here it is IEEE-defined instead.** With
 *       `x` and `y` both zero, `atan2` reports the branch the *signs of those zeros* select — so
 *       `(-0, -0, z)` gives −180° where IRBEM gives 0°. That is the one genuine convention
 *       difference in this file. It is left alone rather than special-cased: `atan2`'s answer makes
 *       @ref sph_to_car → @ref car_to_sph recover the original longitude for every input this
 *       library produces (its polar `x`, `y` are ~1e-17, not zero), which a hard-coded 0° would
 *       not.
 */
inline Position<Frame::SPH> car_to_sph(const Position<Frame::GEO>& car) {
    const double colatitude = std::atan2(std::hypot(car.v[0], car.v[1]), car.v[2]);
    return Position<Frame::SPH>{
        fixarray::vec3d{fixarray::norm(car.v), 90.0 - (colatitude * detail::deg_per_rad),
                        std::atan2(car.v[1], car.v[0]) * detail::deg_per_rad}};
}

/**
 * Radius/latitude/longitude → geodetic: IRBEM's `RLL2GDZ` (`rll_gdz`).
 *
 * RLL is the awkward frame: a *geocentric* radius carried alongside a *geodetic* latitude (IRBEM's
 * `sysaxes` table calls this out — "the latitude is still geodetic latitude and is therefore not
 * interchangeable with SPH"). So the job is not a spherical-to-geodetic conversion; it is to solve
 * for the altitude `h` that puts the point at geocentric radius `r` along the ellipsoid normal at a
 * *given* geodetic latitude. Squaring the forward transform makes that a quadratic in `h`:
 *
 *     r² = (N + h)²cos²φ + (N(1-e²) + h)²sin²φ  =  h² + 2Bh + C,
 *     B = N(1 - e² sin²φ) = aW,      C = N²(cos²φ + (1-e²)² sin²φ),
 *     h = -B + sqrt(B² - C + r²).
 *
 * and the discriminant simplifies exactly, since `B² - C = -N² e⁴ sin²φ cos²φ`:
 *
 *     h = -B + sqrt(r² - N² e⁴ sin²φ cos²φ).
 *
 * The positive root is the one taken: the other places the point on the far side of the axis.
 *
 * @param rll the RLL position — geocentric radius in Earth radii, geodetic latitude in degrees,
 *        east longitude in degrees.
 * @return the geodetic position, with latitude and longitude passed through unchanged (RLL and GDZ
 *         share both by definition) and the solved altitude in km.
 * @complexity O(1) — two transcendentals and two square roots.
 * @alloc none.
 *
 * @note The subtracted term peaks at ~21.4 km, so the discriminant goes negative for points deeper
 *       than about 0.0034 Re from the centre near 45° latitude — there is no ellipsoid normal
 *       through such a point at that latitude. It is clamped to zero, giving `h = -aW`, the
 *       altitude of the foot of the normal. That is a well-defined answer for a degenerate input,
 *       not a silent failure: no radius a caller can physically observe reaches it.
 * @note This is exact where IRBEM's is not. IRBEM returns 14.447685999998612 km for
 *       `rll_gdz(1, 90, 0)`; the closed form gives `Re - b = 14.447685754820668` km. The 2.4e-7 km
 *       gap is IRBEM's own iteration, and it is present at every latitude.
 */
inline Position<Frame::GDZ> rll_to_gdz(const Position<Frame::RLL>& rll) {
    const double lat = rll.latitude() * detail::rad_per_deg;
    const double sin_lat = std::sin(lat);
    const double cos_lat = std::cos(lat);
    const double w = detail::ellipsoid_w(sin_lat);
    const double prime_vertical = detail::wgs84_semi_major_km / w;
    const double radius_km = rll.radius() * detail::earth_radius_km;

    const double offset = prime_vertical * detail::wgs84_eccentricity_squared * sin_lat * cos_lat;
    const double discriminant = std::max(0.0, (radius_km * radius_km) - (offset * offset));
    const double altitude_km =
        std::sqrt(discriminant) - (detail::wgs84_semi_major_km * w);

    return Position<Frame::GDZ>{
        fixarray::vec3d{altitude_km, rll.latitude(), rll.longitude()}};
}

}  // namespace cheatah::space::irbem
