#pragma once

/**
 * @file coords_helio.hpp
 * @brief space.irbem — the heliospheric frames HEE, HAE and HEEQ, and the bridge from GSE.
 *
 * Everything upstream of the magnetosphere is published heliocentrically. Solar-wind monitors (ACE,
 * Wind, DSCOVR), coronal imagers and the interplanetary-field indices that drive every external
 * field model arrive in HEE, HAE or HEEQ; the magnetosphere itself is worked entirely in geocentric
 * frames. Propagating an upstream measurement to the magnetopause therefore crosses this boundary
 * every single time, which makes these six transforms the seam between two whole literatures.
 *
 * It is also the ONE seam in this library where the ORIGIN MOVES. HEE, HAE and HEEQ are centred on
 * the Sun; GSE is centred on the Earth. So GSE↔HEE is a rotation *plus a translation by the
 * Sun–Earth vector* for a position, and a rotation *alone* for a field vector — a field is a
 * direction and a magnitude at a point, and translating it is meaningless. Both objects are three
 * doubles, both come out of the wrong routine looking entirely plausible, and the error is ~23 000
 * Earth radii. That is the classic defect of this file, so the distinction is structural rather
 * than documentary:
 *
 *   - The three purely rotational transforms (HEE↔HAE, HAE↔HEEQ — both frames heliocentric, so no
 *     translation exists to get wrong) are written ONCE over both vector kinds, through
 *     @ref FrameTagged. There is no per-kind body in which a translation could be introduced.
 *   - GSE↔HEE is written out FOUR times — position and field, each direction — so the presence and
 *     the absence of the `+ r₀` sit next to each other on the page and are reviewed as a pair.
 *
 * The epoch-dependent geometry — the Earth's heliocentric ecliptic longitude λ_geo, the Sun–Earth
 * distance r₀, the ascending node Ω of the solar equator and the central-meridian angle θ — is
 * computed ONCE by @ref helio_geometry into a @ref HelioGeometry, which carries the finished
 * rotation matrices. Every transform is then pure arithmetic: no trigonometry, no branch and no
 * allocation per point, which is what an ephemeris of 10⁵ points needs.
 *
 * ## Sources
 *
 * The derivation is Hapgood, M. A. (1992), *Space physics coordinate transformations: a user
 * guide*, Planet. Space Sci. 40, 711–717 (with the 1995 corrigendum, Planet. Space Sci. 43, 1391).
 * The formulation implemented here, including the Eulerian-matrix spelling of each system, is
 * Fränz, M. & Harper, D. (2002), *Heliospheric coordinate systems*, Planet. Space Sci. 50, 217–233
 * (corrected version of 2002-03-12), hereafter **F&H**:
 *
 *   - §3.2.2 "Heliocentric Earth Ecliptic HEE": `T(HAE_D→HEE_D) = E(0, 0, λ_geo)`.
 *   - §3.2.2 "Heliocentric Earth Equatorial HEEQ": `T(HAE_D→HEEQ) = E(Ω, i, θ)`.
 *   - §3.3: `T(HAE_D→GSE) = E(0, 0, λ_geo + 180°)`, whence GSE↔HEE is exactly a half turn about Z.
 *   - eqn. 14: the solar rotation axis, `i = 7.25°`, `Ω = 75.76° + 1.397°·T₀`.
 *   - eqn. 17: `θ = arctan(cos i · tan(λ − Ω))`, the central-meridian angle.
 *   - eqn. 36: λ_geo and r₀ from the mean elements of the Earth–Moon barycentre (F&H Table 4, after
 *     Simon et al. 1994, A&A 282, 663), with the equation-of-centre approximation of the
 *     Astronomical Almanac page C24 — the same approximation Hapgood (1992) used.
 *
 * F&H's Eulerian rotation `E(Ω, θ, φ)` is `<φ, Z>·<θ, X>·<Ω, Z>`, where `<a, Z>` rotates the frame
 * (not the vector) by `a` about Z. Those are the matrices built below.
 *
 * @note **Units.** Positions in the heliospheric frames are in Earth radii, like every other
 *       Cartesian frame in @ref frames.hpp, and like IRBEM's own `sysaxes` table which lists HEE,
 *       HAE and HEEQ as "Re". IRBEM's *dedicated* `GSE2HEE`/`HEE2GSE` entry points document their
 *       heliospheric argument in AU instead — an inconsistency inside IRBEM's own documentation.
 *       Multiply by @ref au_in_earth_radii to reach IRBEM's AU convention at that boundary.
 *
 * @note **Time scale.** @ref helio_geometry takes a Modified Julian Date on the TT scale, which is
 *       what F&H's own worked example uses. Feeding it UTC instead shifts λ_geo by
 *       (TT−UTC)·1.14e-5 °/s ≈ 8e-4 ° in 2026 — an order of magnitude below the 0.01° accuracy of
 *       eqn. 36 itself, so it is a defensible approximation, but it is an approximation.
 *
 * @note **Accuracy.** eqn. 36 is a first-order Keplerian approximation: F&H quote 34 arcsec (0.01°)
 *       for 1950–2050, which at 1 AU is ~3.5 Re of transverse position error in HEE. That, not
 *       arithmetic, is the accuracy of these transforms; see `docs/ERROR_BUDGET.md`.
 */

#include <cmath>
#include <concepts>

#include "cheatah.hpp"
#include "fixarray.hpp"
#include "space/irbem/frames.hpp"

namespace cheatah::space::irbem {

namespace detail::helio {

/// Degrees per radian — the one conversion this file needs, in both directions.
inline constexpr double kDegreesPerRadian = 57.295779513082320876798154814105;

/// The Modified Julian Date of the J2000.0 epoch, 2000-01-01 12:00 TT (JD 2451545.0).
inline constexpr double kMjdJ2000 = 51544.5;

/// Days in the Julian century that F&H's time argument `T₀` is scaled in.
inline constexpr double kDaysPerJulianCentury = 36525.0;

/// Mean longitude of the Earth–Moon barycentre at J2000.0, degrees (F&H Table 4, Simon et al. 1994).
inline constexpr double kEmbMeanLongitudeDeg = 100.4664568;
/// Rate of the EMB mean longitude, degrees per Julian century (F&H Table 4).
inline constexpr double kEmbMeanLongitudeRateDeg = 35999.3728565;
/// Longitude of periapsis of the EMB orbit at J2000.0, degrees (F&H Table 4).
inline constexpr double kEmbPeriapsisDeg = 102.9373481;
/// Rate of the EMB longitude of periapsis, degrees per Julian century (F&H Table 4).
inline constexpr double kEmbPeriapsisRateDeg = 0.3225654;

/// First equation-of-centre coefficient, degrees (F&H eqn. 36; Astronomical Almanac C24).
inline constexpr double kCentreFirstDeg = 1.915;
/// Second equation-of-centre coefficient, degrees (F&H eqn. 36).
inline constexpr double kCentreSecondDeg = 0.020;

/// Mean term of the Sun–Earth distance, AU (F&H eqn. 36).
inline constexpr double kDistanceMeanAu = 1.00014;
/// Coefficient of cos(g) in the Sun–Earth distance, AU (F&H eqn. 36; subtracted).
inline constexpr double kDistanceCosAu = 0.01671;
/// Coefficient of cos(2g) in the Sun–Earth distance, AU (F&H eqn. 36; subtracted).
inline constexpr double kDistanceCos2Au = 0.00014;

/// Longitude of the ascending node of the solar equator on the ecliptic at J2000.0, degrees
/// (F&H eqn. 14; the value traces to Carrington and is kept fixed by convention).
inline constexpr double kSolarNodeDeg = 75.76;
/// Rate of that node, degrees per Julian century (F&H eqn. 14). It absorbs ecliptic precession, so
/// no separate precession matrix is applied on the way to HEEQ.
inline constexpr double kSolarNodeRateDeg = 1.397;

/// Light-aberration correction applied to the Earth's longitude before the central-meridian angle,
/// degrees. F&H §3.2.1: the solar central meridian is the *apparent* sub-Earth point, whose
/// heliocentric longitude is `λ_geo + a` with `a ≈ −20″` (F&H Appendix A.3). Reproducing F&H's own
/// Table 8 requires it: dropping it moves θ by 20″ and the HEEQ position by ~7e-4 Re.
inline constexpr double kAberrationDeg = -20.0 / 3600.0;

/// One astronomical unit in kilometres, IAU 2012 Resolution B2 (a defined constant).
inline constexpr double kAstronomicalUnitKm = 149597870.7;

/// The Earth radius that "Re" means throughout this module, in kilometres: the IGRF/IRBEM
/// geomagnetic reference radius, so that a heliospheric position and a `GEO` position are in the
/// same unit. (F&H's own Table 8 instead uses the 6378.14 km equatorial radius; that choice affects
/// only how a length is *labelled*, never a rotation, and the tests exercise the rotations in F&H's
/// units and the translation in this one.)
inline constexpr double kEarthRadiusKm = 6371.2;

/**
 * Degrees to radians.
 * @param degrees_in the angle in degrees.
 * @return the same angle in radians.
 * @complexity O(1).
 * @alloc none.
 */
constexpr double radians(double degrees_in) { return degrees_in / kDegreesPerRadian; }

/**
 * Radians to degrees.
 * @param radians_in the angle in radians.
 * @return the same angle in degrees.
 * @complexity O(1).
 * @alloc none.
 */
constexpr double degrees(double radians_in) { return radians_in * kDegreesPerRadian; }

/**
 * An angle folded into `[0, 360)`.
 *
 * The mean longitudes grow by 36 000° per century, so by 2100 they are ~7 revolutions from zero.
 * Folding before the trigonometry keeps the argument small — libm's own reduction would be accurate
 * anyway, but the STORED angle is what a caller reads and compares against a published value.
 *
 * @param degrees_in the angle in degrees, of any magnitude and sign.
 * @return the equivalent angle in `[0, 360)`.
 * @complexity O(1).
 * @alloc none.
 */
inline double wrap_360(double degrees_in) {
    const double folded = std::fmod(degrees_in, 360.0);
    return folded < 0.0 ? folded + 360.0 : folded;
}

/**
 * The frame rotation `<a, Z>` — a right-handed rotation of the AXES by @p degrees_in about Z.
 *
 * Rotating the axes is the transpose of rotating the vector, which is why `+sin` sits in the top
 * row: a vector at longitude `a` comes out on the new +X axis.
 *
 * @param degrees_in the rotation angle in degrees.
 * @return the 3×3 rotation matrix, orthonormal with determinant +1.
 * @complexity O(1).
 * @alloc none.
 */
inline fixarray::mat3d rotation_z(double degrees_in) {
    const double angle = radians(degrees_in);
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return fixarray::mat3d{c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0};
}

/**
 * The frame rotation `<a, X>` — a right-handed rotation of the AXES by @p degrees_in about X.
 * @param degrees_in the rotation angle in degrees.
 * @return the 3×3 rotation matrix, orthonormal with determinant +1.
 * @complexity O(1).
 * @alloc none.
 */
inline fixarray::mat3d rotation_x(double degrees_in) {
    const double angle = radians(degrees_in);
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return fixarray::mat3d{1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c};
}

}  // namespace detail::helio

/**
 * One astronomical unit in kilometres.
 * @return 149 597 870.7 km, the defining value of IAU 2012 Resolution B2.
 * @complexity O(1).
 * @alloc none.
 */
constexpr double astronomical_unit_km() { return detail::helio::kAstronomicalUnitKm; }

/**
 * One astronomical unit in Earth radii — the factor between this module's positions and the AU that
 * the heliospheric literature (and IRBEM's own `GSE2HEE` entry point) prefers.
 * @return `astronomical_unit_km() / 6371.2`, ≈ 23 480.33.
 * @complexity O(1).
 * @alloc none.
 */
constexpr double au_in_earth_radii() {
    return detail::helio::kAstronomicalUnitKm / detail::helio::kEarthRadiusKm;
}

/**
 * The inclination of the solar equator to the ecliptic.
 * @return 7.25°, the conventional value of F&H eqn. 14. It is held fixed: later measurements make
 *         the axis direction less well determined, not better, and coordinate work sticks with
 *         Carrington's value so that datasets remain comparable.
 * @complexity O(1).
 * @alloc none.
 */
constexpr double solar_equator_inclination_deg() { return 7.25; }

/**
 * The Sun–Earth geometry at one epoch — everything the six heliospheric transforms need, evaluated
 * once so that no transform pays for trigonometry.
 *
 * A plain aggregate of doubles and two matrices: trivially copyable, nothing allocated, and cheap
 * enough to hold one per ephemeris sample or to rebuild per time step, whichever the caller prefers.
 */
struct HelioGeometry {
    /// The Earth's heliocentric ecliptic longitude λ_geo of date, degrees in `[0, 360)`
    /// (F&H eqn. 36). This is the angle that carries HAE into HEE.
    double earth_longitude_deg = 0.0;

    /// λ_geo corrected for light aberration, degrees in `[0, 360)` — the longitude of the
    /// *apparent* sub-Earth point on the solar disc, which is what defines the central meridian
    /// (F&H §3.2.1).
    double apparent_earth_longitude_deg = 0.0;

    /// The Sun–Earth distance r₀ in astronomical units (F&H eqn. 36).
    double sun_earth_distance_au = 0.0;

    /// The same distance in Earth radii — the translation between the geocentric and heliocentric
    /// origins, in the unit this module's positions use.
    double sun_earth_distance_re = 0.0;

    /// The ascending node Ω of the solar equator on the ecliptic of date, degrees in `[0, 360)`
    /// (F&H eqn. 14).
    double solar_node_deg = 0.0;

    /// The central-meridian angle θ: the heliographic longitude of the sub-Earth point measured
    /// from Ω, degrees in `[0, 360)` (F&H eqn. 17).
    double solar_central_meridian_deg = 0.0;

    /// `<λ_geo, Z>` — the rotation carrying HAE components into HEE (F&H §3.2.2).
    fixarray::mat3d hae_to_hee;

    /// `<θ, Z>·<i, X>·<Ω, Z>` — the rotation carrying HAE components into HEEQ (F&H §3.2.2).
    fixarray::mat3d hae_to_heeq;
};

/**
 * The heliospheric geometry at an epoch — the only place in this file that evaluates a
 * transcendental, and the only place that reads the published ephemeris constants.
 *
 * λ_geo and r₀ come from the Earth–Moon barycentre's mean elements (F&H Table 4) fed through the
 * equation-of-centre approximation of F&H eqn. 36; Ω from eqn. 14; and θ from eqn. 17, evaluated as
 * a two-argument arctangent of `(cos i · sin(λ−Ω), cos(λ−Ω))` so that the quadrant follows from the
 * geometry instead of from a sign convention. F&H's prose ("the quadrant of θ is opposite that of
 * λ−Ω") does not describe this construction; the two-argument form is what reproduces F&H's own
 * worked example (their Table 8 implies θ = 259.899186°, against the 259.89919° they print), so the
 * prose is read as a slip and the arithmetic is trusted.
 *
 * @param mjd_tt the epoch as a Modified Julian Date on the TT scale (`JD_TT − 2400000.5`).
 * @return the geometry, with every angle folded into `[0, 360)` and both rotation matrices built.
 * @complexity O(1) — a fixed count of transcendental evaluations and two 3×3 matrix products,
 *             independent of the epoch.
 * @alloc none; the result is 6 doubles plus 18 more in two matrices, all inline.
 */
inline HelioGeometry helio_geometry(double mjd_tt) {
    namespace h = detail::helio;

    // F&H's time argument: Julian centuries of 36525 days from J2000.0, TT.
    const double t0 = (mjd_tt - h::kMjdJ2000) / h::kDaysPerJulianCentury;

    // F&H Table 4: the EMB's mean longitude and longitude of periapsis, both linear in T0. Their
    // difference is the mean anomaly g. Standing in for the Earth by the barycentre costs ~14
    // arcsec (F&H Table 5) — an order below eqn. 36's own error, and the same choice F&H make.
    const double mean_longitude = h::kEmbMeanLongitudeDeg + h::kEmbMeanLongitudeRateDeg * t0;
    const double periapsis = h::kEmbPeriapsisDeg + h::kEmbPeriapsisRateDeg * t0;
    const double anomaly = h::radians(h::wrap_360(mean_longitude - periapsis));
    const double sin_g = std::sin(anomaly);
    const double cos_g = std::cos(anomaly);
    // sin(2g) and cos(2g) by the double-angle identities rather than two more library calls. The
    // identities are exact in the reals and cost one rounding here, against libm's own ~0.5 ulp —
    // a wash numerically, and two fewer transcendentals per epoch.
    const double sin_2g = 2.0 * sin_g * cos_g;
    const double cos_2g = 1.0 - 2.0 * sin_g * sin_g;

    // F&H eqn. 36 — the first-order equation of centre, and the corresponding radius vector.
    const double longitude =
        h::wrap_360(mean_longitude + h::kCentreFirstDeg * sin_g + h::kCentreSecondDeg * sin_2g);
    const double distance_au =
        h::kDistanceMeanAu - h::kDistanceCosAu * cos_g - h::kDistanceCos2Au * cos_2g;

    // F&H eqn. 14 and eqn. 17: the solar rotation axis, and the Earth's longitude within the solar
    // equator measured from the ascending node.
    const double node = h::wrap_360(h::kSolarNodeDeg + h::kSolarNodeRateDeg * t0);
    const double apparent_longitude = h::wrap_360(longitude + h::kAberrationDeg);
    const double from_node = h::radians(apparent_longitude - node);
    const double inclination = h::radians(solar_equator_inclination_deg());
    const double central_meridian =
        h::wrap_360(h::degrees(std::atan2(std::cos(inclination) * std::sin(from_node),  //
                                          std::cos(from_node))));

    return HelioGeometry{
        .earth_longitude_deg = longitude,
        .apparent_earth_longitude_deg = apparent_longitude,
        .sun_earth_distance_au = distance_au,
        .sun_earth_distance_re = distance_au * au_in_earth_radii(),
        .solar_node_deg = node,
        .solar_central_meridian_deg = central_meridian,
        .hae_to_hee = h::rotation_z(longitude),
        .hae_to_heeq = h::rotation_z(central_meridian) *
                       h::rotation_x(solar_equator_inclination_deg()) * h::rotation_z(node),
    };
}

/**
 * A frame-tagged 3-vector template: @ref Position or @ref FieldVector, and nothing else.
 *
 * The purely rotational transforms are written once over this concept. That is not merely brevity:
 * it means there is no position-specific body in which an origin shift could ever be added to a
 * transform that has no origin shift to make.
 *
 * @tparam V the class template, taking a @ref Frame.
 */
template <template <Frame> class V>
concept FrameTagged = std::same_as<V<Frame::HAE>, Position<Frame::HAE>> ||
                      std::same_as<V<Frame::HAE>, FieldVector<Frame::HAE>>;

/**
 * HAE → HEE: a rotation by the Earth's heliocentric longitude about the ecliptic pole.
 *
 * Both frames are centred on the Sun, so this is a pure rotation for a position and for a field
 * alike — `T(HAE_D→HEE_D) = E(0, 0, λ_geo)`, F&H §3.2.2.
 *
 * @tparam V the vector kind, @ref Position or @ref FieldVector.
 * @param in the HAE components.
 * @param geometry the epoch geometry from @ref helio_geometry.
 * @return the same physical vector, in HEE.
 * @complexity O(1) — nine multiplies and six adds.
 * @alloc none.
 */
template <template <Frame> class V>
    requires FrameTagged<V>
[[nodiscard]] V<Frame::HEE> HAE2HEE(const V<Frame::HAE>& in, const HelioGeometry& geometry) {
    return V<Frame::HEE>{geometry.hae_to_hee * in.v};
}

/**
 * HEE → HAE: the inverse rotation, which for an orthonormal matrix is its transpose.
 *
 * @tparam V the vector kind, @ref Position or @ref FieldVector.
 * @param in the HEE components.
 * @param geometry the epoch geometry from @ref helio_geometry.
 * @return the same physical vector, in HAE.
 * @complexity O(1).
 * @alloc none.
 */
template <template <Frame> class V>
    requires FrameTagged<V>
[[nodiscard]] V<Frame::HAE> HEE2HAE(const V<Frame::HEE>& in, const HelioGeometry& geometry) {
    return V<Frame::HAE>{fixarray::transpose(geometry.hae_to_hee) * in.v};
}

/**
 * HAE → HEEQ: into the frame of the solar equator, with +X on the solar central meridian.
 *
 * `T(HAE_D→HEEQ) = E(Ω, i, θ)` (F&H §3.2.2): swing +X onto the ascending node of the solar equator,
 * tip by the 7.25° inclination, then rotate within the equator until +X lies under the Earth. Both
 * frames are heliocentric, so again there is no translation.
 *
 * @tparam V the vector kind, @ref Position or @ref FieldVector.
 * @param in the HAE components.
 * @param geometry the epoch geometry from @ref helio_geometry.
 * @return the same physical vector, in HEEQ.
 * @complexity O(1).
 * @alloc none.
 */
template <template <Frame> class V>
    requires FrameTagged<V>
[[nodiscard]] V<Frame::HEEQ> HAE2HEEQ(const V<Frame::HAE>& in, const HelioGeometry& geometry) {
    return V<Frame::HEEQ>{geometry.hae_to_heeq * in.v};
}

/**
 * HEEQ → HAE: the inverse of @ref HAE2HEEQ, again by transpose.
 *
 * @tparam V the vector kind, @ref Position or @ref FieldVector.
 * @param in the HEEQ components.
 * @param geometry the epoch geometry from @ref helio_geometry.
 * @return the same physical vector, in HAE.
 * @complexity O(1).
 * @alloc none.
 */
template <template <Frame> class V>
    requires FrameTagged<V>
[[nodiscard]] V<Frame::HAE> HEEQ2HAE(const V<Frame::HEEQ>& in, const HelioGeometry& geometry) {
    return V<Frame::HAE>{fixarray::transpose(geometry.hae_to_heeq) * in.v};
}

// ---- GSE ↔ HEE: the one pair where the origin moves --------------------------------------------
//
// GSE and HEE share the ecliptic pole and differ only in the sense of the Sun–Earth line: +X_GSE
// runs Earth→Sun, +X_HEE runs Sun→Earth. Composing F&H §3.3's `T(HAE→GSE) = E(0, 0, λ_geo + 180°)`
// with §3.2.2's `T(HAE→HEE) = E(0, 0, λ_geo)` leaves exactly `<180°, Z>`: negate X and Y, keep Z,
// with no rounding at all because the half turn is a sign flip rather than a multiply by cos 180°.
//
// The Sun–Earth vector expressed in HEE is `(r₀, 0, 0)`. A POSITION therefore transforms as
// `r_HEE = (r₀ − x_GSE, −y_GSE, z_GSE)`, and the inverse works out to the same expression — the map
// is an involution — while a FIELD VECTOR takes the half turn alone. The four routines are written
// out rather than folded together precisely so that the presence and absence of `r₀` are legible
// side by side.

/**
 * GSE → HEE for a POSITION: half turn about the ecliptic pole, then the shift from the geocentric
 * to the heliocentric origin.
 *
 * @param in the GSE position, in Earth radii.
 * @param geometry the epoch geometry from @ref helio_geometry.
 * @return the HEE position, in Earth radii — the Earth itself (GSE origin) maps to `(r₀, 0, 0)`.
 * @complexity O(1) — two negations and one add.
 * @alloc none.
 */
[[nodiscard]] inline Position<Frame::HEE> GSE2HEE(const Position<Frame::GSE>& in,
                                                  const HelioGeometry& geometry) {
    return Position<Frame::HEE>{
        fixarray::vec3d{geometry.sun_earth_distance_re - in.v[0], -in.v[1], in.v[2]}};
}

/**
 * GSE → HEE for a FIELD VECTOR: the half turn ALONE. A field has a direction and a magnitude but no
 * position, so the change of origin does not touch it.
 *
 * @param in the GSE field, in nanotesla.
 * @param geometry unused — the rotation is a fixed half turn, and the parameter is kept so that the
 *        field and position overloads are called identically and reviewed as a pair.
 * @return the HEE field, in nanotesla.
 * @complexity O(1) — two negations.
 * @alloc none.
 */
[[nodiscard]] inline FieldVector<Frame::HEE> GSE2HEE(const FieldVector<Frame::GSE>& in,
                                                     const HelioGeometry& geometry) {
    static_cast<void>(geometry);
    return FieldVector<Frame::HEE>{fixarray::vec3d{-in.v[0], -in.v[1], in.v[2]}};
}

/**
 * HEE → GSE for a POSITION: the inverse of @ref GSE2HEE, which is the same expression again.
 *
 * @param in the HEE position, in Earth radii.
 * @param geometry the epoch geometry from @ref helio_geometry.
 * @return the GSE position, in Earth radii — the Sun itself (HEE origin) maps to `(r₀, 0, 0)`.
 * @complexity O(1).
 * @alloc none.
 */
[[nodiscard]] inline Position<Frame::GSE> HEE2GSE(const Position<Frame::HEE>& in,
                                                  const HelioGeometry& geometry) {
    return Position<Frame::GSE>{
        fixarray::vec3d{geometry.sun_earth_distance_re - in.v[0], -in.v[1], in.v[2]}};
}

/**
 * HEE → GSE for a FIELD VECTOR: the half turn alone, with no origin shift.
 *
 * @param in the HEE field, in nanotesla.
 * @param geometry unused, as in the forward field overload.
 * @return the GSE field, in nanotesla.
 * @complexity O(1).
 * @alloc none.
 */
[[nodiscard]] inline FieldVector<Frame::GSE> HEE2GSE(const FieldVector<Frame::HEE>& in,
                                                     const HelioGeometry& geometry) {
    static_cast<void>(geometry);
    return FieldVector<Frame::GSE>{fixarray::vec3d{-in.v[0], -in.v[1], in.v[2]}};
}

}  // namespace cheatah::space::irbem
