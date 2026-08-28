#pragma once

/**
 * @file coords_rotations.hpp
 * @brief space.irbem — the Cartesian frame rotations for one epoch, and the transforms built on
 *        them.
 *
 * Every magnetic-field evaluation in this library begins and ends with a frame change. An
 * ephemeris arrives in geographic coordinates, the internal field wants geographic, every external
 * model is defined in solar-magnetospheric, and a drift shell is traced in solar-magnetic — so a
 * single L\* evaluation, which costs ~10⁵ field-model calls, also costs ~10⁵ frame changes. This is
 * the most-executed code in the module, and the foundation everything else stands on.
 *
 * The whole design follows from one observation: the six Cartesian frames here differ only by a
 * rotation, and the rotation depends ONLY on the epoch — not on the point being transformed. It
 * costs a Greenwich sidereal time, a solar ephemeris and half a dozen trigonometric calls; the
 * transform of one position afterwards costs a 3×3 matrix–vector product and nothing else. So
 * @ref Rotations pays the transcendentals once per epoch and stores the matrices, and
 * @ref transform is branch-free arithmetic on a value that fits in registers. Recomputing the
 * angles per point — which is what a naive port does — is roughly two orders of magnitude of
 * wasted work on the hottest path in the library.
 *
 * The second half of the design is that these matrices are ORTHOGONAL. A frame rotation has no
 * scale and no shear, so its inverse is exactly its transpose: a reverse transform is a re-indexed
 * read of the same nine numbers, never a matrix inversion, and never a second stored matrix that
 * could drift out of agreement with the first. Eight matrices therefore cover sixteen directed
 * transforms, and `R · Rᵀ = I` holds to the last bit or two rather than approximately.
 *
 * The frames themselves come from @ref frames.hpp, so the target frame is a template argument and
 * a GEO-to-GSM transform physically cannot be handed a GSM input.
 *
 * ### Sources
 *
 * The transformation chain is Hapgood, *Space physics coordinate transformations: a user guide*,
 * Planet. Space Sci. **40**, 711–717 (1992), with the 1997 corrigendum (Planet. Space Sci. **45**,
 * 1047). Its equation numbers are cited at each step below; Russell, *Cosmic Electrodynamics*
 * **2**, 184 (1971) gives the same chain in a different notation.
 *
 * Greenwich Mean Sidereal Time is the **IAU 1982** expression (Aoki et al., Astron. Astrophys.
 * **105**, 359, 1982), which is what the FK5/mean-equinox convention these frames are defined in
 * calls for; Hapgood's own truncated linear form (his eq. 2) is offered alongside it as
 * @ref GmstModel::Hapgood1992 for bit-comparison against implementations that use it.
 *
 * The geodipole comes from the degree-1 IGRF Gauss coefficients, taken as a parameter — the IGRF
 * table itself lives elsewhere. Its direction follows from the spherical-harmonic definition of
 * the geomagnetic potential (IAGA, *International Geomagnetic Reference Field*), derived in the
 * comment on @ref DipoleCoefficients::axis_geo rather than quoted, so the sign convention is
 * checkable.
 *
 * @note Every angle in the public interface is in DEGREES, matching the papers and IRBEM's own
 *       reporting of the dipole tilt; radians appear only inside the trigonometry.
 * @note The epoch is UT1 expressed as a Julian date. UT1 is what sidereal time is a function of;
 *       feeding UTC instead costs at most |DUT1| < 0.9 s of Earth rotation, i.e. < 0.004° of
 *       longitude, which is below every tolerance in `docs/ERROR_BUDGET.md` but is a real,
 *       deliberate approximation and is the caller's to make.
 */

#include <cmath>
#include <concepts>
#include <cstdint>
#include <numbers>

#include "space/irbem/frames.hpp"

namespace cheatah::space::irbem {

namespace detail {

/// Degrees to radians — the papers are written in degrees, the trigonometry is not.
inline constexpr double kDegToRad = std::numbers::pi / 180.0;

/// Radians to degrees, for reporting angles back in the papers' units.
inline constexpr double kRadToDeg = 180.0 / std::numbers::pi;

/// A quarter turn, in radians — the colatitude/latitude complement in the GEO→MAG construction.
inline constexpr double kHalfPi = std::numbers::pi / 2.0;

/// Seconds of time per degree of rotation: a full turn is 86400 s of sidereal measure over 360°.
inline constexpr double kSecondsPerDegree = 240.0;

/// The Julian date of the J2000.0 epoch, 2000 January 1, 12:00 — the origin of every series here.
inline constexpr double kJdJ2000 = 2451545.0;

/// Days in a Julian century, the unit of the time argument of every series here.
inline constexpr double kDaysPerJulianCentury = 36525.0;

/**
 * Reduce an angle to `[0, 360)`.
 * @param degrees the angle, any magnitude or sign.
 * @return the equivalent angle in `[0, 360)`.
 * @complexity O(1).
 * @alloc none.
 */
inline double wrap_degrees(double degrees) {
    const double wrapped = std::fmod(degrees, 360.0);
    return wrapped < 0.0 ? wrapped + 360.0 : wrapped;
}

/// Hapgood's time argument: the epoch split at the preceding midnight (his §3, eqs. 2 and 5).
struct HapgoodEpoch {
    double centuries;  ///< Julian centuries from J2000.0 to the midnight preceding the epoch.
    double ut_hours;   ///< Universal time since that midnight, in hours.
};

/**
 * Split an epoch the way Hapgood's series are parameterized.
 *
 * Hapgood writes the sidereal time and the solar ephemeris as a function of `T0`, measured in
 * centuries to the PRECEDING MIDNIGHT, plus the hours of UT since it — not as a function of the
 * instant. The split is not cosmetic: the linear coefficients differ between the two terms (36000.770
 * per century versus 15.04107 per hour), so evaluating his formula with a single continuous argument
 * changes the answer by about half a degree.
 *
 * @param jd_ut1 the epoch as a UT1 Julian date.
 * @return the century count to the preceding midnight and the hours of UT since it.
 * @complexity O(1).
 * @alloc none.
 */
inline HapgoodEpoch hapgood_epoch(double jd_ut1) {
    // Modified Julian date runs from midnight, so its integer part IS the preceding midnight.
    const double mjd = jd_ut1 - 2400000.5;
    const double midnight = std::floor(mjd);
    return HapgoodEpoch{(midnight - 51544.5) / kDaysPerJulianCentury, (mjd - midnight) * 24.0};
}

/**
 * The frame rotation about the X axis — a passive rotation, carrying components from the old frame
 * to a frame turned by @p radians about X in the right-handed sense.
 * @param radians the rotation angle.
 * @return the 3×3 rotation.
 * @complexity O(1), two transcendental calls.
 * @alloc none.
 */
inline fixarray::mat3d rot_x(double radians) {
    const double c = std::cos(radians);
    const double s = std::sin(radians);
    return fixarray::mat3d{1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c};
}

/**
 * The frame rotation about the Y axis, in the same passive, right-handed sense as @ref rot_x.
 *
 * @note Hapgood's `<θ, Y>` is the TRANSPOSE of this: his Y matrix is not the cyclic continuation of
 *       his X and Z matrices. That is why the tilt and pole rotations below carry the opposite sign
 *       to the angle in his eqs. (11) and (12). The compositions agree; only the spelling differs,
 *       and the defining properties are asserted in the tests rather than trusted to the sign.
 * @param radians the rotation angle.
 * @return the 3×3 rotation.
 * @complexity O(1), two transcendental calls.
 * @alloc none.
 */
inline fixarray::mat3d rot_y(double radians) {
    const double c = std::cos(radians);
    const double s = std::sin(radians);
    return fixarray::mat3d{c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c};
}

/**
 * The frame rotation about the Z axis, in the same passive, right-handed sense as @ref rot_x.
 * @param radians the rotation angle.
 * @return the 3×3 rotation.
 * @complexity O(1), two transcendental calls.
 * @alloc none.
 */
inline fixarray::mat3d rot_z(double radians) {
    const double c = std::cos(radians);
    const double s = std::sin(radians);
    return fixarray::mat3d{c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0};
}

}  // namespace detail

/**
 * Which published series to use for Greenwich Mean Sidereal Time.
 *
 * They differ by roughly two arcseconds near J2000 — irrelevant to any tolerance in the error
 * budget, but not irrelevant when reproducing another implementation bit for bit, which is the only
 * reason the choice is exposed.
 */
enum class GmstModel : std::uint8_t {
    /// The IAU 1982 series (Aoki et al. 1982) — the default, and the FK5-consistent one.
    Iau1982,
    /// Hapgood 1992 eq. (2), a linear truncation kept for comparison with codes that use it.
    Hapgood1992,
};

/**
 * Greenwich Mean Sidereal Time from the IAU 1982 series.
 *
 * The series (Aoki et al., Astron. Astrophys. 105, 359, 1982) is conventionally written for 0ʰ UT1
 * and then advanced by the elapsed day; the equivalent whole-instant form used here is
 *
 *     GMST [s] = 67310.54841 + (876600ʰ·3600 + 8640184.812866)·T + 0.093104·T² − 6.2e-6·T³
 *
 * with `T` in Julian centuries of UT1 from J2000.0. The leading 67310.54841 s is the familiar
 * 24110.54841 s at 0ʰ plus the half day between midnight and the J2000.0 noon epoch.
 *
 * The large linear coefficient is evaluated exactly rather than as written: `876600·3600 = 3155760000`
 * is precisely `36525 · 86400`, so that term is exactly `86400 · d` for `d` days since J2000.0, and
 * since the result is reduced modulo a day only its FRACTIONAL part matters. Dropping the whole days
 * before multiplying keeps the sum from carrying an integer of order 10⁹ that is about to be
 * discarded. Measured against an extended-precision evaluation of the same series over ±100 years,
 * that is worth a factor of ~300: 1.5e-11° against 4.6e-9°, pinned by
 * `IrbemGmst.TheFractionalDayFoldingIsWorthAFactorOfHundreds`.
 *
 * Neither figure is anywhere near mattering physically — both are microarcseconds, and the
 * coordinate-transform line of `docs/ERROR_BUDGET.md` is 1e-10 RELATIVE. It is done because it is
 * free, not because the naive form would break anything.
 *
 * @param jd_ut1 the epoch as a UT1 Julian date.
 * @return Greenwich Mean Sidereal Time in degrees, reduced to `[0, 360)`.
 * @complexity O(1).
 * @alloc none.
 */
inline double gmst_iau1982_degrees(double jd_ut1) {
    const double days = jd_ut1 - detail::kJdJ2000;
    const double t = days / detail::kDaysPerJulianCentury;
    const double whole_turns_removed = 86400.0 * (days - std::floor(days));
    const double seconds = 67310.54841 + whole_turns_removed + (8640184.812866 * t) +
                           (0.093104 * t * t) - (6.2e-6 * t * t * t);
    return detail::wrap_degrees(seconds / detail::kSecondsPerDegree);
}

/**
 * Greenwich Mean Sidereal Time from Hapgood 1992 eq. (2):
 * `θ = 100.461 + 36000.770·T0 + 15.04107·H` degrees, with `T0` and `H` as in @ref detail::hapgood_epoch.
 *
 * This is a linear truncation of the IAU series and drops its quadratic term, so it drifts slowly
 * away from @ref gmst_iau1982_degrees — about 2 arcseconds near J2000, growing to a few tens of
 * arcseconds a century out. It exists so a result can be reproduced against codes that use it.
 *
 * @param jd_ut1 the epoch as a UT1 Julian date.
 * @return Greenwich Mean Sidereal Time in degrees, reduced to `[0, 360)`.
 * @complexity O(1).
 * @alloc none.
 */
inline double gmst_hapgood_degrees(double jd_ut1) {
    const detail::HapgoodEpoch epoch = detail::hapgood_epoch(jd_ut1);
    return detail::wrap_degrees(100.461 + (36000.770 * epoch.centuries) + (15.04107 * epoch.ut_hours));
}

/**
 * Greenwich Mean Sidereal Time, from whichever published series @p model names.
 * @param jd_ut1 the epoch as a UT1 Julian date.
 * @param model which series to evaluate.
 * @return Greenwich Mean Sidereal Time in degrees, reduced to `[0, 360)`.
 * @complexity O(1).
 * @alloc none.
 */
inline double gmst_degrees(double jd_ut1, GmstModel model) {
    return model == GmstModel::Hapgood1992 ? gmst_hapgood_degrees(jd_ut1)
                                           : gmst_iau1982_degrees(jd_ut1);
}

/**
 * The Sun's apparent position, to the accuracy the GSE frame definition needs.
 *
 * Only the ecliptic longitude and the obliquity enter the transformations — GSE's X axis is the
 * Earth–Sun DIRECTION, so the Sun's distance is irrelevant and its ecliptic latitude is zero by
 * construction of the ecliptic.
 */
struct SolarEphemeris {
    double mean_anomaly_deg;         ///< The Sun's mean anomaly `M`, in degrees.
    double mean_longitude_deg;       ///< The Sun's mean longitude `Λ`, in degrees.
    double ecliptic_longitude_deg;   ///< The Sun's apparent ecliptic longitude `λ☉`, in degrees.
    double obliquity_deg;            ///< The obliquity of the ecliptic `ε`, in degrees.
};

/**
 * The solar ephemeris of Hapgood 1992 eq. (5), itself the *Almanac for Computers* low-precision
 * series:
 *
 *     M   = 357.528 + 35999.050·T0 + 0.04107·H
 *     Λ   = 280.460 + 36000.772·T0 + 0.04107·H
 *     λ☉  = Λ + (1.915 − 0.0048·T0)·sin M + 0.020·sin 2M
 *     ε   = 23.439 − 0.013·T0
 *
 * all in degrees, with `T0` and `H` as in @ref detail::hapgood_epoch. The `sin M` term is the
 * equation of centre truncated at second order in the eccentricity; the series is accurate to about
 * 0.01° in λ☉ over 1950–2050, which propagates to 0.01° in the GSE/GSM axes — far inside the
 * budget for a frame whose defining direction is the Sun.
 *
 * @param jd_ut1 the epoch as a UT1 Julian date.
 * @return the mean anomaly, mean longitude, apparent ecliptic longitude and obliquity, in degrees.
 *         The two longitudes are reduced to `[0, 360)`; the obliquity is not, being small.
 * @complexity O(1).
 * @alloc none.
 */
inline SolarEphemeris solar_ephemeris(double jd_ut1) {
    const detail::HapgoodEpoch epoch = detail::hapgood_epoch(jd_ut1);
    const double mean_anomaly =
        detail::wrap_degrees(357.528 + (35999.050 * epoch.centuries) + (0.04107 * epoch.ut_hours));
    const double mean_longitude =
        detail::wrap_degrees(280.460 + (36000.772 * epoch.centuries) + (0.04107 * epoch.ut_hours));
    const double m = mean_anomaly * detail::kDegToRad;
    const double ecliptic_longitude = detail::wrap_degrees(
        mean_longitude + ((1.915 - (0.0048 * epoch.centuries)) * std::sin(m)) +
        (0.020 * std::sin(2.0 * m)));
    return SolarEphemeris{mean_anomaly, mean_longitude, ecliptic_longitude,
                          23.439 - (0.013 * epoch.centuries)};
}

/**
 * The unit vector from the Earth to the Sun, in GEI.
 *
 * The Sun sits on the ecliptic by definition, so its ecliptic latitude is zero and the standard
 * ecliptic-to-equatorial rotation by the obliquity reduces to
 * `(cos λ☉, cos ε · sin λ☉, sin ε · sin λ☉)`.
 *
 * @param sun the solar ephemeris for the epoch, from @ref solar_ephemeris.
 * @return the Earth-to-Sun unit vector in GEI (true equator and equinox of date).
 * @complexity O(1), three transcendental calls.
 * @alloc none.
 */
inline fixarray::vec3d sun_direction_gei(const SolarEphemeris& sun) {
    const double lambda = sun.ecliptic_longitude_deg * detail::kDegToRad;
    const double eps = sun.obliquity_deg * detail::kDegToRad;
    const double sin_lambda = std::sin(lambda);
    return fixarray::vec3d{std::cos(lambda), std::cos(eps) * sin_lambda, std::sin(eps) * sin_lambda};
}

/**
 * The degree-1 IGRF Gauss coefficients — the centred geomagnetic dipole for one epoch.
 *
 * The IGRF table itself is owned elsewhere; this is the three numbers the frame definitions
 * actually need, in nanotesla, at the epoch of interest (interpolated or extrapolated by whoever
 * owns the table).
 */
struct DipoleCoefficients {
    double g10;  ///< The axial dipole coefficient `g₁⁰`, in nT. Negative for the present-day Earth.
    double g11;  ///< The equatorial dipole coefficient `g₁¹`, in nT.
    double h11;  ///< The equatorial dipole coefficient `h₁¹`, in nT.

    /**
     * The dipole moment magnitude `B₀ = √(g₁⁰² + g₁¹² + h₁¹²)`, the field strength the centred
     * dipole would produce at the magnetic equator on the reference sphere.
     * @return `B₀` in nanotesla; ~30 000 nT for the present-day Earth.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] double moment_nt() const {
        return std::sqrt((g10 * g10) + (g11 * g11) + (h11 * h11));
    }

    /**
     * The dipole axis in GEO, as a unit vector pointing toward the NORTH geomagnetic pole — the +Z
     * axis of both @ref Frame::MAG and @ref Frame::SM.
     *
     * Derivation, so the sign is checkable rather than folklore. The geomagnetic scalar potential's
     * degree-1 part is
     *
     *     V = a·(a/r)²·[ g₁⁰·cos θ + (g₁¹·cos φ + h₁¹·sin φ)·sin θ ]
     *
     * and substituting `cos θ = z/r`, `sin θ cos φ = x/r`, `sin θ sin φ = y/r` collapses it to
     * `V = a³ (m⃗ · r⃗)/r³` with `m⃗ = (g₁¹, h₁¹, g₁⁰)`. That `m⃗` is the physical dipole moment, and
     * because `g₁⁰` is negative it points into the southern hemisphere — which is exactly why the
     * Earth's "north magnetic pole" is a magnetic south pole. The geomagnetic frames orient their
     * +Z along the pole in the NORTHERN hemisphere, so the axis is `−m⃗`, normalized.
     *
     * @return the unit vector `(−g₁¹, −h₁¹, −g₁⁰)/B₀` in geographic Cartesian coordinates.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] fixarray::vec3d axis_geo() const {
        const double b0 = moment_nt();
        return fixarray::vec3d{-g11 / b0, -h11 / b0, -g10 / b0};
    }

    /**
     * The GEOCENTRIC geographic latitude of the north geomagnetic pole.
     *
     * @note Published pole positions are usually GEODETIC. At these latitudes the two differ by
     *       about 0.066°, which is exactly the gap between the value returned here and the
     *       commonly quoted one (80.589° geocentric against 80.65° geodetic for IGRF-13 2020,
     *       80.313° against 80.37° for DGRF 2015). Nothing in this header wants a geodetic angle —
     *       the frames are geocentric — so no conversion is done, but the difference is real and
     *       must not be mistaken for an error in the coefficients.
     * @return degrees, north positive; 80.589° for the IGRF-13 2020 epoch.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] double north_pole_latitude_deg() const {
        return std::asin(axis_geo()[2]) * detail::kRadToDeg;
    }

    /**
     * The geographic east longitude of the north geomagnetic pole.
     * @note Undefined for an exactly axial dipole (`g₁¹ = h₁¹ = 0`), where the pole sits on the
     *       rotation axis and has no longitude at all; `atan2` of two signed zeros then returns
     *       −180°. No IGRF epoch has both coefficients exactly zero, so no branch is spent on it.
     * @return degrees in `(-180, 180]`, east positive; −72.680° for the IGRF-13 2020 epoch.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] double north_pole_longitude_deg() const {
        const fixarray::vec3d axis = axis_geo();
        return std::atan2(axis[1], axis[0]) * detail::kRadToDeg;
    }
};

/**
 * Every rotation among the Cartesian frames, for ONE epoch and one dipole.
 *
 * This is the object the whole header exists to produce. Building it costs the sidereal time, the
 * solar ephemeris and about a dozen trigonometric calls; transforming a position afterwards costs a
 * 3×3 matrix–vector product. Build one per epoch, then transform the whole ephemeris against it —
 * never rebuild it per point.
 *
 * Only eight matrices are stored, because these are orthogonal: the inverse of each is its
 * transpose, which @ref rotation_matrix takes rather than storing a second copy that could drift.
 * Sixteen directed transforms therefore come out of eight stored ones, and `R·Rᵀ = I` holds to
 * roundoff by construction rather than by arithmetic luck.
 *
 * The construction chain, with Hapgood 1992's equation numbers:
 *
 * | matrix | Hapgood | definition |
 * |---|---|---|
 * | `gei_to_geo` | T1, eq. 2 | rotate about Z by the Greenwich sidereal time |
 * | `gei_to_gse` | T2, eq. 3 | rotate about X by the obliquity, then about Z by the solar longitude |
 * | `gse_to_gsm` | T3, eq. 6 | rotate about X until the dipole has no GSE-Y component |
 * | `gsm_to_sm`  | T4, eq. 11 | rotate about Y by the dipole tilt |
 * | `geo_to_mag` | T5, eq. 12 | rotate about Z to the pole's longitude, then about Y to its latitude |
 *
 * The three composites are products of those: `geo_to_gse = T2·T1ᵀ`, `geo_to_gsm = T3·geo_to_gse`,
 * `geo_to_sm = T4·geo_to_gsm`. They are stored, not recomputed, so that GEO — the frame ephemerides
 * arrive in — reaches every other frame in a single product.
 *
 * @note The pairs that are NOT stored (GEI↔GSM, GEI↔SM, GEI↔MAG, GSE↔SM, GSE↔MAG, GSM↔MAG,
 *       SM↔MAG) are a compile error rather than a silent chain through an intermediate frame.
 *       Compose them explicitly through GEO if they are needed; making that visible is the point.
 */
struct Rotations {
    double jd_ut1;                        ///< The epoch these matrices were built for, as a UT1 Julian date.
    double gmst_deg;                      ///< Greenwich Mean Sidereal Time at the epoch, in degrees.
    double sun_ecliptic_longitude_deg;    ///< The Sun's apparent ecliptic longitude `λ☉`, in degrees.
    double obliquity_deg;                 ///< The obliquity of the ecliptic `ε`, in degrees.
    double dipole_tilt_deg;               ///< The geodipole tilt `μ`: the angle of the dipole axis from
                                          ///< the GSM Z axis, positive when the north dipole leans
                                          ///< sunward. Ranges about ±35° over a year.
    double gsm_dipole_angle_deg;          ///< Hapgood's `ψ`: the angle of the dipole's projection in the
                                          ///< GSE Y–Z plane, measured from Z toward Y. This is the angle
                                          ///< `gse_to_gsm` undoes.
    fixarray::vec3d dipole_geo;           ///< The dipole axis (toward the north geomagnetic pole), unit, in GEO.
    fixarray::vec3d dipole_gse;           ///< The same axis in GSE — the vector `ψ` and `μ` are read off.
    fixarray::mat3d gei_to_geo;           ///< Hapgood T1: inertial to Earth-fixed.
    fixarray::mat3d gei_to_gse;           ///< Hapgood T2: inertial to solar-ecliptic.
    fixarray::mat3d gse_to_gsm;           ///< Hapgood T3: solar-ecliptic to solar-magnetospheric.
    fixarray::mat3d gsm_to_sm;            ///< Hapgood T4: solar-magnetospheric to solar-magnetic.
    fixarray::mat3d geo_to_mag;           ///< Hapgood T5: Earth-fixed to centred-dipole geomagnetic.
    fixarray::mat3d geo_to_gse;           ///< `T2·T1ᵀ`, precomputed: the ephemeris frame to solar-ecliptic.
    fixarray::mat3d geo_to_gsm;           ///< `T3·T2·T1ᵀ`, precomputed: the ephemeris frame to the frame
                                          ///< every external field model is defined in.
    fixarray::mat3d geo_to_sm;            ///< `T4·T3·T2·T1ᵀ`, precomputed: the ephemeris frame to the frame
                                          ///< drift shells are traced in.

    /**
     * Build the rotations for one epoch and one dipole.
     *
     * This is where every transcendental in the module is spent. The steps, in order:
     * T1 from the sidereal time; T2 from the solar ephemeris; the dipole axis carried from GEO into
     * GSE through `T2·T1ᵀ`; `ψ` and T3 as the rotation about X that kills the dipole's GSE-Y
     * component (Hapgood eqs. 9–10, which is precisely the statement that GSM's Z axis lies in the
     * plane of the dipole and the Sun); `μ` and T4 as the rotation about Y that carries the
     * remaining tilt onto Z (his eq. 11, so SM's Z axis IS the dipole axis); and T5 from the pole's
     * own latitude and longitude.
     *
     * @param jd_ut1 the epoch as a UT1 Julian date.
     * @param dipole the degree-1 IGRF coefficients for the epoch.
     * @param model which sidereal-time series to use; the IAU 1982 one unless a bit-comparison
     *        against a Hapgood-based implementation is wanted.
     * @return the fully populated rotation set.
     * @complexity O(1) — a fixed ~14 transcendental calls and four 3×3 products, independent of how
     *             many points will later be transformed with it.
     * @alloc none. Every intermediate is an inline-storage `fixarray` value.
     */
    [[nodiscard]] static Rotations at(double jd_ut1, const DipoleCoefficients& dipole,
                                      GmstModel model = GmstModel::Iau1982) {
        const double gmst = gmst_degrees(jd_ut1, model);
        const SolarEphemeris sun = solar_ephemeris(jd_ut1);

        const fixarray::mat3d t1 = detail::rot_z(gmst * detail::kDegToRad);
        const fixarray::mat3d t2 = detail::rot_z(sun.ecliptic_longitude_deg * detail::kDegToRad) *
                                   detail::rot_x(sun.obliquity_deg * detail::kDegToRad);
        const fixarray::mat3d geo_to_gse = t2 * fixarray::transpose(t1);

        const fixarray::vec3d dipole_geo = dipole.axis_geo();
        const fixarray::vec3d dipole_gse = geo_to_gse * dipole_geo;

        // Hapgood eq. (10). ψ is the dipole's clock angle in the GSE Y–Z plane; undoing it about X
        // leaves the dipole with zero Y, which is the definition of GSM.
        const double psi = std::atan2(dipole_gse[1], dipole_gse[2]);
        const fixarray::mat3d t3 = detail::rot_x(-psi);
        const fixarray::mat3d geo_to_gsm = t3 * geo_to_gse;

        // Hapgood eq. (11). X is untouched by T3, and the Y–Z magnitude is preserved by it, so the
        // tilt can be read straight off the GSE components. hypot, not sqrt of a sum of squares: the
        // dipole is a unit vector, but the same expression is reused when it is not.
        const double mu = std::atan2(dipole_gse[0], std::hypot(dipole_gse[1], dipole_gse[2]));
        const fixarray::mat3d t4 = detail::rot_y(mu);
        const fixarray::mat3d geo_to_sm = t4 * geo_to_gsm;

        // Hapgood eq. (12): swing the pole's longitude onto the X–Z plane, then tip its latitude up
        // to +Z. The resulting Y axis is (ẑ_GEO × dipole) normalized, which is MAG's definition.
        const double pole_lat = std::asin(dipole_geo[2]);
        const double pole_lon = std::atan2(dipole_geo[1], dipole_geo[0]);
        const fixarray::mat3d t5 = detail::rot_y(detail::kHalfPi - pole_lat) * detail::rot_z(pole_lon);

        return Rotations{jd_ut1,
                         gmst,
                         sun.ecliptic_longitude_deg,
                         sun.obliquity_deg,
                         mu * detail::kRadToDeg,
                         psi * detail::kRadToDeg,
                         dipole_geo,
                         dipole_gse,
                         t1,
                         t2,
                         t3,
                         t4,
                         t5,
                         geo_to_gse,
                         geo_to_gsm,
                         geo_to_sm};
    }
};

namespace detail {

/**
 * The compile-time table of which stored matrix carries frame @p From to frame @p To.
 *
 * A specialization exists for exactly the eight directions @ref Rotations stores, and names the
 * member rather than wrapping it in an accessor — so the dispatch is a pointer-to-member constant
 * resolved at compile time, with no function to call and nothing to inline away.
 *
 * @tparam To the destination frame.
 * @tparam From the source frame.
 */
template <Frame To, Frame From>
struct StoredRotation {};

/// GEI to GEO is Hapgood's T1, stored directly.
template <>
struct StoredRotation<Frame::GEO, Frame::GEI> {
    /// The member of @ref Rotations holding this rotation.
    static constexpr fixarray::mat3d Rotations::*member = &Rotations::gei_to_geo;
};

/// GEI to GSE is Hapgood's T2, stored directly.
template <>
struct StoredRotation<Frame::GSE, Frame::GEI> {
    /// The member of @ref Rotations holding this rotation.
    static constexpr fixarray::mat3d Rotations::*member = &Rotations::gei_to_gse;
};

/// GSE to GSM is Hapgood's T3, stored directly.
template <>
struct StoredRotation<Frame::GSM, Frame::GSE> {
    /// The member of @ref Rotations holding this rotation.
    static constexpr fixarray::mat3d Rotations::*member = &Rotations::gse_to_gsm;
};

/// GSM to SM is Hapgood's T4, stored directly.
template <>
struct StoredRotation<Frame::SM, Frame::GSM> {
    /// The member of @ref Rotations holding this rotation.
    static constexpr fixarray::mat3d Rotations::*member = &Rotations::gsm_to_sm;
};

/// GEO to MAG is Hapgood's T5, stored directly.
template <>
struct StoredRotation<Frame::MAG, Frame::GEO> {
    /// The member of @ref Rotations holding this rotation.
    static constexpr fixarray::mat3d Rotations::*member = &Rotations::geo_to_mag;
};

/// GEO to GSE is the precomputed composite `T2·T1ᵀ`.
template <>
struct StoredRotation<Frame::GSE, Frame::GEO> {
    /// The member of @ref Rotations holding this rotation.
    static constexpr fixarray::mat3d Rotations::*member = &Rotations::geo_to_gse;
};

/// GEO to GSM is the precomputed composite `T3·T2·T1ᵀ`.
template <>
struct StoredRotation<Frame::GSM, Frame::GEO> {
    /// The member of @ref Rotations holding this rotation.
    static constexpr fixarray::mat3d Rotations::*member = &Rotations::geo_to_gsm;
};

/// GEO to SM is the precomputed composite `T4·T3·T2·T1ᵀ`.
template <>
struct StoredRotation<Frame::SM, Frame::GEO> {
    /// The member of @ref Rotations holding this rotation.
    static constexpr fixarray::mat3d Rotations::*member = &Rotations::geo_to_sm;
};

/// Whether @ref Rotations stores the @p From to @p To rotation in that direction.
template <Frame To, Frame From>
concept StoredDirection = requires { StoredRotation<To, From>::member; };

}  // namespace detail

/**
 * The rotation carrying frame @p From to frame @p To at this epoch.
 *
 * Three constrained overloads share this name and this documentation, because they are three cases
 * of one operation and a caller never chooses between them:
 *
 * - @p To equal to @p From is the identity. It exists so that generic code parameterized on a
 *   target frame compiles for every frame including its own, instead of needing a special case at
 *   every call site.
 * - A pair @ref Rotations stores in that direction is a read of nine already-computed numbers.
 * - The remaining pairs are stored the other way round, and are returned as the TRANSPOSE. A frame
 *   rotation is orthogonal — no scale, no shear — so its inverse is exactly its transpose. That is
 *   why the reverse direction is not stored: the transpose is a re-indexed read that cannot
 *   disagree with the forward matrix, whereas a second stored copy could, and a general inversion
 *   would cost a determinant and nine divisions to reproduce numbers already held exactly.
 *
 * Any other pair fails @ref RotationAvailable and does not compile.
 *
 * @tparam To the destination frame.
 * @tparam From the source frame.
 * @param rotations the epoch's rotations. Unread by the identity case, which keeps the argument
 *        anyway so the identity is not a different call shape from every other transform.
 * @return the 3×3 rotation from @p From to @p To.
 * @complexity O(1).
 * @alloc none.
 */
template <Frame To, Frame From>
    requires CartesianFrame<To> && CartesianFrame<From> && (To == From)
[[nodiscard]] fixarray::mat3d rotation_matrix([[maybe_unused]] const Rotations& rotations) {
    return fixarray::mat3d::identity();
}

/// The stored-direction case of @ref rotation_matrix; see it for the full description.

template <Frame To, Frame From>
    requires CartesianFrame<To> && CartesianFrame<From> && (To != From) &&
             detail::StoredDirection<To, From>
[[nodiscard]] fixarray::mat3d rotation_matrix(const Rotations& rotations) {
    return rotations.*detail::StoredRotation<To, From>::member;
}

/// The transposed-direction case of @ref rotation_matrix; see it for the full description.

template <Frame To, Frame From>
    requires CartesianFrame<To> && CartesianFrame<From> && (To != From) &&
             (!detail::StoredDirection<To, From>) && detail::StoredDirection<From, To>
[[nodiscard]] fixarray::mat3d rotation_matrix(const Rotations& rotations) {
    return fixarray::transpose(rotations.*detail::StoredRotation<From, To>::member);
}

/**
 * Whether a rotation exists between two frames at all: both Cartesian, and either the same frame or
 * one of the eight pairs @ref Rotations stores, in either direction.
 *
 * The `transform` function template is constrained on this rather than left to fail inside its
 * body, so an unsupported
 * pair is reported as an unsatisfied constraint at the call site — and so a `requires` test for
 * transformability answers correctly, which a body-level failure would not.
 *
 * @tparam To the destination frame.
 * @tparam From the source frame.
 */
template <Frame To, Frame From>
concept RotationAvailable =
    CartesianFrame<To> && CartesianFrame<From> &&
    (To == From || detail::StoredDirection<To, From> || detail::StoredDirection<From, To>);

/// The frame-tagged vector templates @ref transform moves between frames: a position or a field.
/// Both are one `vec3d` and nothing else, and both rotate identically — a rotation has no
/// translation, so a position transforms exactly as a direction does.
template <template <Frame> class V>
concept FrameTaggedVector = std::same_as<V<Frame::GEO>, Position<Frame::GEO>> ||
                            std::same_as<V<Frame::GEO>, FieldVector<Frame::GEO>>;

/**
 * Transform a frame-tagged position or field vector into frame @p To.
 *
 * The whole cost of a transform, once @p rotations exists: one 3×3 matrix–vector product, no
 * branches, no allocation, and the frame checked entirely at compile time. Call it as
 * `transform<Frame::GSM>(p, rotations)` — the source frame and whether this is a position or a
 * field are deduced from the argument, so a GEO-to-GSM transform simply cannot be handed a GSM
 * input.
 *
 * Supported pairs are the eight @ref Rotations stores and their eight inverses; anything else is a
 * compile error, and must be composed through GEO deliberately.
 *
 * @tparam To the destination frame.
 * @tparam V the tagged vector template — @ref Position or @ref FieldVector; deduced.
 * @tparam From the source frame; deduced.
 * @param value the position or field vector to transform.
 * @param rotations the epoch's rotations, built once by @ref Rotations::at.
 * @return the same physical quantity, expressed in frame @p To.
 * @complexity O(1): nine multiplies and six adds.
 * @alloc none.
 */
template <Frame To, template <Frame> class V, Frame From>
    requires FrameTaggedVector<V> && RotationAvailable<To, From>
[[nodiscard]] V<To> transform(V<From> value, const Rotations& rotations) {
    return V<To>{rotation_matrix<To, From>(rotations) * value.v};
}

}  // namespace cheatah::space::irbem
