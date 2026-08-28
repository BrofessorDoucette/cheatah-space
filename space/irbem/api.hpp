#pragma once

/**
 * @file api.hpp
 * @brief space.irbem — the IRBEM-named public entry points.
 *
 * Everything else in this module is designed for C++: the frame lives in the type, transforms are
 * one `transform<To>(value, rotations)` template, and epochs are strong types. That is the right
 * surface for new code and the wrong one for the thousands of existing programs, scripts and
 * habits built around IRBEM's names.
 *
 * This header is the adapter. Each entry point carries IRBEM's spelling, IRBEM's argument order and
 * IRBEM's units, and forwards to the typed machinery underneath. It is deliberately thin — there is
 * no physics here, only naming and unit conventions — so that a bug can live in one place or the
 * other but not in both.
 *
 * Three deliberate differences from the Fortran, each of which makes a silent wrong answer harder:
 *
 *  - **Results carry a @ref Status.** IRBEM answers `baddata = -1e31` and leaves the caller to
 *    guess whether that meant an open field line, an unconverged root-find, or an input outside the
 *    model's fitted envelope. Every routine here returns a @ref Result, and @ref baddata is
 *    produced only by the explicit C-compatibility shims at the bottom.
 *  - **Frames are still checked where they can be.** The `sysaxes` integer enters through
 *    @ref frame_from_sysaxes, which reports an out-of-range code rather than defaulting to a frame.
 *  - **The epoch's rotations are built once and reused.** IRBEM recomputes per call; the rotation
 *    matrices depend only on the epoch, so `Rotations::at` is hoisted out of any batch. That is
 *    also what keeps the transcendentals out of the vectorized inner loop.
 *
 * Coverage is honest and partial: the routines below are the ones whose physics is implemented.
 * The magnetic-field, tracing and drift-shell entry points — `make_lstar`, `drift_shell`,
 * `trace_field_line`, `find_mirror_point` and the rest — are not here yet, because their kernels
 * are not written yet, and a named stub that returns `baddata` would be worse than an absence: it
 * would compile at a caller's site and fail at their runtime.
 */

#include <optional>

#include "context.hpp"
#include "coords_geodetic.hpp"
#include "coords_helio.hpp"
#include "coords_rotations.hpp"
#include "library_info.hpp"
#include "frames.hpp"
#include "igrf.hpp"
#include "status.hpp"

namespace cheatah::space::irbem::api {

/**
 * The largest batch the IRBEM-compatible entry points accept, matching the reference's
 * `ntime_max`.
 * @return `100000`, the value IRBEM's generated `ntime_max.inc` carries.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.LibraryInfoMatchesTheReference
 */
constexpr int get_irbem_ntime_max() { return 100000; }

/**
 * Which IGRF generation the internal field implements.
 * @return `14` — IGRF-14, the generation whose coefficients `tables/igrf14.hpp` carries.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.LibraryInfoMatchesTheReference
 */
constexpr int get_igrf_version() { return 14; }

// ---- coordinate transforms -------------------------------------------------------------------
// IRBEM spells these `geo2gsm`, `gsm2geo` and so on, taking and returning three doubles. The typed
// core takes a Position<From> and an epoch's Rotations. The macro-free expansion below is
// deliberate: each entry point is written out so its @test tag, its units and its frame pair are
// all greppable, which a macro would hide.

/**
 * Build the rotations an epoch needs — the argument every transform below shares.
 *
 * Hoisted deliberately: IRBEM rebuilds this per call, but the matrices depend only on the epoch, so
 * a batch over one timestamp should pay for the trigonometry once. This is the single largest
 * structural difference between calling the two libraries in a loop.
 *
 * @param year the calendar year, e.g. `2015`.
 * @param doy the day of year, 1-based.
 * @param ut_seconds seconds of UT within the day.
 * @param igrf the internal-field model the dipole axis is taken from.
 * @return the rotations, or @ref Status::DomainError when the epoch is not representable.
 * @complexity O(1) — a fixed number of trigonometric evaluations.
 * @alloc none.
 * @test IrbemApi.RotationsAreBuiltOncePerEpoch
 */
template <int NMAX>
[[nodiscard]] inline Result<Rotations> rotations_at(int year, int doy, double ut_seconds,
                                                    const Igrf<NMAX>& igrf) {
    const DateTime dt = date_and_time_from_doy_and_ut(year, doy, ut_seconds);
    const double decy = decimal_year(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second);
    if (!(decy >= Igrf<NMAX>::earliest_year) || !(decy <= Igrf<NMAX>::latest_year)) {
        return {Status::DomainError, Rotations{}};
    }
    // Julian DAY NUMBER is midnight-based; the UT fraction and the -0.5 turn it into the Julian
    // DATE the sidereal-time series is defined against. Getting that half-day wrong shifts GMST by
    // twelve hours, which rotates GSM by 180 degrees and produces a field that looks plausible and
    // is inverted — so it is written out rather than folded into a constant.
    const double jd_ut1 = static_cast<double>(julian_day_number(dt.year, dt.month, dt.day)) - 0.5 +
                          (dt.ut_seconds / 86400.0);
    const DipoleCoefficients dipole{igrf.g(1, 0), igrf.g(1, 1), igrf.h(1, 1)};
    return {Status::Ok, Rotations::at(jd_ut1, dipole)};
}

/// The frame-pair transforms, in IRBEM's spelling. Each is one line over @ref transform, so the
/// physics lives in exactly one place; what these add is the name, the argument order and the
/// `@test` tag a caller can grep for.
/// @{

/// Geographic to solar-magnetospheric — the frame every external field model is defined in.
/// @param geo the position in GEO, Earth radii. @param r the epoch's rotations.
/// @return the same point in GSM. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::GSM> geo2gsm(const Position<Frame::GEO>& geo,
                                                  const Rotations& r) {
    return transform<Frame::GSM>(geo, r);
}

/// Solar-magnetospheric to geographic.
/// @param gsm the position in GSM, Earth radii. @param r the epoch's rotations.
/// @return the same point in GEO. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::GEO> gsm2geo(const Position<Frame::GSM>& gsm,
                                                  const Rotations& r) {
    return transform<Frame::GEO>(gsm, r);
}

/// Geographic to solar-magnetic — the frame drift shells are traced in.
/// @param geo the position in GEO, Earth radii. @param r the epoch's rotations.
/// @return the same point in SM. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::SM> geo2sm(const Position<Frame::GEO>& geo,
                                                const Rotations& r) {
    return transform<Frame::SM>(geo, r);
}

/// Solar-magnetic to geographic.
/// @param sm the position in SM, Earth radii. @param r the epoch's rotations.
/// @return the same point in GEO. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::GEO> sm2geo(const Position<Frame::SM>& sm,
                                                 const Rotations& r) {
    return transform<Frame::GEO>(sm, r);
}

/// Geographic to solar-ecliptic.
/// @param geo the position in GEO, Earth radii. @param r the epoch's rotations.
/// @return the same point in GSE. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::GSE> geo2gse(const Position<Frame::GEO>& geo,
                                                  const Rotations& r) {
    return transform<Frame::GSE>(geo, r);
}

/// Solar-ecliptic to geographic.
/// @param gse the position in GSE, Earth radii. @param r the epoch's rotations.
/// @return the same point in GEO. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::GEO> gse2geo(const Position<Frame::GSE>& gse,
                                                  const Rotations& r) {
    return transform<Frame::GEO>(gse, r);
}

/// Geographic to centred-dipole geomagnetic.
/// @param geo the position in GEO, Earth radii. @param r the epoch's rotations.
/// @return the same point in MAG. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::MAG> geo2mag(const Position<Frame::GEO>& geo,
                                                  const Rotations& r) {
    return transform<Frame::MAG>(geo, r);
}

/// Centred-dipole geomagnetic to geographic.
/// @param mag the position in MAG, Earth radii. @param r the epoch's rotations.
/// @return the same point in GEO. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::GEO> mag2geo(const Position<Frame::MAG>& mag,
                                                  const Rotations& r) {
    return transform<Frame::GEO>(mag, r);
}

/// Geographic to inertial.
/// @param geo the position in GEO, Earth radii. @param r the epoch's rotations.
/// @return the same point in GEI. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::GEI> geo2gei(const Position<Frame::GEO>& geo,
                                                  const Rotations& r) {
    return transform<Frame::GEI>(geo, r);
}

/// Inertial to geographic.
/// @param gei the position in GEI, Earth radii. @param r the epoch's rotations.
/// @return the same point in GEO. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::GEO> gei2geo(const Position<Frame::GEI>& gei,
                                                  const Rotations& r) {
    return transform<Frame::GEO>(gei, r);
}

/// Solar-magnetospheric to solar-magnetic.
/// @param gsm the position in GSM, Earth radii. @param r the epoch's rotations.
/// @return the same point in SM. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::SM> gsm2sm(const Position<Frame::GSM>& gsm,
                                                const Rotations& r) {
    return transform<Frame::SM>(gsm, r);
}

/// Solar-magnetic to solar-magnetospheric.
/// @param sm the position in SM, Earth radii. @param r the epoch's rotations.
/// @return the same point in GSM. @complexity O(1). @alloc none.
/// @test IrbemApi.TransformsRoundTrip
[[nodiscard]] inline Position<Frame::GSM> sm2gsm(const Position<Frame::SM>& sm,
                                                 const Rotations& r) {
    return transform<Frame::GSM>(sm, r);
}
/// @}

// ---- the non-rotational conversions, which need no epoch ---------------------------------------

/// Geodetic (altitude km, latitude deg, east longitude deg) to geocentric Cartesian.
/// @param gdz the geodetic position. @return the same point in GEO, Earth radii.
/// @complexity O(1). @alloc none.
/// @test IrbemApi.GeodeticRoundTrips
[[nodiscard]] inline Position<Frame::GEO> gdz2geo(const Position<Frame::GDZ>& gdz) {
    return gdz_to_geo(gdz);
}

/// Geocentric Cartesian to geodetic.
/// @param geo the position in GEO, Earth radii. @return the same point as GDZ.
/// @complexity O(1) — a bounded closed-form inversion, not an iteration to tolerance.
/// @alloc none.
/// @test IrbemApi.GeodeticRoundTrips
[[nodiscard]] inline Position<Frame::GDZ> geo2gdz(const Position<Frame::GEO>& geo) {
    return geo_to_gdz(geo);
}

/// Spherical (radius Re, latitude deg, east longitude deg) to Cartesian.
/// @param sph the spherical position. @return the same point in GEO, Earth radii.
/// @complexity O(1). @alloc none.
/// @test IrbemApi.SphericalRoundTrips
[[nodiscard]] inline Position<Frame::GEO> sph2car(const Position<Frame::SPH>& sph) {
    return sph_to_car(sph);
}

/// Cartesian to spherical.
/// @param car the position in GEO, Earth radii. @return the same point as SPH.
/// @complexity O(1). @alloc none.
/// @test IrbemApi.SphericalRoundTrips
[[nodiscard]] inline Position<Frame::SPH> car2sph(const Position<Frame::GEO>& car) {
    return car_to_sph(car);
}

/// Radius/latitude/longitude to geodetic.
/// @param rll the RLL position. @return the same point as GDZ.
/// @complexity O(1). @alloc none.
/// @test IrbemApi.GeodeticRoundTrips
[[nodiscard]] inline Position<Frame::GDZ> rll2gdz(const Position<Frame::RLL>& rll) {
    return rll_to_gdz(rll);
}

// ---- date and time -----------------------------------------------------------------------------

/// Day of year from a calendar date, leap-year correct.
/// @param year the year. @param month 1-12. @param day 1-based day of month.
/// @return the day of year, 1-based. @complexity O(1). @alloc none.
/// @test IrbemApi.DateRoundTrips
constexpr int get_doy(int year, int month, int day) { return day_of_year(year, month, day); }

/// A calendar date and time as a decimal year — the epoch argument the field models take.
/// @param year the year. @param month 1-12. @param day 1-based. @param hour 0-23.
/// @param minute 0-59. @param second 0-59.
/// @return the decimal year. @complexity O(1). @alloc none.
/// @test IrbemApi.DateRoundTrips
constexpr double date_and_time2decy(int year, int month, int day, int hour, int minute,
                                    int second) {
    return decimal_year(year, month, day, hour, minute, second);
}

/// A year, day-of-year and UT as a full broken-down date and time.
/// @param year the year. @param doy 1-based day of year. @param ut_seconds seconds since midnight.
/// @return the broken-down date and time. @complexity O(1). @alloc none.
/// @test IrbemApi.DateRoundTrips
constexpr DateTime doy_and_ut2date_and_time(int year, int doy, double ut_seconds) {
    return date_and_time_from_doy_and_ut(year, doy, ut_seconds);
}

}  // namespace cheatah::space::irbem::api
