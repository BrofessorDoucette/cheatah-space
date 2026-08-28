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
 * no physics here, only naming, unit conventions and the runtime-to-typed frame dispatch — so that
 * a bug can live in one place or the other but not in both.
 *
 * Four deliberate differences from the Fortran, each of which makes a silent wrong answer harder:
 *
 *  - **Results carry a @ref Status.** IRBEM answers `baddata = -1e31` and leaves the caller to
 *    guess whether that meant an open field line, an unconverged root-find, or an input outside the
 *    model's fitted envelope. Every routine here returns a @ref Result. The one place @ref baddata
 *    survives is inside @ref MagneticCoordinates, whose six fields are IRBEM's six output ARRAYS
 *    and have to carry a per-field "not computed" marker of their own — the record's status lives
 *    beside it in a parallel span, exactly as IRBEM's does, and the sentinel is what fills a field
 *    that was never reached.
 *  - **Frames are still checked where they can be.** The `sysaxes` integer enters through
 *    @ref frame_from_sysaxes, which reports an out-of-range code rather than defaulting to a frame.
 *    Measured, the reference prints `sysaxesOUT out of range !` on stdout and fills the output with
 *    `baddata`; @ref coord_trans_vec returns @ref Status::DomainError and touches nothing.
 *  - **The epoch's rotations are built once and reused.** IRBEM recomputes per call; the rotation
 *    matrices depend only on the epoch, so `Rotations::at` is hoisted out of any batch. That is
 *    also what keeps the transcendentals out of the vectorized inner loop.
 *  - **The heliospheric frames stay in Earth radii.** This one is a genuine divergence and is
 *    called out rather than buried. IRBEM's `sysaxes` table lists HEE, HAE and HEEQ as "Re", but
 *    its `GSE2HEE` entry point returns AU — measured, `gse2hee1_` maps GSE `(1, 0, 0)` Re to HEE
 *    `(1.01655079, 0, 0)` at 2015-180 12:00 UT, which is `r₀ − 1 Re` expressed in AU. The two
 *    conventions disagree inside IRBEM's own documentation, and @ref Position carries Earth radii
 *    for every other Cartesian frame in this module, so Earth radii is what these return.
 *    Multiply by @ref au_in_earth_radii to reach IRBEM's AU spelling at that one boundary.
 *
 * ## The batch forms are the point, not a convenience
 *
 * Every routine here that takes more than one point has an `ntime` form, because a per-point loop
 * cannot be accelerated: the parallelism in this library is entirely ACROSS points, and it only
 * becomes available if the caller hands over the whole batch at once. Which lane a batch actually
 * runs on is measured, not assumed, and the answer differs by routine:
 *
 * | batch lane | arithmetic intensity | today |
 * |---|---|---|
 * | @ref coord_trans_vec, @ref transform_vec, @ref get_mlt_vec | ~0.5 flops/byte | host, always |
 * | @ref make_lstar_vec, @ref make_lstar_shell_splitting | ~9 400 flops/byte | device above ~512 points, 18.6× |
 *
 * A frame rotation is nine multiplies over 48 bytes moved. Measured on an RTX 3070 Ti, the module's
 * own streaming dipole kernel at that intensity LOSES to the host, 0.69×; the trace kernel wins
 * 48.9× at 65 536 lines. So the coordinate lanes here return `false` for "the device serviced this
 * call" and will keep doing so: the shape exists so the answer is askable, and so a future device
 * lane needs no signature change, not because a device is coming for these. The L\* lanes return
 * the real answer — `true` when the device ran the traces — because for them it is the whole point:
 * one L\* point is only `Nder = 25` field lines, an order below the ~512-line crossover, so a
 * single `make_lstar` runs on the host however much hardware is present and only `ntime × Nder`
 * reaches the device. Note that 18.6× is the WHOLE-PIPELINE figure through @ref make_lstar_vec,
 * not the 48.9× of the bare trace kernel: an L\* also spends a serial root-find and a flux
 * quadrature, and quoting the kernel's number for the entry point would be quoting the best stage
 * as if it were the total.
 *
 * ## What is here, and what is not
 *
 * Coverage is honest and partial. Implemented: the whole coordinate-transformation group
 * (`COORD_TRANS_VEC`, the eleven geographic pairs, `SPH2CAR`/`CAR2SPH`, `RLL2GDZ`, and the six
 * heliospheric pairs), `GET_MLT`, `MAKE_LSTAR`, `MAKE_LSTAR_SHELL_SPLITTING` and `LSTAR_PHI` at
 * `kext = 0`, the date-and-time group (`JULDAY`, `CALDAT`, `GET_DOY`,
 * `DECY2DATE_AND_TIME`, `DATE_AND_TIME2DECY`, `DOY_AND_UT2DATE_AND_TIME`) and the library-info
 * group (`GET_IRBEM_NTIME_MAX`, `GET_IGRF_VERSION`, `IRBEM_FORTRAN_VERSION`,
 * `IRBEM_FORTRAN_RELEASE`).
 *
 * Absent, because the physics underneath them is not written yet, and a named stub returning
 * `baddata` would be strictly worse than an absence — it would compile at a caller's site and fail
 * at their runtime:
 *
 *  - **The rest of the drift-shell family** — `LANDI2LSTAR`, `LANDI2LSTAR_SHELL_SPLITTING`,
 *    `EMPIRICALLSTAR`, `DRIFT_SHELL`, `DRIFT_BOUNCE_ORBIT`. The first two invert `(L_m, I)` to
 *    `L*` without a starting position, `EMPIRICALLSTAR` is a fitted formula rather than a
 *    computation, and the last two return the traced shell itself. `driftshell.hpp` computes `L*`
 *    from a POSITION and returns scalars, so none of these five is a wrapper of what exists.
 *  - **The field-line points of interest** — `FIND_MIRROR_POINT`, `FIND_MAGEQUATOR`,
 *    `FIND_FOOT_POINT`, `TRACE_FIELD_LINE`, `TRACE_FIELD_LINE_TOWARD_EARTH`. @ref lstar.hpp traces
 *    without storing a path, by design; these return the path, which is a different (bandwidth-
 *    bound) routine rather than a wrapper of the one that exists.
 *  - **The field evaluators** — `GET_FIELD_MULTI`, `GET_BDERIVS`, `COMPUTE_GRAD_CURV_CURL`,
 *    `GET_HEMI_MULTI`. These are defined for every `kext`, and only `kext = 0` (IGRF internal) is
 *    implemented; the external models (Tsyganenko and the rest) are absent, so exposing the routine
 *    would mean a name that answers for one of fifteen keys.
 *  - **The radiation-belt, atmosphere and orbit groups** — `FLY_IN_NASA_AEAP`, `GET_AE8_AP8_FLUX`,
 *    `FLY_IN_AFRL_CRRES`, `GET_CRRES_FLUX`, `FLY_IN_IGE`, `FLY_IN_MEO_GNSS`, `SHIELDOSE2`,
 *    `MSIS86`, `MSISE90`, `NRLMSISE00`, `SGP4_TLE`, `SGP4_ELE`, `RV2COE`. No part of any of these
 *    exists in this module.
 *
 * ## One convention difference worth knowing before you diff against the reference
 *
 * The reference samples the IGRF secular variation ONCE PER YEAR, at the year's MIDPOINT — and it
 * says so: `docs/source/api/general_information.rst` documents `options(2) = 0` as "initialize IGRF
 * field once per year (year.5)", with any `options(2) = n` meaning an `n`-day update cadence
 * instead. Confirmed by measurement, since only the measurement settles what the entry points that
 * take NO options do: GEO→MAG is a pure dipole rotation, and sweeping day-of-year through 2015
 * against `coord_trans_vec1_` gives a residual of 5.2e-4 Re at doy 1 falling linearly through
 * 8.5e-6 at doy 180 and climbing symmetrically to 5.1e-4 at doy 364. Feed this module an
 * `Igrf<>::at(year + 0.5)` and the same comparison collapses to 3.3e-16 — bit level.
 *
 * Two consequences, and they differ by entry point. `coord_trans_vec1_` and `get_mlt1_` carry no
 * `options` argument at all, so year.5 is FORCED there and a differential test that uses the exact
 * date is measuring that convention rather than either library's arithmetic. `make_lstar1_` DOES
 * carry one: at `options(2) = 1` the reference updates daily, and its `Blocal` then reproduces this
 * module's exact-date `Igrf<13>` to 3.7e-7 relative rather than the 2.2e-6 it manages at
 * `options(2) = 0` (measured, GEO `(4, 0.7, -0.4)`, 2015-180 12:00 UT). This module keeps the exact
 * date throughout, because a dipole axis a third of a year stale is an error, not a feature.
 *
 * The second half of "matched options" is the TRUNCATION: the reference's internal field stops at
 * degree 10, and this module defaults to IGRF-14's published degree 13. On `Blocal` — a pure field
 * evaluation with no algorithm in it — the two choices are worth two and a half orders: `Igrf<10>`
 * at mid-year reproduces `make_lstar1_`'s `Blocal` to 6e-16 relative, `Igrf<13>` at mid-year to
 * 2e-9, and `Igrf<13>` at the exact date only to 3e-4. Above ~3 Re the degree 11–13 terms are
 * negligible and the epoch dominates; near the surface it is the other way round.
 *
 * At matched epochs, what remains is model difference with a name attached to each leg. The figures
 * below are the worst relative componentwise gap over a 296-point × 8-epoch sweep (radii 1.02–10
 * Re, latitudes ±89°, longitudes ±175°, 1965–2029) rather than over one point at one longitude,
 * because a single point hides two of them entirely: GEO exact and SPH 2.4e-16 once the longitude
 * conventions are matched (see @ref detail::sph_in_irbem_longitude — before that fold the raw gap
 * is a full 360° of longitude for every point with a negative GEO y); MAG 1.6e-15; RLL 4.2e-8 and
 * GDZ 6.8e-5 km of altitude off the polar axis (the reference's iterative geodetic latitude against
 * this module's closed form — and the gap grows towards the surface, so the 3e-10 a single
 * geosynchronous point shows is not the number to quote); GEI 1.9e-4 (sidereal time); GSM and SM
 * 2.5e-4 and GSE 6.3e-4 (the solar ephemeris).
 *
 * ON THE POLAR AXIS THE REFERENCE IS WRONG AND THIS MODULE IS NOT, so that leg is a divergence
 * rather than a residual. Measured: `coord_trans_vec1_(1 → 0)` at GEO `(0, 0, 1)` returns a
 * geodetic altitude of exactly 0 km and at `(0, 0, 1.5)` returns 3178.376157 km. At the pole the
 * geodetic height is `|z| − b` with `b = a(1 − f) = 6356.752314245` km and the reference's own
 * `Re = 6371.2` km (measured: it reports −6.937 km at GEO `(1, 0, 0)`), i.e. 14.447685755 km and
 * 3200.047685755 km — which is what this module returns, exactly. The reference is 14.4 km low at
 * `r = 1` and 21.7 km low at `r = 1.5`; its `(0, 0, 0)` answer is a latitude of 180°. RLL and SPH
 * are unaffected: both carry the radius rather than an altitude, and both agree at the pole.
 *
 * Two `sysaxes` regions are also deliberately unreachable through @ref coord_trans_vec, for
 * different reasons. Codes 9–11 (HEE, HAE, HEEQ) name frames this module DOES implement, but they
 * need a @ref HelioGeometry — a different epoch object from @ref Rotations — so they are reached
 * through the named pairs below instead. Codes 12–14 (TOD, J2000, TEME) have no @ref Frame here at
 * all: IRBEM documents TOD as identical to its GEI, while J2000 and TEME need precession and
 * nutation series this module does not carry.
 */

#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

#include "context.hpp"
#include "coords_geodetic.hpp"
#include "coords_helio.hpp"
#include "coords_rotations.hpp"
#include "driftshell.hpp"
#include "library_info.hpp"
#include "frames.hpp"
#include "igrf.hpp"
#include "status.hpp"

namespace cheatah::space::irbem::api {

// ---- library information -----------------------------------------------------------------------
// Four constants describing THIS implementation. Each forwards to @ref library_info.hpp rather than
// repeating the literal: the number and the prose explaining where it came from belong together,
// and two copies of `100000` in one module is exactly the drift this module is built to avoid.

/**
 * The largest batch the IRBEM-compatible entry points accept, matching the reference's
 * `ntime_max`.
 * @return `100000`, the value IRBEM's generated `ntime_max.inc` carries; confirmed by calling
 *         `get_irbem_ntime_max1_` on the shipped shared library.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.LibraryInfoMatchesTheReference
 */
constexpr int get_irbem_ntime_max() { return max_batch_times(); }

/**
 * Which IGRF generation the internal field implements.
 * @return `14` — IGRF-14, the generation whose coefficients `tables/igrf14.hpp` carries. The
 *         shipped reference returns the same, measured.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.LibraryInfoMatchesTheReference
 */
constexpr int get_igrf_version() { return igrf_generation(); }

/**
 * Implementation version — the shim standing in for IRBEM's `IRBEM_FORTRAN_VERSION`.
 * @return a monotonically increasing integer identifying THIS C++ implementation. It is NOT an
 *         IRBEM revision: see @ref implementation_version for why returning a plausible-looking one
 *         would be worse than useless to a caller feature-detecting against it.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.LibraryInfoMatchesTheReference
 */
constexpr int irbem_fortran_version() { return implementation_version(); }

/**
 * Implementation release tag — the shim standing in for `IRBEM_FORTRAN_RELEASE`.
 * @return a static, never-empty view naming this implementation in words, so a log line carrying it
 *         cannot be mistaken for the reference's own tag (measured: a git short hash, `e7cecb0`,
 *         space-padded into an 80-character buffer).
 * @complexity O(1).
 * @alloc none — the view refers to a string literal with static storage duration.
 * @test IrbemApi.LibraryInfoMatchesTheReference
 */
constexpr std::string_view irbem_fortran_release() { return implementation_release(); }

// ---- coordinate transforms -----------------------------------------------------------------
// IRBEM spells these `geo2gsm`, `gsm2geo` and so on, taking and returning three doubles. The typed
// core takes a Position<From> and an epoch's Rotations. The macro-free expansion below is
// deliberate: each entry point is written out so its test tag, its units and its frame pair are
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

/// Cartesian to spherical. East longitude comes back in `(-180, 180]`, this module's convention
/// everywhere; @ref coord_trans with `sysaxes_out = 7` folds the SAME answer into the reference's
/// `[0, 360)` instead. Both name one meridian; see @ref detail::sph_in_irbem_longitude.
/// @param car the position in GEO, Earth radii. @return the same point as SPH.
/// @complexity O(1). @alloc none.
/// @test IrbemApi.SphericalRoundTrips
/// @test IrbemApi.CoordTransUsesTheReferenceLongitudeConventions
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

// ---- the heliospheric pairs --------------------------------------------------------------------
// The seam between the geocentric and the heliocentric literatures. See @ref coords_helio.hpp for
// why GSE↔HEE is written out per vector kind while the other four are written once: GSE and HEE do
// not share an origin, so a POSITION picks up the Sun–Earth vector and a FIELD must not, and the
// ~23 000 Re error that confusion produces is the classic defect of that file.
//
// Units are Earth radii on BOTH sides — see the file brief's fourth bullet, and multiply by
// @ref au_in_earth_radii for IRBEM's AU spelling of the heliospheric half.
/// @{

/// Heliocentric Aries Ecliptic to Heliocentric Earth Ecliptic — a rotation by the Earth's
/// heliocentric longitude, with no origin shift (both frames are centred on the Sun).
/// @tparam V the vector kind, @ref Position or @ref FieldVector; deduced.
/// @param hae the HAE components. @param g the epoch geometry from @ref helio_geometry.
/// @return the same physical vector in HAE. @complexity O(1). @alloc none.
/// @test IrbemApi.HeliosphericPairsRoundTrip
template <template <Frame> class V>
    requires FrameTagged<V>
[[nodiscard]] inline V<Frame::HEE> hae2hee(const V<Frame::HAE>& hae, const HelioGeometry& g) {
    return HAE2HEE<V>(hae, g);
}

/// Heliocentric Earth Ecliptic to Heliocentric Aries Ecliptic.
/// @tparam V the vector kind, @ref Position or @ref FieldVector; deduced.
/// @param hee the HEE components. @param g the epoch geometry from @ref helio_geometry.
/// @return the same physical vector in HAE. @complexity O(1). @alloc none.
/// @test IrbemApi.HeliosphericPairsRoundTrip
template <template <Frame> class V>
    requires FrameTagged<V>
[[nodiscard]] inline V<Frame::HAE> hee2hae(const V<Frame::HEE>& hee, const HelioGeometry& g) {
    return HEE2HAE<V>(hee, g);
}

/// Heliocentric Aries Ecliptic to Heliocentric Earth Equatorial — into the frame of the solar
/// equator, with +X on the solar central meridian.
/// @tparam V the vector kind, @ref Position or @ref FieldVector; deduced.
/// @param hae the HAE components. @param g the epoch geometry from @ref helio_geometry.
/// @return the same physical vector in HEEQ. @complexity O(1). @alloc none.
/// @test IrbemApi.HeliosphericPairsRoundTrip
template <template <Frame> class V>
    requires FrameTagged<V>
[[nodiscard]] inline V<Frame::HEEQ> hae2heeq(const V<Frame::HAE>& hae, const HelioGeometry& g) {
    return HAE2HEEQ<V>(hae, g);
}

/// Heliocentric Earth Equatorial to Heliocentric Aries Ecliptic.
/// @tparam V the vector kind, @ref Position or @ref FieldVector; deduced.
/// @param heeq the HEEQ components. @param g the epoch geometry from @ref helio_geometry.
/// @return the same physical vector in HAE. @complexity O(1). @alloc none.
/// @test IrbemApi.HeliosphericPairsRoundTrip
template <template <Frame> class V>
    requires FrameTagged<V>
[[nodiscard]] inline V<Frame::HAE> heeq2hae(const V<Frame::HEEQ>& heeq, const HelioGeometry& g) {
    return HEEQ2HAE<V>(heeq, g);
}

/// Solar-ecliptic to Heliocentric Earth Ecliptic, for a POSITION: the half turn about the ecliptic
/// pole AND the shift from the geocentric to the heliocentric origin.
/// @param gse the GSE position, Earth radii. @param g the epoch geometry.
/// @return the HEE position, Earth radii. @complexity O(1). @alloc none.
/// @test IrbemApi.HeliosphericOriginShiftAppliesToPositionsOnly
[[nodiscard]] inline Position<Frame::HEE> gse2hee(const Position<Frame::GSE>& gse,
                                                  const HelioGeometry& g) {
    return GSE2HEE(gse, g);
}

/// Solar-ecliptic to Heliocentric Earth Ecliptic, for a FIELD VECTOR: the half turn ALONE.
/// @param gse the GSE field, nT. @param g the epoch geometry; unread, kept so the two overloads are
///        called identically and reviewed as a pair.
/// @return the HEE field, nT. @complexity O(1). @alloc none.
/// @test IrbemApi.HeliosphericOriginShiftAppliesToPositionsOnly
[[nodiscard]] inline FieldVector<Frame::HEE> gse2hee(const FieldVector<Frame::GSE>& gse,
                                                     const HelioGeometry& g) {
    return GSE2HEE(gse, g);
}

/// Heliocentric Earth Ecliptic to solar-ecliptic, for a POSITION.
/// @param hee the HEE position, Earth radii. @param g the epoch geometry.
/// @return the GSE position, Earth radii. @complexity O(1). @alloc none.
/// @test IrbemApi.HeliosphericOriginShiftAppliesToPositionsOnly
[[nodiscard]] inline Position<Frame::GSE> hee2gse(const Position<Frame::HEE>& hee,
                                                  const HelioGeometry& g) {
    return HEE2GSE(hee, g);
}

/// Heliocentric Earth Ecliptic to solar-ecliptic, for a FIELD VECTOR: the half turn alone.
/// @param hee the HEE field, nT. @param g the epoch geometry; unread, as in the forward overload.
/// @return the GSE field, nT. @complexity O(1). @alloc none.
/// @test IrbemApi.HeliosphericOriginShiftAppliesToPositionsOnly
[[nodiscard]] inline FieldVector<Frame::GSE> hee2gse(const FieldVector<Frame::HEE>& hee,
                                                     const HelioGeometry& g) {
    return HEE2GSE(hee, g);
}
/// @}

// ---- magnetic local time -----------------------------------------------------------------------

/// Hours of magnetic local time per radian of SM longitude: a full turn is 24 h.
inline constexpr double hours_per_radian = 12.0 / std::numbers::pi;

/**
 * Magnetic local time at a geographic position — IRBEM's `GET_MLT`.
 *
 * MLT is the SM azimuth read as a clock: `MLT = 12 + λ_SM / 15°`, so local magnetic noon (the SM
 * +X half-plane, which by construction contains the Sun-facing side of the dipole meridian) is
 * 12 h and midnight is 0 h. Nothing but the SM rotation enters it, which is why this routine takes
 * the epoch's @ref Rotations and no field model.
 *
 * Verified as a black box against the shipped reference: over forty points spanning 1965–2029,
 * `get_mlt1_` agrees with this expression to 9.6e-5 h — 0.34 s of local time — once the dipole
 * epochs are matched (see the file brief). That is the measurement that fixed the definition; no
 * Fortran was read.
 *
 * @param geo the position in GEO, Earth radii.
 * @param r the epoch's rotations, built once by @ref rotations_at.
 * @return MLT in hours, folded into `[0, 24)`. @ref Status::DomainError for a non-finite input
 *         component, and for a point ON the dipole axis, where the SM meridian — and therefore MLT
 *         — is not defined; the value returned there is the 12.0 that `atan2(0, 0)` produces, so a
 *         caller that ignores the status still gets a number in range rather than a NaN.
 * @complexity O(1) — one 3×3 product and one `atan2`.
 * @alloc none.
 * @test IrbemApi.MltIsTheSolarMagneticClockAngle
 * @test IrbemApi.MltMatchesTheReference
 */
[[nodiscard]] inline Result<double> get_mlt(const Position<Frame::GEO>& geo, const Rotations& r) {
    // The INPUT is what is checked, not the SM components derived from it: a rotation is
    // orthonormal, so it maps finite to finite, and a NaN anywhere in the input poisons every SM
    // component through the matrix product (`0 * NaN` is NaN) — checking the output would collapse
    // three distinguishable caller errors into one unreachable test.
    if (!std::isfinite(geo.v[0]) || !std::isfinite(geo.v[1]) || !std::isfinite(geo.v[2])) {
        return {Status::DomainError, 0.0};
    }
    const Position<Frame::SM> sm = transform<Frame::SM>(geo, r);
    if (sm.v[0] == 0.0 && sm.v[1] == 0.0) return {Status::DomainError, 12.0};
    // `atan2` has range [-pi, pi] (it reaches -pi exactly, for a negative-zero ordinate — and a
    // negative zero DOES survive the rotation, since the product sums three negative zeros). So the
    // raw hour angle is monotone in the SM azimuth over [0, 24], and only the ENDPOINT needs
    // folding: midnight arrives from one side as 24 and has to leave as 0.
    //
    // There is deliberately no LOWER fold, and that rests on an exact fact rather than on luck:
    // `(-pi) * (12/pi)` in IEEE double is exactly -12, so the -pi end evaluates to exactly +0.0.
    // Both operands are fixed doubles and the multiply is correctly rounded, so this holds on every
    // IEEE-754 platform; IrbemApi.MltIsTheSolarMagneticClockAngle asserts it directly. A branch that
    // cannot run is a branch that cannot be tested, so it is not written.
    const double raw = 12.0 + (std::atan2(sm.v[1], sm.v[0]) * hours_per_radian);
    return {Status::Ok, raw >= 24.0 ? raw - 24.0 : raw};
}

/**
 * Magnetic local time for a whole ephemeris — the `ntime` form of @ref get_mlt.
 *
 * @param geo the positions in GEO, Earth radii.
 * @param rotations either ONE rotation set shared by every point, or one per point. The first is
 *        the common case and is why this form exists: a satellite pass is one epoch's worth of
 *        geometry over thousands of samples, and rebuilding it per sample is ~14 transcendental
 *        calls against the nine multiplies of the transform itself.
 * @param mlt receives one MLT per input, hours in `[0, 24)`; same length as @p geo.
 * @param statuses receives each point's status; same length as @p geo.
 * @return @ref Status::Ok when every point converted, the first failing point's status otherwise,
 *         and @ref Status::DomainError for a length mismatch or a @p rotations span that is neither
 *         1 nor `geo.size()` long. The value is whether a device serviced the call: always `false`,
 *         and see the file brief's batch-lane table for why that is a measurement rather than a
 *         missing feature.
 * @complexity O(n).
 * @alloc none — everything is written through the caller's spans.
 * @test IrbemApi.MltBatchAgreesWithTheScalarLane
 * @test IrbemApi.BatchRejectsMismatchedSpans
 */
[[nodiscard]] inline Result<bool> get_mlt_vec(std::span<const Position<Frame::GEO>> geo,
                                              std::span<const Rotations> rotations,
                                              std::span<double> mlt, std::span<Status> statuses) {
    const std::size_t n = geo.size();
    if (mlt.size() != n || statuses.size() != n) return {Status::DomainError, false};
    if (rotations.size() != 1 && rotations.size() != n) return {Status::DomainError, false};
    const bool shared = rotations.size() == 1;
    Status worst = Status::Ok;
    for (std::size_t i = 0; i < n; ++i) {
        const Result<double> one = get_mlt(geo[i], rotations[shared ? 0 : i]);
        mlt[i] = one.value;
        statuses[i] = one.status;
        worst = first_failure(worst, one.status);
    }
    return {worst, false};
}

// ---- the runtime frame-pair dispatcher ---------------------------------------------------------

namespace detail {

/**
 * Carry a raw component triple from frame @p f into GEO — the first half of @ref coord_trans.
 *
 * Every `sysaxes` frame has a DIRECT relationship with GEO: the five Cartesian ones each have a
 * rotation @ref Rotations stores (or its exact transpose), and the three angular ones are closed
 * forms in @ref coords_geodetic.hpp. So the dispatcher is a hub-and-spoke through GEO rather than a
 * 9×9 table — 18 legs instead of 81, and every pair composed of two legs that are each tested on
 * their own.
 *
 * @param f the source frame.
 * @param v the components, in whatever units @p f means; see @ref kind_of.
 * @param r the epoch's rotations. Unread for the angular frames, whose conversions carry no epoch.
 * @return the same point in GEO, Earth radii, or @ref Status::DomainError for a heliospheric frame
 *         — those need a @ref HelioGeometry rather than a @ref Rotations, and are unreachable from
 *         @ref coord_trans, which admits only the `sysaxes` codes @ref frame_from_sysaxes accepts.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.CoordTransHubRejectsHeliosphericFrames
 */
[[nodiscard]] inline Result<fixarray::vec3d> to_geo(Frame f, const fixarray::vec3d& v,
                                                    const Rotations& r) {
    switch (f) {
        case Frame::GDZ: return {Status::Ok, gdz_to_geo(Position<Frame::GDZ>{v}).v};
        case Frame::GEO: return {Status::Ok, v};
        case Frame::GSM: return {Status::Ok, transform<Frame::GEO>(Position<Frame::GSM>{v}, r).v};
        case Frame::GSE: return {Status::Ok, transform<Frame::GEO>(Position<Frame::GSE>{v}, r).v};
        case Frame::SM: return {Status::Ok, transform<Frame::GEO>(Position<Frame::SM>{v}, r).v};
        case Frame::GEI: return {Status::Ok, transform<Frame::GEO>(Position<Frame::GEI>{v}, r).v};
        case Frame::MAG: return {Status::Ok, transform<Frame::GEO>(Position<Frame::MAG>{v}, r).v};
        case Frame::SPH: return {Status::Ok, sph_to_car(Position<Frame::SPH>{v}).v};
        case Frame::RLL: return {Status::Ok, gdz_to_geo(rll_to_gdz(Position<Frame::RLL>{v})).v};
        case Frame::HEE:
        case Frame::HAE:
        case Frame::HEEQ: break;
    }
    return {Status::DomainError, fixarray::vec3d{}};
}

/**
 * GEO to RLL — the one leg of the hub that @ref coords_geodetic.hpp does not already name.
 *
 * It supplies `rll_to_gdz` and not its converse, because the converse is a composition rather than
 * a solve. RLL is `(geocentric radius, GEODETIC latitude, east longitude)`, so the radius is
 * `|x_GEO|` directly and the other two components are GDZ's, unchanged. That `|x_GEO|` comes out
 * exact rather than approximate is what `rll_to_gdz` solves its quadratic FOR — it picks the
 * altitude that puts the point at the requested geocentric radius — so the two legs invert each
 * other by construction, which `IrbemApi.CoordTransCoversEveryFramePair` checks.
 *
 * @param v the point in GEO, Earth radii.
 * @return the RLL components: radius in Earth radii, geodetic latitude and east longitude in
 *         degrees.
 * @complexity O(1) — one closed-form geodetic inversion and one norm.
 * @alloc none.
 * @test IrbemApi.CoordTransCoversEveryFramePair
 */
[[nodiscard]] inline fixarray::vec3d geo_to_rll(const fixarray::vec3d& v) {
    const Position<Frame::GDZ> gdz = geo_to_gdz(Position<Frame::GEO>{v});
    return fixarray::vec3d{fixarray::norm(v), gdz.latitude(), gdz.longitude()};
}

/**
 * GEO to SPH in the REFERENCE's longitude convention — east longitude folded into `[0, 360)`.
 *
 * This is the one place @ref coord_trans deliberately departs from the typed @ref car2sph beside
 * it, and it is a convention, not arithmetic: the two answers name the same meridian. It exists
 * because the reference's own two angular conventions disagree with each other, which is not
 * something a caller can be expected to guess.
 *
 * Measured on `coord_trans_vec1_` at 2015-180 12:00 UT: GEO `(1.5, -1.5, 4)` comes back as
 * `(4.5276925690687087, 62.06164727039765, 315)` for `sysaxesOUT = 7` (SPH) but
 * `(4.527692569…, 62.09681281…, -45)` for `sysaxesOUT = 8` (RLL) and the same `-45` for
 * `sysaxesOUT = 0` (GDZ). So the reference wraps SPH into `[0, 360)` and leaves GDZ and RLL in
 * `(-180, 180]`. Every point with a negative GEO y — half the sky — differs by exactly 360°
 * between the two spellings, which is why a golden set drawn from a single positive-longitude
 * point cannot see it. `IrbemApi.CoordTransUsesTheReferenceLongitudeConventions` pins both halves.
 *
 * The fold is `< 0`, not `<= 0`: the reference returns `-0` for GEO `(1, -0, 0)`, measured, and a
 * negative zero is not less than zero, so the same expression reproduces that too.
 *
 * @param v the point in GEO, Earth radii.
 * @return the SPH components: radius in Earth radii, geocentric latitude in degrees, east
 *         longitude in degrees folded into `[0, 360)`.
 * @complexity O(1) — one Cartesian-to-spherical conversion and one compare.
 * @alloc none.
 * @test IrbemApi.CoordTransUsesTheReferenceLongitudeConventions
 */
[[nodiscard]] inline fixarray::vec3d sph_in_irbem_longitude(const fixarray::vec3d& v) {
    fixarray::vec3d sph = car_to_sph(Position<Frame::GEO>{v}).v;
    if (sph[2] < 0.0) sph[2] += 360.0;
    return sph;
}

/**
 * Carry a GEO point out into frame @p f — the second half of @ref coord_trans, and the exact
 * inverse of @ref to_geo leg for leg.
 *
 * @param f the destination frame.
 * @param v the point in GEO, Earth radii.
 * @param r the epoch's rotations. Unread for the angular frames.
 * @return the components in @p f, or @ref Status::DomainError for a heliospheric frame, as in
 *         @ref to_geo. The SPH leg goes through @ref sph_in_irbem_longitude rather than
 *         @ref car_to_sph directly, so its longitude carries the reference's `[0, 360)` convention;
 *         @ref to_geo needs no counterpart, since @ref sph_to_car is periodic and accepts either
 *         spelling.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.CoordTransHubRejectsHeliosphericFrames
 * @test IrbemApi.CoordTransUsesTheReferenceLongitudeConventions
 */
[[nodiscard]] inline Result<fixarray::vec3d> from_geo(Frame f, const fixarray::vec3d& v,
                                                      const Rotations& r) {
    switch (f) {
        case Frame::GDZ: return {Status::Ok, geo_to_gdz(Position<Frame::GEO>{v}).v};
        case Frame::GEO: return {Status::Ok, v};
        case Frame::GSM: return {Status::Ok, transform<Frame::GSM>(Position<Frame::GEO>{v}, r).v};
        case Frame::GSE: return {Status::Ok, transform<Frame::GSE>(Position<Frame::GEO>{v}, r).v};
        case Frame::SM: return {Status::Ok, transform<Frame::SM>(Position<Frame::GEO>{v}, r).v};
        case Frame::GEI: return {Status::Ok, transform<Frame::GEI>(Position<Frame::GEO>{v}, r).v};
        case Frame::MAG: return {Status::Ok, transform<Frame::MAG>(Position<Frame::GEO>{v}, r).v};
        case Frame::SPH: return {Status::Ok, sph_in_irbem_longitude(v)};
        case Frame::RLL: return {Status::Ok, geo_to_rll(v)};
        case Frame::HEE:
        case Frame::HAE:
        case Frame::HEEQ: break;
    }
    return {Status::DomainError, fixarray::vec3d{}};
}

}  // namespace detail

/**
 * One point through the runtime frame-pair dispatcher — the scalar half of `COORD_TRANS_VEC`.
 *
 * This is the ONE routine in the module where the frames are runtime values rather than template
 * arguments, which is exactly why it is written as a hub through GEO over @ref frame_from_sysaxes
 * rather than as a switch that picks a default. An unrecognised code is reported; it never silently
 * becomes GEO. (The reference prints `sysaxesOUT out of range !` to stdout and returns `baddata` —
 * measured — which a caller redirecting stdout will never see.)
 *
 * @param sysaxes_in the input frame's IRBEM code, `0..8`.
 * @param sysaxes_out the output frame's IRBEM code, `0..8`.
 * @param in the input components, in whatever units @p sysaxes_in means; see @ref kind_of.
 * @param r the epoch's rotations, built once by @ref rotations_at.
 * @return the components in the output frame. @ref Status::DomainError when either code is outside
 *         `0..8` (see the file brief on codes 9–14) or when any input component is not finite.
 * @note SPH output carries east longitude in `[0, 360)`, which is the reference's convention for
 *       `sysaxes = 7` and NOT the `(-180, 180]` the typed @ref car2sph returns — see
 *       @ref detail::sph_in_irbem_longitude for the measurement. GDZ and RLL output keep
 *       `(-180, 180]`, because that is what the reference returns for codes 0 and 8. On the polar
 *       axis this module's GDZ altitude is right and the reference's is not; the file brief has
 *       the numbers.
 * @complexity O(1) — at most two 3×3 products plus one closed-form geodetic conversion per side.
 * @alloc none.
 * @test IrbemApi.CoordTransCoversEveryFramePair
 * @test IrbemApi.CoordTransReportsBadSysaxes
 * @test IrbemApi.CoordTransMatchesTheReference
 * @test IrbemApi.CoordTransUsesTheReferenceLongitudeConventions
 * @test IrbemApi.CoordTransIsRightWhereTheReferenceIsWrongOnThePolarAxis
 */
[[nodiscard]] inline Result<fixarray::vec3d> coord_trans(int sysaxes_in, int sysaxes_out,
                                                         const fixarray::vec3d& in,
                                                         const Rotations& r) {
    const std::optional<Frame> from = frame_from_sysaxes(sysaxes_in);
    const std::optional<Frame> to = frame_from_sysaxes(sysaxes_out);
    if (!from.has_value() || !to.has_value()) return {Status::DomainError, fixarray::vec3d{}};
    if (!std::isfinite(in[0]) || !std::isfinite(in[1]) || !std::isfinite(in[2])) {
        return {Status::DomainError, fixarray::vec3d{}};
    }
    // Same frame in and out is the identity, not a round trip through GEO. Routing it through the
    // hub would cost two conversions and return something a bit or two away from what came in; a
    // caller normalising a mixed-frame ephemeris to a common frame passes many such points.
    if (*from == *to) return {Status::Ok, in};
    const Result<fixarray::vec3d> geo = detail::to_geo(*from, in, r);
    if (!geo.ok()) return geo;  // unreachable via frame_from_sysaxes; kept so the hub is total
    return detail::from_geo(*to, geo.value, r);
}

/**
 * `COORD_TRANS_VEC` — a whole ephemeris through one frame pair.
 *
 * @param sysaxes_in the input frame's IRBEM code, `0..8`.
 * @param sysaxes_out the output frame's IRBEM code, `0..8`.
 * @param in the input components, one triple per point.
 * @param out receives the converted components; same length as @p in.
 * @param rotations either ONE rotation set shared by every point, or one per point — IRBEM takes a
 *        per-point `iyr`/`idoy`/`secs`, and both shapes are real: an ephemeris resampled to a
 *        single epoch wants the first, a multi-year survey the second.
 * @param statuses receives each point's status; same length as @p in.
 * @return @ref Status::Ok when every point converted, the first failing point's status otherwise,
 *         and @ref Status::DomainError for a bad frame code or a length mismatch — in which case
 *         nothing is written at all, since the fault is in the call and not in the data. The value
 *         is whether a device serviced the call: always `false`; see the file brief's batch table.
 * @complexity O(n).
 * @alloc none.
 * @test IrbemApi.CoordTransVecAgreesWithTheScalarLane
 * @test IrbemApi.CoordTransVecTakesPerPointEpochs
 * @test IrbemApi.BatchRejectsMismatchedSpans
 */
[[nodiscard]] inline Result<bool> coord_trans_vec(int sysaxes_in, int sysaxes_out,
                                                  std::span<const fixarray::vec3d> in,
                                                  std::span<fixarray::vec3d> out,
                                                  std::span<const Rotations> rotations,
                                                  std::span<Status> statuses) {
    const std::size_t n = in.size();
    if (out.size() != n || statuses.size() != n) return {Status::DomainError, false};
    if (rotations.size() != 1 && rotations.size() != n) return {Status::DomainError, false};
    if (!frame_from_sysaxes(sysaxes_in).has_value() ||
        !frame_from_sysaxes(sysaxes_out).has_value()) {
        return {Status::DomainError, false};
    }
    const bool shared = rotations.size() == 1;
    Status worst = Status::Ok;
    for (std::size_t i = 0; i < n; ++i) {
        const Result<fixarray::vec3d> one =
            coord_trans(sysaxes_in, sysaxes_out, in[i], rotations[shared ? 0 : i]);
        out[i] = one.value;
        statuses[i] = one.status;
        worst = first_failure(worst, one.status);
    }
    return {worst, false};
}

/**
 * The typed `ntime` transform: a whole ephemeris through ONE compile-time frame pair.
 *
 * @ref coord_trans_vec is the runtime-code lane and pays a switch per point per side; this is the
 * lane for code that knows its frames, and it hoists the 3×3 out of the loop when the epoch is
 * shared, leaving nine multiplies and six adds per point with no branch and no dispatch. Both
 * exist because both callers exist — a generic converter driven by a configuration file, and an
 * inner loop that has known it was going GEO→GSM since it was written.
 *
 * @tparam To the destination frame.
 * @tparam V the tagged vector template — @ref Position or @ref FieldVector; deduced.
 * @tparam From the source frame; deduced.
 * @param in the input values.
 * @param out receives the transformed values; same length as @p in.
 * @param rotations either ONE rotation set shared by every point, or one per point.
 * @return @ref Status::Ok, or @ref Status::DomainError for a length mismatch or a @p rotations span
 *         that is neither 1 nor `in.size()` long. The value is whether a device serviced the call:
 *         always `false`; see the file brief's batch table for the measurement behind that.
 * @complexity O(n) — one 3×3 product per point, and ONE matrix build per call when the epoch is
 *             shared rather than one per point.
 * @alloc none.
 * @test IrbemApi.TypedBatchAgreesWithTheScalarTransform
 * @test IrbemApi.BatchRejectsMismatchedSpans
 */
template <Frame To, template <Frame> class V, Frame From>
    requires FrameTaggedVector<V> && RotationAvailable<To, From>
[[nodiscard]] inline Result<bool> transform_vec(std::span<const V<From>> in, std::span<V<To>> out,
                                                std::span<const Rotations> rotations) {
    const std::size_t n = in.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (rotations.size() != 1 && rotations.size() != n) return {Status::DomainError, false};
    if (rotations.size() == 1) {
        const fixarray::mat3d m = rotation_matrix<To, From>(rotations[0]);
        for (std::size_t i = 0; i < n; ++i) out[i] = V<To>{m * in[i].v};
        return {Status::Ok, false};
    }
    for (std::size_t i = 0; i < n; ++i) out[i] = transform<To>(in[i], rotations[i]);
    return {Status::Ok, false};
}

// ---- the MAKE_LSTAR family -----------------------------------------------------------------
// The six numbers IRBEM's flagship routine returns, over `driftshell.hpp`'s `L*` and `lstar.hpp`'s
// invariants. Everything here is `kext = 0` — the internal field alone — because that is what
// `driftshell.hpp` takes. `ext_t89.hpp` exists in this module, but nothing wires an external field
// into the drift-shell root-find yet, so a caller who asks for one is TOLD, not quietly given the
// internal-field answer under an external field's name. That distinction matters more here than
// anywhere else in this header: an internal-only L* at geosynchronous during a storm is not a
// slightly worse T89 L*, it is a different number.

/**
 * The six per-point outputs of IRBEM's `MAKE_LSTAR`, in one record.
 *
 * IRBEM returns them as six parallel arrays. They are one array of records here because the six are
 * computed together from one drift shell, and six spans a caller has to keep in step is six chances
 * to hand one of them the wrong length.
 */
struct MagneticCoordinates {
    double lm = baddata;      ///< McIlwain's L, Earth radii — IRBEM's `Lm`.
    double lstar = baddata;   ///< Roederer's L\*, Earth radii — IRBEM's `Lstar`.
    double blocal = baddata;  ///< `|B|` at the point, nT — IRBEM's `Blocal`.
    double bmin = baddata;    ///< The minimum `|B|` on the field line, nT — IRBEM's `Bmin`.
    double xj = baddata;      ///< The second invariant `I`, Earth radii — IRBEM's `XJ`.
    double mlt = baddata;     ///< Magnetic local time, hours — IRBEM's `MLT`.

    /// Two records are equal when every field is.
    /// @param a the left record. @param b the right record.
    /// @return whether all six fields match exactly.
    /// @complexity O(1). @alloc none.
    /// @test IrbemApi.MakeLstarVecAgreesWithTheScalarLane
    friend constexpr bool operator==(const MagneticCoordinates& a, const MagneticCoordinates& b) {
        return a.lm == b.lm && a.lstar == b.lstar && a.blocal == b.blocal && a.bmin == b.bmin &&
               a.xj == b.xj && a.mlt == b.mlt;
    }
};

namespace detail {

/**
 * Whether this module can answer for an IRBEM `kext` key at all.
 *
 * @param kext the external-field key.
 * @return @ref Status::Ok for @ref ExternalModel::None, @ref Status::ParametersMissing for any
 *         other recognised model — this module HAS no external field in the drift-shell path, so
 *         the parameters the model would need are missing in the strongest possible sense — and
 *         @ref Status::DomainError for a key outside `0..14`.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.MakeLstarReportsAnUnsupportedExternalModel
 */
[[nodiscard]] constexpr Status external_model_supported(ExternalModel kext) {
    if (!is_recognised(kext)) return Status::DomainError;
    return kext == ExternalModel::None ? Status::Ok : Status::ParametersMissing;
}

/**
 * Fill one @ref MagneticCoordinates from a computed shell and its point.
 *
 * Split out because @ref make_lstar_vec and @ref make_lstar_shell_splitting index their shells
 * differently — one shell per point against one per (point, pitch angle) — and the assembly of the
 * six outputs from a @ref DriftShell must not be written twice and drift.
 *
 * @param shell the drift shell for this point.
 * @param start the point, GEO, Earth radii.
 * @param rotations the epoch's rotations, for the MLT.
 * @return the six outputs; `mlt` is @ref baddata when the point sits on the dipole axis, where MLT
 *         is undefined, and everything else is still filled.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.MakeLstarMatchesTheOracle
 */
[[nodiscard]] inline MagneticCoordinates coordinates_of(const DriftShell& shell,
                                                        const Position<Frame::GEO>& start,
                                                        const Rotations& rotations) {
    const Result<double> mlt = get_mlt(start, rotations);
    return MagneticCoordinates{shell.lm,   shell.lstar, shell.b_local,
                               shell.b_min, shell.invariant_i, to_baddata(mlt)};
}

}  // namespace detail

/**
 * `MAKE_LSTAR` for a whole ephemeris — **this is the routine to call**, and the reason the batch
 * forms exist at all.
 *
 * One `L*` point is `Nder = 25` independent root-finds of a few field lines each. That is an order
 * of magnitude below the ~512-line crossover `gpu/dispatch.hpp` measures, so a single point runs on
 * the host no matter what hardware is present. Hand over `ntime` points and every stage of the
 * root-find becomes `ntime × Nder` wide, which is where the trace kernel's measured 48.9× lands. A
 * loop calling @ref make_lstar per point cannot be accelerated — that is a property of the problem,
 * not of this implementation.
 *
 * Measured through THIS entry point on an RTX 3070 Ti against its own fp64 host lane
 * (`-O3 -march=native`, IGRF-14, points spread over L = 3…5, `Nder = 25`):
 *
 * | points | host | device | |
 * |---|---|---|---|
 * | 64 | 13.1 ms/point | 20.2 ms/point | device **loses** 0.65× |
 * | 512 | 13.8 ms/point | 0.74 ms/point | device wins **18.6×** |
 *
 * The crossover is the one `gpu/dispatch.hpp` measures, reached here in POINTS rather than lines
 * because each point contributes `Nder` of them. Re-measured on an RTX 3070 Ti with the device flag
 * confirmed `true`: 20.16 ms/point at n = 64 against the host's 13.13 (device loses 0.65×) and
 * 0.739 ms/point at n = 512 against 13.77 (device wins 18.6×). Both lanes are bit-identical
 * run-to-run, and the device lane also serves the call on llvmpipe (8.27 ms/point) and on an Intel
 * iGPU (14.42), so the shape is not NVIDIA-only.
 *
 * The device lane traces in fp32, and over 512 points spread across L = 3…5 the two lanes' `L*`
 * differ by at most 4.2e-3 absolute and `L_m` by 4.0e-3 — inside `docs/ERROR_BUDGET.md`'s 0.01
 * absolute on `L*`, which is what makes the fp32 trace admissible, but ABOVE its 1e-3 on `L_m`, so
 * a caller who needs `L_m` to that budget must pin the host lane with
 * `CHEATAH_SPACE_IRBEM_NO_GPU=1`. (An earlier revision of this comment quoted 2.8e-3 and 2.4e-3;
 * those are one point set's numbers, not the routine's.)
 *
 * ONE CONSEQUENCE THAT BITES: which lane runs depends on the BATCH SIZE, so `make_lstar_vec` over
 * `n` points and `n` separate @ref make_lstar calls are NOT bit-identical on a machine with a
 * device. `n = 1` is 25 field lines and stays on the host; `n = 6` is 150, which crosses
 * `irbem_igrf_f32`'s measured crossover of 128, and the flux quadrature then runs in fp32 on the
 * device. The two answers agree to the budget above, not to the last bit. `==` between the two
 * lanes is only meaningful with the device pinned off, which is what
 * `IrbemApi.MakeLstarVecAgreesWithTheScalarLane` does.
 *
 * The epoch is ONE @ref Rotations for the whole batch, not one per point as IRBEM takes. That is
 * deliberate: the drift shell is organised about the dipole axis, so points at different epochs
 * belong to different shells and cannot share a dispatch. Splitting the batch internally would hand
 * back the speed the batch form exists to provide, silently. Group by epoch in the caller, which is
 * what a caller must do anyway to get the throughput.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param rotations the epoch's rotations, shared by every point.
 * @param starts the points, GEO, Earth radii.
 * @param out receives one record per point; same length as @p starts.
 * @param statuses receives each point's status; same length as @p starts.
 * @param kext IRBEM's external-field key. Only @ref ExternalModel::None is supported; see
 *        @ref detail::external_model_supported and the section comment above for why an
 *        unsupported key is refused rather than silently answered with the internal field.
 * @param pitch_angle_deg the local pitch angle; 90° is IRBEM's `MAKE_LSTAR` convention.
 * @param opt the drift-shell resolution; @ref DriftShellOptions::from_irbem translates IRBEM's
 *        `options(3)`/`options(4)` so a differential comparison is at matched resolution.
 * @return the statuses @ref make_lstar_batch reports, or @ref Status::ParametersMissing /
 *         @ref Status::DomainError for an unsupported or unrecognised @p kext, in which case
 *         nothing is computed and every output stays at @ref baddata. The value is `true` when the
 *         device serviced the traces.
 * @complexity O(points × Nder × (trials + iterations)) traces; concurrent on the device.
 * @alloc as @ref make_lstar_batch — O(rounds) vectors for the root-find, nothing per trace.
 * @test IrbemApi.MakeLstarVecAgreesWithTheScalarLane
 * @test IrbemApi.MakeLstarVecDeviceLaneAgreesWithTheHostLane
 * @test IrbemApi.MakeLstarReportsAnUnsupportedExternalModel
 * @test IrbemApi.MakeLstarMatchesTheOracle
 */
template <int NMAX>
[[nodiscard]] inline Result<bool> make_lstar_vec(const Igrf<NMAX>& model, const Rotations& rotations,
                                                 std::span<const Position<Frame::GEO>> starts,
                                                 std::span<MagneticCoordinates> out,
                                                 std::span<Status> statuses,
                                                 ExternalModel kext = ExternalModel::None,
                                                 double pitch_angle_deg = 90.0,
                                                 const DriftShellOptions& opt = {}) {
    const std::size_t n = starts.size();
    if (out.size() != n || statuses.size() != n) return {Status::DomainError, false};
    const Status supported = detail::external_model_supported(kext);
    if (!is_ok(supported)) {
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = MagneticCoordinates{};
            statuses[i] = supported;
        }
        return {supported, false};
    }
    if (n == 0) return {Status::Ok, false};

    std::vector<DriftShell> shells(n);
    std::vector<double> pitch(n, pitch_angle_deg);
    const Result<bool> r =
        make_lstar_batch(model, rotations, starts, pitch, shells, statuses, opt);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = detail::coordinates_of(shells[i], starts[i], rotations);
    }
    return r;
}

/**
 * `MAKE_LSTAR` for one point — the reference lane, not the fast one.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param rotations the epoch's rotations.
 * @param start the point, GEO, Earth radii.
 * @param kext IRBEM's external-field key; only @ref ExternalModel::None is supported.
 * @param pitch_angle_deg the local pitch angle; 90° is IRBEM's `MAKE_LSTAR` convention.
 * @param opt the drift-shell resolution.
 * @return the six outputs and the point's status. Every output is @ref baddata when @p kext is not
 *         supported; otherwise the record is filled even on a failed shell, because `L_m`, `I` and
 *         the fields are what make the failure diagnosable.
 * @complexity O(Nder × (trials + iterations)) traces — ~10⁵ field evaluations, ~15 ms on the host.
 * @alloc as @ref make_lstar_vec at `n = 1`.
 * @test IrbemApi.MakeLstarMatchesTheOracle
 * @test IrbemApi.MakeLstarVecAgreesWithTheScalarLane
 */
template <int NMAX>
[[nodiscard]] inline Result<MagneticCoordinates> make_lstar(const Igrf<NMAX>& model,
                                                            const Rotations& rotations,
                                                            const Position<Frame::GEO>& start,
                                                            ExternalModel kext = ExternalModel::None,
                                                            double pitch_angle_deg = 90.0,
                                                            const DriftShellOptions& opt = {}) {
    const std::array<Position<Frame::GEO>, 1> starts{start};
    std::array<MagneticCoordinates, 1> out{};
    std::array<Status, 1> st{Status::Ok};
    const Result<bool> r = make_lstar_vec(model, rotations, starts, out, st, kext, pitch_angle_deg, opt);
    return {r.status == Status::DomainError ? Status::DomainError : st[0], out[0]};
}

/**
 * `MAKE_LSTAR_SHELL_SPLITTING` — the same six outputs for every (point, pitch angle) pair.
 *
 * Shell splitting is the phenomenon that particles of different pitch angle at the SAME point drift
 * on DIFFERENT shells once the field is not a dipole, so `L*` is a function of `α` and not just of
 * position. Computationally it is the batch form's best case: `ntime × Nipa` points that all share
 * one epoch, so the whole grid is one dispatch of `ntime × Nipa × Nder` field lines.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param rotations the epoch's rotations, shared by every point.
 * @param starts the points, GEO, Earth radii.
 * @param pitch_angles_deg the local pitch angles, applied to every point — IRBEM's `alpha`.
 * @param out receives `starts.size() × pitch_angles_deg.size()` records, POINT-MAJOR: the record
 *        for point `i` and pitch angle `k` is at `i * pitch_angles_deg.size() + k`. IRBEM's
 *        Fortran array is `(Nipa, ntime)`, which is the same layout read the other way round.
 * @param statuses receives one status per record; same length as @p out.
 * @param kext IRBEM's external-field key; only @ref ExternalModel::None is supported.
 * @param opt the drift-shell resolution.
 * @return the statuses @ref make_lstar_batch reports, or @ref Status::ParametersMissing /
 *         @ref Status::DomainError for an unsupported @p kext. @ref Status::DomainError also for a
 *         length mismatch. The value is `true` when the device serviced the traces.
 * @complexity O(points × angles × Nder × (trials + iterations)) traces; concurrent on the device.
 * @alloc one point vector and one pitch vector of `points × angles`, plus what
 *        @ref make_lstar_batch allocates; nothing per trace.
 * @test IrbemApi.ShellSplittingVariesLstarWithPitchAngle
 * @test IrbemApi.MakeLstarReportsAnUnsupportedExternalModel
 */
template <int NMAX>
[[nodiscard]] inline Result<bool> make_lstar_shell_splitting(
    const Igrf<NMAX>& model, const Rotations& rotations,
    std::span<const Position<Frame::GEO>> starts, std::span<const double> pitch_angles_deg,
    std::span<MagneticCoordinates> out, std::span<Status> statuses,
    ExternalModel kext = ExternalModel::None, const DriftShellOptions& opt = {}) {
    const std::size_t points = starts.size();
    const std::size_t angles = pitch_angles_deg.size();
    const std::size_t total = points * angles;
    if (out.size() != total || statuses.size() != total) return {Status::DomainError, false};
    const Status supported = detail::external_model_supported(kext);
    if (!is_ok(supported)) {
        for (std::size_t i = 0; i < total; ++i) {
            out[i] = MagneticCoordinates{};
            statuses[i] = supported;
        }
        return {supported, false};
    }
    if (total == 0) return {Status::Ok, false};

    // The grid is flattened into ONE batch rather than looped over pitch angles: `ntime × Nipa`
    // shells at one epoch are exactly the shape the device lane wants, and a loop over angles would
    // submit `Nipa` batches of `ntime × Nder` instead of one of `ntime × Nipa × Nder`.
    std::vector<Position<Frame::GEO>> grid;
    std::vector<double> pitch;
    grid.reserve(total);
    pitch.reserve(total);
    for (std::size_t i = 0; i < points; ++i) {
        for (std::size_t k = 0; k < angles; ++k) {
            grid.push_back(starts[i]);
            pitch.push_back(pitch_angles_deg[k]);
        }
    }
    std::vector<DriftShell> shells(total);
    const Result<bool> r = make_lstar_batch(model, rotations, grid, pitch, shells, statuses, opt);
    for (std::size_t i = 0; i < total; ++i) {
        out[i] = detail::coordinates_of(shells[i], grid[i], rotations);
    }
    return r;
}

/**
 * `LSTAR_PHI` — convert between Roederer's `L*` and the third invariant `Φ`.
 *
 * `Φ = 2π k₀ / L*`, the relation `driftshell.hpp` computes `L*` with, and the one the reference
 * documents for `LSTAR_PHI` as `Φ = 2π B₀/L*`. The reference takes `iyear` and `idoy`, and MEASURED
 * it uses that epoch's own dipole moment rather than a frozen one: `lstar_phi1_(whichinv = 1)` at
 * `L* = 4` returns 48599.478358059241 nT·Re² in 1965 and 46590.805305975191 in 2029, whose implied
 * `B₀` — 30939.3888 and 29660.6279 nT — is @ref dipole_moment of `Igrf<>::at(year + 0.5)` to the
 * last bit. So this routine reproduces the reference EXACTLY, 0 to 1.5e-16 relative over 1965,
 * 1990, 2015 and 2029 in both directions, and there is no convention to reconcile. Taking `k₀` from
 * the epoch is therefore agreement, not the divergence an earlier revision of this comment claimed.
 * @ref mcilwain_l documents why a frozen moment would be wrong: it has fallen several percent over
 * the last century, and freezing it shows up as a constant offset at every shell rather than as an
 * obvious error.
 *
 * @param whichinv IRBEM's direction code: `1` converts `L*` to `Φ`, `2` converts `Φ` to `L*`.
 * @param value the `L*` in Earth radii (`whichinv = 1`) or the `Φ` in nT·Re² (`whichinv = 2`).
 * @param dipole_moment_nt the epoch's dipole moment `k₀` in nT; use @ref dipole_moment.
 * @return the converted quantity, or @ref Status::DomainError for a `whichinv` outside `1..2`, a
 *         non-positive or non-finite @p value, or a non-positive moment. The relation is a
 *         reciprocal, so zero has no image in either direction and is refused rather than turned
 *         into an infinity.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.LstarPhiInvertsItself
 * @test IrbemApi.LstarPhiMatchesTheReference
 */
[[nodiscard]] inline Result<double> lstar_phi(int whichinv, double value,
                                              double dipole_moment_nt) {
    if (whichinv < 1 || whichinv > 2) return {Status::DomainError, 0.0};
    if (!(value > 0.0) || !std::isfinite(value)) return {Status::DomainError, 0.0};
    if (!(dipole_moment_nt > 0.0) || !std::isfinite(dipole_moment_nt)) {
        return {Status::DomainError, 0.0};
    }
    // The map is its own inverse up to the constant, so ONE expression serves both directions —
    // which is also why the two cannot disagree.
    return {Status::Ok, 2.0 * std::numbers::pi * dipole_moment_nt / value};
}

// ---- date and time -----------------------------------------------------------------------------

/// Day of year from a calendar date, leap-year correct.
/// @param year the year. @param month 1-12. @param day 1-based day of month.
/// @return the day of year, 1-based. @complexity O(1). @alloc none.
/// @test IrbemApi.DateRoundTrips
constexpr int get_doy(int year, int month, int day) { return day_of_year(year, month, day); }

/// The Julian Day Number of a calendar date — IRBEM's `JULDAY`.
/// @param year the year. @param month 1-12. @param day 1-based day of month.
/// @return the Julian Day Number, which begins at NOON of the given date.
/// @complexity O(1). @alloc none.
/// @test IrbemApi.DateRoundTrips
constexpr std::int64_t julday(int year, int month, int day) {
    return julian_day_number(year, month, day);
}

/// The calendar date of a Julian Day Number — IRBEM's `CALDAT`, the inverse of @ref julday.
/// @param jdn the Julian Day Number.
/// @return the year, month and day. @complexity O(1). @alloc none.
/// @test IrbemApi.DateRoundTrips
constexpr CalendarDate caldat(std::int64_t jdn) { return calendar_date(jdn); }

/// A calendar date and time as a decimal year — the epoch argument the field models take.
/// @param year the year. @param month 1-12. @param day 1-based. @param hour 0-23.
/// @param minute 0-59. @param second 0-59.
/// @return the decimal year. @complexity O(1). @alloc none.
/// @test IrbemApi.DateRoundTrips
constexpr double date_and_time2decy(int year, int month, int day, int hour, int minute,
                                    int second) {
    return decimal_year(year, month, day, hour, minute, second);
}

/// A decimal year back to a full broken-down date and time — IRBEM's `DECY2DATE_AND_TIME`.
/// @param decy the decimal year, where `yyyy.0` is January 1st at 00:00.
/// @return the broken-down date and time. @complexity O(1). @alloc none.
/// @test IrbemApi.DateRoundTrips
constexpr DateTime decy2date_and_time(double decy) {
    return date_and_time_from_decimal_year(decy);
}

/// A year, day-of-year and UT as a full broken-down date and time.
/// @param year the year. @param doy 1-based day of year. @param ut_seconds seconds since midnight.
/// @return the broken-down date and time. @complexity O(1). @alloc none.
/// @test IrbemApi.DateRoundTrips
constexpr DateTime doy_and_ut2date_and_time(int year, int doy, double ut_seconds) {
    return date_and_time_from_doy_and_ut(year, doy, ut_seconds);
}

}  // namespace cheatah::space::irbem::api
