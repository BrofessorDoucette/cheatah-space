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
 *    model's fitted envelope. Every routine here returns a @ref Result, and @ref baddata is
 *    produced only by the explicit C-compatibility shims at the bottom.
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
 * | `trace_invariant_batch` (@ref lstar.hpp) | ~9 400 flops/byte | device above 512 lines, 48.9× |
 *
 * A frame rotation is nine multiplies over 48 bytes moved. Measured on an RTX 3070 Ti, the module's
 * own streaming dipole kernel at that intensity LOSES to the host, 0.69×; the trace kernel wins
 * 48.9× at 65 536 lines. So the coordinate lanes here return `false` for "the device serviced this
 * call" and will keep doing so: the shape exists so the answer is askable, and so a future device
 * lane needs no signature change, not because a device is coming for these.
 *
 * ## What is here, and what is not
 *
 * Coverage is honest and partial. Implemented: the whole coordinate-transformation group
 * (`COORD_TRANS_VEC`, the eleven geographic pairs, `SPH2CAR`/`CAR2SPH`, `RLL2GDZ`, and the six
 * heliospheric pairs), `GET_MLT`, the date-and-time group (`JULDAY`, `CALDAT`, `GET_DOY`,
 * `DECY2DATE_AND_TIME`, `DATE_AND_TIME2DECY`, `DOY_AND_UT2DATE_AND_TIME`) and the library-info
 * group (`GET_IRBEM_NTIME_MAX`, `GET_IGRF_VERSION`, `IRBEM_FORTRAN_VERSION`,
 * `IRBEM_FORTRAN_RELEASE`).
 *
 * Absent, because the physics underneath them is not written yet, and a named stub returning
 * `baddata` would be strictly worse than an absence — it would compile at a caller's site and fail
 * at their runtime:
 *
 *  - **The drift-shell family** — `MAKE_LSTAR`, `MAKE_LSTAR_SHELL_SPLITTING`, `LANDI2LSTAR`,
 *    `LANDI2LSTAR_SHELL_SPLITTING`, `EMPIRICALLSTAR`, `LSTAR_PHI`, `DRIFT_SHELL`,
 *    `DRIFT_BOUNCE_ORBIT`. These need `driftshell.hpp`, which does not exist in this tree at the
 *    time of writing; @ref lstar.hpp supplies `B_min`, `I` and `L_m` but not Roederer's `L*`.
 *    When `driftshell.hpp` lands, the `MAKE_LSTAR` forwarders belong here and nowhere else.
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

#include "context.hpp"
#include "coords_geodetic.hpp"
#include "coords_helio.hpp"
#include "coords_rotations.hpp"
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
 * Verified as a black box against the shipped reference: `get_mlt1_` agrees with this expression to
 * 4e-5 h at 2015-180 12:00 UT and to 1.2e-3 h at 2020-090, the residual being the two libraries'
 * different solar-ephemeris and dipole-axis series rather than a different definition. That is the
 * measurement that fixed the definition; no Fortran was read.
 *
 * @param geo the position in GEO, Earth radii.
 * @param r the epoch's rotations, built once by @ref rotations_at.
 * @return MLT in hours, folded into `[0, 24)`. @ref Status::DomainError for a non-finite input, and
 *         for a point ON the dipole axis, where the SM meridian — and therefore MLT — is not
 *         defined; the value returned there is the 12.0 that `atan2(0, 0)` produces, so a caller
 *         that ignores the status still gets a number in range rather than a NaN.
 * @complexity O(1) — one 3×3 product, one `atan2` and one `fmod`.
 * @alloc none.
 * @test IrbemApi.MltIsTheSolarMagneticClockAngle
 * @test IrbemApi.MltMatchesTheReference
 */
[[nodiscard]] inline Result<double> get_mlt(const Position<Frame::GEO>& geo, const Rotations& r) {
    const Position<Frame::SM> sm = transform<Frame::SM>(geo, r);
    if (!std::isfinite(sm.v[0]) || !std::isfinite(sm.v[1])) return {Status::DomainError, 0.0};
    if (sm.v[0] == 0.0 && sm.v[1] == 0.0) return {Status::DomainError, 12.0};
    // fmod rather than a single subtract: atan2 returns (-pi, pi], so the raw hour angle lands in
    // (0, 24], and the endpoint has to come back to 0 rather than out at 24.
    const double raw = 12.0 + (std::atan2(sm.v[1], sm.v[0]) * hours_per_radian);
    const double folded = std::fmod(raw, 24.0);
    return {Status::Ok, folded < 0.0 ? folded + 24.0 : folded};
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
 * Carry a GEO point out into frame @p f — the second half of @ref coord_trans, and the exact
 * inverse of @ref to_geo leg for leg.
 *
 * The RLL leg is the one that is not simply a named inverse: @ref coords_geodetic.hpp supplies
 * `rll_to_gdz` and not its converse, because the converse is a composition rather than a solve.
 * RLL is `(geocentric radius, GEODETIC latitude, east longitude)`, so the radius is `|x_GEO|`
 * directly and the other two components are GDZ's, unchanged. That `|x_GEO|` is exact rather than
 * approximate is what `rll_to_gdz` solves its quadratic FOR — it picks the altitude that puts the
 * point at the requested geocentric radius — so the two legs invert each other by construction.
 *
 * @param f the destination frame.
 * @param v the point in GEO, Earth radii.
 * @param r the epoch's rotations. Unread for the angular frames.
 * @return the components in @p f, or @ref Status::DomainError for a heliospheric frame, as in
 *         @ref to_geo.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemApi.CoordTransHubRejectsHeliosphericFrames
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
        case Frame::SPH: return {Status::Ok, car_to_sph(Position<Frame::GEO>{v}).v};
        case Frame::RLL: {
            const Position<Frame::GDZ> gdz = geo_to_gdz(Position<Frame::GEO>{v});
            return {Status::Ok,
                    fixarray::vec3d{fixarray::norm(v), gdz.latitude(), gdz.longitude()}};
        }
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
 * @complexity O(1) — at most two 3×3 products plus one closed-form geodetic conversion per side.
 * @alloc none.
 * @test IrbemApi.CoordTransCoversEveryFramePair
 * @test IrbemApi.CoordTransReportsBadSysaxes
 * @test IrbemApi.CoordTransMatchesTheReference
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
