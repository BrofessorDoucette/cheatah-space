#pragma once

/**
 * @file ext_mead.hpp
 * @brief space.irbem — Mead & Fairfield (1975), the simplest EXTERNAL magnetospheric field model
 *        in the IRBEM set (`kext = 1`), and the one whose published form the oracle reproduces
 *        exactly.
 *
 * Where @ref ext_t89.hpp is four current systems with a warped sheet, this is a **polynomial**: each
 * GSM-derived component of the external field is a quadratic in position and linear in the
 * geodipole tilt, fitted by least squares to twelve satellite-years of magnetometer data sorted into
 * four Kp groups. There is no loop, no square root, no exponential and no data-dependent branch —
 * about 50 flops per point — which makes it both the cheapest model of the family and the easiest
 * to verify, since every derivative is a coefficient.
 *
 * ## The model, and where each piece comes from
 *
 * Mead & Fairfield, *A quantitative magnetospheric model derived from spacecraft magnetometer
 * data*, J. Geophys. Res. **80**(4):523-534 (1975). The paper fits the difference between the
 * measured field and the internal reference field as a power series in **solar magnetic**
 * coordinates — the frame in which the dipole axis is `z` — with the positions **aberrated** by the
 * mean solar-wind aberration angle of 4 degrees (the apparent flow direction is shifted from the
 * Sun-Earth line by the Earth's orbital motion, `atan(30/430)`), and with the tilt angle `psi`
 * entering **in degrees**. Terms odd in `y` are excluded from `B_x` and `B_z`, and even ones from
 * `B_y`, by the noon-midnight symmetry of the fit. In this file's notation, with `(x_m, y_m, z)`
 * the aberrated solar-magnetic position and `psi` the tilt in degrees:
 *
 *     B_x = a_1 z + a_2 x_m z + psi (a_3 + a_4 x_m + a_5 x_m^2 + a_6 y_m^2 + a_7 z^2)
 *     B_y = b_1 y_m z + psi (b_2 y_m + b_3 x_m y_m)
 *     B_z = c_1 + c_2 x_m + c_3 x_m^2 + c_4 y_m^2 + c_5 z^2 + psi (c_6 z + c_7 x_m z)
 *
 * — seventeen coefficients per Kp group, four groups: `Kp = {0, 0+}`, `{1-, 1, 1+, 2-}`,
 * `{2, 2+, 3-}` and `Kp >= 3` (the paper's Tables 2 and 3; @ref mead_coefficient_sets here). The
 * paper's units are gammas (nanotesla) with distances in Earth radii, which is what this file
 * uses too.
 *
 * ## Two rotations, one of which is not undone — MEASURED, NOT ASSUMED
 *
 * This header takes GSM in and gives GSM out. It rotates GSM to SM about `y` by the tilt, exactly as
 * @ref t89_components does; then rotates the SM **position** about the dipole axis by the 4-degree
 * aberration to get `(x_m, y_m)`; evaluates the three polynomials; and rotates the resulting SM
 * **components** back to GSM by the tilt alone. The aberration is applied to where the field is
 * evaluated and NOT to the components that come out. That is what the paper's fit defines and what
 * the IRBEM oracle computes — it is not a guess, it is the result of the provenance experiment
 * below, and it matters at the 0.5% level (the `y`-odd terms the aberration induces in `B_x` and
 * `B_z` are exactly `tan(4 deg)` times their `x` partners).
 *
 * ## Kp is a BIN, not a number
 *
 * Exactly as T89: four coefficient sets, no interpolation between them. @ref mead_kp_bin does the
 * binning and its thresholds in IRBEM's Kp x 10 scaling are 4, 20 and 30 — between `0+` and `1-`,
 * between `2-` and `2`, and between `3-` and `3`. A caller who blends two sets is inventing a model
 * that was never fitted.
 *
 * ## What this implementation is, and what IRBEM's `kext = 1` is — ORACLE PARITY
 *
 * Both the functional form AND the coefficients are published, and the differential harness
 * [`tools/oracle/mead_diff.cpp`](../../../tools/oracle/mead_diff.cpp) confirms that IRBEM's
 * `kext = 1` is that model and nothing else, by measurement rather than by reading its source:
 *
 * 1. With the oracle's external field isolated as `kext = 1` minus `kext = 0` (the IGRF term
 *    cancels exactly) and a FREE 20-term quadratic-in-position, linear-in-tilt basis fitted per
 *    component, the fit in GSM leaves a 19% residual in its worst component; the same fit in SM
 *    with the tilt angle leaves **7.4e-15** — roundoff. A fit in SM with `sin(psi)` in place of
 *    `psi` leaves 5.0e-3, so the tilt enters as the angle, and no `psi^2` term is needed. That
 *    settles the frame and the tilt variable without an assumption.
 * 2. The ratio of the fitted `y z` to `x z` coefficients of `B_x` is `-0.0699268 = -tan(4.000 deg)`
 *    in every bin, and the fitted `x^2` and `y^2` tilt terms of `B_y` are equal and opposite — the
 *    signature of a rotated `x y` term. Refitting on the seventeen-term basis above with the
 *    position rotated by 4 degrees recovers **every one of the 68 coefficients as a three- or
 *    four-significant-figure decimal to 1e-9 relative** (e.g. `c_1 = -9.41000000042`), i.e. the
 *    table as printed.
 * 3. With those printed values hard-coded — @ref mead_coefficient_sets — the implementation agrees
 *    with the oracle over 1800 scattered points per bin at 1.2-16.7 R_E across six epochs of
 *    tilt (1989, 2003, 2015 x 3, 2022): RMS relative deviation 2.2e-9, 2.2e-9, 2.0e-9, 2.0e-9 for
 *    bins 1-4, worst single-point relative deviation 6.0e-8, 5.0e-8, 7.8e-8, 2.8e-8 — inside the
 *    1e-6 parity budget by more than an order of magnitude; at the three tilts +0.0002, +25.64
 *    and -30.42 degrees separately, worst 2.2e-9, 8.2e-9 and 5.0e-9. The 1e-9 floor is the
 *    oracle's own rounded internal constants, not chased: chasing it would mean reading the
 *    Fortran. @ref IrbemMead.MatchesTheIrbemOracleToParity asserts the budget in the suite.
 *
 * Two checks on THIS side need no oracle. The paper's fit was constrained to be divergence-free,
 * and in the unaberrated frame `div B = 0` for a polynomial of this shape is three linear identities
 * per bin — `a_2 + b_1 + 2 c_5 = 0`, `a_4 + b_2 + c_6 = 0`, `2 a_5 + b_3 + c_7 = 0` — which the
 * printed coefficients satisfy to 1e-4, the rounding of a table given to four decimals
 * (@ref IrbemMead.PublishedCoefficientsAreDivergenceFree, a transcription check with no way to pass
 * by accident). And because every component is a polynomial of degree two, a second-order central
 * difference of it is EXACT rather than `O(h^2)`-accurate: @ref IrbemMead.DivergenceVanishesEverywhere
 * measures the stencil residual at `h = 0.5` and `h = 5e-4` and finds the same number at both —
 * 7.0e-2 nT/R_E for the printed table over a 14 x 9 x 8 R_E box at tilts to 31.5 degrees, and
 * 5.9e-14 (h = 0.5) / 6.4e-11 (h = 5e-4, the roundoff of a 30 nT field differenced over 1e-3 R_E)
 * once the identities are enforced exactly. The 7.0e-2 is dominated not by the table's rounding but
 * by the aberration: rotating the position and not the components leaves `psi y_m sin(4 deg)
 * (2 a_6 - b_3)`, a property of the model as defined and as the oracle computes it, and a fact
 * about it rather than a defect in this file. A wrong derivative anywhere in the evaluator makes
 * the residual grow with the coordinate, which that test also checks by breaking one coefficient
 * and watching the stencil return exactly the predicted `2 delta z_SM`.
 *
 * @note Nothing on a hot path allocates, and nothing but the device lane can throw. The evaluator,
 *       the scalar entry points and the fp32 host lane touch nothing but the caller's spans and one
 *       stack parameter block — @ref IrbemMead.NothingOnTheHeapInTheHotPath counts. The one
 *       exception is @ref mead_field_batch's DEVICE lane, which stages `3N` floats each way and
 *       forwards whatever `gpu::dispatch_batch` throws; its `@alloc` says so.
 */

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

#include "context.hpp"
#include "coords_rotations.hpp"
#include "frames.hpp"
#include "igrf.hpp"
#include "status.hpp"

// The device lane is opt-in by include path, exactly as ext_t89.hpp's is. Without
// cheatah-gpu-linalg every routine here still compiles and runs on the host.
#if __has_include("cheatah_gpu_linalg/context.hpp")
#  include "gpu/dispatch.hpp"
/// 1 when this translation unit can reach the Mead-Fairfield device kernel; 0 when host-only.
#  define CHEATAH_SPACE_IRBEM_MEAD_GPU 1
#else
/// 1 when this translation unit can reach the Mead-Fairfield device kernel; 0 when host-only.
#  define CHEATAH_SPACE_IRBEM_MEAD_GPU 0
#endif

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// The one parameter the paper FIXED: the solar-wind aberration
// -------------------------------------------------------------------------------------------

/// The mean solar-wind aberration angle, degrees, by which the solar-magnetic position is rotated
/// about the dipole axis before the polynomials are evaluated. Mead & Fairfield (1975) prepared
/// their data in aberrated coordinates with this mean value; the oracle experiment in the file
/// brief measures it back as `atan(0.0699268) = 4.000 deg`.
inline constexpr double mead_aberration_deg = 4.0;

/// `sin(4 deg)`, to the last bit of binary64 — asserted against `std::sin` by
/// @ref IrbemMead.AberrationConstantsAreExactTrigonometry, since `std::sin` is not `constexpr`.
inline constexpr double mead_aberration_sin = 0.0697564737441253;

/// `cos(4 deg)`, likewise.
inline constexpr double mead_aberration_cos = 0.9975640502598242;

/// Degrees per radian, `180 / pi`: the tilt enters the model in degrees.
inline constexpr double mead_deg_per_rad = 180.0 / std::numbers::pi;

// -------------------------------------------------------------------------------------------
// The fitted parameters — one set per Kp bin
// -------------------------------------------------------------------------------------------

/// How many coefficients one Kp group carries: seven for `B_x`, three for `B_y`, seven for `B_z`.
inline constexpr std::size_t mead_coefficient_count = 17;

/**
 * One Kp group's coefficients — Mead & Fairfield (1975), one column of Tables 2-3.
 *
 * Named rather than indexed so that the evaluator reads as the three equations in the file brief
 * and each divergence identity names the three coefficients it relates. Templated on the scalar
 * type for the same reason @ref T89Parameters is: the device kernel receives these as `float`, and
 * the host lane it is verified against must round them FIRST and then evaluate, or a disagreement
 * between the lanes cannot be attributed. @ref mead_parameters is the converter.
 *
 * Units: nanotesla, with positions in Earth radii and the tilt in degrees — each coefficient's unit
 * is whatever makes its term come out in nT.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` for the device-mirroring one.
 * @test IrbemMead.PublishedCoefficientsAreDivergenceFree
 */
template <std::floating_point T>
struct MeadParameters {
    T bx_z;    ///< `a_1`: the `z` term of `B_x`.
    T bx_xz;   ///< `a_2`: the `x_m z` term of `B_x`.
    T bx_t;    ///< `a_3`: the tilt-proportional constant of `B_x`.
    T bx_tx;   ///< `a_4`: the `psi x_m` term of `B_x`.
    T bx_txx;  ///< `a_5`: the `psi x_m^2` term of `B_x`.
    T bx_tyy;  ///< `a_6`: the `psi y_m^2` term of `B_x`.
    T bx_tzz;  ///< `a_7`: the `psi z^2` term of `B_x`.
    T by_yz;   ///< `b_1`: the `y_m z` term of `B_y`.
    T by_ty;   ///< `b_2`: the `psi y_m` term of `B_y`.
    T by_txy;  ///< `b_3`: the `psi x_m y_m` term of `B_y`.
    T bz_1;    ///< `c_1`: the constant of `B_z` — the ring-current-plus-tail depression at the origin.
    T bz_x;    ///< `c_2`: the `x_m` term of `B_z` — dayside compression, nightside stretching.
    T bz_xx;   ///< `c_3`: the `x_m^2` term of `B_z`.
    T bz_yy;   ///< `c_4`: the `y_m^2` term of `B_z`.
    T bz_zz;   ///< `c_5`: the `z^2` term of `B_z`.
    T bz_tz;   ///< `c_6`: the `psi z` term of `B_z`.
    T bz_txz;  ///< `c_7`: the `psi x_m z` term of `B_z`.
};

/// The fp64 reference spelling of a coefficient set.
using MeadCoefficients = MeadParameters<double>;

/// How many Kp bins the model has — four, all published.
inline constexpr std::size_t mead_bin_count = 4;

/**
 * Mead & Fairfield (1975) Tables 2-3: the fitted coefficients, indexed by `bin - 1`.
 *
 * The columns are, in order, `Kp = {0, 0+}`, `{1-, 1, 1+, 2-}`, `{2, 2+, 3-}` and `Kp >= 3`.
 *
 * Every number here is a printed table value and every one is checked two ways: the three
 * divergence identities per bin (@ref IrbemMead.PublishedCoefficientsAreDivergenceFree) and the
 * exact recovery from the IRBEM oracle described in the file brief, which reproduces each of these
 * 68 decimals to 1e-9 relative. The monotone deepening of `c_1` from -9.41 nT to -22.9 nT across
 * the columns is the paper's central result and @ref IrbemMead.PublishedCoefficientsAreOrderedByDisturbance
 * pins it.
 *
 * @test IrbemMead.PublishedCoefficientsAreDivergenceFree
 * @test IrbemMead.PublishedCoefficientsAreOrderedByDisturbance
 */
inline constexpr std::array<MeadCoefficients, mead_bin_count> mead_coefficient_sets{{
    // Kp = 0, 0+
    {1.793, -0.0579, 0.298, -0.0257, -0.0003, -0.00147, 0.00105,
     -0.1011, -0.0198, 0.00009,
     -9.41, 1.507, 0.1316, 0.0836, 0.0795, 0.0455, 0.00051},
    // Kp = 1-, 1, 1+, 2-
    {2.179, -0.0703, 0.302, -0.0299, -0.00062, -0.00122, 0.00095,
     -0.1184, -0.0257, -0.00028,
     -11.96, 1.787, 0.1588, 0.0977, 0.0943, 0.0557, 0.00153},
    // Kp = 2, 2+, 3-
    {3.316, -0.0639, 0.430, -0.0325, -0.00044, -0.00127, 0.00045,
     -0.1654, -0.0308, 0.00022,
     -19.88, 2.023, 0.2272, 0.1323, 0.1146, 0.0633, 0.00067},
    // Kp >= 3
    {3.948, -0.0291, 0.517, -0.0386, -0.00104, -0.00129, -0.00114,
     -0.191, -0.0350, 0.00023,
     -22.90, 2.270, 0.2650, 0.1554, 0.1100, 0.0736, 0.00185},
}};

/**
 * The Kp bin Mead & Fairfield (1975) uses, from Kp in IRBEM's OMNI2 scaling.
 *
 * Four groups: `{0, 0+}`, `{1-, 1, 1+, 2-}`, `{2, 2+, 3-}` and `Kp >= 3`. In the Kp x 10 scaling
 * the attainable values are 0, 3, 7, 10, 13, 17, 20, 23, 27, 30, ..., so the group edges are
 * `4` (between `0+ = 3` and `1- = 7`), `20` (`2- = 17` from `2 = 20`) and `30` (`3- = 27` from
 * `3 = 30`). The edges are inclusive on the upper side — `20` is bin 3 and `30` is bin 4 — which
 * the oracle switch points in `tools/oracle/mead_diff.cpp` confirm. A negative or NaN Kp (which
 * @ref check_validity reports separately) yields bin 1 and anything at or above `30` yields bin 4.
 *
 * @param kp_times_ten Kp in IRBEM's slot-1 scaling, i.e. Kp x 10, nominally 0..90.
 * @return the bin, `1..4`.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemMead.KpBinsFollowThePublishedGroups
 */
[[nodiscard]] inline int mead_kp_bin(double kp_times_ten) {
    if (!(kp_times_ten >= 4.0)) return 1;  // the negation also catches NaN and every negative Kp
    if (kp_times_ten < 20.0) return 2;
    if (kp_times_ten < 30.0) return 3;
    return 4;                              // and +infinity
}

/**
 * The coefficients for @p bin, in the scalar type the caller's lane evaluates in.
 *
 * The `float` instantiation is how the host lane is made to run the same arithmetic as the device
 * kernel: rounding the coefficients on the way in, rather than evaluating in `double` and rounding
 * the answer, is what makes a host-vs-device disagreement attributable to the DEVICE.
 *
 * @tparam T the scalar type; `double` or `float`.
 * @param bin the Kp bin, `1..4`. Values outside that range are clamped, because this function is
 *        called after @ref mead_field has already decided what status the answer carries and must
 *        not be able to index out of the table.
 * @return the coefficient set, by value, converted to @p T.
 * @complexity O(1) — 17 conversions, all folded for a compile-time @p bin.
 * @alloc none; the returned object is inline storage.
 * @test IrbemMead.ParametersRoundTripThroughFloat
 * @test IrbemMead.ParametersClampAnOutOfRangeBin
 */
template <std::floating_point T>
[[nodiscard]] constexpr MeadParameters<T> mead_parameters(int bin) {
    int index = bin - 1;
    if (index < 0) index = 0;
    if (index >= static_cast<int>(mead_bin_count)) index = static_cast<int>(mead_bin_count) - 1;
    const MeadCoefficients& s = mead_coefficient_sets[static_cast<std::size_t>(index)];
    return MeadParameters<T>{
        static_cast<T>(s.bx_z),  static_cast<T>(s.bx_xz), static_cast<T>(s.bx_t),
        static_cast<T>(s.bx_tx), static_cast<T>(s.bx_txx), static_cast<T>(s.bx_tyy),
        static_cast<T>(s.bx_tzz), static_cast<T>(s.by_yz), static_cast<T>(s.by_ty),
        static_cast<T>(s.by_txy), static_cast<T>(s.bz_1),  static_cast<T>(s.bz_x),
        static_cast<T>(s.bz_xx), static_cast<T>(s.bz_yy),  static_cast<T>(s.bz_zz),
        static_cast<T>(s.bz_tz), static_cast<T>(s.bz_txz)};
}

// -------------------------------------------------------------------------------------------
// The evaluator
// -------------------------------------------------------------------------------------------

/**
 * The Mead-Fairfield external field at one GSM point, as three components in nanotesla.
 *
 * This is the whole model. The order is: GSM to SM by the tilt (a rotation about `y`, the same one
 * @ref t89_components makes and for the same reason — the paper's frame has the dipole axis as
 * `z`); the SM position to its aberrated twin `(x_m, y_m)` by a rotation about `z` through 4
 * degrees; the three polynomials of the file brief; and the SM components back to GSM. The
 * aberration rotates the POSITION only — see the brief for the measurement that says so.
 *
 * The tilt arrives three times over: as its sine and cosine for the frame rotation, and as the
 * angle in degrees for the polynomials. All three are properties of the epoch, not of the point,
 * so they are paid once per timestamp by the caller rather than recomputed here — there is no
 * trigonometry in this function at all. The three MUST describe the same angle; the entry points
 * below guarantee it by deriving all three from one `tilt_rad`.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device kernel.
 * @param p the coefficients for the Kp bin; see @ref mead_parameters.
 * @param sin_tilt `sin(psi)`, the dipole tilt's sine.
 * @param cos_tilt `cos(psi)`.
 * @param tilt_deg `psi` in DEGREES — the paper's unit for the tilt, and the unit the coefficients
 *        are printed in.
 * @param x the GSM x coordinate, R_E.
 * @param y the GSM y coordinate, R_E.
 * @param z the GSM z coordinate, R_E.
 * @return `{B_x, B_y, B_z}` in GSM, nanotesla.
 * @complexity O(1) — about 50 flops; no loop, no branch, no transcendental.
 * @alloc none.
 * @test IrbemMead.DivergenceVanishesEverywhere
 * @test IrbemMead.ZeroTiltIsMirrorSymmetricAboutTheEquator
 * @test IrbemMead.TheTiltDependenceIsExactlyLinear
 * @test IrbemMead.TheAberrationRotatesTheNoonMidnightPlaneByFourDegrees
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, 3> mead_components(const MeadParameters<T>& p, T sin_tilt,
                                                      T cos_tilt, T tilt_deg, T x, T y, T z) {
    // ---- GSM -> SM. A rotation about y by the tilt; SM's z axis IS the dipole axis. ----------
    const T xs = (x * cos_tilt) - (z * sin_tilt);
    const T ys = y;
    const T zs = (x * sin_tilt) + (z * cos_tilt);

    // ---- the aberrated position: SM rotated about the dipole axis by 4 degrees ---------------
    const T ca = static_cast<T>(mead_aberration_cos);
    const T sa = static_cast<T>(mead_aberration_sin);
    const T xm = (xs * ca) - (ys * sa);
    const T ym = (xs * sa) + (ys * ca);

    // ---- the three polynomials, in SM components ---------------------------------------------
    const T xm2 = xm * xm;
    const T ym2 = ym * ym;
    const T zs2 = zs * zs;
    const T bx_sm = (p.bx_z * zs) + (p.bx_xz * xm * zs) +
                    (tilt_deg * (p.bx_t + (p.bx_tx * xm) + (p.bx_txx * xm2) + (p.bx_tyy * ym2) +
                                 (p.bx_tzz * zs2)));
    const T by_sm = (p.by_yz * ym * zs) + (tilt_deg * ((p.by_ty * ym) + (p.by_txy * xm * ym)));
    const T bz_sm = p.bz_1 + (p.bz_x * xm) + (p.bz_xx * xm2) + (p.bz_yy * ym2) + (p.bz_zz * zs2) +
                    (tilt_deg * ((p.bz_tz * zs) + (p.bz_txz * xm * zs)));

    // ---- SM -> GSM. The inverse of the tilt rotation; a transpose, not a solve. --------------
    return {(bx_sm * cos_tilt) + (bz_sm * sin_tilt), by_sm,
            (-bx_sm * sin_tilt) + (bz_sm * cos_tilt)};
}

/**
 * The Mead-Fairfield external field at one GSM point, in `double` — the reference lane.
 *
 * @param p the position, GSM, in Earth radii.
 * @param sin_tilt `sin(psi)`; @ref HotState::sin_tilt holds it, precomputed per epoch.
 * @param cos_tilt `cos(psi)`.
 * @param tilt_deg `psi` in degrees; `HotState::tilt_rad * mead_deg_per_rad`.
 * @param bin the Kp bin, `1..4`; out-of-range values are clamped by @ref mead_parameters.
 * @return the external field at @p p, GSM, in nanotesla.
 * @complexity O(1); see @ref mead_components.
 * @alloc none.
 * @test IrbemMead.ReferenceLaneMatchesTheComponentForm
 */
[[nodiscard]] inline FieldVector<Frame::GSM> mead_field_at(Position<Frame::GSM> p, double sin_tilt,
                                                           double cos_tilt, double tilt_deg,
                                                           int bin) {
    const std::array<double, 3> b = mead_components<double>(
        mead_parameters<double>(bin), sin_tilt, cos_tilt, tilt_deg, p.v[0], p.v[1], p.v[2]);
    return FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
}

/**
 * The Mead-Fairfield external field, with the model's own verdict on whether to believe it here.
 *
 * The value is **always** returned, including when the status is @ref Status::OutOfValidityRange —
 * `status.hpp`'s standing rule. What is refused outright is input that has no meaning: a
 * non-finite coordinate, tilt or Kp, a point inside the Earth, or a tilt beyond a right angle.
 * The last is not a singularity as it is in @ref t89_field (there is no `tan(psi)` here, and
 * `psi = 90 deg` itself has a value) but a definition: the geodipole tilt is the angle between
 * the dipole axis and the GSM `z` axis and lies in `[-90, 90]` degrees, the same bound
 * @ref make_field_context enforces. Within it, and with a finite radius, the three quadratics
 * cannot overflow — every squared coordinate is bounded by the finite `r^2` and every coefficient
 * by 30 — so no check on the OUTPUT is needed, and none is made.
 *
 * What IS checked against the published envelope, both through `status.hpp` so the rules live in
 * one place: Kp against `0 <= Kp <= 9` (in IRBEM's Kp x 10 scaling) and the position against the
 * paper's `r <= 17 R_E`. The Kp check happens before the binning, so a caller who passes Kp = 12
 * gets the `Kp >= 3` set AND is told the model was never fitted there.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians; positive when the north dipole leans sunward.
 * @param kp_times_ten Kp in IRBEM's `maginput` slot-1 scaling, i.e. Kp x 10, nominally 0..90.
 * @return the field and its caveat. @ref Status::DomainError (with a zero field) for a non-finite
 *         input, a radius inside the Earth, or `|psi| > pi/2`;
 *         @ref Status::OutOfValidityRange for a Kp or a radius outside the published envelope, with
 *         the field still computed; otherwise @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemMead.OutOfRangeKpIsReportedButStillEvaluated
 * @test IrbemMead.TheRadialEnvelopeIsCheckedFromBothSides
 * @test IrbemMead.NonFiniteInputIsADomainError
 * @test IrbemMead.ATiltBeyondARightAngleIsADomainError
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> mead_field(Position<Frame::GSM> p,
                                                                double tilt_rad,
                                                                double kp_times_ten) {
    const FieldVector<Frame::GSM> zero{};
    if (!std::isfinite(p.v[0]) || !std::isfinite(p.v[1]) || !std::isfinite(p.v[2]) ||
        !std::isfinite(tilt_rad) || !std::isfinite(kp_times_ten)) {
        return {Status::DomainError, zero};
    }
    if (std::fabs(tilt_rad) > max_tilt_rad) return {Status::DomainError, zero};

    const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
    const Status where = check_position(ExternalModel::MeadFairfield1975, r, p.v[0]);
    if (where == Status::DomainError) return {Status::DomainError, zero};

    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Kp)] = kp_times_ten;
    const Status drives = check_validity(ExternalModel::MeadFairfield1975, drivers);

    const FieldVector<Frame::GSM> b =
        mead_field_at(p, std::sin(tilt_rad), std::cos(tilt_rad), tilt_rad * mead_deg_per_rad,
                      mead_kp_bin(kp_times_ten));
    return {first_failure(drives, where), b};
}

/**
 * The Mead-Fairfield external field for a whole epoch's worth of state — the production entry.
 *
 * Reads the tilt and Kp straight out of @ref HotState. The sine and cosine the context precomputed
 * are what the three-argument overload would compute again from the angle, so this overload defers
 * to it rather than duplicating the finiteness and envelope logic: the trigonometry is two calls
 * per point, which is nothing beside the checks, and the batch lane below is the one that pays for
 * neither.
 *
 * @param p the position, GSM, in Earth radii.
 * @param ctx the epoch's context; its `hot()` block carries `psi` and Kp.
 * @return the field and its caveat, exactly as the three-argument overload.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemMead.ContextOverloadAgreesWithTheExplicitOne
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> mead_field(Position<Frame::GSM> p,
                                                                const FieldContext& ctx) {
    return mead_field(p, ctx.hot().tilt_rad, ctx.hot().kp);
}

// -------------------------------------------------------------------------------------------
// The batch lanes
// -------------------------------------------------------------------------------------------

/**
 * The Mead-Fairfield field over a whole batch, on the CPU, in `float`.
 *
 * The host twin of `irbem_mead_f32`: the same expressions, in the same order, in the same
 * precision, from coefficients rounded to `float` FIRST — what makes a disagreement between the
 * two lanes attributable to the device. It is also the lane a machine with no GPU actually runs.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, GSM, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, nanotesla; overwritten in full.
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`.
 * @param tilt_deg `psi` in degrees.
 * @param bin the Kp bin, `1..4`.
 * @return `false` when @p pos is not a whole number of points or @p out is a different length, in
 *         which case nothing is written; `true` otherwise.
 * @complexity O(N).
 * @alloc none: the loop is over caller-provided spans and one stack parameter set.
 * @test IrbemMead.HostFloatLaneTracksTheReferenceLane
 * @test IrbemMead.HostFloatLaneRejectsMismatchedSpans
 */
[[nodiscard]] inline bool mead_field_host(std::span<const float> pos, std::span<float> out,
                                          float sin_tilt, float cos_tilt, float tilt_deg,
                                          int bin) {
    if (pos.size() % 3 != 0 || out.size() != pos.size()) return false;
    const MeadParameters<float> p = mead_parameters<float>(bin);
    const std::size_t n = pos.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const std::array<float, 3> b =
            mead_components<float>(p, sin_tilt, cos_tilt, tilt_deg, pos[(3 * i) + 0],
                                   pos[(3 * i) + 1], pos[(3 * i) + 2]);
        out[(3 * i) + 0] = b[0];
        out[(3 * i) + 1] = b[1];
        out[(3 * i) + 2] = b[2];
    }
    return true;
}

/// How many `float` scalars the device kernel's parameter buffer holds: `sin(psi)`, `cos(psi)`,
/// `psi` in degrees, then the seventeen coefficients. Asserted against the kernel registry.
inline constexpr std::size_t mead_param_count = 3 + mead_coefficient_count;

/**
 * Pack the epoch's tilt and a Kp bin's coefficients into the kernel's parameter buffer.
 *
 * The layout is the kernel's ABI and is stated in exactly two places — here and the comment above
 * `irbem_mead_f32` in `irbem.slang`: `[0] sin psi, [1] cos psi, [2] psi in degrees, [3..9] a_1..a_7,
 * [10..12] b_1..b_3, [13..19] c_1..c_7`. A test evaluates both lanes on the same points, which is
 * what actually keeps the two statements in step.
 *
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`.
 * @param tilt_deg `psi` in degrees.
 * @param bin the Kp bin, `1..4`.
 * @return the parameter block, `mead_param_count` floats, by value.
 * @complexity O(1).
 * @alloc none — the block is the returned object's own inline array.
 * @test IrbemMead.ParameterBlockCarriesTheTiltThenTheCoefficients
 */
[[nodiscard]] inline std::array<float, mead_param_count> mead_param_block(float sin_tilt,
                                                                          float cos_tilt,
                                                                          float tilt_deg, int bin) {
    const MeadParameters<float> p = mead_parameters<float>(bin);
    return {sin_tilt, cos_tilt, tilt_deg,
            p.bx_z,   p.bx_xz,  p.bx_t,  p.bx_tx, p.bx_txx, p.bx_tyy, p.bx_tzz,
            p.by_yz,  p.by_ty,  p.by_txy,
            p.bz_1,   p.bz_x,   p.bz_xx, p.bz_yy, p.bz_zz,  p.bz_tz,  p.bz_txz};
}

/**
 * The batch's position caveat, accumulated one point at a time.
 *
 * The same device as @ref T89PositionFold and for the same reason: one @ref Status for N points can
 * only be the worst of them, and computing that honestly must not cost a second pass over 100 MB
 * of positions. The Mead-Fairfield envelope is a radius band, monotone, so the smallest and largest
 * `r^2` decide the batch and no per-point `sqrt` is paid.
 *
 * @test IrbemMead.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
struct MeadPositionFold {
    /// The smallest `r^2` seen, R_E^2; `+inf` until the first point.
    double r2_lo = std::numeric_limits<double>::infinity();
    /// The largest `r^2` seen, R_E^2; zero until the first point.
    double r2_hi = 0.0;
    /// The smallest GSM `x` seen, R_E; `+inf` until the first point. The envelope publishes no
    /// `x` bound, but @ref check_position reads one, so it is folded for the same two lookups.
    double x_lo = std::numeric_limits<double>::infinity();
    /// False once any point has had a non-finite coordinate — a NaN radius compares false against
    /// everything and would otherwise slip through both extremes.
    bool finite = true;

    /**
     * Fold one position in.
     * @param p the position, GSM, in Earth radii.
     * @complexity O(1) — one fused radius, three comparisons, no `sqrt`.
     * @alloc none.
     * @test IrbemMead.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    constexpr void add(const Position<Frame::GSM>& p) {
        const double r2 = (p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]);
        finite = finite && std::isfinite(r2);
        r2_lo = r2 < r2_lo ? r2 : r2_lo;
        r2_hi = r2 > r2_hi ? r2 : r2_hi;
        x_lo = p.v[0] < x_lo ? p.v[0] : x_lo;
    }

    /**
     * What the batch's positions say about the model's envelope.
     * @return @ref Status::DomainError when any point is not finite or is inside the Earth,
     *         @ref Status::OutOfValidityRange when any point is beyond the published `r <= 17 R_E`,
     *         otherwise @ref Status::Ok.
     * @complexity O(1) — two square roots and two envelope lookups for the whole batch.
     * @alloc none.
     * @test IrbemMead.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    [[nodiscard]] Status verdict() const {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double lo = finite ? std::sqrt(r2_lo) : nan;
        const double hi = finite ? std::sqrt(r2_hi) : nan;
        return first_failure(check_position(ExternalModel::MeadFairfield1975, lo, x_lo),
                             check_position(ExternalModel::MeadFairfield1975, hi, x_lo));
    }
};

/**
 * The Mead-Fairfield field over a whole batch of GSM points, on the device when that is worth it.
 *
 * The shape mirrors @ref t89_field_batch exactly — device above the registry's measured crossover,
 * the fp64 host loop otherwise, one folded @ref Status for the batch, and a returned value that
 * says which lane actually served the call. What differs is the arithmetic intensity: this model
 * is ~50 flops for 24 bytes in and 12 out, about 1.4 flops/byte, which is within a factor of three
 * of the streaming dipole kernel that LOSES on this seam and an order of magnitude below T89's
 * ~11. The `irbem_mead_f32` row of `gpu/dispatch.hpp` carries the measurement and the verdict.
 *
 * **The batch reports the same caveats the scalar entry point does, folded over the whole batch.**
 * If any point is beyond the published `r <= 17 R_E` the batch says @ref Status::OutOfValidityRange
 * and is still computed in full; if any point is inside the Earth or not finite it says
 * @ref Status::DomainError and every output is zeroed.
 *
 * @param points the positions, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians, `|psi| <= pi/2` as for @ref mead_field.
 * @param kp_times_ten Kp in IRBEM's slot-1 scaling; binned by @ref mead_kp_bin.
 * @param out receives one field vector per input, GSM, nanotesla; same length as @p points.
 * @return @ref Status::DomainError on a length mismatch, a non-finite tilt or Kp, a tilt beyond a
 *         right angle, or a point that is not finite or is inside the Earth, and then every output
 *         is zeroed;
 *         @ref Status::OutOfValidityRange when Kp or any point's radius is outside the published
 *         range, with every point still computed; otherwise @ref Status::Ok. The value is `true`
 *         exactly when the device lane serviced the call.
 * @complexity O(N); on the device those N run concurrently over `ceil(N/256)` workgroups.
 * @alloc the device lane stages positions and results into two `std::vector<float>` of `3N`; the
 *        host lane allocates nothing.
 * @test IrbemMead.BatchAgreesWithTheReferenceLane
 * @test IrbemMead.BatchRejectsMismatchedSpans
 * @test IrbemMead.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
[[nodiscard]] inline Result<bool> mead_field_batch(
    std::span<const Position<Frame::GSM>> points, double tilt_rad, double kp_times_ten,
    std::span<FieldVector<Frame::GSM>> out) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(tilt_rad) || !std::isfinite(kp_times_ten)) return {Status::DomainError, false};
    if (std::fabs(tilt_rad) > max_tilt_rad) return {Status::DomainError, false};

    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Kp)] = kp_times_ten;
    const Status drives = check_validity(ExternalModel::MeadFairfield1975, drivers);
    if (n == 0) return {drives, false};

    const int bin = mead_kp_bin(kp_times_ten);
    const double sin_tilt = std::sin(tilt_rad);
    const double cos_tilt = std::cos(tilt_rad);
    const double tilt_deg = tilt_rad * mead_deg_per_rad;
    MeadPositionFold fold;

#if CHEATAH_SPACE_IRBEM_MEAD_GPU
    if (gpu::prefer_gpu("irbem_mead_f32", n) &&
        std::filesystem::exists(gpu::shader_path("irbem_mead_f32"))) {
        std::vector<float> pos(3 * n);
        std::vector<float> raw(3 * n);
        for (std::size_t i = 0; i < n; ++i) {
            fold.add(points[i]);
            pos[(3 * i) + 0] = static_cast<float>(points[i].v[0]);
            pos[(3 * i) + 1] = static_cast<float>(points[i].v[1]);
            pos[(3 * i) + 2] = static_cast<float>(points[i].v[2]);
        }
        const Status where = fold.verdict();
        if (where == Status::DomainError) {
            for (std::size_t i = 0; i < n; ++i) out[i] = FieldVector<Frame::GSM>{};
            return {Status::DomainError, false};
        }
        const std::array<float, mead_param_count> block =
            mead_param_block(static_cast<float>(sin_tilt), static_cast<float>(cos_tilt),
                             static_cast<float>(tilt_deg), bin);
        gpu::dispatch_batch("irbem_mead_f32", pos, raw, std::span<const float>(block));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = FieldVector<Frame::GSM>{
                fixarray::vec3d{raw[(3 * i) + 0], raw[(3 * i) + 1], raw[(3 * i) + 2]}};
        }
        return {first_failure(drives, where), true};
    }
#endif

    const MeadParameters<double> p = mead_parameters<double>(bin);
    for (std::size_t i = 0; i < n; ++i) {
        fold.add(points[i]);
        const std::array<double, 3> b = mead_components<double>(
            p, sin_tilt, cos_tilt, tilt_deg, points[i].v[0], points[i].v[1], points[i].v[2]);
        out[i] = FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
    }
    const Status where = fold.verdict();
    if (where == Status::DomainError) {
        for (std::size_t i = 0; i < n; ++i) out[i] = FieldVector<Frame::GSM>{};
        return {Status::DomainError, false};
    }
    return {first_failure(drives, where), false};
}

// -------------------------------------------------------------------------------------------
// The superposition with IGRF — what a tracer follows
// -------------------------------------------------------------------------------------------

/**
 * IGRF plus Mead & Fairfield (1975), as a single field.
 *
 * The same shape as @ref TotalFieldT89, for the same reason: a particle experiences one field, and
 * @ref trace_invariant and @ref make_lstar follow whatever satisfies @ref GeoFieldModel without
 * knowing there are two of them. IGRF is defined in GEO and the external model in GSM; the rotation
 * between them depends only on the epoch and is built once in @ref Rotations.
 *
 * There is no device trace lane for this model — traces through it run on the host through the
 * generic @ref trace_invariant, which is where a model this cheap belongs: the external term adds
 * ~50 flops to each of ~600 IGRF evaluations of a trace, so a dedicated kernel would buy almost
 * nothing over `irbem_trace_i_f32`'s measured curve.
 *
 * @tparam NMAX the internal field's truncation degree. 10 reproduces IRBEM's own choice.
 * @test IrbemMead.TotalFieldSuperposesInternalAndExternal
 */
template <int NMAX = 10>
class TotalFieldMead {
  public:
    /// The internal part's truncation degree — what generic staging and buffer sizing read.
    static constexpr int degree = NMAX;

    /**
     * @param internal the internal field, already built for the epoch.
     * @param rotations the epoch's frame rotations — built once, reused for every point.
     * @param kp_times_ten Kp in IRBEM's slot-1 scaling (Kp x 10, nominally 0..90). Binned by
     *        @ref mead_kp_bin: values inside one group give identical fields by construction.
     * @complexity O(1) — three pointers and a double.
     * @alloc none.
     * @test IrbemMead.TotalFieldSuperposesInternalAndExternal
     */
    constexpr TotalFieldMead(const Igrf<NMAX>& internal, const Rotations& rotations,
                             double kp_times_ten)
        : internal_(&internal), rotations_(&rotations), kp_times_ten_(kp_times_ten) {}

    /**
     * The total field at a geographic point.
     *
     * @param p the position, GEO, Earth radii.
     * @return `B_internal + B_external`, in GEO, nT. When the external model refuses the point —
     *         a non-finite input or a point inside the Earth — the INTERNAL field is returned
     *         alone rather than a zero or a NaN, exactly as @ref TotalFieldT89 does; an
     *         out-of-validity external field is still added, because extrapolating the fit is the
     *         caller's decision and @ref external_status is how they learn it was made.
     * @complexity One IGRF evaluation, one Mead-Fairfield evaluation, two 3x3 rotations.
     * @alloc none.
     * @test IrbemMead.TotalFieldSuperposesInternalAndExternal
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::GEO>& p) const {
        const FieldVector<Frame::GEO> b_int = internal_->evaluate(p);
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        const double tilt_rad = rotations_->dipole_tilt_deg * (std::numbers::pi / 180.0);
        const Result<FieldVector<Frame::GSM>> b_ext = mead_field(p_gsm, tilt_rad, kp_times_ten_);
        if (b_ext.status == Status::DomainError) return b_int;
        const FieldVector<Frame::GEO> b_ext_geo = transform<Frame::GEO>(b_ext.value, *rotations_);
        return FieldVector<Frame::GEO>{b_int.v + b_ext_geo.v};
    }

    /**
     * Whether the external model answered at @p p, and if not, why.
     * @param p the position, GEO, Earth radii.
     * @return the external model's status; @ref Status::Ok when it contributed without caveat.
     * @complexity One Mead-Fairfield evaluation and one rotation.
     * @alloc none.
     * @test IrbemMead.TotalFieldReportsWhenTheExternalModelDeclines
     */
    [[nodiscard]] Status external_status(const Position<Frame::GEO>& p) const {
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        const double tilt_rad = rotations_->dipole_tilt_deg * (std::numbers::pi / 180.0);
        return mead_field(p_gsm, tilt_rad, kp_times_ten_).status;
    }

    /// The activity level this field was built for, in IRBEM's Kp x 10 scaling.
    /// @return the value passed to the constructor.
    /// @complexity O(1). @alloc none.
    /// @test IrbemMead.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double kp_times_ten() const { return kp_times_ten_; }

    /// The epoch's frame rotations.
    /// @return the rotations this field was built with.
    /// @complexity O(1). @alloc none.
    /// @test IrbemMead.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Rotations& rotations() const { return *rotations_; }

    /// The internal part's Gauss coefficient `g(n, m)`, in nT — the INTERNAL field's, because a
    /// superposition with a non-potential external term has no expansion of its own; see
    /// @ref TotalFieldT89::g for the full argument.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemMead.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double g(int n, int m) const { return internal_->g(n, m); }

    /// The internal part's `h(n, m)`, in nT — see @ref g for why this forwards.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemMead.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double h(int n, int m) const { return internal_->h(n, m); }

    /// The internal field alone — what `dipole_moment` and any staging need.
    /// @return the internal model.
    /// @complexity O(1). @alloc none.
    /// @test IrbemMead.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Igrf<NMAX>& internal() const { return *internal_; }

  private:
    const Igrf<NMAX>* internal_;
    const Rotations* rotations_;
    double kp_times_ten_;
};

}  // namespace cheatah::space::irbem
