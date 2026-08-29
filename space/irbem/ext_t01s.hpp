#pragma once

/**
 * @file ext_t01s.hpp
 * @brief space.irbem — the storm-time external field: Tsyganenko, Singer & Kasper (2003),
 *        IRBEM's `kext = 10` ("T01 storm"), the member of this family built FOR the main phase.
 *
 * T89 (`ext_t89.hpp`) is the quiet-to-moderate workhorse and is binned in Kp. This model is the
 * other end of the family: fitted to 37 major storms of 1996-2000 (Dst minima to -300 nT and
 * beyond) and driven CONTINUOUSLY by six quantities — Dst, the solar-wind dynamic pressure, the
 * IMF By and Bz, and two solar-wind history integrals G2 and G3 — so that the ring current can
 * deepen, the tail current can move earthward, and the dawn-dusk asymmetry of the partial ring
 * current can grow, all smoothly, as a storm develops. Every one of those responses is measured
 * against the oracle in the brief below, and @ref t01s_field is where they enter a drift-shell
 * trace.
 *
 * ## The paper, and what it does and does not publish
 *
 * Tsyganenko, Singer & Kasper, *Storm-time distortion of the inner magnetosphere: How severe can it
 * get?*, J. Geophys. Res. **108**(A5), 1209, doi:10.1029/2002JA009808 (2003). The model is the
 * Tsyganenko (2002a) structure — *A model of the near magnetosphere with a dawn-dusk asymmetry:
 * 1. Mathematical structure*, JGR **107**(A8), 1179 — re-fitted to storm data with a storm-time
 * driver set. The paper publishes the MODULE STRUCTURE (the magnetopause/Chapman-Ferraro
 * shielding, the cross-tail current sheet, the axisymmetric and the partial ring currents, the
 * Birkeland current systems, each driven by its own variable) and the DRIVERS (the IRBEM `kext`
 * table's wording: "uses Dst, Pdyn, By, Bz, G2, G3 ... there is no upper or lower limit for those
 * inputs ... valid for xGSM >= -15 Re"; `G2 = <a Vsw Bs>` with `a = 0.005`, `G3 = <Vsw Dsw Bs / 2000>`,
 * `Bs = max(-Bz, 0)`). It does NOT publish the numerical coefficient set: several hundred numbers,
 * most of them shielding-field expansions, that exist only in the author's code distribution — and
 * that distribution is GPL-3.0, which this MIT clean room may not read.
 *
 * **So the provenance verdict is: published form, unpublished coefficients.** The T89 lesson applied
 * in advance (see `ext_t89.hpp`'s brief): there is no table to transcribe, so there is nothing to
 * reach oracle parity WITH, and pretending otherwise would be the silent disagreement this module
 * refuses to ship. What this header does instead, in the open:
 *
 * 1. It carries the published structure as a sum of **108 exactly divergence-free spatial modes**
 *    (@ref t01s_mode_fields), each the curl of an analytic vector potential or the gradient of an
 *    analytic harmonic function, grouped by the module they stand for:
 *    - the **cross-tail and ring currents** as thickened, warped, hinged current discs (the
 *      construction of Tsyganenko 1989, eqs. 7-17: the three disc potentials `P1 = 1/(S+u)`,
 *      `P2 = -1/(S(S+u))`, `P3 = 1/S^3`, thickened by `|z| -> sqrt(z^2 + D^2)`, warped by the
 *      eq. (11) sheet surface, truncated sunward by the eq. (12)-(13) factor `W`), extended here
 *      with **azimuthal weights** `h(x, y) in {1, y, x, x^2 - y^2, x y}` so that the partial
 *      ring current's `sin phi` and `cos phi` structure — the 2002a paper's dawn-dusk asymmetry —
 *      and its second harmonic are representable;
 *    - the **partial ring current's closure** as thin equatorial sheets of RADIAL current,
 *      `B = curl(Lambda z_hat)`, `Lambda = G(rho) H(z) q(phi)` — inward at dusk, outward at dawn
 *      — which the oracle probe below shows to be the dominant storm-time feature at 4 R_E;
 *    - the **magnetopause shielding** as box harmonics `exp(x/s) cos(y/p) sin(z/r)` and their
 *      tilt-odd and dawn-dusk-odd partners (the general form of Tsyganenko 1995, JGR 100, 5599,
 *      eq. 24, used again in 2002a §2);
 *    - the **interconnection field** as the two uniform penetration fields `y_hat` and `z_hat`.
 * 2. Each mode's amplitude is a **linear combination of eleven smooth driver features**
 *    (@ref t01s_features): `1, Dst, sqrt(Pdyn), Pdyn, G2, G3, G3^2/100, By, Bs, Bz,
 *    Dst sqrt(Pdyn)/10`. Linear in smooth features is what makes the model continuous in every
 *    driver by construction — the property @ref IrbemT01s.FieldIsContinuousInEveryDriver pins
 *    down, and the property a Kp-binned model like T89 does not have.
 * 3. The 108 x 11 coefficients (@ref t01s_coefficients) were determined by **weighted least
 *    squares against the IRBEM oracle's `kext = 10` minus `kext = 0` field, run as a black box**
 *    (`tools/oracle/t01s_diff.cpp`), over 71 514 samples: 48 driver states (the four corpus regimes,
 *    the four real storm events, and 40 storm-shaped random states) x three dipole tilts
 *    (+0.0, +25.6, -30.4 degrees) x a grid inside the Shue et al. (1998) magnetopause with
 *    `x >= -15 R_E` and `r >= 2.2 R_E`. Fitting a black box's OUTPUT is observation, not
 *    transcription; no source of the oracle or of Tsyganenko's distribution was read.
 *
 * ## The gap, measured
 *
 * Measured by [`tools/oracle/t01s_diff.cpp`](../../../tools/oracle/t01s_diff.cpp) against the
 * `-O2` oracle on a HOLDOUT set — two fresh epochs (tilts -23.2 and +19.3 degrees), a jittered
 * off-grid point set, and driver states never fitted — with the external field isolated as
 * `kext = 10` minus `kext = 0` so IGRF cancels exactly, and the tilt taken from the oracle itself:
 *
 * | regime (belts, 2.2-10 R_E, in the magnetopause) | N | RMS abs | p99 abs | RMS rel |
 * |---|---|---|---|---|
 * | quiet (Dst -8, Pdyn 1.8)          | 760 | 4.18 nT | 9.49 nT | 19.3% |
 * | moderate (Dst -42, G3 20)         | 746 | 23.5 | 69.4 | 44.6% |
 * | storm (Dst -150, Pdyn 9, G3 80)   | 640 | 66.9 | 191 | 42.5% |
 * | extreme (Dst -350, Pdyn 28, G3 180) | 592 | 103 | 324 | 33.5% |
 * | Halloween 2003 (Dst -383)         | 598 | 116 | 359 | 35.6% |
 * | March 1989 (Dst -589)             | 592 | 159 | 388 | 38.6% |
 * | St Patrick's 2015 (Dst -223)      | 626 | 89.0 | 292 | 40.5% |
 * | Starlink 2022 (Dst -66)           | 662 | 52.9 | 158 | 48.7% |
 *
 * Whole holdout set, 22 036 points to `x = -15`: RMS 57.8 nT, 38.2% relative; training set 41.4%.
 * Training and holdout agree, so the number is a model-family gap and not an overfit.
 *
 * **Where the gap is, and why more modes do not close it.** With every one of the 108 amplitudes
 * free for ONE driver state (no driver parameterization at all), the residual floor is 15.4%
 * quiet, 37.6% storm, 29.0% extreme, 24.0% for March 1989 — within a few points of the calibrated
 * model's numbers above, so the driver features are not the limit; the SPATIAL basis is. Probing
 * the oracle directly says what it lacks: at storm time the oracle's field at `(0, ±4, z)` carries
 * a sheet of RADIAL current of half-thickness ~0.5 R_E across which `B_x` jumps by ~180 nT, and at
 * `r = 4`, 40 degrees latitude, `B_x` rises from 9 to 154 nT within 30 degrees of azimuth at dawn
 * and dusk — the partial ring current's field-aligned closure. Adding the radial-sheet family took
 * the storm floor from 52% to 38%; the remainder is the Birkeland sheets, which are current-carrying
 * structures off the equator that neither a disc nor a harmonic can represent. Those are the
 * Tsyganenko (2002a) conical-harmonic and deformed-FAC modules whose coefficients are exactly the
 * unpublished part. The gap is stated here rather than hidden: a caller doing storm-time L\*
 * through this model gets the published storm response — smooth deepening with Dst, sunward
 * compression with Pdyn, southward-Bz and G2/G3 loading of the tail — with a field-level error of
 * 35-45% RMS in the main phase against IRBEM's implementation, and 19% in quiet time.
 *
 * ## Two independent checks that THIS side is right regardless
 *
 * Neither can pass by accident. @ref IrbemT01s.DivergenceVanishesEverywhere verifies that every
 * analytic derivative in the disc, sheet and harmonic families is the derivative it claims to be:
 * `|div B|` falls as `h^2` from ~1e-6 to ~1e-10 relative as the stencil step goes from 1e-2 to 1e-4
 * R_E, which is the signature of an exactly divergence-free field sampled by a second-order stencil
 * and of nothing else — and the same test, fed a deliberately non-solenoidal field, fails. And
 * @ref IrbemT01s.OracleSamplesLandInsideTheDocumentedGap pins a dozen oracle values, storm and
 * quiet, noon and midnight, so the calibration cannot silently drift.
 *
 * @note Nothing on a hot path allocates, and nothing but the device lane can throw. The evaluator,
 *       the scalar entry points and the fp32 host lane touch nothing but the caller's spans and
 *       stack arrays; the amplitude block is 108 doubles on the stack. The one exception is
 *       @ref t01s_field_batch's DEVICE lane, which stages `3N` floats each way and forwards whatever
 *       `gpu::dispatch_batch` throws; its `@alloc` says so.
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
/// 1 when this translation unit can reach the T01S device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_T01S_GPU 1
#else
/// 1 when this translation unit can reach the T01S device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_T01S_GPU 0
#endif

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// The drivers
// -------------------------------------------------------------------------------------------

/**
 * The six drivers IRBEM's `kext = 10` reads, in the paper's units.
 *
 * Gathered into one aggregate rather than passed as six doubles so that a call site cannot swap
 * G2 and G3, which are both dimensionless, both non-negative and both "about 10 in a storm".
 *
 * @test IrbemT01s.DriversAreReadFromTheMaginputSlots
 */
struct T01sDrivers {
    double dst;     ///< Dst, nT (`maginput` slot 2).
    double pdyn;    ///< Solar-wind dynamic pressure, nPa (slot 5).
    double by_imf;  ///< IMF By, GSM, nT (slot 6).
    double bz_imf;  ///< IMF Bz, GSM, nT (slot 7).
    double g2;      ///< `G2 = <a Vsw Bs>`, `a = 0.005`, an hourly average (slot 9).
    double g3;      ///< `G3 = <Vsw Dsw Bs / 2000>`, an hourly average (slot 10).
};

/**
 * The six drivers, read out of a full `maginput` vector.
 *
 * @param d the whole 25-slot vector, in IRBEM's order; the other nineteen slots are never read,
 *        which is what makes a vector full of fill values in the slots this model ignores still
 *        a valid input for it.
 * @return the six the model uses.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01s.DriversAreReadFromTheMaginputSlots
 */
[[nodiscard]] constexpr T01sDrivers t01s_drivers(const DriverSet& d) {
    return T01sDrivers{d[static_cast<std::size_t>(Driver::Dst)],
                       d[static_cast<std::size_t>(Driver::Pdyn)],
                       d[static_cast<std::size_t>(Driver::ByIMF)],
                       d[static_cast<std::size_t>(Driver::BzIMF)],
                       d[static_cast<std::size_t>(Driver::G2)],
                       d[static_cast<std::size_t>(Driver::G3)]};
}

// -------------------------------------------------------------------------------------------
// The geometry — fixed, not fitted per driver state
// -------------------------------------------------------------------------------------------

/**
 * The current-sheet geometry shared by every disc and sheet mode, held fixed across all driver
 * states.
 *
 * These are the non-linear parameters of the Tsyganenko (1989) sheet construction — the hinge, the
 * warp, the truncation — chosen once from that paper's fitted range (its `R_c`, `G`, `x_0`, `D_y`
 * columns of Table 1 bracket the values here) and then HELD, so that every driver dependence of the
 * model is carried by the linear amplitudes and the field stays linear in the calibrated
 * coefficients. Moving the sheet with the drivers is what the 2003 paper does for its tail inner
 * edge; here the same effect is spanned by discs at several radial scales instead, which is what
 * keeps the amplitude solve linear.
 *
 * @test IrbemT01s.GeometryIsWithinThePublishedSheetRange
 */
struct T01sGeometry {
    double delta_y;   ///< `delta`, R_E^-1 — the rate at which the tail sheet thickens towards
                      ///< its flanks, T89 eq. (13); applied to the truncated modes only.
    double r_hinge;   ///< `R_c`, R_E — the hinging distance of the sheet surface, T89 eq. (11).
    double warp_g;    ///< `G`, R_E — the amplitude of the sheet's dawn-dusk bending, eq. (11).
    double l_y;       ///< `L_y`, R_E — the dawn-dusk scale of that bending, eq. (11).
    double hinge2;    ///< The `16` in eq. (11)'s square root: `(4 R_E)^2`.
    double x0_w;      ///< `x_0`, R_E — where the tail truncation factor `W` transitions, eq. (13).
    double dx_w;      ///< `D_x`, R_E — the sunward scale of `W`, eq. (13).
    double dy_w;      ///< `D_y`, R_E — the dawn-dusk scale of `W`, eq. (13).
};

/// The geometry the coefficient table was calibrated with. Changing any value here invalidates
/// @ref t01s_coefficients; the harness re-derives both together.
inline constexpr T01sGeometry t01s_geometry{
    /* delta_y */ 0.01, /* r_hinge */ 8.0, /* warp_g */ 4.0, /* l_y */ 10.0,
    /* hinge2 */ 16.0,  /* x0_w */ 6.0,    /* dx_w */ 13.0,  /* dy_w */ 20.0};

// -------------------------------------------------------------------------------------------
// The mode tables
// -------------------------------------------------------------------------------------------

/**
 * One current-carrying disc mode.
 *
 * The three disc potentials are Tsyganenko (1989) eqs. (7)-(9): with `xi = sqrt(z_r^2 + D^2)`,
 * `u = a + xi`, `S = sqrt(rho^2 + u^2)`, family 1 is `P1 = 1/(S + u)` (the paper's `A^(1)`, the
 * field of a thin disc), family 2 is `P2 = dP1/da = -1/(S (S + u))` and family 3 is
 * `P3 = dP2/da = 1/S^3` — the only one with a finite magnetic moment, hence the ring current.
 * `P` is the potential over radius: `A = P (-y, x, 0)`. The azimuthal weight multiplies `P` INSIDE
 * the potential, so the field stays exactly divergence-free whatever the weight.
 *
 * @test IrbemT01s.EveryModeIsDivergenceFree
 */
struct T01sDiscMode {
    int family;       ///< 1: `P = 1/(S+u)`; 2: `P = -1/(S(S+u))`; 3: `P = 1/S^3`.
    double a;         ///< The radial scale `a`, R_E.
    double d;         ///< The half-thickness `D`, R_E (at midnight, for truncated modes).
    bool truncated;   ///< Multiply by the tail truncation factor `W(x, y)` and thicken flankward.
    int azimuth;      ///< Homogeneous weight `h(x, y)`: 0: `1`; 1: `y`; 2: `x`; 3: `x^2 - y^2`;
                      ///< 4: `x y`. Degrees 0, 1, 1, 2, 2.
};

/**
 * One equatorial radial-current sheet mode — the partial ring current's closure system.
 *
 * The curl of `Lambda z_hat` with `Lambda = G(rho) H(z_r) q`, `G = (rho^2 + a^2)^(-3/2)`,
 * `H = z_r / sqrt(z_r^2 + d^2)` and `q` an azimuthal weight. A field with no `z` component in the
 * sheet's own frame, that flips sign across the sheet — which is the field of a thin sheet of
 * current flowing RADIALLY in the equatorial plane, inward on one side and outward on the other
 * for `q = y`. Exactly divergence-free because it is a curl, whatever `G`, `H` and `q` are.
 *
 * @test IrbemT01s.EveryModeIsDivergenceFree
 */
struct T01sSheetMode {
    double a;      ///< Radial scale, R_E.
    double d;      ///< Half-thickness, R_E.
    int azimuth;   ///< 1: `q = y` (sin phi); 2: `q = x` (cos phi); 3: `x^2 - y^2`; 4: `x y`.
};

/**
 * One harmonic mode — the gradient of a box harmonic, curl-free and divergence-free.
 *
 * `gamma = exp(x / s) f(y / p) g(z / r)` with `1/s^2 = 1/p^2 + 1/r^2`, which is what makes
 * `Laplacian(gamma) = 0` exactly and therefore `div grad gamma = 0`. Tsyganenko (1995) eq. (24)
 * and 2002a §2 use this family for the magnetopause shielding of every module; here it also
 * stands in for the shielding of the modules this header does not carry.
 *
 * @test IrbemT01s.EveryModeIsDivergenceFree
 */
struct T01sHarmonicMode {
    int kind;   ///< 0: `cos(y/p) sin(z/r)` (tilt-even); 1: `sin(psi) cos(y/p) cos(z/r)`
                ///< (tilt-odd); 2: `sin(y/p) sin(z/r)` (dawn-dusk-odd); 3: `sin(psi) sin(y/p)
                ///< cos(z/r)` (both).
    double p;   ///< Dawn-dusk scale, R_E.
    double r;   ///< North-south scale, R_E.
};

/// The disc modes: the tail (truncated, flank-thickening) at three radial scales, and the ring
/// current at four scales with five azimuthal weights each in two families, plus four thin rings
/// where the storm-time ring current concentrates.
inline constexpr std::array<T01sDiscMode, 58> t01s_disc_modes{{
    // tail: truncated, warped, thickening towards the flanks
    {1, 12.0, 2.0, true, 0}, {2, 12.0, 2.0, true, 0}, {3, 12.0, 2.0, true, 0},
    {1, 12.0, 2.0, true, 1}, {1, 12.0, 2.0, true, 2},
    {1, 6.0, 1.5, true, 0},  {3, 6.0, 1.5, true, 0},
    {1, 20.0, 3.0, true, 0}, {2, 20.0, 3.0, true, 0},
    // ring currents at four radial scales, five azimuthal weights each
    {3, 2.5, 0.8, false, 0}, {3, 2.5, 0.8, false, 1}, {3, 2.5, 0.8, false, 2}, {3, 2.5, 0.8, false, 3}, {3, 2.5, 0.8, false, 4},
    {3, 3.5, 1.0, false, 0}, {3, 3.5, 1.0, false, 1}, {3, 3.5, 1.0, false, 2}, {3, 3.5, 1.0, false, 3}, {3, 3.5, 1.0, false, 4},
    {3, 5.0, 1.5, false, 0}, {3, 5.0, 1.5, false, 1}, {3, 5.0, 1.5, false, 2}, {3, 5.0, 1.5, false, 3}, {3, 5.0, 1.5, false, 4},
    {3, 7.0, 2.0, false, 0}, {3, 7.0, 2.0, false, 1}, {3, 7.0, 2.0, false, 2}, {3, 7.0, 2.0, false, 3}, {3, 7.0, 2.0, false, 4},
    // the P1 family at the same scales (a different radial profile)
    {1, 2.5, 0.8, false, 0}, {1, 2.5, 0.8, false, 1}, {1, 2.5, 0.8, false, 2}, {1, 2.5, 0.8, false, 3}, {1, 2.5, 0.8, false, 4},
    {1, 3.5, 1.0, false, 0}, {1, 3.5, 1.0, false, 1}, {1, 3.5, 1.0, false, 2}, {1, 3.5, 1.0, false, 3}, {1, 3.5, 1.0, false, 4},
    {1, 5.0, 1.5, false, 0}, {1, 5.0, 1.5, false, 1}, {1, 5.0, 1.5, false, 2}, {1, 5.0, 1.5, false, 3}, {1, 5.0, 1.5, false, 4},
    {1, 7.0, 2.0, false, 0}, {1, 7.0, 2.0, false, 1},
    // thin rings where the storm-time ring current concentrates
    {3, 3.0, 0.4, false, 0}, {3, 3.0, 0.4, false, 1}, {3, 3.0, 0.4, false, 2},
    {3, 3.5, 0.4, false, 0}, {3, 3.5, 0.4, false, 1}, {3, 3.5, 0.4, false, 2},
    {3, 4.0, 0.4, false, 0}, {3, 4.0, 0.4, false, 1}, {3, 4.0, 0.4, false, 2},
    {3, 4.5, 0.4, false, 0}, {3, 4.5, 0.4, false, 1}, {3, 4.5, 0.4, false, 2},
}};

/// The radial-current sheet modes: five radial scales at two thicknesses for the `sin phi` and
/// `cos phi` closure patterns, and two scales for the second harmonic.
inline constexpr std::array<T01sSheetMode, 24> t01s_sheet_modes{{
    {2.0, 0.4, 1}, {2.0, 0.4, 2}, {3.0, 0.4, 1}, {3.0, 0.4, 2}, {4.0, 0.4, 1}, {4.0, 0.4, 2},
    {5.0, 0.4, 1}, {5.0, 0.4, 2}, {7.0, 0.4, 1}, {7.0, 0.4, 2},
    {2.0, 1.2, 1}, {2.0, 1.2, 2}, {3.0, 1.2, 1}, {3.0, 1.2, 2}, {4.0, 1.2, 1}, {4.0, 1.2, 2},
    {5.0, 1.2, 1}, {5.0, 1.2, 2}, {7.0, 1.2, 1}, {7.0, 1.2, 2},
    {3.0, 0.6, 3}, {3.0, 0.6, 4}, {5.0, 0.6, 3}, {5.0, 0.6, 4},
}};

/// The harmonic modes: eight tilt-even, seven tilt-odd, six dawn-dusk-odd and three doubly odd
/// box harmonics over dawn-dusk scales of 6, 12 and 24 R_E and north-south scales of 4, 8 and
/// 16 R_E.
inline constexpr std::array<T01sHarmonicMode, 24> t01s_harmonic_modes{{
    {0, 6.0, 4.0},  {0, 6.0, 8.0},  {0, 12.0, 4.0}, {0, 12.0, 8.0}, {0, 12.0, 16.0}, {0, 24.0, 8.0}, {0, 24.0, 16.0}, {0, 24.0, 4.0},
    {1, 6.0, 4.0},  {1, 6.0, 8.0},  {1, 12.0, 4.0}, {1, 12.0, 8.0}, {1, 24.0, 8.0}, {1, 24.0, 16.0}, {1, 24.0, 4.0},
    {2, 6.0, 4.0},  {2, 6.0, 8.0},  {2, 12.0, 4.0}, {2, 12.0, 8.0}, {2, 24.0, 8.0}, {2, 24.0, 16.0},
    {3, 6.0, 4.0},  {3, 12.0, 8.0}, {3, 24.0, 8.0},
}};

/// How many spatial modes the model sums: the discs, the sheets, the harmonics, and the two
/// uniform penetration fields `y_hat` and `z_hat`.
inline constexpr std::size_t t01s_mode_count =
    t01s_disc_modes.size() + t01s_sheet_modes.size() + t01s_harmonic_modes.size() + 2;

/// How many smooth driver features each mode's amplitude is linear in; see @ref t01s_features.
inline constexpr std::size_t t01s_feature_count = 11;

/// The sunward limit past which the model is refused outright, R_E. The shielding harmonics carry
/// `exp(x / s)` with `s` as small as 3.3 R_E, so a position far enough sunward — nothing a trace can
/// reach, but exactly what a kilometres-for-Earth-radii unit confusion produces — overflows them:
/// past ~300 R_E in `float`. The magnetopause is never beyond 15 R_E, so 100 R_E refuses nothing a
/// caller could mean and keeps every lane, including the fp32 device lane, finite.
inline constexpr double t01s_max_x_gsm = 100.0;

// -------------------------------------------------------------------------------------------
// The calibrated coefficients (mode-major: coefficient[mode][feature])
// -------------------------------------------------------------------------------------------

/**
 * The amplitude coefficients: row `m` is mode `m`'s weight on each of the eleven features.
 *
 * **These are NOT published values.** The 2003 paper prints no coefficient table (see the file
 * brief); every number here was determined by `tools/oracle/t01s_diff.cpp`, a weighted linear
 * least-squares fit of THIS header's mode fields to the IRBEM oracle's `kext = 10` external field
 * observed as a black box, with a ridge of 1e-3 on the column-scaled normal matrix so that nearly
 * collinear modes cannot cancel through amplitudes the fp32 device lane could not resolve (the
 * largest storm-regime amplitude is 1.4e5 against ~1e6 without the ridge, for 1.5 points of RMS).
 * The mode order is @ref t01s_disc_modes, then @ref t01s_sheet_modes, then
 * @ref t01s_harmonic_modes, then `y_hat`, then `z_hat`; the feature order is @ref t01s_features'.
 *
 * @test IrbemT01s.OracleSamplesLandInsideTheDocumentedGap
 * @test IrbemT01s.AmplitudesAreLinearInTheFeatures
 */
inline constexpr std::array<std::array<double, t01s_feature_count>, t01s_mode_count>
    t01s_coefficients{{
        {{97.6708514, -0.732333147, 44.9929728, 6.06786867, 4.4249999, 2.22121672, -0.489490974, 2.55065805, 3.19426232, -3.24016116, -1.68363252}},
        {{4207.19121, -26.9944784, 1370.7791, 171.500812, 157.259972, 21.6972394, 42.42706, 48.0572867, 95.5096816, -26.3758417, -9.08848693}},
        {{4073.5204, 233.894807, -7680.04963, -2762.00889, -651.552945, 52.6328269, -393.743966, 927.315779, 81.6857558, 239.682755, 618.89943}},
        {{-0.691537011, -0.0430497432, 0.428375621, 0.191808951, 0.0129309468, 0.0368833847, -0.0175473012, 0.00128833988, -0.0240992386, 0.0121188165, -0.0105536297}},
        {{-10.6424101, 0.126144517, -3.33926621, 0.58668451, -0.0152186758, -0.183743897, 0.389065415, 0.253352537, 0.493767063, -0.257867519, -0.239392645}},
        {{1.84826428, 0.0372866731, -1.65045817, -1.38512974, 0.61084521, 1.07209648, -0.555768844, 0.83167943, 0.220930564, -0.786259296, -0.297875175}},
        {{-4453.96995, 62.0070466, 430.656591, 449.874126, 80.404217, 95.1444568, 177.417118, -196.300316, -24.9793229, 113.312287, 49.0500923}},
        {{500.874155, -2.31837692, 162.203527, 11.785956, 16.7479575, 5.52594489, -0.692435085, 11.3339521, 11.3639861, -6.08850528, -0.32011862}},
        {{8705.83154, -6.29876817, 1042.62416, -182.110038, 145.316897, -21.9029391, 1.7520994, 175.835248, 142.75779, 34.0590083, 120.878649}},
        {{460.974909, -21.2240248, 176.065628, -24.9909932, 29.0628206, 14.3808669, -6.55923228, 21.9849518, 5.94877992, -8.95549713, -26.0036574}},
        {{52.8668382, 0.493126236, 19.8102078, 6.61031811, 2.80321773, 1.49050759, 0.411677735, 0.976326696, 2.45489985, -2.79725326, 2.32881035}},
        {{-45.9794328, 3.44950175, -27.9465949, -15.6444222, -8.6343553, -22.6941381, 2.39186854, -4.18609969, -0.321872941, 4.68125638, -10.0828948}},
        {{-2.15707395, -0.127180783, 0.127413465, -1.23779908, -0.227485817, 1.63532087, -0.697495608, -0.266679355, 0.227496489, 0.134014956, 1.10044711}},
        {{-1.37830613, -0.128545813, 2.4470819, 0.444942615, 0.0819906378, 0.078125668, -0.0123063055, 0.01258612, -0.0181816674, -0.0671567051, 0.147359784}},
        {{-665.273771, 7.2723971, 509.72493, 149.604622, -11.8069067, -2.90198628, -36.0331547, 22.2849215, -4.98742382, -6.28482335, 63.9578935}},
        {{11.3640017, 5.93172228, -88.7932904, -17.1330647, -5.24403018, -3.91938947, -0.372926961, -2.15461219, -4.60559636, 5.27643448, 6.04520258}},
        {{328.193608, -9.46696072, 286.820179, 80.3910004, 49.7816529, 21.872909, 6.7781898, 34.7439474, 24.8177248, -29.0324003, -9.3931906}},
        {{-13.9838763, -0.335649775, -4.13994662, 0.619190539, 0.81048176, 1.6127092, 0.644147605, 0.772571776, 0.0739437225, -0.186953817, -0.395236549}},
        {{1.80653421, 0.157298486, -1.40248223, -0.626633114, -0.100658531, -0.0685902617, -0.0308396939, 0.0565069058, -0.0100335221, 0.0736490675, 0.429248906}},
        {{-2064.58207, 61.6589708, -431.377306, -9.14358858, -46.8983791, -30.5386574, -10.5646132, -113.566343, -68.4588449, 80.9769164, 131.314372}},
        {{127.51501, 3.88401617, -156.0961, -26.2472657, -3.83318292, -4.0404287, 2.05023796, -3.69577428, 1.7307879, 1.99228967, -13.6832767}},
        {{-1250.22395, 15.2248869, -496.682703, -69.4969963, -52.2492606, -45.8812573, -38.0001064, -32.0746961, -73.2900602, 42.3472738, 39.6363723}},
        {{14.3345426, 0.0458360005, -7.53387962, 1.97004566, -1.01441188, -1.84547613, -1.40962998, 0.588043285, -0.865559696, 0.142095241, -0.141144839}},
        {{-1.46221573, 0.241159324, -8.88715732, -1.69435057, -0.295174709, -0.201019742, 0.0818999069, -0.293842019, -0.196090963, 0.352083301, -0.496570784}},
        {{8299.38668, 39.863299, 281.256926, -124.323778, 192.353554, 27.7422955, 77.7876255, 94.2813639, 341.354249, -125.210395, -47.4278866}},
        {{256.162606, -4.65663728, 5.18559225, 31.9985891, 7.41369122, 4.63529622, 0.0882760651, 2.68466727, 12.1672616, -7.49224455, -27.4193696}},
        {{2164.29682, -6.21360965, 330.836972, 114.185045, 65.1210746, -50.929778, 25.8367749, -4.10723921, 69.8750731, -33.9147691, -82.4966662}},
        {{68.9730944, -0.0483236005, -2.57013421, 6.70222813, 2.91204562, -1.92022349, -1.61267948, -1.8944935, 0.205700076, -1.57657828, -0.443660195}},
        {{17.6898907, -0.145469858, 1.93134254, 3.22319788, 0.554971891, 0.180282142, 0.133407394, 0.185504723, 0.90018183, -0.538136835, -2.39428857}},
        {{61.6535356, -0.266323797, 3.81223977, 0.663764079, 2.44706368, 0.365767971, 0.89510671, 1.73143002, 2.11951553, -1.44749472, -2.66321384}},
        {{-0.554132062, -0.00575276933, 0.224170251, 0.114225004, -0.00177262854, 0.00972949031, -0.00936360735, -0.000584903684, -0.0175475897, 0.0064383524, 0.0366321028}},
        {{10.4323754, -0.0546124526, 0.438256838, -0.435171913, 0.192407798, 0.0925783261, 0.0356874747, 0.172447322, 0.196978927, 0.0165752228, -0.0233610969}},
        {{-0.326982716, -7.17009549e-06, -0.14879063, -0.00213323054, 0.00119745554, 0.00529189939, 0.0145271859, -0.00197374935, -0.00292205222, 0.00424275897, -0.00400278732}},
        {{-0.0163259587, -0.000386133059, 0.0048098748, 0.00250773779, 0.000478343541, 0.000290822595, -0.000510301972, 0.000497737907, -0.000595561368, 8.37928407e-05, 0.00216261797}},
        {{36.5871287, 0.265822688, -8.72115369, -1.15343722, 0.638830692, -0.74575104, 0.771743481, 1.00849473, 1.44705157, -1.06642928, -2.65436317}},
        {{-0.573634846, -0.00207792787, 0.0989746586, 0.0657166408, -0.0077357802, 0.00392574542, -0.00938385998, -0.00365653221, -0.022824926, 0.0124584046, 0.0401651228}},
        {{6.62379419, 0.000409757747, -0.985042703, -0.573082487, -0.0698287997, -0.0494687889, 0.0116233027, 0.0369968005, 0.117899153, 0.0571421206, -0.0176896497}},
        {{-0.0974825146, -7.8540468e-05, -0.0945817887, -0.00208346221, -0.00322891722, -3.41329785e-05, 0.00854948527, -0.00252910318, 0.00211856543, 0.00317282684, -0.00128551908}},
        {{-0.021801414, -0.000459868591, 0.00767858724, 0.00277526629, 1.03917153e-05, 0.000425641822, -0.000571026705, 0.000380975908, -0.000735970698, 0.000122234485, 0.00293220798}},
        {{-53.960343, 1.6316959, -44.9104852, -6.71659467, -4.81882138, -3.71284607, 0.310200418, -1.45669502, -1.4013715, 0.739970565, -2.09244601}},
        {{-0.358255733, 0.00463853451, -0.129391869, -0.0402713006, -0.013496373, -0.00675612117, -0.00652977725, -0.00397544746, -0.0199016776, 0.0144220811, 0.040306123}},
        {{-8.53280034, 0.156549006, -5.45730296, -0.886712781, -0.879937712, -0.428059547, -0.0582082606, -0.471651817, -0.261531196, 0.247109514, 0.0144051915}},
        {{0.532999271, 0.0014595907, 0.0529929528, -0.00446735422, -0.017891048, -0.0150576302, -0.00877868277, -0.000963033034, 0.0140461684, 0.00224279913, 0.0102714729}},
        {{-0.0495697344, -0.000733451039, 0.0100693719, 0.00317024081, -0.000997927317, 0.000688793005, -0.000571253588, -0.000360109401, -0.0017510244, 0.000705872242, 0.00405116771}},
        {{-193.176022, 3.42057029, -99.7864853, -16.2562946, -12.1163513, -7.52810426, -0.228672489, -5.16179285, -6.71439864, 4.39331598, -0.93508789}},
        {{0.550678935, 0.00187740601, -0.0665594905, -0.087416298, 0.0144399607, -0.00267670411, 0.0110871301, 0.0163243663, 0.032777774, -0.0250256966, 0.01300094}},
        {{-142.014534, 1.42960388, -201.233274, -85.1230443, -30.8842219, -14.0465549, 3.27424687, -20.4267356, -7.78741922, 5.78623523, -15.6276497}},
        {{-30.6369398, -1.13817356, -0.811389777, -0.169607005, -0.187725463, 0.517010769, -0.508000116, -0.477462159, -0.926048276, -0.203851968, 0.490704742}},
        {{-175.263791, 8.32267599, -158.417941, -61.3433856, -36.968163, -24.9611771, -8.03811171, -14.3638976, -7.88934135, 10.6655675, 5.76009243}},
        {{-14.5632474, -0.223756996, 139.906917, 16.6136751, -15.369202, -3.3033107, -1.32242671, 1.47103002, 12.5337814, -12.9703407, -1.34844796}},
        {{-22.8161443, 0.117673962, -8.28246978, -0.0550618805, -1.03112625, -0.276812407, -0.655359458, -0.113089643, -1.64243916, 0.903341624, 2.02618168}},
        {{46.0238053, -1.03247358, 53.3010628, 1.28528254, 0.95589098, 7.7129662, 0.243762408, 3.24099827, 6.35801226, -6.67783509, 1.07258735}},
        {{-220.757596, 0.601531387, 159.299535, 20.1331386, 1.4433446, 9.94463595, -1.94016216, 2.04362338, -2.02621495, 0.624223384, 20.7184522}},
        {{2.1895905, -0.0220743271, -2.3616949, 4.23488886, 0.133154711, 0.0760974817, -0.0425477008, 0.48169314, -0.450484896, 0.235159504, 0.0466013472}},
        {{-67.3702734, -5.55135298, 118.372597, 31.2725499, 15.9272668, 27.4135085, 0.736790677, 5.14122361, -0.675743542, -6.33049576, 5.64512711}},
        {{-237.505867, -2.86972373, 56.7857578, -33.2375141, 40.2539155, 34.5667639, 5.46415768, -1.2833648, -25.1152345, 26.1805796, 36.5102306}},
        {{47.6058341, -1.37990957, 21.3138808, 13.673366, 3.27779728, 1.59484917, 1.09137974, 1.63536938, 2.7652016, -2.34145965, -4.49713377}},
        {{-205.551716, -9.57872795, 167.945622, 58.3190342, 27.0982504, 43.1096975, -0.342789904, 0.782167647, -11.7396534, -1.78936058, 9.78014121}},
        {{102.455216, -17.8554242, 223.475468, 48.7597718, 55.2110555, 90.5131974, 7.74234185, 15.4816798, 9.40714123, -12.8341566, 20.6813208}},
        {{48.4087705, 2.36968221, -42.2166756, -7.69345468, -0.757956488, -1.10759131, 0.460216172, -0.859777726, 0.424379737, 0.450926179, -0.376210051}},
        {{-1821.41651, -7.53271216, 2.86458845, 95.1275785, -25.1904556, 76.0944249, 12.2298355, -22.1885793, -82.3889056, 68.914216, 39.9741147}},
        {{121.515563, 4.47745849, -52.634555, -10.4270448, 0.0887140111, -1.39346283, 1.18185809, -0.755725497, 3.30133512, -1.45042036, 1.81085017}},
        {{1110.09902, 2.8775261, -99.7560191, -34.4699784, -11.8942065, -77.3509544, -47.0205881, 32.0388844, 28.3393378, -11.3354338, 5.21189799}},
        {{54.3233916, -1.93622396, 91.5246651, 35.9053219, 8.04718349, 5.28540062, 0.435294691, 3.21685346, 3.81752108, -5.47734349, 1.58138408}},
        {{3955.77569, 33.7141087, -1047.04875, -467.767838, 20.3275567, -272.339668, -92.7042202, 109.683234, 147.697927, -111.638998, -25.5219824}},
        {{-266.914816, -1.02749127, 38.6005844, 27.4535685, -4.31456479, 1.79193802, -7.51115253, 0.789230983, -15.4311682, 7.74310574, 17.7817911}},
        {{-777.494595, 93.9963343, -63.1580231, 279.304029, -113.877003, -478.298856, 347.99832, -208.77513, 0.826811509, 44.5663535, -434.832626}},
        {{-102.135724, -5.22499608, -105.09591, 7.09921883, 1.78169311, 6.36586278, 1.87244346, 1.15691218, -3.75507855, -2.30937433, -14.5629376}},
        {{460.331742, 14.205242, -405.047675, -95.3487776, -47.8929466, -71.4426917, -35.7064787, -11.681722, 1.04434429, 5.45180224, -4.25531853}},
        {{-33.6105358, -2.53882449, 35.5125228, 13.4128246, 1.3965978, 1.63155021, -0.560028822, 0.993256784, -0.390915166, -0.344505068, -1.76520383}},
        {{237.541788, 3.53496454, 59.9633592, 175.985958, -6.63411024, -24.5621624, 37.3528369, -11.3781725, 28.1027846, -14.2646258, -65.0803451}},
        {{-48.811414, 0.0342615774, -49.279182, -16.0952653, -5.19135695, -2.60296886, -1.08546297, -1.21070579, -1.7945508, 2.38021159, 0.795783358}},
        {{752.243443, -2.7568422, -174.385298, 85.4519376, 39.8300943, 16.9396954, 38.9471775, -21.997966, 51.4987068, -46.1891411, -36.3595557}},
        {{-55.1811304, -6.41044606, 48.4442047, 13.7390433, 3.04860765, 2.73499831, 1.36379683, 0.940436748, 1.00013343, -1.35912433, -7.95496734}},
        {{-663.689781, 8.83549266, -1510.729, -505.134274, 21.2198615, 99.961542, -45.834243, 8.33798354, -52.3073461, 4.84574456, 169.583759}},
        {{-83.8263686, -3.51402955, -11.85708, 1.62933582, -0.898186854, -0.305104187, 0.210343985, -1.96648648, -5.88071309, 4.89325171, -4.7703634}},
        {{-5155.40532, -151.170614, 5457.78446, 1246.76235, 15.0124896, 722.970225, 69.1135705, 2.56605133, -260.155878, 170.795887, 232.201759}},
        {{1164.20982, 38.1144802, -450.782577, -109.950576, -7.10028113, -23.346259, 9.64823453, -7.80181242, 28.8411122, -9.65816123, 0.834791238}},
        {{0.54707965, -0.0709196309, -0.192144999, 0.0196082094, -0.00689000506, 0.0248082556, -0.00148281498, 0.000694914763, 0.00716256822, -0.00634705674, -0.180139141}},
        {{-87.6061948, -2.71426841, 53.1962007, 2.10596333, 4.48178758, 10.4895869, -6.01774846, 2.85095472, -1.37388303, -0.694338982, 7.10616855}},
        {{-0.41645912, 0.172037694, 0.17251899, -0.0783718327, 0.0171566314, -0.0776208192, 0.011048104, -0.0113161853, -0.0700584199, 0.0184645225, 0.800193249}},
        {{345.844815, 8.93497443, -185.994326, 5.37188876, -12.7594352, -40.7150069, 24.0296112, -8.68589431, 3.04118532, 1.62696754, -25.7745376}},
        {{6.48436911, -0.156996234, -0.397951031, -0.813001939, -0.135163453, 0.000878424928, -0.16075951, -0.2752249, -0.00417722185, 0.176967707, 0.322650959}},
        {{5.21740849, 0.426956323, 2.42854096, -1.43314202, 0.285392665, -1.30591589, 1.91297168, 0.62034503, 0.607449311, -0.73394075, -1.31583111}},
        {{2.27334133, 0.0953290889, -1.16409092, -0.774708152, 0.195310239, 0.146930227, 0.0752642374, 0.156568733, 0.394514836, -0.344025557, -0.138963289}},
        {{-82.1744672, 0.571223123, -5.072884, -4.94160552, -1.48213827, -0.490255722, -0.34967087, -0.375770688, -2.53381759, 0.993676929, 1.33718068}},
        {{-164.162193, 0.615947514, -15.2397881, -22.9790621, -1.72271531, -0.415383265, -3.82613476, -0.267309007, -7.80203169, 4.61128978, 7.63272567}},
        {{-16.9756106, -0.9299415, 66.0757491, 16.2500322, 0.471561971, 1.6639687, -0.0911557654, -0.221257173, 0.120432273, -0.11960142, -1.33728468}},
        {{283.676312, -4.88280584, 258.132979, 50.7853839, 11.6723713, 8.83248426, 0.486103838, 3.68559882, 11.4030513, -7.13696584, -3.85607313}},
        {{-4.45090928, 0.117031783, -0.306644149, 0.254976286, -0.0272959345, 0.00652623198, -0.127179767, 0.0847144472, -0.213640933, 0.135307189, -0.29317723}},
        {{21.9378148, 0.147185941, -4.52052253, -4.06831449, 0.294363642, -0.543301993, 0.670139278, -0.0238058747, 1.14446601, -0.0443480631, -0.423261787}},
        {{-33.2153461, 0.229620227, -5.00148163, -0.188261157, -0.910943782, -0.45552498, -0.22100783, -0.032788878, -1.94796441, 1.01378362, 0.946739014}},
        {{-48.4744713, 0.670578757, -16.4495417, -3.86517962, -2.37043156, -0.931711293, -0.925283884, -1.29223589, -3.06275762, 1.49913295, 1.28908777}},
        {{29.9403606, -2.19896747, 124.45185, 33.2970441, 4.42516659, 5.00020427, 1.34053627, 1.21589448, 5.76364103, -5.44396314, -3.96182685}},
        {{2.32273873, -0.417114302, 11.9502214, -30.151817, 4.63545716, 2.5274857, -4.81833905, 4.17535806, -2.84051837, 1.97379869, 8.97388665}},
        {{-42.4818537, 1.20717193, 96.8399605, 3.94149055, -6.11891038, -5.66720897, 4.55871089, -7.13579081, 0.940474954, 5.02912492, -6.69411876}},
        {{9.48729742, -0.697392861, 28.619451, 8.05842436, 1.69152097, 1.87153199, 0.391213915, 1.11492383, 1.27804287, -1.55734374, -1.60792066}},
        {{0.482682179, 0.0060278457, -0.32464389, -0.103807962, 9.29087109e-05, -0.0126299062, 0.0139294892, -0.00754978483, 0.0306106936, 0.0105205534, -0.0723434375}},
        {{-0.429457583, -0.0131063521, -0.00157655041, 0.335217974, 0.00827991995, 0.0131539878, 0.0228756451, -0.00955415981, -0.00255471639, 0.0151219701, -0.179581524}},
        {{-0.394103911, -0.00632833818, -0.00936696837, 0.293050002, -0.0137998326, 0.0291475757, -0.0333945202, -0.000120462745, -0.0440169399, 0.0112289492, 0.186059808}},
        {{1.09846865, 0.108321438, -0.440890733, -0.965631233, -0.0351673028, -0.0994728398, -0.126133317, 0.0656639244, -0.0234581512, -0.021075252, 0.543719345}},
        {{0.970192274, -0.0641654618, 3.06073864, -0.362014743, 0.22077352, 0.0944052699, 0.190665561, 0.0112356839, 0.118780685, -0.218334621, 0.208109671}},
        {{-21.5741787, -0.549216449, 10.6709415, 0.437347953, -0.295358257, -0.20724512, 0.174945983, -0.0279813967, -0.307349295, 0.183762395, -0.35615989}},
        {{-0.649155464, -0.0117079422, 0.194089695, 0.181221767, 0.00644539366, 0.000541448024, -0.00254999349, 0.0179939059, -0.0408618152, -0.0224096185, 0.104522447}},
        {{-0.42329603, -0.0860646243, 0.815194889, 0.110776566, -0.0233709284, -0.0171275786, -0.0566446021, 0.0314165493, -0.00236967602, 0.0297248453, -0.137056832}},
        {{1.87574892, 0.200894756, -1.71753726, -1.23912526, 0.00273312444, -0.0566700291, 0.179239976, -0.153596574, 0.313406796, 0.147176675, 0.0978955879}},
        {{-0.111433133, 0.00232698482, 0.0498289618, 0.0433590339, 0.0189465093, 0.00432069397, -0.00620642995, 0.485450712, -0.430108081, -0.437859948, 0.00364042899}},
        {{-2.10842126, 0.0866885224, -4.843243, -1.65762886, -0.187689661, -0.188182851, 0.0716535365, -0.152473548, -0.480513168, 0.337533748, 0.111083926}},
    }};

// -------------------------------------------------------------------------------------------
// Driver features and amplitudes
// -------------------------------------------------------------------------------------------

/**
 * The eleven smooth driver features every mode amplitude is linear in.
 *
 * Chosen from the driver responses measured on the oracle before anything was fitted (the probe
 * is in the harness's history; the numbers are at 6.6 R_E noon/midnight): the G2 response is
 * linear to four digits, so G2 enters linearly; the G3 response is strongly super-linear (the tail
 * current moves earthward and the midnight `B_z` changes SIGN between G3 = 40 and 100), so G3 enters
 * with its square; the Pdyn response is sub-linear, so `sqrt(Pdyn)` and `Pdyn` both enter; the
 * Dst response is linear with a pressure-dependent slope, so `Dst` and `Dst sqrt(Pdyn)` enter (the
 * Burton et al. 1975 pressure correction of Dst is of that form); By enters linearly, as a
 * near-uniform penetration field of `0.32 By`; Bz acts through its SOUTHWARD part — the response
 * per nT is 0.64 southward and 0.06 northward — so the smooth rectifier
 * `Bs = (sqrt(Bz^2 + 1) - Bz) / 2` enters beside `Bz` itself.
 *
 * `sqrt(Pdyn)` is evaluated as `sqrt(max(Pdyn, 0))`: a negative pressure is unphysical, is not a
 * value any model was fitted to, and must not turn a whole batch into NaN.
 *
 * @param d the drivers.
 * @return `{1, Dst, sqrt(Pdyn), Pdyn, G2, G3, G3^2/100, By, Bs, Bz, Dst sqrt(Pdyn)/10}`. The two
 *         scalings keep the coefficient table's entries within a few orders of magnitude of each
 *         other; they are part of the feature definition, not of the table.
 * @complexity O(1) — two square roots.
 * @alloc none.
 * @test IrbemT01s.FeaturesAreSmoothAndRectifyBzSouthward
 */
[[nodiscard]] inline std::array<double, t01s_feature_count> t01s_features(const T01sDrivers& d) {
    const double sqrtp = std::sqrt(d.pdyn > 0.0 ? d.pdyn : 0.0);
    const double bs = 0.5 * (std::sqrt((d.bz_imf * d.bz_imf) + 1.0) - d.bz_imf);
    return {1.0,
            d.dst,
            sqrtp,
            d.pdyn,
            d.g2,
            d.g3,
            d.g3 * d.g3 * 0.01,
            d.by_imf,
            bs,
            d.bz_imf,
            d.dst * sqrtp * 0.1};
}

/**
 * The 108 mode amplitudes for a driver state.
 *
 * Always accumulated in `double` and then converted, because the coefficient table is `double`
 * and the conversion to `float` is the one place the device lane's inputs are rounded: the host
 * fp32 lane and the kernel receive the SAME rounded amplitudes, which is what makes a disagreement
 * between them attributable to the device.
 *
 * @tparam T the scalar type the amplitudes are wanted in; `double` for the reference lane, `float`
 *         for the device-mirroring one.
 * @param d the drivers.
 * @return one amplitude per mode, in @ref t01s_coefficients' row order.
 * @complexity O(108 x 11) — 1188 multiply-adds, comparable to one evaluation of the modes.
 * @alloc none.
 * @test IrbemT01s.AmplitudesAreLinearInTheFeatures
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, t01s_mode_count> t01s_amplitudes(const T01sDrivers& d) {
    const std::array<double, t01s_feature_count> f = t01s_features(d);
    std::array<T, t01s_mode_count> a{};
    for (std::size_t m = 0; m < t01s_mode_count; ++m) {
        double s = 0.0;
        for (std::size_t k = 0; k < t01s_feature_count; ++k) s += t01s_coefficients[m][k] * f[k];
        a[m] = static_cast<T>(s);
    }
    return a;
}

// -------------------------------------------------------------------------------------------
// The spatial modes
// -------------------------------------------------------------------------------------------

/**
 * Rotate an SM-frame vector to GSM: a rotation about `y` by the dipole tilt.
 *
 * @tparam T the scalar type.
 * @param bx_sm the SM x component. @param by_sm the SM y component. @param bz_sm the SM z
 *        component. @param sin_tilt `sin(psi)`. @param cos_tilt `cos(psi)`.
 * @return the same vector in GSM components.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01s.ZeroTiltIsMirrorSymmetricAboutTheEquator
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, 3> rotate_sm_to_gsm(T bx_sm, T by_sm, T bz_sm, T sin_tilt,
                                                       T cos_tilt) {
    return {(bx_sm * cos_tilt) + (bz_sm * sin_tilt), by_sm,
            (-bx_sm * sin_tilt) + (bz_sm * cos_tilt)};
}

/**
 * The 108 divergence-free mode fields at one GSM point, before any amplitude is applied.
 *
 * This is the spatial half of the model, and the half the divergence test exercises. Three
 * families, each exactly solenoidal by construction:
 *
 * **Discs** (in SM, where the dipole axis is `z`, then rotated). With `P(rho^2, xi)` one of the
 * three disc potentials, `Q = W P h` the potential over radius (`W` the truncation, `h` the
 * azimuthal weight) and `A = Q (-y, x, 0)`, the curl gives `B_x = -x dQ/dz`, `B_y = -y dQ/dz`,
 * `B_z = 2Q + x dQ/dx + y dQ/dy`, where every derivative is TOTAL — through `rho`, through the
 * sheet surface `z_s(x, y)`, through the thickness `D(y)` and through `W(x, y)`. With
 * `G_z = -W P_xi z_r / xi`, `base = W (2P + 2 rho^2 P_rho2) + (x W_x + y W_y) P + W P_xi D y D_y / xi`
 * and `n` the degree of `h`, that is `B_x = x h G_z`, `B_y = y h G_z`,
 * `B_z = h (base + n W P) + B_x z_s,x + B_y z_s,y` — the T89 eqs. (14)-(17) form generalised to a
 * weighted potential. Dropping any one derivative breaks `div B = 0`, which is why the finite-
 * difference divergence test is the main test this file carries.
 *
 * **Sheets** (in SM): `B = curl(Lambda z_hat) = (dLambda/dy, -dLambda/dx, 0)`, with the sheet
 * warp entering through `H(z_r)` and its derivative.
 *
 * **Harmonics** (in GSM directly): `B = grad gamma` for the four box-harmonic kinds; the
 * tilt-odd kinds carry `sin(psi)` as a factor.
 *
 * The geometry is a parameter rather than the global so that the calibration harness can refine
 * it; every production caller passes @ref t01s_geometry.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device kernel.
 * @param geo the sheet geometry; @ref t01s_geometry in production.
 * @param sin_tilt `sin(psi)`; a property of the epoch, paid once per timestamp.
 * @param cos_tilt `cos(psi)`; must be non-zero — the sheet surface carries `tan(psi)` — which
 *        @ref t01s_check enforces before anything gets here.
 * @param x the GSM x coordinate, R_E. @param y the GSM y, R_E. @param z the GSM z, R_E.
 * @param out receives the 108 fields, `{B_x, B_y, B_z}` each, in GSM, in the table order.
 * @complexity O(1) — ~3 700 flops, 24 `exp`, 24 `sin`/`cos` pairs and ~130 square roots, with no
 *             data-dependent branch: the table loops have constant trip counts and the `family`
 *             and `azimuth` selects are on constants.
 * @alloc none.
 * @test IrbemT01s.EveryModeIsDivergenceFree
 * @test IrbemT01s.DivergenceVanishesEverywhere
 * @test IrbemT01s.ZeroTiltIsMirrorSymmetricAboutTheEquator
 */
template <std::floating_point T>
inline void t01s_mode_fields(const T01sGeometry& geo, T sin_tilt, T cos_tilt, T x, T y, T z,
                             std::array<std::array<T, 3>, t01s_mode_count>& out) {
    const T half = static_cast<T>(0.5);
    const T one = static_cast<T>(1);
    const T two = static_cast<T>(2);
    const T three = static_cast<T>(3);

    const T tan_tilt = sin_tilt / cos_tilt;

    // ---- GSM -> SM for the current-carrying modes -------------------------------------------
    const T xs = (x * cos_tilt) - (z * sin_tilt);
    const T ys = y;
    const T zs = (x * sin_tilt) + (z * cos_tilt);

    // ---- the warped, hinged current sheet (T89 eq. 11) and its slopes -----------------------
    const T hinge = xs + static_cast<T>(geo.r_hinge);
    const T hinge_root = std::sqrt((hinge * hinge) + static_cast<T>(geo.hinge2));
    const T y2 = ys * ys;
    const T y4 = y2 * y2;
    const T ly2 = static_cast<T>(geo.l_y) * static_cast<T>(geo.l_y);
    const T ly4 = ly2 * ly2;
    const T bend_den = y4 + ly4;
    const T warp_g = static_cast<T>(geo.warp_g);
    const T sheet = (half * tan_tilt * (hinge - hinge_root)) - (warp_g * sin_tilt * y4 / bend_den);
    const T dsheet_dx = half * tan_tilt * (one - (hinge / hinge_root));
    const T dsheet_dy =
        -warp_g * sin_tilt * static_cast<T>(4) * ys * y2 * ly4 / (bend_den * bend_den);

    const T zr = zs - sheet;
    const T rho2 = (xs * xs) + (ys * ys);

    // ---- the tail truncation factor W(x, y) (T89 eqs. 12-13) and x W_x + y W_y ---------------
    const T xw = xs - static_cast<T>(geo.x0_w);
    const T dxw2 = static_cast<T>(geo.dx_w) * static_cast<T>(geo.dx_w);
    const T w_root = std::sqrt((xw * xw) + dxw2);
    const T w_x = half * (one - (xw / w_root));
    const T dyw2 = static_cast<T>(geo.dy_w) * static_cast<T>(geo.dy_w);
    const T w_y = one / (one + (y2 / dyw2));
    const T w_trunc = w_x * w_y;
    const T dw_dx = w_y * (-half * dxw2 / (w_root * w_root * w_root));
    const T dw_dy = w_x * (-two * ys / dyw2) * w_y * w_y;
    const T grad_w_trunc = (xs * dw_dx) + (ys * dw_dy);

    // ---- the disc modes -----------------------------------------------------------------------
    const T delta_y = static_cast<T>(geo.delta_y);
    std::size_t slot = 0;
    for (const T01sDiscMode& m : t01s_disc_modes) {
        const T w = m.truncated ? w_trunc : one;
        const T grad_w = m.truncated ? grad_w_trunc : static_cast<T>(0);
        // Truncated (tail) modes thicken towards the flanks; ring modes have a fixed thickness.
        const T d = static_cast<T>(m.d) + (m.truncated ? delta_y * y2 : static_cast<T>(0));
        const T y_dd_dy = m.truncated ? (ys * two * delta_y * ys) : static_cast<T>(0);
        const T xi = std::sqrt((zr * zr) + (d * d));
        const T u = static_cast<T>(m.a) + xi;
        const T s2 = rho2 + (u * u);
        const T s = std::sqrt(s2);
        T p;        // P
        T p_xi;     // dP/dxi
        T p_base;   // 2P + 2 rho^2 dP/drho^2
        if (m.family == 1) {
            const T su = s + u;
            p = one / su;
            p_xi = -one / (s * su);
            p_base = one / s;
        } else if (m.family == 2) {
            const T su = s + u;
            const T s3 = s2 * s;
            p = -one / (s * su);
            p_xi = one / s3;
            p_base = -u / s3;
        } else {
            const T s5 = s2 * s2 * s;
            p = one / (s2 * s);
            p_xi = -three * u / s5;
            p_base = ((two * u * u) - rho2) / s5;
        }
        const T g_z = -w * p_xi * zr / xi;
        const T base = (w * p_base) + (grad_w * p) + (w * p_xi * d * y_dd_dy / xi);
        T h;       // the homogeneous azimuthal weight
        T degree;  // its degree, so x h_x + y h_y = degree * h
        if (m.azimuth == 0) {
            h = one;
            degree = static_cast<T>(0);
        } else if (m.azimuth == 1) {
            h = ys;
            degree = one;
        } else if (m.azimuth == 2) {
            h = xs;
            degree = one;
        } else if (m.azimuth == 3) {
            h = (xs * xs) - (ys * ys);
            degree = two;
        } else {
            h = xs * ys;
            degree = two;
        }
        const T bx = xs * h * g_z;
        const T by = ys * h * g_z;
        T bz = h * (base + (degree * w * p));
        bz += (bx * dsheet_dx) + (by * dsheet_dy);
        out[slot++] = rotate_sm_to_gsm(bx, by, bz, sin_tilt, cos_tilt);
    }

    // ---- the radial-current sheet modes, in SM ---------------------------------------------------
    const T zr2 = zr * zr;
    for (const T01sSheetMode& m : t01s_sheet_modes) {
        const T a2 = static_cast<T>(m.a) * static_cast<T>(m.a);
        const T d2 = static_cast<T>(m.d) * static_cast<T>(m.d);
        const T ra = rho2 + a2;
        const T ra_root = std::sqrt(ra);
        const T g = one / (ra * ra_root);                 // G(rho)
        const T gp = -three / (ra * ra * ra_root);        // (dG/drho) / rho
        const T hz_root = std::sqrt(zr2 + d2);
        const T hz = zr / hz_root;                        // H(z_r)
        const T hzp = d2 / (hz_root * hz_root * hz_root); // dH/dz_r
        T q;
        T dq_dx;
        T dq_dy;
        if (m.azimuth == 1) {
            q = ys;
            dq_dx = static_cast<T>(0);
            dq_dy = one;
        } else if (m.azimuth == 2) {
            q = xs;
            dq_dx = one;
            dq_dy = static_cast<T>(0);
        } else if (m.azimuth == 3) {
            q = (xs * xs) - (ys * ys);
            dq_dx = two * xs;
            dq_dy = -two * ys;
        } else {
            q = xs * ys;
            dq_dx = ys;
            dq_dy = xs;
        }
        const T dl_dx = (q * ((gp * xs * hz) - (g * hzp * dsheet_dx))) + (g * hz * dq_dx);
        const T dl_dy = (q * ((gp * ys * hz) - (g * hzp * dsheet_dy))) + (g * hz * dq_dy);
        out[slot++] = rotate_sm_to_gsm(dl_dy, -dl_dx, static_cast<T>(0), sin_tilt, cos_tilt);
    }

    // ---- the harmonic modes, in GSM -------------------------------------------------------------
    for (const T01sHarmonicMode& m : t01s_harmonic_modes) {
        const T p = static_cast<T>(m.p);
        const T r = static_cast<T>(m.r);
        const T sc = one / std::sqrt((one / (p * p)) + (one / (r * r)));
        const T e = std::exp(x / sc);
        const T sy = std::sin(y / p);
        const T cy = std::cos(y / p);
        const T sz = std::sin(z / r);
        const T cz = std::cos(z / r);
        std::array<T, 3> b{};
        if (m.kind == 0) {
            b = {e * cy * sz / sc, -e * sy * sz / p, e * cy * cz / r};
        } else if (m.kind == 1) {
            b = {sin_tilt * e * cy * cz / sc, -sin_tilt * e * sy * cz / p,
                 -sin_tilt * e * cy * sz / r};
        } else if (m.kind == 2) {
            b = {e * sy * sz / sc, e * cy * sz / p, e * sy * cz / r};
        } else {
            b = {sin_tilt * e * sy * cz / sc, sin_tilt * e * cy * cz / p,
                 -sin_tilt * e * sy * sz / r};
        }
        out[slot++] = b;
    }

    // ---- the uniform penetration fields ---------------------------------------------------------
    out[slot++] = {static_cast<T>(0), one, static_cast<T>(0)};
    out[slot++] = {static_cast<T>(0), static_cast<T>(0), one};
}

/**
 * The T01S external field at one GSM point, as three components in nanotesla, from precomputed
 * amplitudes.
 *
 * The whole model is this: the 108 mode fields of @ref t01s_mode_fields, each weighted by its
 * amplitude. Written against an amplitude array rather than the drivers so that a caller
 * evaluating many points at one driver state — a batch, a trace — pays for @ref t01s_amplitudes
 * once.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device kernel.
 * @param amp the mode amplitudes; see @ref t01s_amplitudes.
 * @param sin_tilt `sin(psi)`. @param cos_tilt `cos(psi)`, non-zero.
 * @param x the GSM x coordinate, R_E. @param y the GSM y, R_E. @param z the GSM z, R_E.
 * @return `{B_x, B_y, B_z}` in GSM, nanotesla.
 * @complexity O(1); one @ref t01s_mode_fields plus 324 multiply-adds.
 * @alloc none — the mode block is a stack array.
 * @test IrbemT01s.ReferenceLaneMatchesTheComponentForm
 * @test IrbemT01s.HostFloatLaneTracksTheReferenceLane
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, 3> t01s_components(
    const std::array<T, t01s_mode_count>& amp, T sin_tilt, T cos_tilt, T x, T y, T z) {
    std::array<std::array<T, 3>, t01s_mode_count> modes{};
    t01s_mode_fields<T>(t01s_geometry, sin_tilt, cos_tilt, x, y, z, modes);
    T bx = static_cast<T>(0);
    T by = static_cast<T>(0);
    T bz = static_cast<T>(0);
    for (std::size_t m = 0; m < t01s_mode_count; ++m) {
        bx += amp[m] * modes[m][0];
        by += amp[m] * modes[m][1];
        bz += amp[m] * modes[m][2];
    }
    return {bx, by, bz};
}

// -------------------------------------------------------------------------------------------
// The scalar entry points
// -------------------------------------------------------------------------------------------

/**
 * The T01S external field at one GSM point, in `double`, from precomputed amplitudes — the
 * reference lane.
 *
 * @param p the position, GSM, in Earth radii.
 * @param sin_tilt `sin(psi)`; @ref HotState::sin_tilt holds it, precomputed per epoch.
 * @param cos_tilt `cos(psi)`; must be non-zero, which @ref t01s_check enforces.
 * @param amp the mode amplitudes for the epoch's drivers; see @ref t01s_amplitudes.
 * @return the external field at @p p, GSM, in nanotesla.
 * @complexity O(1); see @ref t01s_components.
 * @alloc none.
 * @test IrbemT01s.ReferenceLaneMatchesTheComponentForm
 */
[[nodiscard]] inline FieldVector<Frame::GSM> t01s_field_at(
    Position<Frame::GSM> p, double sin_tilt, double cos_tilt,
    const std::array<double, t01s_mode_count>& amp) {
    const std::array<double, 3> b =
        t01s_components<double>(amp, sin_tilt, cos_tilt, p.v[0], p.v[1], p.v[2]);
    return FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
}

/**
 * The model's verdict on a point and a driver state, without evaluating the field.
 *
 * Two kinds of refusal, kept apart because a caller acts on them differently. Arithmetic that has
 * no answer is @ref Status::DomainError: a non-finite input, a radius inside the Earth, a tilt of
 * exactly ±90 degrees (the sheet surface carries `tan(psi)`), or a position sunward of
 * @ref t01s_max_x_gsm, where the shielding exponentials overflow the fp32 lane. Everything else is
 * a value with, at most, a caveat — and for THIS model the caveat is narrow by publication: the
 * IRBEM table states "no upper or lower limit" for its six drivers and one spatial bound,
 * `x_GSM >= -15 R_E`. So @ref Status::OutOfValidityRange fires for the tail beyond `x = -15` and
 * for NOTHING else: a Dst of -1000 nT is `Ok` here because the model's publication says so, not
 * because it is a good idea. Both rules live in `status.hpp`'s envelope table; this function only
 * asks it.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians.
 * @param d the drivers.
 * @return @ref Status::DomainError, @ref Status::OutOfValidityRange or @ref Status::Ok as above.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01s.OutOfRangePositionIsReportedButStillEvaluated
 * @test IrbemT01s.DriversHaveNoPublishedLimits
 * @test IrbemT01s.NonFiniteInputIsADomainError
 * @test IrbemT01s.RightAngleTiltIsADomainError
 * @test IrbemT01s.FarSunwardExtrapolationIsADomainErrorNotANaN
 */
[[nodiscard]] inline Status t01s_check(Position<Frame::GSM> p, double tilt_rad,
                                       const T01sDrivers& d) {
    if (!std::isfinite(p.v[0]) || !std::isfinite(p.v[1]) || !std::isfinite(p.v[2]) ||
        !std::isfinite(tilt_rad)) {
        return Status::DomainError;
    }
    if (!(std::fabs(tilt_rad) < max_tilt_rad)) return Status::DomainError;
    if (p.v[0] > t01s_max_x_gsm) return Status::DomainError;

    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Dst)] = d.dst;
    drivers[static_cast<std::size_t>(Driver::Pdyn)] = d.pdyn;
    drivers[static_cast<std::size_t>(Driver::ByIMF)] = d.by_imf;
    drivers[static_cast<std::size_t>(Driver::BzIMF)] = d.bz_imf;
    drivers[static_cast<std::size_t>(Driver::G2)] = d.g2;
    drivers[static_cast<std::size_t>(Driver::G3)] = d.g3;
    const Status drives = check_validity(ExternalModel::Tsyganenko2001Storm, drivers);
    if (drives == Status::DomainError) return Status::DomainError;

    const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
    const Status where = check_position(ExternalModel::Tsyganenko2001Storm, r, p.v[0]);
    return first_failure(where, drives);
}

/**
 * The T01S external field, with the model's own verdict on whether it should be believed here.
 *
 * The value is **always** returned, including when the status is @ref Status::OutOfValidityRange
 * — `status.hpp`'s standing rule. What is refused outright, with a zero field, is what
 * @ref t01s_check refuses, plus one check on the OUTPUT: a non-finite component, which no input
 * @ref t01s_check accepts can produce but which this function will not let past regardless, because
 * a NaN that escapes here surfaces a hundred RK4 steps later with nothing left pointing at its
 * cause.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians; positive when the north dipole leans sunward.
 * @param d the drivers.
 * @return the field and its caveat. @ref Status::DomainError (with a zero field) per
 *         @ref t01s_check or for a non-finite result; @ref Status::OutOfValidityRange for a
 *         position tailward of `x = -15 R_E`, with the field still computed; otherwise
 *         @ref Status::Ok.
 * @complexity O(1) — one amplitude solve and one evaluation, ~5 000 flops.
 * @alloc none.
 * @test IrbemT01s.OutOfRangePositionIsReportedButStillEvaluated
 * @test IrbemT01s.NonFiniteInputIsADomainError
 * @test IrbemT01s.FarSunwardExtrapolationIsADomainErrorNotANaN
 * @test IrbemT01s.OracleSamplesLandInsideTheDocumentedGap
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> t01s_field(Position<Frame::GSM> p,
                                                                double tilt_rad,
                                                                const T01sDrivers& d) {
    const FieldVector<Frame::GSM> zero{};
    const Status verdict = t01s_check(p, tilt_rad, d);
    if (verdict == Status::DomainError) return {Status::DomainError, zero};
    const FieldVector<Frame::GSM> b =
        t01s_field_at(p, std::sin(tilt_rad), std::cos(tilt_rad), t01s_amplitudes<double>(d));
    if (!std::isfinite(b.v[0]) || !std::isfinite(b.v[1]) || !std::isfinite(b.v[2])) {
        return {Status::DomainError, zero};
    }
    return {verdict, b};
}

/**
 * The T01S external field for a whole epoch's worth of state — the production entry point.
 *
 * Dst, Pdyn, By and Bz come from the @ref HotState cache line; G2 and G3 are not on it (it was
 * laid out for the pre-2001 models) and are read from the cold driver vector, which costs this
 * overload one cache miss per call that the batch lanes do not pay.
 *
 * @param p the position, GSM, in Earth radii.
 * @param ctx the epoch's context.
 * @return the field and its caveat, exactly as the three-argument overload.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01s.ContextOverloadAgreesWithTheExplicitOne
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> t01s_field(Position<Frame::GSM> p,
                                                                const FieldContext& ctx) {
    return t01s_field(p, ctx.hot().tilt_rad, t01s_drivers(ctx.drivers()));
}

// -------------------------------------------------------------------------------------------
// The batch lanes
// -------------------------------------------------------------------------------------------

/**
 * The T01S field over a whole batch, on the CPU, in `float`.
 *
 * The host twin of `irbem_t01s_f32`: the same expressions, in the same order, in the same
 * precision, from the same `float`-rounded amplitudes. That is what makes a disagreement between
 * the two lanes attributable to the device rather than to the arithmetic having been written
 * differently on the two sides. It is also the lane a machine with no GPU actually runs.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, GSM, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, nanotesla; overwritten in full.
 * @param sin_tilt `sin(psi)`. @param cos_tilt `cos(psi)`; must be non-zero.
 * @param amp the mode amplitudes, already rounded to `float`; see @ref t01s_amplitudes.
 * @return `false` when @p pos is not a whole number of points or @p out is a different length, in
 *         which case nothing is written; `true` otherwise.
 * @complexity O(N).
 * @alloc none.
 * @test IrbemT01s.HostFloatLaneTracksTheReferenceLane
 * @test IrbemT01s.HostFloatLaneRejectsMismatchedSpans
 */
[[nodiscard]] inline bool t01s_field_host(std::span<const float> pos, std::span<float> out,
                                          float sin_tilt, float cos_tilt,
                                          const std::array<float, t01s_mode_count>& amp) {
    if (pos.size() % 3 != 0 || out.size() != pos.size()) return false;
    const std::size_t n = pos.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const std::array<float, 3> b = t01s_components<float>(
            amp, sin_tilt, cos_tilt, pos[(3 * i) + 0], pos[(3 * i) + 1], pos[(3 * i) + 2]);
        out[(3 * i) + 0] = b[0];
        out[(3 * i) + 1] = b[1];
        out[(3 * i) + 2] = b[2];
    }
    return true;
}

/// How many `float` scalars the device kernel's parameter buffer holds: `sin(psi)`, `cos(psi)`,
/// then the 108 amplitudes. Asserted against the kernel registry by the GPU suite.
inline constexpr std::size_t t01s_param_count = 2 + t01s_mode_count;

/**
 * Pack the epoch's tilt and the driver state's amplitudes into the kernel's parameter buffer.
 *
 * The layout is the kernel's ABI and is stated in exactly two places — here and the comment above
 * `irbem_t01s_f32` in `irbem.slang`. A test evaluates both lanes on the same points, which is what
 * actually keeps the two statements in step.
 *
 * @param sin_tilt `sin(psi)`. @param cos_tilt `cos(psi)`.
 * @param d the drivers.
 * @return the parameter block, `t01s_param_count` floats, by value.
 * @complexity O(1) — one amplitude solve.
 * @alloc none — the block is the returned object's own inline array.
 * @test IrbemT01s.ParameterBlockCarriesTheTiltThenTheAmplitudes
 */
[[nodiscard]] inline std::array<float, t01s_param_count> t01s_param_block(float sin_tilt,
                                                                          float cos_tilt,
                                                                          const T01sDrivers& d) {
    const std::array<float, t01s_mode_count> amp = t01s_amplitudes<float>(d);
    std::array<float, t01s_param_count> block{};
    block[0] = sin_tilt;
    block[1] = cos_tilt;
    for (std::size_t m = 0; m < t01s_mode_count; ++m) block[2 + m] = amp[m];
    return block;
}

/**
 * The batch's position caveat, accumulated one point at a time.
 *
 * A batch returns ONE @ref Status for N points, so it can only be the worst of them, and computing
 * that must not cost a second pass over 100 MB of positions. This model's envelope is monotone in
 * three quantities — the smallest radius (inside the Earth is a domain error), the smallest `x`
 * (tailward of -15 R_E is out of validity) and the largest `x` (sunward of @ref t01s_max_x_gsm is
 * a domain error) — so three running extremes decide the whole batch, and the radius is compared
 * as a SQUARE so there is no per-point `sqrt`.
 *
 * @test IrbemT01s.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
struct T01sPositionFold {
    /// The smallest `r^2` seen, R_E^2; `+inf` until the first point.
    double r2_lo = std::numeric_limits<double>::infinity();
    /// The smallest GSM `x` seen, R_E; `+inf` until the first point.
    double x_lo = std::numeric_limits<double>::infinity();
    /// The largest GSM `x` seen, R_E; `-inf` until the first point.
    double x_hi = -std::numeric_limits<double>::infinity();
    /// False once any point has had a non-finite coordinate. Tracked separately because a NaN
    /// compares false against everything and would otherwise slip through every extreme.
    bool finite = true;

    /**
     * Fold one position in.
     * @param p the position, GSM, in Earth radii.
     * @complexity O(1) — one fused radius and three compares.
     * @alloc none.
     * @test IrbemT01s.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    constexpr void add(const Position<Frame::GSM>& p) {
        const double r2 = (p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]);
        finite = finite && std::isfinite(r2);
        r2_lo = r2 < r2_lo ? r2 : r2_lo;
        x_lo = p.v[0] < x_lo ? p.v[0] : x_lo;
        x_hi = p.v[0] > x_hi ? p.v[0] : x_hi;
    }

    /**
     * What the batch's positions say about the model's envelope.
     * @return @ref Status::DomainError when any point is not finite, is inside the Earth, or is
     *         sunward of @ref t01s_max_x_gsm; @ref Status::OutOfValidityRange when any point is
     *         tailward of `x = -15 R_E`; otherwise @ref Status::Ok.
     * @complexity O(1) — one square root and one envelope lookup for the whole batch.
     * @alloc none.
     * @test IrbemT01s.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    [[nodiscard]] Status verdict() const {
        if (!finite || x_hi > t01s_max_x_gsm) return Status::DomainError;
        return check_position(ExternalModel::Tsyganenko2001Storm, std::sqrt(r2_lo), x_lo);
    }
};

/**
 * The T01S field over a whole batch of GSM points, on the device when that is worth it.
 *
 * **This is the routine to call for more than a handful of points.** @ref t01s_field is the
 * reference lane: it is what the batch is verified against, and what runs when there is no device
 * or the batch is too small to pay for one.
 *
 * T01S is ~3 700 flops for 24 bytes in and 12 out — ~100 flops/byte, nine times T89's, which
 * already wins 15x on this seam. The measured crossover is in the `irbem_t01s_f32` registry row.
 *
 * The batch reports the same caveats the scalar entry point does, folded over the whole batch:
 * any point tailward of `x = -15` makes the batch @ref Status::OutOfValidityRange (still computed
 * in full); any point that is not finite, inside the Earth or sunward of @ref t01s_max_x_gsm makes
 * it @ref Status::DomainError, and then every output is zeroed, which is exactly what
 * @ref t01s_field does with its one point.
 *
 * @param points the positions, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians.
 * @param d the drivers.
 * @param out receives one field vector per input, GSM, nanotesla; same length as @p points.
 * @return @ref Status::DomainError on a length mismatch, a non-finite driver, a tilt at which the
 *         model has no value, or a refused point, with every output zeroed;
 *         @ref Status::OutOfValidityRange when any point is tailward of -15 R_E, with every point
 *         still computed; otherwise @ref Status::Ok. The value is `true` exactly when the device
 *         lane serviced the call.
 * @complexity O(N); on the device those N run concurrently over `ceil(N/256)` workgroups.
 * @alloc the device lane stages positions and results into two `std::vector<float>` of `3N`; the
 *        host lane allocates nothing.
 * @test IrbemT01s.BatchAgreesWithTheReferenceLane
 * @test IrbemT01s.BatchRejectsMismatchedSpans
 * @test IrbemT01s.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 * @test IrbemT01sGpu.DeviceLaneAgreesWithTheHostLane
 */
[[nodiscard]] inline Result<bool> t01s_field_batch(std::span<const Position<Frame::GSM>> points,
                                                   double tilt_rad, const T01sDrivers& d,
                                                   std::span<FieldVector<Frame::GSM>> out) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(tilt_rad)) return {Status::DomainError, false};
    if (!(std::fabs(tilt_rad) < max_tilt_rad)) return {Status::DomainError, false};
    if (!std::isfinite(d.dst) || !std::isfinite(d.pdyn) || !std::isfinite(d.by_imf) ||
        !std::isfinite(d.bz_imf) || !std::isfinite(d.g2) || !std::isfinite(d.g3)) {
        return {Status::DomainError, false};
    }
    if (n == 0) return {Status::Ok, false};

    const double sin_tilt = std::sin(tilt_rad);
    const double cos_tilt = std::cos(tilt_rad);
    T01sPositionFold fold;

#if CHEATAH_SPACE_IRBEM_T01S_GPU
    if (gpu::prefer_gpu("irbem_t01s_f32", n) &&
        std::filesystem::exists(gpu::shader_path("irbem_t01s_f32"))) {
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
        const std::array<float, t01s_param_count> block =
            t01s_param_block(static_cast<float>(sin_tilt), static_cast<float>(cos_tilt), d);
        gpu::dispatch_batch("irbem_t01s_f32", pos, raw, std::span<const float>(block));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = FieldVector<Frame::GSM>{
                fixarray::vec3d{raw[(3 * i) + 0], raw[(3 * i) + 1], raw[(3 * i) + 2]}};
        }
        return {where, true};
    }
#endif

    const std::array<double, t01s_mode_count> amp = t01s_amplitudes<double>(d);
    for (std::size_t i = 0; i < n; ++i) {
        fold.add(points[i]);
        const std::array<double, 3> b = t01s_components<double>(
            amp, sin_tilt, cos_tilt, points[i].v[0], points[i].v[1], points[i].v[2]);
        out[i] = FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
    }
    const Status where = fold.verdict();
    if (where == Status::DomainError) {
        for (std::size_t i = 0; i < n; ++i) out[i] = FieldVector<Frame::GSM>{};
        return {Status::DomainError, false};
    }
    return {where, false};
}

// -------------------------------------------------------------------------------------------
// The total field
// -------------------------------------------------------------------------------------------

/**
 * IGRF plus T01S, as a single field — the storm-time counterpart of @ref TotalFieldT89.
 *
 * Satisfies @ref GeoFieldModel, so it drops into @ref trace_invariant, @ref make_lstar and
 * everything above them: a drift shell traced through this type is a drift shell in the
 * storm-time magnetosphere the drivers describe. The amplitudes are solved ONCE, at construction,
 * because a trace evaluates the field thousands of times per line at one driver state and the
 * 1 188-term amplitude solve would otherwise be paid at every RK4 stage.
 *
 * There is no device trace kernel for this type yet; @ref trace_invariant runs it on the host.
 *
 * @tparam NMAX the internal field's truncation degree; 10 reproduces IRBEM's own choice.
 * @test IrbemT01s.TotalFieldSuperposesInternalAndExternal
 * @test IrbemT01s.TotalFieldReportsWhenTheExternalModelDeclines
 * @test IrbemT01s.TotalFieldTracesAStormTimeFieldLine
 */
template <int NMAX = 10>
class TotalFieldT01S {
  public:
    /// The internal part's truncation degree — what generic staging and buffer sizing read.
    static constexpr int degree = NMAX;

    /**
     * @param internal the internal field, already built for the epoch.
     * @param rotations the epoch's frame rotations — built once, reused for every point.
     * @param drivers the storm-time driver state the whole field is evaluated at.
     */
    TotalFieldT01S(const Igrf<NMAX>& internal, const Rotations& rotations,
                   const T01sDrivers& drivers)
        : internal_(&internal),
          rotations_(&rotations),
          drivers_(drivers),
          tilt_rad_(rotations.dipole_tilt_deg * (std::numbers::pi / 180.0)),
          sin_tilt_(std::sin(tilt_rad_)),
          cos_tilt_(std::cos(tilt_rad_)),
          amp_(t01s_amplitudes<double>(drivers)) {}

    /**
     * The total field at a geographic point.
     *
     * @param p the position, GEO, Earth radii.
     * @return `B_internal + B_external`, in GEO, nT. When the external model refuses the point —
     *         a domain error — the INTERNAL field is returned alone rather than a zero or a NaN,
     *         exactly as @ref TotalFieldT89::evaluate does and for the same reason; an
     *         out-of-validity point still gets the external contribution, which is the
     *         `status.hpp` rule. Callers who need to know ask @ref external_status.
     * @complexity One IGRF evaluation, one T01S evaluation (no amplitude solve), two rotations.
     * @alloc none.
     * @test IrbemT01s.TotalFieldSuperposesInternalAndExternal
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::GEO>& p) const {
        const FieldVector<Frame::GEO> b_int = internal_->evaluate(p);
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        if (t01s_check(p_gsm, tilt_rad_, drivers_) == Status::DomainError) return b_int;
        const FieldVector<Frame::GSM> b_ext = t01s_field_at(p_gsm, sin_tilt_, cos_tilt_, amp_);
        if (!std::isfinite(b_ext.v[0]) || !std::isfinite(b_ext.v[1]) ||
            !std::isfinite(b_ext.v[2])) {
            return b_int;
        }
        const FieldVector<Frame::GEO> b_ext_geo = transform<Frame::GEO>(b_ext, *rotations_);
        return FieldVector<Frame::GEO>{b_int.v + b_ext_geo.v};
    }

    /**
     * Whether the external model answered at @p p, and if not, why.
     * @param p the position, GEO, Earth radii.
     * @return the external model's status; @ref Status::Ok when it contributed without caveat.
     * @complexity One T01S evaluation and one rotation.
     * @alloc none.
     * @test IrbemT01s.TotalFieldReportsWhenTheExternalModelDeclines
     */
    [[nodiscard]] Status external_status(const Position<Frame::GEO>& p) const {
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        return t01s_field(p_gsm, tilt_rad_, drivers_).status;
    }

    /// The driver state this field was built for.
    /// @return the value passed to the constructor.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01s.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const T01sDrivers& drivers() const { return drivers_; }

    /// The epoch's frame rotations.
    /// @return the rotations this field was built with.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01s.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Rotations& rotations() const { return *rotations_; }

    /// The internal part's Gauss coefficient `g(n, m)`, in nT — forwarded, for the reason
    /// @ref TotalFieldT89::g states: the only thing a caller can mean by it.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01s.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double g(int n, int m) const { return internal_->g(n, m); }

    /// The internal part's `h(n, m)`, in nT — see @ref g.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01s.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double h(int n, int m) const { return internal_->h(n, m); }

    /// The internal field alone.
    /// @return the internal model.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01s.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Igrf<NMAX>& internal() const { return *internal_; }

  private:
    const Igrf<NMAX>* internal_;
    const Rotations* rotations_;
    T01sDrivers drivers_;
    double tilt_rad_;
    double sin_tilt_;
    double cos_tilt_;
    std::array<double, t01s_mode_count> amp_;
};

}  // namespace cheatah::space::irbem
