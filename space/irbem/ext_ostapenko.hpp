#pragma once

/**
 * @file ext_ostapenko.hpp
 * @brief space.irbem — Ostapenko & Maltsev (1997), IRBEM's `kext = 8`: the external field as a
 *        fourth-order polynomial whose 17 amplitudes are linear regressions on Dst, Pdyn, Kp and
 *        IMF Bz.
 *
 * Where @ref t89_field is a physical construction — current sheets with a shape, a thickness and a
 * truncation — this model is a **statistical** one. Ostapenko & Maltsev took 14 073 vector
 * measurements from the Fairfield et al. (1994) database inside `3 <= r`, `rho <= 10 R_E`,
 * `|z| <= 7 R_E`, subtracted the internal field, and fitted what was left with the most general
 * divergence-free fourth-order polynomial that has the magnetosphere's symmetries: six terms even
 * about the noon-midnight meridian and axially symmetric, six carrying the day-night asymmetry,
 * and five proportional to the sine of the dipole tilt. Each of those 17 amplitudes is a linear
 * function of four activity parameters. That is the whole model: 85 numbers, no transcendental in
 * the evaluator beyond the tilt's sine and cosine, and — unlike T89 — **continuous** in every
 * driver, because the regression is linear in each.
 *
 * ## The model, and where each piece comes from
 *
 * Ostapenko & Maltsev, *Relation of the magnetic field in the magnetosphere to the geomagnetic and
 * solar wind activity*, J. Geophys. Res. **102**(A8):17467-17473 (1997). Its eq. (2):
 *
 *     B_ext = sum_{i=1..17} sum_{k=0..4} a_ik * A~_k * b_i
 *
 * with `b_i` the 17 harmonics of the paper's Table 1, `a_ik` the 85 relation coefficients of its
 * Table 4, and `A~_k` the four activity parameters NORMALIZED by eq. (3), `A~ = (A - <A>) / sigma_A`,
 * with `A~_0 = 1` and the means and dispersions of Table 2 (@ref om97_normalization_published).
 * Positions enter the polynomials in units of **10 R_E** (the paper: "we have used the normalized
 * coordinates r~ = r / 10 R_E so that the coefficients are expressed in nanoteslas"), in **solar
 * magnetic** cylindrical coordinates with `z` along the dipole axis and `phi = 0` at noon.
 *
 * ## How Table 1 is evaluated here — Cartesian, not cylindrical
 *
 * The table is printed in `(b_rho, b_phi, b_z)` with `cos phi` and `sin phi` factors, which is
 * undefined on the dipole axis (`rho = 0`, where `phi` is not a number) and needs a `sqrt` and a
 * division per point everywhere else. Every one of the 17 is nonetheless a POLYNOMIAL in Cartesian
 * coordinates — the curl-free ones are gradients of solid harmonics and the rest are low-order
 * monomials — so @ref om97_basis evaluates the Cartesian form: `cos phi = x / rho`,
 * `sin phi = y / rho` substituted and the `rho` cancelled by hand. The two agree to roundoff at
 * every off-axis point (@ref IrbemOm97.CartesianFormsMatchThePrintedCylindricalTable is a
 * transcription check with no way to pass by accident), the Cartesian form is total on the axis,
 * and there is no data-dependent branch anywhere in the evaluator, which is what the device lane
 * wants.
 *
 * **One printed row is corrected, and the paper itself is the authority for the correction.** Its
 * text states that the curl-free harmonics are "expressed in terms of the associated Legendre
 * functions in the Schmidt normalization". Rows 7, 8 and 16 are exactly that — the gradients of
 * `r^n S_n^1(cos theta) cos phi` for `n = 2, 4, 3`, whose leading factors `sqrt(3)`, `sqrt(10)` and
 * `sqrt(3/2)` are the Schmidt factors `sqrt(2 (n-1)! / (n+1)!)` times the Legendre leading
 * coefficient. Row 17 is `n = 5`, for which the same rule gives `sqrt(15) / 8`; the table prints
 * `1 / (8 sqrt(15))`, a factor of **15** smaller — the reciprocal of the Schmidt factor where its
 * square root should be. This header uses the normalization the paper's text prescribes and its
 * other rows obey. @ref IrbemOm97.CurlFreeHarmonicsAreSchmidtNormalizedSolidHarmonics verifies all
 * eight curl-free rows against an independent Legendre construction, and the black-box measurement
 * below confirms that the fitted `a_17k` belong to the Schmidt-normalized row: with the row as
 * printed the reference implementation's harmonic-17 amplitude comes out 15x the tabulated one.
 *
 * ## What this implementation is, and what IRBEM's `kext = 8` is — MEASURED, NOT ASSUMED
 *
 * Both the functional form and every coefficient are published, so the target is agreement with
 * IRBEM's `kext = 8` at the internal-field standard, and the differential harness measures how
 * close the PUBLISHED numbers get and exactly what accounts for the rest.
 *
 * Measured by [`tools/oracle/ostapenko_diff.cpp`](../../../tools/oracle/ostapenko_diff.cpp)
 * against the `-O2` oracle, with the external field isolated as `kext = 8` minus `kext = 0` so
 * the internal IGRF term cancels exactly and the tilt taken from the oracle itself:
 *
 * - **The functional form is the oracle's, to roundoff.** Regressing the oracle's external field
 *   onto this header's 17 basis fields — the Cartesian Table 1, in SM, with the Schmidt row 17 —
 *   leaves an RMS residual of 2e-12 nT at every tilt and every driver set tried (signal RMS
 *   ~22 nT). Evaluating the same basis in GSM instead of SM, or with the tilt sign flipped, leaves
 *   11 and 16 nT. There is no structural difference to document, and no free-refit floor.
 * - **The oracle is linear in every driver to 1e-12**, which is the regression form and nothing
 *   else.
 * - **The difference is the rounding of the printed tables, and nothing else.** Column by column,
 *   the oracle's driver sensitivities are the tabulated `a_ik` scaled by one factor per driver —
 *   0.9889 for Dst, 1.0092 for Pdyn, 0.9601 for Kp, 0.9986 for Bz — i.e. the reference
 *   implementation carries the dispersions `sigma_A` to more digits than the paper's two
 *   (25 / 1.9 / 1.3 / 3.7), and likewise the means. With Table 2 as printed the two differ by an
 *   RMS of 0.37-0.46 nT at quiet and moderate activity and 0.82-1.9 nT at storm and extreme
 *   (0.9-1.8% of the external field, every tilt). Letting only the eight normalization scalars
 *   float, with the published `a_ik` held, explains the oracle's 85 measured intercepts and slopes
 *   to an RMS of 0.0026 nT — which is the rounding floor of a two-decimal table (0.0029) — and
 *   drops the field deviation to 0.012-0.13 nT RMS, 3-6e-4 relative. That floor is the precision
 *   of Table 4 as printed: the oracle evidently carries the `a_ik` to more digits as well, and a
 *   clean room that transcribes the paper cannot go below the paper's own rounding. The recovered
 *   scalars are @ref om97_normalization_measured, offered for callers who need IRBEM agreement and
 *   labelled for what they are: a black-box measurement, not a citation. The published constants
 *   are the default because they are the ones a reader can check against the paper.
 *
 * Two independent checks confirm THIS side regardless of any oracle:
 * @ref IrbemOm97.DivergenceVanishesEverywhere — every harmonic is divergence-free by construction
 * and a second-order stencil residual falls as `h^2` — and the transcription checks above.
 *
 * ## Validity, from the paper rather than the IRBEM table
 *
 * IRBEM's `kext` table names this model's four drivers and publishes no ranges, so `status.hpp`'s
 * row 8 is all-unbounded and @ref check_validity never reports it. The paper itself is more
 * specific, and this header applies what the paper says (@ref om97_fitted_region): the fit covers
 * `3 R_E <= r`, `rho_SM <= 10 R_E`, `|z_SM| <= 7 R_E`; the abstract describes the result as valid
 * "for weak and moderate geomagnetic activity"; and the discussion states that "at intense storms
 * (Dst < -200 nT) our model becomes invalid because it yields a negative total field in the
 * equatorial plane at r > 5 R_E". Outside any of those the value is still returned, with
 * @ref Status::OutOfValidityRange — a polynomial extrapolates without complaint, which is precisely
 * why the caller has to be told.
 *
 * @note Nothing on a hot path allocates, and nothing but the device lane can throw. Measured under
 *       the suite's allocation counter (@ref IrbemOm97.NothingOnTheHeapInTheHotPath): scalar calls
 *       and host batches add nothing to the process's allocation count. The one exception is
 *       @ref om97_field_batch's DEVICE lane, which stages `3N` floats each way and forwards whatever
 *       `gpu::dispatch_batch` throws; its `@alloc` and the dispatch header say so.
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
#  include <filesystem>
#  include "gpu/dispatch.hpp"
/// 1 when this translation unit can reach the OM97 device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_OM97_GPU 1
#else
/// 1 when this translation unit can reach the OM97 device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_OM97_GPU 0
#endif

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// The published tables
// -------------------------------------------------------------------------------------------

/// How many harmonics Table 1 has: six axially symmetric, six day-night asymmetric, five tilt.
inline constexpr std::size_t om97_harmonic_count = 17;

/// How many regressors each amplitude has: the constant term and the four activity parameters,
/// in the paper's column order `{a_i0, a_iDst, a_ip, a_iKp, a_iIMFz}`.
inline constexpr std::size_t om97_regressor_count = 5;

/// The length that positions are divided by before they enter the polynomials, in Earth radii.
/// The paper's `r~ = r / 10 R_E`, chosen so that the coefficients come out in nanotesla.
inline constexpr double om97_length_scale_re = 10.0;

/**
 * The means and dispersions that turn a raw driver into the paper's normalized parameter,
 * `A~ = (A - <A>) / sigma_A` — eq. (3), with the four rows of Table 2.
 *
 * A struct rather than two arrays so that a caller who wants a different set — see
 * @ref om97_normalization_measured — passes one object and cannot pair a mean with the wrong
 * dispersion. Kp here is the paper's **numeric** Kp (`0..9`, with `3-` as 2.7), not the
 * `Kp x 10` of IRBEM's driver vector; the entry points do that division.
 *
 * @test IrbemOm97.PublishedNormalizationIsTable2
 */
struct Om97Normalization {
    /// `<Dst>`, nT.
    double dst_mean;
    /// `sigma_Dst`, nT.
    double dst_sigma;
    /// `<p>`, the solar wind dynamic pressure mean, nPa.
    double pdyn_mean;
    /// `sigma_p`, nPa.
    double pdyn_sigma;
    /// `<Kp>`, numeric.
    double kp_mean;
    /// `sigma_Kp`, numeric.
    double kp_sigma;
    /// `<IMF Bz>`, nT.
    double bz_mean;
    /// `sigma_Bz`, nT.
    double bz_sigma;
};

/**
 * Table 2 of the paper, exactly as printed: `Dst -17 / 25`, `p 2.2 / 1.9`, `Kp 2.3 / 1.3`,
 * `IMFz 0.0 / 3.7`. This is the default normalization everywhere in this header, because it is the
 * one a reader can check against the source.
 *
 * @test IrbemOm97.PublishedNormalizationIsTable2
 */
inline constexpr Om97Normalization om97_normalization_published{
    /* dst */ -17.0, 25.0, /* pdyn */ 2.2, 1.9, /* kp */ 2.3, 1.3, /* bz */ 0.0, 3.7};

/**
 * The normalization IRBEM's `kext = 8` is MEASURED to use — a black-box recovery, not a citation.
 *
 * The paper prints Table 2 to two significant figures; the reference implementation evidently
 * carries the dataset's means and dispersions to more. These eight values are what
 * [`tools/oracle/ostapenko_diff.cpp`](../../../tools/oracle/ostapenko_diff.cpp) recovers by
 * least squares from the oracle's response — its 17 amplitude intercepts and 68 driver slopes,
 * each measured to 1e-12 — under the model "the published `a_ik`, with the eight normalization
 * scalars free". The fit's residual is what says the model is right: it is at the rounding of the
 * two-decimal `a_ik` table and not above it. The recovered digits are quoted to the precision that
 * fit supports and no further.
 *
 * Use this when the requirement is agreement with IRBEM; use @ref om97_normalization_published
 * when the requirement is the paper. The difference is 1-2% of the external field.
 *
 * @test IrbemOm97.MeasuredNormalizationIsCloseToThePublishedOne
 */
inline constexpr Om97Normalization om97_normalization_measured{
    /* dst */ -16.94, 25.28, /* pdyn */ 2.28, 1.883, /* kp */ 2.31, 1.354, /* bz */ 0.02, 3.705};

/// One harmonic's five relation coefficients, `{a_i0, a_iDst, a_ip, a_iKp, a_iIMFz}`, in nT.
using Om97Row = std::array<double, om97_regressor_count>;

/**
 * Table 4 of the paper: the 85 relation coefficients `a_ik`, in nanotesla, row `i - 1` for
 * harmonic `i`. Columns are the paper's: the constant term, then the Dst, p, Kp and IMF Bz
 * regressors on the NORMALIZED parameters.
 *
 * Transcribed from the paper and then checked two ways that cannot both pass by accident: the
 * paper's own Figure 5 profiles (the near-Earth `B_z` depression deepens with negative Dst and
 * compresses with pressure — @ref IrbemOm97.EquatorialBzFollowsThePapersProfiles), and the
 * black-box differential against IRBEM described in the file brief, which recovers every column
 * to its printed precision once the normalization is accounted for.
 *
 * @test IrbemOm97.RelationCoefficientsAreTable4
 */
inline constexpr std::array<Om97Row, om97_harmonic_count> om97_relation_coefficients{{
    //   a_i0     a_iDst   a_ip     a_iKp    a_iIMFz
    {-43.39, 20.36, -0.75, -3.29, -1.63},     // 1  : uniform B_z
    {40.31, -13.73, 6.45, 3.41, 2.89},        // 2  : grad(r^3 P_3)
    {7.05, -3.29, 2.29, 0.11, 0.80},          // 3  : grad(r^5 P_5)
    {133.47, -61.19, 25.06, 2.58, 11.37},     // 4  : B_z = rho^2       (ring current)
    {-36.97, 24.56, -12.77, 4.13, -4.27},     // 5  : B_z = rho^4       (ring current)
    {-130.10, 46.96, -25.64, -14.84, -9.64},  // 6  : (-2 rho z^3, 0, z^4)
    {23.63, 4.91, 7.94, 4.32, 2.69},          // 7  : grad(r^2 S_2^1 cos phi)
    {-1.71, -0.74, 0.29, -0.96, -0.51},       // 8  : grad(r^4 S_4^1 cos phi)
    {1.27, -15.85, 7.17, -1.62, -10.93},      // 9  : uniform B_x = z
    {-21.44, -10.04, -7.19, -9.93, -5.02},    // 10 : B_z = rho^3 cos phi
    {4.49, 3.18, -9.65, -0.69, 3.44},         // 11 : B_x = z^3
    {23.31, -6.88, 12.32, 3.59, -6.12},       // 12 : (0, rho^2 z sin phi, -rho z^2 cos phi / 2)
    {23.18, -2.01, 9.06, -0.19, 0.15},        // 13 : sin(psi) grad(r^2 P_2)
    {-1.97, -1.38, 2.97, -3.41, 1.41},        // 14 : sin(psi) grad(r^4 P_4)
    {12.97, -4.30, 6.49, -0.13, -2.11},       // 15 : sin(psi) uniform B_x
    {5.72, -3.19, 2.22, 0.35, 0.09},          // 16 : sin(psi) grad(r^3 S_3^1 cos phi)
    {4.91, -1.08, 3.57, -1.00, -0.62},        // 17 : sin(psi) grad(r^5 S_5^1 cos phi)
}};

/**
 * Where the paper says the fit holds, and the storm level at which it says the fit fails.
 *
 * The spatial box is the abstract's data domain — "inner boundary r = 3 R_E and outer boundary
 * (x^2 + y^2)^(1/2) = 10 R_E, |z| <= 7 R_E" — in SM coordinates, because that is the frame the
 * data were binned in. The Dst floor is the discussion's "at intense storms (Dst < -200 nT) our
 * model becomes invalid". None of this is in IRBEM's `kext` table, which is why `status.hpp`'s row
 * for this model is unbounded and this header carries its own.
 *
 * @test IrbemOm97.FittedRegionIsThePapers
 */
struct Om97FittedRegion {
    /// Inner boundary of the data, geocentric radius, R_E.
    double r_min;
    /// Outer boundary of the data, cylindrical radius in SM, R_E.
    double rho_sm_max;
    /// Half-height of the data domain, `|z_SM|`, R_E.
    double abs_z_sm_max;
    /// Below this Dst the paper declares the model invalid, nT. Inclusive: `-200` itself is inside.
    double dst_min;
};

/**
 * The fitted region, as the paper states it.
 *
 * @test IrbemOm97.FittedRegionIsThePapers
 */
inline constexpr Om97FittedRegion om97_fitted_region{3.0, 10.0, 7.0, -200.0};

/**
 * The four drivers this model reads, in the units IRBEM's driver vector carries them.
 *
 * A struct rather than four doubles so that a call site cannot transpose Dst and Bz — both are
 * nanotesla, both are usually negative when interesting, and a compiler cannot tell them apart.
 *
 * @test IrbemOm97.ContextOverloadAgreesWithTheExplicitOne
 */
struct Om97Drivers {
    /// Dst, nT — @ref Driver::Dst.
    double dst;
    /// Solar wind dynamic pressure, nPa — @ref Driver::Pdyn.
    double pdyn;
    /// Kp in IRBEM's slot-1 scaling, **Kp x 10** — @ref Driver::Kp. Divided by ten before it meets
    /// the paper's numeric-Kp normalization.
    double kp_times_ten;
    /// IMF Bz in GSM, nT — @ref Driver::BzIMF.
    double bz_imf;
};

// -------------------------------------------------------------------------------------------
// The amplitudes — eq. (2)'s inner sum, done once per driver set
// -------------------------------------------------------------------------------------------

/**
 * The 17 harmonic amplitudes for one driver set: `A_i = sum_k a_ik A~_k`, eq. (2)'s inner sum.
 *
 * Done once per epoch rather than per point, because the drivers are a property of the epoch —
 * the same economy @ref HotState makes with the tilt's trigonometry. The result is what the
 * evaluator, the host batch lane and the device kernel all consume, so the three lanes cannot
 * disagree about what the drivers meant.
 *
 * Computed in `double` and rounded to @p T at the end — the fp32 lane wants the amplitudes the
 * kernel receives, and the kernel receives them as floats made from the best available doubles.
 *
 * @tparam T the scalar type of the result; `double` for the reference lane, `float` for the
 *         device-mirroring one.
 * @param d the drivers, in IRBEM's units; Kp is divided by ten here.
 * @param norm the normalization; the published Table 2 by default.
 * @return the amplitudes, nanotesla, element `i - 1` for harmonic `i`.
 * @complexity O(1) — 85 multiply-adds.
 * @alloc none; the returned object is inline storage.
 * @test IrbemOm97.AmplitudesAtTheMeansAreTheConstantColumn
 * @test IrbemOm97.AmplitudesAreLinearInEveryDriver
 */
template <std::floating_point T>
[[nodiscard]] constexpr std::array<T, om97_harmonic_count> om97_amplitudes(
    const Om97Drivers& d, const Om97Normalization& norm = om97_normalization_published) {
    const double n_dst = (d.dst - norm.dst_mean) / norm.dst_sigma;
    const double n_p = (d.pdyn - norm.pdyn_mean) / norm.pdyn_sigma;
    const double n_kp = ((d.kp_times_ten / 10.0) - norm.kp_mean) / norm.kp_sigma;
    const double n_bz = (d.bz_imf - norm.bz_mean) / norm.bz_sigma;
    std::array<T, om97_harmonic_count> out{};
    for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
        const Om97Row& a = om97_relation_coefficients[i];
        out[i] = static_cast<T>(a[0] + (a[1] * n_dst) + (a[2] * n_p) + (a[3] * n_kp) +
                                (a[4] * n_bz));
    }
    return out;
}

// -------------------------------------------------------------------------------------------
// The evaluator
// -------------------------------------------------------------------------------------------

/// One point's worth of basis: the 17 harmonics' three SM Cartesian components each.
template <std::floating_point T>
using Om97Basis = std::array<std::array<T, 3>, om97_harmonic_count>;

/**
 * Table 1 at one point — the 17 harmonics, each as its three SM Cartesian components.
 *
 * Every row is the printed cylindrical entry with `cos phi = x / rho`, `sin phi = y / rho`
 * substituted and `rho` cancelled, which turns each into a polynomial: `B_x = b_rho cos phi -
 * b_phi sin phi`, `B_y = b_rho sin phi + b_phi cos phi`. For a row of the form
 * `(f cos phi, g sin phi, h cos phi)` that is `B_x = (f x^2 - g y^2) / rho^2`,
 * `B_y = (f + g) x y / rho^2`, `B_z = h x / rho`, and in every row the numerators carry the
 * `rho^2` or `rho` needed to cancel. The tilt rows carry their own `sin psi`, as printed.
 *
 * The four irrational leading factors are the Schmidt-normalized Legendre coefficients:
 * `sqrt(3)` (row 7, `n = 2`), `sqrt(10)` (row 8, `n = 4`), `sqrt(3/2)` (row 16, `n = 3`) and
 * `sqrt(15) / 8` (row 17, `n = 5` — see the file brief for why the printed `1 / (8 sqrt(15))` is
 * a misprint of that).
 *
 * @tparam T the scalar type.
 * @param sin_tilt `sin(psi)`; multiplies rows 13-17.
 * @param x SM x, in units of 10 R_E.
 * @param y SM y, in units of 10 R_E.
 * @param z SM z, in units of 10 R_E.
 * @return the 17 harmonics at the point, SM Cartesian, dimensionless (they multiply amplitudes in
 *         nT).
 * @complexity O(1) — about 150 flops, no branch, no transcendental.
 * @alloc none; the returned object is inline storage.
 * @test IrbemOm97.CartesianFormsMatchThePrintedCylindricalTable
 * @test IrbemOm97.CurlFreeHarmonicsAreSchmidtNormalizedSolidHarmonics
 * @test IrbemOm97.EveryHarmonicIsDivergenceFree
 */
template <std::floating_point T>
[[nodiscard]] constexpr Om97Basis<T> om97_basis(T sin_tilt, T x, T y, T z) {
    // The irrational leading factors, to double precision: sqrt(3), sqrt(10), sqrt(3/2),
    // sqrt(15) / 8 and sqrt(15) / 2. Literals rather than std::sqrt so this stays constexpr.
    const T sqrt3 = std::numbers::sqrt3_v<T>;
    const T sqrt10 = static_cast<T>(3.1622776601683795);
    const T sqrt3_2 = static_cast<T>(1.224744871391589);
    const T sqrt15_8 = static_cast<T>(0.48412291827592713);
    const T sqrt15_2 = static_cast<T>(1.9364916731037085);
    const T sqrt15 = static_cast<T>(3.872983346207417);
    const T zero = static_cast<T>(0);
    const T one = static_cast<T>(1);
    const T two = static_cast<T>(2);
    const T three = static_cast<T>(3);
    const T four = static_cast<T>(4);
    const T five = static_cast<T>(5);
    const T half = static_cast<T>(0.5);
    const T quarter = static_cast<T>(0.25);

    const T x2 = x * x;
    const T y2 = y * y;
    const T z2 = z * z;
    const T rho2 = x2 + y2;
    const T rho4 = rho2 * rho2;
    const T z3 = z2 * z;
    const T z4 = z2 * z2;
    const T xy = x * y;
    const T xz = x * z;
    const T yz = y * z;
    const T sp = sin_tilt;

    Om97Basis<T> b{};
    // ---- the six axially symmetric harmonics ---------------------------------------------
    // 1: uniform B_z.
    b[0] = {zero, zero, one};
    // 2: grad(z^3 - 3 rho^2 z / 2), the r^3 P_3 solid harmonic.
    b[1] = {-three * xz, -three * yz, three * (z2 - (half * rho2))};
    // 3: grad(r^5 P_5).
    const T s3 = (two * z2) - (static_cast<T>(1.5) * rho2);
    b[2] = {-five * xz * s3, -five * yz * s3,
            five * (z4 - (three * rho2 * z2) + (static_cast<T>(0.375) * rho4))};
    // 4, 5: the ring-current rows, B_z = rho^2 and rho^4.
    b[3] = {zero, zero, rho2};
    b[4] = {zero, zero, rho4};
    // 6: (-2 rho z^3, 0, z^4) — divergence-free, not curl-free.
    b[5] = {-two * x * z3, -two * y * z3, z4};

    // ---- the six day-night asymmetric harmonics --------------------------------------------
    // 7: grad(sqrt(3) rho z cos phi) = sqrt(3) (z, 0, x).
    b[6] = {sqrt3 * z, zero, sqrt3 * x};
    // 8: grad(sqrt(10) rho z (z^2 - 3 rho^2 / 4) cos phi).
    b[7] = {sqrt10 * z * (z2 - (quarter * ((static_cast<T>(9) * x2) + (three * y2)))),
            -static_cast<T>(1.5) * sqrt10 * xy * z, three * sqrt10 * x * (z2 - (quarter * rho2))};
    // 9: B_x = z.
    b[8] = {z, zero, zero};
    // 10: B_z = rho^3 cos phi = rho^2 x.
    b[9] = {zero, zero, rho2 * x};
    // 11: B_x = z^3.
    b[10] = {z3, zero, zero};
    // 12: (0, rho^2 z sin phi, -rho z^2 cos phi / 2) = (-y^2 z, x y z, -x z^2 / 2).
    b[11] = {-y2 * z, xy * z, -half * x * z2};

    // ---- the five tilt harmonics, each carrying sin(psi) as printed --------------------------
    // 13: sin(psi) grad(z^2 - rho^2 / 2), the r^2 P_2 solid harmonic.
    b[12] = {-x * sp, -y * sp, two * z * sp};
    // 14: sin(psi) grad(r^4 P_4).
    const T s14 = (half * rho2) - (two * z2);
    b[13] = {three * x * s14 * sp, three * y * s14 * sp, two * z * ((two * z2) - (three * rho2)) * sp};
    // 15: sin(psi) uniform B_x.
    b[14] = {sp, zero, zero};
    // 16: sin(psi) grad(sqrt(3/2) rho (2 z^2 - rho^2 / 2) cos phi).
    b[15] = {sqrt3_2 * ((two * z2) - (half * ((three * x2) + y2))) * sp, -sqrt3_2 * xy * sp,
             four * sqrt3_2 * xz * sp};
    // 17: sin(psi) grad((sqrt(15) / 8) rho (8 z^4 - 12 rho^2 z^2 + rho^4) cos phi).
    b[16] = {sqrt15_8 * ((static_cast<T>(8) * z4) - (static_cast<T>(12) * rho2 * z2) + rho4 -
                         (static_cast<T>(24) * x2 * z2) + (four * rho2 * x2)) * sp,
             sqrt15_2 * xy * (rho2 - (static_cast<T>(6) * z2)) * sp,
             sqrt15 * xz * ((four * z2) - (three * rho2)) * sp};
    return b;
}

/**
 * The OM97 external field at one GSM point, as three components in nanotesla.
 *
 * The whole model in one straight line: rotate GSM to SM about `y` by the tilt (SM's `z` IS the
 * dipole axis, which is the frame the harmonics are defined in), divide by 10 R_E, evaluate the 17
 * harmonics, weight them by the amplitudes, rotate back. No loop with a data-dependent trip count,
 * no branch, no transcendental: the tilt arrives as its sine and cosine.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device kernel.
 * @param amp the 17 amplitudes for the epoch's drivers; see @ref om97_amplitudes.
 * @param sin_tilt `sin(psi)`, the dipole tilt's sine; positive when the north dipole leans sunward.
 * @param cos_tilt `cos(psi)`.
 * @param x the GSM x coordinate, R_E.
 * @param y the GSM y coordinate, R_E.
 * @param z the GSM z coordinate, R_E.
 * @return `{B_x, B_y, B_z}` in GSM, nanotesla.
 * @complexity O(1) — about 250 flops, no branch, no transcendental.
 * @alloc none.
 * @test IrbemOm97.DivergenceVanishesEverywhere
 * @test IrbemOm97.ZeroTiltIsMirrorSymmetricAboutTheEquator
 * @test IrbemOm97.DawnDuskSymmetryHoldsAtEveryTilt
 */
template <std::floating_point T>
[[nodiscard]] constexpr std::array<T, 3> om97_components(
    const std::array<T, om97_harmonic_count>& amp, T sin_tilt, T cos_tilt, T x, T y, T z) {
    const T scale = static_cast<T>(1.0 / om97_length_scale_re);
    // ---- GSM -> SM: one rotation about y by the tilt, then the paper's 10 R_E normalization ----
    const T xs = ((x * cos_tilt) - (z * sin_tilt)) * scale;
    const T ys = y * scale;
    const T zs = ((x * sin_tilt) + (z * cos_tilt)) * scale;

    const Om97Basis<T> b = om97_basis<T>(sin_tilt, xs, ys, zs);
    T bx = static_cast<T>(0);
    T by = static_cast<T>(0);
    T bz = static_cast<T>(0);
    for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
        bx += amp[i] * b[i][0];
        by += amp[i] * b[i][1];
        bz += amp[i] * b[i][2];
    }
    // ---- SM -> GSM: the transpose of the rotation above ----------------------------------
    return {(bx * cos_tilt) + (bz * sin_tilt), by, (-bx * sin_tilt) + (bz * cos_tilt)};
}

/**
 * The OM97 external field at one GSM point, in `double` — the reference lane.
 *
 * @param p the position, GSM, in Earth radii.
 * @param sin_tilt `sin(psi)`; @ref HotState::sin_tilt holds it, precomputed per epoch.
 * @param cos_tilt `cos(psi)`.
 * @param amp the epoch's amplitudes, from @ref om97_amplitudes.
 * @return the external field at @p p, GSM, in nanotesla.
 * @complexity O(1); see @ref om97_components.
 * @alloc none.
 * @test IrbemOm97.ReferenceLaneMatchesTheComponentForm
 */
[[nodiscard]] inline FieldVector<Frame::GSM> om97_field_at(
    Position<Frame::GSM> p, double sin_tilt, double cos_tilt,
    const std::array<double, om97_harmonic_count>& amp) {
    const std::array<double, 3> b =
        om97_components<double>(amp, sin_tilt, cos_tilt, p.v[0], p.v[1], p.v[2]);
    return FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
}

// -------------------------------------------------------------------------------------------
// Validity
// -------------------------------------------------------------------------------------------

/**
 * Whether a point is inside the region the paper fitted, from its geocentric radius and its SM
 * cylindrical coordinates.
 *
 * Separate from `status.hpp`'s @ref check_position because the box is stated in SM, which that
 * function does not know about, and because IRBEM's table — which that function reads — publishes
 * no limit for this model at all. The bounds are closed: a point exactly on the boundary is inside.
 *
 * @param r the geocentric radius, R_E.
 * @param rho_sm the cylindrical radius in SM, R_E.
 * @param abs_z_sm `|z_SM|`, R_E.
 * @return @ref Status::OutOfValidityRange outside the paper's box, otherwise @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOm97.PositionValidityIsCheckedFromBothSides
 */
[[nodiscard]] constexpr Status om97_check_fitted_region(double r, double rho_sm, double abs_z_sm) {
    if (r < om97_fitted_region.r_min || rho_sm > om97_fitted_region.rho_sm_max ||
        abs_z_sm > om97_fitted_region.abs_z_sm_max) {
        return Status::OutOfValidityRange;
    }
    return Status::Ok;
}

/**
 * Whether the drivers are ones the paper says the model can be believed at.
 *
 * Two layers, in `status.hpp`'s order. First @ref check_validity for this model's row, which
 * catches a non-finite driver among the four this model reads (and, the IRBEM table being
 * unbounded for `kext = 8`, nothing else). Then the paper's own caveat, which the table does not
 * carry: `Dst < -200 nT` is declared invalid. The interval is closed, so `-200` itself is inside.
 *
 * @param d the drivers.
 * @return @ref Status::DomainError for a non-finite driver, @ref Status::OutOfValidityRange below
 *         the paper's Dst floor, otherwise @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOm97.DriverValidityIsCheckedFromBothSides
 */
[[nodiscard]] inline Status om97_check_drivers(const Om97Drivers& d) {
    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Kp)] = d.kp_times_ten;
    drivers[static_cast<std::size_t>(Driver::Dst)] = d.dst;
    drivers[static_cast<std::size_t>(Driver::Pdyn)] = d.pdyn;
    drivers[static_cast<std::size_t>(Driver::BzIMF)] = d.bz_imf;
    const Status table = check_validity(ExternalModel::OstapenkoMaltsev1997, drivers);
    if (table != Status::Ok) return table;
    if (d.dst < om97_fitted_region.dst_min) return Status::OutOfValidityRange;
    return Status::Ok;
}

/**
 * The OM97 external field, with the model's own verdict on whether it should be believed here.
 *
 * The value is **always** returned, including when the status is @ref Status::OutOfValidityRange —
 * `status.hpp`'s standing rule. What is refused outright is arithmetic that has no answer: a
 * non-finite input, a point inside the Earth, a tilt that is not an angle to an axis, or an
 * extrapolation so far that the fourth-order polynomial overflows (a unit-confusion bug's
 * signature — kilometres where Earth radii were meant — and the one way this evaluator can produce
 * a non-finite number from finite input). Unlike @ref t89_field there is no refusal at
 * `|psi| = pi/2`: the model carries `sin psi` and the rotation, never `tan psi`, so a right-angle
 * tilt is merely unphysical, not undefined.
 *
 * What is checked against the paper: the position against the fitted box
 * (@ref om97_check_fitted_region) and the drivers against its Dst floor (@ref om97_check_drivers).
 * Both report, neither suppresses.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians; positive when the north dipole leans sunward.
 * @param d the drivers, IRBEM units (Kp x 10).
 * @param norm the normalization; the published Table 2 by default.
 * @return the field and its caveat. @ref Status::DomainError (with a zero field) for a non-finite
 *         input, a radius inside the Earth, `|psi| > pi/2`, or a non-finite result;
 *         @ref Status::OutOfValidityRange for a position outside the fitted box or a Dst below
 *         -200 nT, with the field still computed; otherwise @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOm97.PositionValidityIsCheckedFromBothSides
 * @test IrbemOm97.DriverValidityIsCheckedFromBothSides
 * @test IrbemOm97.NonFiniteInputIsADomainError
 * @test IrbemOm97.RightAngleTiltIsNotRefusedButBeyondItIs
 * @test IrbemOm97.AnOverflowingExtrapolationIsADomainErrorNotANaN
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> om97_field(
    Position<Frame::GSM> p, double tilt_rad, const Om97Drivers& d,
    const Om97Normalization& norm = om97_normalization_published) {
    const FieldVector<Frame::GSM> zero{};
    if (!std::isfinite(p.v[0]) || !std::isfinite(p.v[1]) || !std::isfinite(p.v[2]) ||
        !std::isfinite(tilt_rad)) {
        return {Status::DomainError, zero};
    }
    if (std::fabs(tilt_rad) > max_tilt_rad) return {Status::DomainError, zero};
    const Status drives = om97_check_drivers(d);
    if (drives == Status::DomainError) return {Status::DomainError, zero};

    const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
    const Status where = check_position(ExternalModel::OstapenkoMaltsev1997, r, p.v[0]);
    if (where == Status::DomainError) return {Status::DomainError, zero};

    const double sin_tilt = std::sin(tilt_rad);
    const double cos_tilt = std::cos(tilt_rad);
    const double xs = (p.v[0] * cos_tilt) - (p.v[2] * sin_tilt);
    const double zs = (p.v[0] * sin_tilt) + (p.v[2] * cos_tilt);
    const Status fitted =
        om97_check_fitted_region(r, std::sqrt((xs * xs) + (p.v[1] * p.v[1])), std::fabs(zs));

    const FieldVector<Frame::GSM> b =
        om97_field_at(p, sin_tilt, cos_tilt, om97_amplitudes<double>(d, norm));
    if (!std::isfinite(b.v[0]) || !std::isfinite(b.v[1]) || !std::isfinite(b.v[2])) {
        return {Status::DomainError, zero};
    }
    return {first_failure(drives, first_failure(where, fitted)), b};
}

/**
 * The OM97 external field for a whole epoch's worth of state — the production entry point.
 *
 * Reads the tilt and the four drivers straight out of @ref HotState. The amplitudes are still
 * recomputed per call (85 multiply-adds); a caller evaluating a whole batch at one epoch should use
 * @ref om97_field_batch, which computes them once.
 *
 * @param p the position, GSM, in Earth radii.
 * @param ctx the epoch's context; its `hot()` block carries the tilt, Kp, Dst, Pdyn and Bz.
 * @param norm the normalization; the published Table 2 by default.
 * @return the field and its caveat, exactly as the explicit overload.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOm97.ContextOverloadAgreesWithTheExplicitOne
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> om97_field(
    Position<Frame::GSM> p, const FieldContext& ctx,
    const Om97Normalization& norm = om97_normalization_published) {
    const HotState& h = ctx.hot();
    return om97_field(p, h.tilt_rad, Om97Drivers{h.dst, h.pdyn, h.kp, h.bz_imf}, norm);
}

// -------------------------------------------------------------------------------------------
// The batch lanes
// -------------------------------------------------------------------------------------------

/**
 * The OM97 field over a whole batch, on the CPU, in `float`.
 *
 * The host twin of `irbem_om97_f32`: the same expressions, in the same order, in the same
 * precision, from the same float amplitudes. That is what makes a disagreement between the two
 * lanes attributable to the device rather than to the arithmetic having been written differently
 * on the two sides. It is also the lane a machine with no GPU actually runs for `float` batches.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, GSM, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, nanotesla; overwritten in full.
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`.
 * @param amp the epoch's amplitudes, already rounded to `float`.
 * @return `false` when @p pos is not a whole number of points or @p out is a different length, in
 *         which case nothing is written; `true` otherwise.
 * @complexity O(N).
 * @alloc none.
 * @test IrbemOm97.HostFloatLaneTracksTheReferenceLane
 * @test IrbemOm97.HostFloatLaneRejectsMismatchedSpans
 */
[[nodiscard]] inline bool om97_field_host(std::span<const float> pos, std::span<float> out,
                                          float sin_tilt, float cos_tilt,
                                          const std::array<float, om97_harmonic_count>& amp) {
    if (pos.size() % 3 != 0 || out.size() != pos.size()) return false;
    const std::size_t n = pos.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const std::array<float, 3> b = om97_components<float>(
            amp, sin_tilt, cos_tilt, pos[(3 * i) + 0], pos[(3 * i) + 1], pos[(3 * i) + 2]);
        out[(3 * i) + 0] = b[0];
        out[(3 * i) + 1] = b[1];
        out[(3 * i) + 2] = b[2];
    }
    return true;
}

/// How many `float` scalars the device kernel's parameter buffer holds: `sin(psi)`, `cos(psi)`,
/// then the 17 amplitudes. Asserted against the kernel registry by the suite.
inline constexpr std::size_t om97_param_count = 2 + om97_harmonic_count;

/**
 * Pack the epoch's tilt and amplitudes into the kernel's parameter buffer.
 *
 * The layout is the kernel's ABI and is stated in exactly two places — here and the comment above
 * `irbem_om97_f32` in `irbem.slang`. A test evaluates both lanes on the same points, which is what
 * actually keeps the two statements in step.
 *
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`.
 * @param amp the epoch's amplitudes, in `float`.
 * @return the parameter block, @ref om97_param_count floats, by value.
 * @complexity O(1).
 * @alloc none — the block is the returned object's own inline array.
 * @test IrbemOm97.ParameterBlockCarriesTheTiltThenTheAmplitudes
 */
[[nodiscard]] constexpr std::array<float, om97_param_count> om97_param_block(
    float sin_tilt, float cos_tilt, const std::array<float, om97_harmonic_count>& amp) {
    std::array<float, om97_param_count> block{};
    block[0] = sin_tilt;
    block[1] = cos_tilt;
    for (std::size_t k = 0; k < om97_harmonic_count; ++k) block[2 + k] = amp[k];
    return block;
}

/**
 * The batch's position caveat, accumulated one point at a time.
 *
 * One @ref Status for N points can only be the worst of them, and computing it must not cost a
 * second pass over the positions. The fitted box is three monotone bounds — a minimum radius, a
 * maximum SM cylindrical radius, a maximum `|z_SM|` — so the extremes decide the batch, and they
 * are folded as SQUARES (radius and cylindrical radius) so there is no per-point `sqrt`. The SM
 * `z` is one fused multiply-add from the tilt's sine and cosine, and `rho_SM^2 = r^2 - z_SM^2`
 * because a rotation preserves `r`.
 *
 * @test IrbemOm97.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
struct Om97PositionFold {
    /// `sin(psi)` for the batch's epoch; zero tilt until the batch sets it.
    double sin_tilt = 0.0;
    /// `cos(psi)` for the batch's epoch; zero tilt until the batch sets it.
    double cos_tilt = 1.0;
    /// The smallest `r^2` seen, R_E^2; `+inf` until the first point.
    double r2_lo = std::numeric_limits<double>::infinity();
    /// The largest `rho_SM^2` seen, R_E^2; zero until the first point.
    double rho2_hi = 0.0;
    /// The largest `|z_SM|` seen, R_E; zero until the first point.
    double abs_z_hi = 0.0;
    /// False once any point has had a non-finite coordinate. Tracked separately because a NaN
    /// compares false against everything and would otherwise slip through every extreme.
    bool finite = true;

    /**
     * Fold one position in.
     * @param p the position, GSM, in Earth radii.
     * @complexity O(1) — one fused radius, one rotation component, four comparisons, no `sqrt`.
     * @alloc none.
     * @test IrbemOm97.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    constexpr void add(const Position<Frame::GSM>& p) {
        const double r2 = (p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]);
        const double zs = (p.v[0] * sin_tilt) + (p.v[2] * cos_tilt);
        const double rho2 = r2 - (zs * zs);
        finite = finite && std::isfinite(r2);
        r2_lo = r2 < r2_lo ? r2 : r2_lo;
        rho2_hi = rho2 > rho2_hi ? rho2 : rho2_hi;
        const double az = zs < 0.0 ? -zs : zs;
        abs_z_hi = az > abs_z_hi ? az : abs_z_hi;
    }

    /**
     * What the batch's positions say about the model's envelope.
     * @return @ref Status::DomainError when any point is not finite or is inside the Earth,
     *         @ref Status::OutOfValidityRange when any point is outside the paper's fitted box,
     *         otherwise @ref Status::Ok.
     * @complexity O(1) — two square roots and two envelope checks for the whole batch.
     * @alloc none.
     * @test IrbemOm97.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    [[nodiscard]] Status verdict() const {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double r_lo = finite ? std::sqrt(r2_lo) : nan;
        const Status where = check_position(ExternalModel::OstapenkoMaltsev1997, r_lo, 0.0);
        if (where != Status::Ok) return where;
        return om97_check_fitted_region(r_lo, std::sqrt(rho2_hi), abs_z_hi);
    }
};

/**
 * The OM97 field over a whole batch of GSM points, on the device when that is worth it.
 *
 * **This is the routine to call for more than a handful of points.** @ref om97_field is the
 * reference lane: it is what the batch is verified against, and what runs when there is no device
 * or the batch is too small to pay for one.
 *
 * The amplitudes are computed ONCE for the batch, which is the whole economy of the epoch model:
 * per point the kernel does ~250 flops of polynomial over 24 bytes in and 12 out, a little under
 * T89's ~11 flops/byte and well above the streaming dipole's 0.5. Its measured crossover is in the
 * `irbem_om97_f32` row of `gpu/dispatch.hpp`.
 *
 * **The batch reports the same caveats the scalar entry point does, folded over the whole batch.**
 * A DOMAIN-ERROR batch is zeroed in full and never reaches the device; an out-of-validity batch is
 * computed in full.
 *
 * @param points the positions, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians.
 * @param d the drivers, IRBEM units (Kp x 10).
 * @param out receives one field vector per input, GSM, nanotesla; same length as @p points.
 * @param norm the normalization; the published Table 2 by default.
 * @return @ref Status::DomainError on a length mismatch, a non-finite tilt or driver, `|psi| >
 *         pi/2`, or a point that is not finite or is inside the Earth, and then every output is
 *         zeroed; @ref Status::OutOfValidityRange when Dst is below the paper's floor or any point
 *         is outside the fitted box, with every point still computed; otherwise @ref Status::Ok.
 *         The value is `true` exactly when the device lane serviced the call.
 * @complexity O(N); on the device those N run concurrently over `ceil(N/256)` workgroups.
 * @alloc the device lane stages positions and results into two `std::vector<float>` of `3N`; the
 *        host lane allocates nothing.
 * @test IrbemOm97.BatchAgreesWithTheReferenceLane
 * @test IrbemOm97.BatchRejectsMismatchedSpans
 * @test IrbemOm97.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
[[nodiscard]] inline Result<bool> om97_field_batch(
    std::span<const Position<Frame::GSM>> points, double tilt_rad, const Om97Drivers& d,
    std::span<FieldVector<Frame::GSM>> out,
    const Om97Normalization& norm = om97_normalization_published) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(tilt_rad) || std::fabs(tilt_rad) > max_tilt_rad) {
        return {Status::DomainError, false};
    }
    const Status drives = om97_check_drivers(d);
    if (drives == Status::DomainError) return {Status::DomainError, false};
    if (n == 0) return {drives, false};

    const double sin_tilt = std::sin(tilt_rad);
    const double cos_tilt = std::cos(tilt_rad);
    Om97PositionFold fold{sin_tilt, cos_tilt};

#if CHEATAH_SPACE_IRBEM_OM97_GPU
    if (gpu::prefer_gpu("irbem_om97_f32", n) &&
        std::filesystem::exists(gpu::shader_path("irbem_om97_f32"))) {
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
        const std::array<float, om97_param_count> block = om97_param_block(
            static_cast<float>(sin_tilt), static_cast<float>(cos_tilt),
            om97_amplitudes<float>(d, norm));
        gpu::dispatch_batch("irbem_om97_f32", pos, raw, std::span<const float>(block));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = FieldVector<Frame::GSM>{
                fixarray::vec3d{raw[(3 * i) + 0], raw[(3 * i) + 1], raw[(3 * i) + 2]}};
        }
        return {first_failure(drives, where), true};
    }
#endif

    const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(d, norm);
    for (std::size_t i = 0; i < n; ++i) {
        fold.add(points[i]);
        const std::array<double, 3> b = om97_components<double>(
            amp, sin_tilt, cos_tilt, points[i].v[0], points[i].v[1], points[i].v[2]);
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
// The total field — IGRF plus OM97, as one GeoFieldModel
// -------------------------------------------------------------------------------------------

/**
 * IGRF plus Ostapenko & Maltsev (1997), as a single field.
 *
 * Satisfies @ref GeoFieldModel, so it drops into @ref trace_invariant and @ref make_lstar exactly
 * as @ref TotalFieldT89 does; the frames are handled the same way (IGRF in GEO, the external model
 * in GSM, the epoch's @ref Rotations built once). Every evaluation goes through @ref om97_field so
 * that the validity verdict @ref external_status reports is the one the value was computed under;
 * the 85 multiply-adds that rebuilds the amplitudes are noise beside the IGRF evaluation.
 *
 * @tparam NMAX the internal field's truncation degree. 10 reproduces IRBEM's own choice.
 * @test IrbemOm97.TotalFieldSuperposesInternalAndExternal
 */
template <int NMAX = 10>
class TotalFieldOm97 {
  public:
    /// The internal part's truncation degree — what generic staging and buffer sizing read.
    static constexpr int degree = NMAX;

    /**
     * @param internal the internal field, already built for the epoch.
     * @param rotations the epoch's frame rotations — built once, reused for every point.
     * @param drivers the epoch's drivers, IRBEM units (Kp x 10).
     * @param norm the normalization; the published Table 2 by default.
     * @complexity O(1) — pointers, a copy of the drivers, one degrees-to-radians.
     * @alloc none.
     * @test IrbemOm97.TotalFieldSuperposesInternalAndExternal
     */
    constexpr TotalFieldOm97(const Igrf<NMAX>& internal, const Rotations& rotations,
                             const Om97Drivers& drivers,
                             const Om97Normalization& norm = om97_normalization_published)
        : internal_(&internal),
          rotations_(&rotations),
          drivers_(drivers),
          norm_(norm),
          tilt_rad_(rotations.dipole_tilt_deg * (std::numbers::pi / 180.0)) {}

    /**
     * The total field at a geographic point.
     *
     * @param p the position, GEO, Earth radii.
     * @return `B_internal + B_external`, in GEO, nT. When the external model refuses the point
     *         (@ref Status::DomainError) the INTERNAL field is returned alone rather than a zero
     *         or a NaN, for the reason @ref TotalFieldT89::evaluate gives; out-of-validity
     *         extrapolations are included, since a polynomial has a value there and the trace
     *         needs a continuous field.
     * @complexity One IGRF evaluation, one OM97 evaluation, two 3x3 rotations.
     * @alloc none.
     * @test IrbemOm97.TotalFieldSuperposesInternalAndExternal
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::GEO>& p) const {
        const FieldVector<Frame::GEO> b_int = internal_->evaluate(p);
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        const Result<FieldVector<Frame::GSM>> b_ext =
            om97_field(p_gsm, tilt_rad_, drivers_, norm_);
        if (b_ext.status == Status::DomainError) return b_int;
        const FieldVector<Frame::GEO> b_ext_geo = transform<Frame::GEO>(b_ext.value, *rotations_);
        return FieldVector<Frame::GEO>{b_int.v + b_ext_geo.v};
    }

    /**
     * Whether the external model answered at @p p, and if not, why.
     * @param p the position, GEO, Earth radii.
     * @return the external model's status; @ref Status::Ok when it contributed without caveat.
     * @complexity One OM97 evaluation and one rotation.
     * @alloc none.
     * @test IrbemOm97.TotalFieldReportsWhenTheExternalModelDeclines
     */
    [[nodiscard]] Status external_status(const Position<Frame::GEO>& p) const {
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        return om97_field(p_gsm, tilt_rad_, drivers_, norm_).status;
    }

    /// The drivers this field was built for.
    /// @return the value passed to the constructor.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOm97.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Om97Drivers& drivers() const { return drivers_; }

    /// The epoch's frame rotations.
    /// @return the rotations this field was built with.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOm97.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Rotations& rotations() const { return *rotations_; }

    /// The internal part's Gauss coefficient `g(n, m)`, in nT — forwarded for the reason
    /// @ref TotalFieldT89::g gives: the dipole moment and the trace step sizing are questions
    /// about the internal field.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOm97.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double g(int n, int m) const { return internal_->g(n, m); }

    /// The internal part's `h(n, m)`, in nT — see @ref g.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOm97.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double h(int n, int m) const { return internal_->h(n, m); }

    /// The internal field alone.
    /// @return the internal model.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOm97.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Igrf<NMAX>& internal() const { return *internal_; }

  private:
    const Igrf<NMAX>* internal_;
    const Rotations* rotations_;
    Om97Drivers drivers_;
    Om97Normalization norm_;
    double tilt_rad_;
};

}  // namespace cheatah::space::irbem
