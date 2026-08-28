#pragma once

/**
 * @file ext_t89.hpp
 * @brief space.irbem — Tsyganenko (1989), the first EXTERNAL magnetospheric field model, and the
 *        one every later one is validated through.
 *
 * Everything this module has built so far describes the field the Earth makes. That field is not
 * what a particle at L = 6 actually feels. Above about 4 R_E the magnetosphere's own currents — the
 * ring current, the cross-tail current sheet, the currents closing on the magnetopause — contribute
 * tens of nanotesla against an internal field of a few hundred, and at geosynchronous during a storm
 * they dominate the shape of the field line entirely. A drift shell traced through IGRF alone is not
 * a drift shell; it is a dipole exercise. @ref t89_field_at is where the real magnetosphere enters.
 *
 * ## Why T89 is the first one
 *
 * It is **pure straight-line closed-form arithmetic**: no loops, no data-dependent branches, no
 * iteration, roughly 400 flops with exactly one `exp` and thirteen square roots — and no trigonometry
 * at all in the inner evaluator, because the tilt's sine and cosine are properties of the epoch and
 * are paid once per timestamp (@ref HotState precomputes them).
 * Its successors (T96, T01, TS05) are fixed-trip-count harmonic expansions an order of magnitude
 * heavier. That makes T89 both the easiest to verify — every term can be differentiated by hand and
 * checked against a finite difference — and by far the best GPU kernel of the family: ~400 flops
 * over 24 bytes in and 12 out is ~11 flops/byte, the same regime as IGRF (~20), which wins 8.96x on
 * this seam, and nothing like the streaming dipole (0.5), which loses. Measured, it wins 15.1x at
 * 2^22 points and crosses over at ~1300 — see the `irbem_t89_f32` row of `gpu/dispatch.hpp`.
 *
 * ## The model, and where each equation comes from
 *
 * Tsyganenko, *A magnetospheric magnetic field model with a warped tail current sheet*,
 * Planet. Space Sci. **37**(1):5-20 (1989). Every formula below carries its equation number from
 * that paper, and the whole construction is one idea applied four times:
 *
 * 1. Solve for the vector potential of an infinitely thin axially symmetric current disc. With
 *    `A = {0, A(rho, z), 0}` in cylindrical coordinates the paper's eqs. (1)-(5) invert to the
 *    strikingly simple `A^(1) = rho / (S + a + |z|)`, `S = sqrt(rho^2 + (a + |z|)^2)` — eq. (7),
 *    written there in the equivalent form `rho^-1 {[(a+|z|)^2 + rho^2]^(1/2) - (a+|z|)}`.
 * 2. Differentiate with respect to the scale length `a` to get a family with progressively steeper
 *    radial fall-off: `A^(2) = dA^(1)/da` (eq. 8) and `A^(3) = dA^(2)/da` (eq. 9). `A^(1)` and
 *    `A^(2)` model the tail current, `A^(3)` — the only one with a finite magnetic moment — the ring
 *    current.
 * 3. Give the sheet a **thickness** by replacing `|z|` with `xi = sqrt(z^2 + D^2)` and a **shape**
 *    by replacing `z` with `z - z_s(x, y, psi)`, the warped surface of eq. (11). Both substitutions
 *    are made in the POTENTIAL, not in the field, which is why the result stays divergence-free
 *    however the sheet is bent — the property @ref IrbemT89.DivergenceVanishesEverywhere pins down.
 * 4. Truncate it to the magnetotail by multiplying the potential by `W(x, y)` (eq. 12/13), which is
 *    ~1 in the central tail and falls to zero towards the subsolar magnetopause and the flanks.
 *
 * The remaining two contributions are the currents that close the system: a pair of planar sheets at
 * `z = ∓R_T` carrying the return current across the high-latitude magnetopause (eqs. 18-19, an
 * `A^(1)` with `a = D = 0`), and a low-order polynomial in `y` and `z` with an `exp(x / dx)` radial
 * envelope standing in for the Chapman-Ferraro and Birkeland systems (eq. 20).
 *
 * ## Two coordinate systems, on purpose
 *
 * The paper is explicit and it matters: **eqs. (12)-(17) are evaluated in SM** — solar magnetic,
 * where the dipole axis IS `z` — because the current sheet is defined by its departure from the
 * dipole equatorial plane; **eqs. (18)-(20) are evaluated in GSM**, because the closure sheets are
 * parallel to the GSM equator. This header takes GSM in and gives GSM out, and does the one rotation
 * about `y` by the dipole tilt internally (see @ref t89_components). Evaluating the tail term in GSM
 * instead is the obvious-looking simplification and puts the current sheet in the wrong place by
 * `x tan(psi)` — 4 R_E at 20 R_E down-tail at a 12-degree tilt.
 *
 * ## Kp is a BIN, not a number
 *
 * T89 is not a continuous function of Kp. Tsyganenko sorted the data into intervals and fitted a
 * separate coefficient set to each (Table 1). @ref t89_kp_bin in `status.hpp` does the binning;
 * @ref t89_coefficient_sets holds the sets. A caller who interpolates between two of them is
 * inventing a model that was never fitted.
 *
 * ## What this implementation is, and what IRBEM's `kext = 4` is — MEASURED, NOT ASSUMED
 *
 * This is the model **as published in 1989**: the functional form of eqs. (11)-(20) and the
 * coefficients of the paper's Table 1. IRBEM's `kext = 4` is labelled "Tsyganenko [1989c]" and
 * evaluates the **1992 revision**, which the community documentation describes (Oulu space-physics
 * reference, magbase.rssi.ru/REFMAN/SPPHTEXT/bmodels.html) as differing in two ways: ISEE-1/2 data
 * were added to the original IMP–HEOS fit, and **two terms were added to the tail field modes,
 * modulating the tail current by the geodipole tilt** — a STRUCTURAL change to the equations, not
 * a re-fit. That description matches what the free-refit experiment below measures independently:
 * the revision was distributed only as code (GPL-3.0, which this MIT clean-room implementation may
 * not read) and its equations were never published, so the published 1989 form is what a clean
 * room can implement, and the differential harness reports the model-family gap instead of hiding
 * it.
 *
 * Measured by [`tools/oracle/t89_diff.cpp`](../../../tools/oracle/t89_diff.cpp) against the `-O2`
 * oracle, at three dipole tilts (+0.0002, +25.64, -30.42 degrees), with the external field isolated
 * as `kext = 4` minus `kext = 0` so the internal IGRF term cancels exactly and the tilt taken from
 * the oracle itself so a frame difference cannot masquerade as a model difference:
 *
 * | Kp bin | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
 * |---|---|---|---|---|---|---|---|
 * | belts, 3-10 R_E: RMS abs | 2.18 nT | 2.50 | 3.90 | 4.58 | 7.01 | 8.96 | 24.9 |
 * | belts, 3-10 R_E: RMS rel | 13.4% | 12.5% | 16.2% | 15.8% | 19.8% | 20.4% | 42.0% |
 * | out to 35 R_E: RMS rel | 14.1% | 14.0% | 16.2% | 19.4% | 29.2% | 33.6% | 58.7% |
 *
 * Bin 7 is the worst by a wide margin and that is exactly the documented gap: IRBEM's revision
 * splits `Kp >= 5-` into two bins and this one cannot, because the seventh set was never published.
 *
 * Carried through the whole drift-shell chain (`make_lstar` over @ref TotalFieldT89, matched
 * `options(3,4) = 9`, L = 3..6.6), the same model-family gap bounds L\*: worst `|dL*|` vs the
 * oracle at `kext = 4` is 0.028–0.088 across bins 1–6 and 0.123 at bin 7 — against 0.0066 for the
 * internal field alone, where both sides evaluate the SAME model. The order-of-magnitude contrast
 * is the point: the numerics agree, the model families differ, and the gap sits exactly where the
 * field-level table above predicts it.
 *
 * The same harness runs the experiment that says the rest of the gap is **structural, not a
 * re-fit**: letting all 19 linear coefficients of the published form float freely against the
 * oracle leaves an RMS residual of 0.63 nT (5.7% of the external field), and letting all 9
 * non-linear parameters float too only reaches 0.44 nT (4.0%). A re-fit of the same equations would
 * drive that to roundoff. It does not go there, so IRBEM's "T89c" is not the published functional
 * form with different numbers in it.
 *
 * Two independent checks confirm that THIS side of the comparison is right regardless, and neither
 * can pass by accident: @ref IrbemT89.PublishedCoefficientsAreDivergenceFree verifies all 24 of the
 * divergence identities the paper says its eq. (20) coefficients satisfy (they hold to 1.3e-3, the
 * rounding of a four-significant-figure table), and @ref IrbemT89.DivergenceVanishesEverywhere
 * verifies that every analytic derivative in the tail, ring-current and closure terms is the
 * derivative it claims to be — with eq. (20) switched off, `|div B|` falls as `h^2` from 9.0e-04 to
 * 9.0e-08 nT/R_E as the difference step goes from 1e-2 to 1e-4 R_E, which is the signature of an
 * exactly divergence-free field sampled by a second-order stencil and of nothing else.
 *
 * @note Nothing on a hot path allocates, and nothing but the device lane can throw. The evaluator,
 *       the scalar entry points and the fp32 host lane touch nothing but the caller's spans and one
 *       stack parameter block — measured under `valgrind --tool=memcheck`: 100 000 scalar
 *       `t89_field` calls plus ten 4096-point host batches plus fifty `t89_field_batch` calls come
 *       to "total heap usage: 6 allocs", all six the harness's own set-up vectors. The one
 *       exception is @ref t89_field_batch's DEVICE lane, which stages `3N` floats each way and
 *       forwards whatever `gpu::dispatch_batch` throws when a device that answered
 *       `gpu::available()` then fails to launch; its `@alloc` and the dispatch header say so.
 */

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "context.hpp"
#include "frames.hpp"
#include "status.hpp"

// The device lane is opt-in by include path, exactly as lstar.hpp's is. Without
// cheatah-gpu-linalg every routine here still compiles and runs on the host.
#if __has_include("cheatah_gpu_linalg/context.hpp")
#  include "gpu/dispatch.hpp"
/// 1 when this translation unit can reach the T89 device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_T89_GPU 1
#else
/// 1 when this translation unit can reach the T89 device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_T89_GPU 0
#endif

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// The parameters the paper FIXED, as opposed to the ones it fitted
// -------------------------------------------------------------------------------------------

/**
 * The T89 parameters that are the same in every Kp bin.
 *
 * Tsyganenko (1989) §3 and §4 name these explicitly as held constant during the fit, "from a priori
 * considerations and preliminary test runs": the §3 group because the data cannot separate them from
 * the free parameters they trade against (the paper's own example is that `gamma_t` and `gamma_rc`
 * "yield nearly the same effects" on the dayside, so one of the two had to be pinned), the §4 group
 * because the closure sheets contribute too little of the total field to constrain their geometry.
 *
 * They are gathered into one `constexpr` object rather than scattered as literals so that the
 * evaluator below reads as the paper's equations and every magic number has a name and a citation.
 *
 * @test IrbemT89.FixedParametersAreThePublishedOnes
 */
struct T89FixedParameters {
    /// `L_y`, R_E — the dawn-dusk scale of the current sheet's transverse bending, eq. (11).
    /// Fixed at 10 R_E "in accordance with results of Fairfield (1980) and Gosling et al. (1986)".
    double l_y;
    /// `D_x`, R_E — the sunward scale over which the truncation factor `W` falls off, eq. (13).
    double d_x;
    /// `L_RC`, R_E — the scale distance of the ring-current thickness profile `h_RC`, eq. (13).
    double l_rc;
    /// `L_T`, R_E — the scale distance of the tail thickness profile `h_T`, eq. (13).
    double l_t;
    /// `gamma_T`, R_E — how much thicker the tail sheet is on the dayside than the nightside,
    /// eq. (13). The one of the two thickness increments the paper pinned; `gamma_RC` stayed free.
    double gamma_t;
    /// `delta`, R_E^-1 — the rate at which the tail sheet thickens towards its flanks, eq. (13).
    double delta;
    /// `gamma_1`, R_E — the amplitude of the far-tail thickening `h_1` that removes the spurious
    /// `B_z` reversals the paper reports beyond `x ~ -15 R_E`, eq. (13).
    double gamma_1;
    /// `R_T`, R_E — the distance of the two planar closure current sheets from the GSM equatorial
    /// plane, §4. They sit at `z = ∓R_T`, i.e. outside the modelling region.
    double r_t;
    /// `x_0c`, R_E — where the closure sheets' own truncation factor `W_c` falls off, eq. (19).
    double x_0c;
    /// `L_xc^2`, R_E^2 — the squared sunward scale of `W_c`, eq. (19). Squared because that is how
    /// the paper states it and how it enters: never as `L_xc`.
    double l_xc2;
    /// `D_yc^2`, R_E^2 — the squared dawn-dusk scale of `W_c`, eq. (19). Squared, as above.
    double d_yc2;
    /// The `16` inside eq. (11)'s square root: `(4 R_E)^2`, the scale over which the sheet's hinge
    /// is rounded off. Named rather than inlined so the equation reads as printed.
    double hinge_scale2;
    /// The `16` in eq. (13)'s `h_1`, R_E — the tailward distance at which the far-tail thickening is
    /// half-grown.
    double h1_offset;
    /// The `36` in eq. (13)'s `h_1`: `(6 R_E)^2`, the squared scale over which it grows.
    double h1_scale2;
};

/**
 * The fixed parameters, exactly as Tsyganenko (1989) §3 and §4 state them.
 *
 * @test IrbemT89.FixedParametersAreThePublishedOnes
 */
inline constexpr T89FixedParameters t89_fixed{
    /* l_y */ 10.0,     /* d_x */ 13.0,   /* l_rc */ 5.0,       /* l_t */ 6.3,
    /* gamma_t */ 4.0,  /* delta */ 0.01, /* gamma_1 */ 1.0,    /* r_t */ 30.0,
    /* x_0c */ 4.0,     /* l_xc2 */ 50.0, /* d_yc2 */ 20.0,     /* hinge_scale2 */ 16.0,
    /* h1_offset */ 16.0, /* h1_scale2 */ 36.0};

// -------------------------------------------------------------------------------------------
// The fitted parameters — one set per Kp bin
// -------------------------------------------------------------------------------------------

/// How many coefficients eq. (12), (18) and (20) between them carry: `C_1 .. C_19`.
inline constexpr std::size_t t89_linear_count = 19;

/**
 * One Kp bin's fitted parameters — Tsyganenko (1989) Table 1, one column.
 *
 * Templated on the scalar type for one reason and it is not generality: the device kernel receives
 * these as `float`, so the host lane that the kernel is verified against must round them to `float`
 * FIRST and then evaluate, or the two lanes are not running the same arithmetic and a disagreement
 * cannot be attributed. @ref t89_parameters is the converter, and the fp64 and fp32 instantiations
 * are the only two that exist.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` for the device-mirroring one.
 * @test IrbemT89.PublishedCoefficientsAreDivergenceFree
 */
template <std::floating_point T>
struct T89Parameters {
    /// The nineteen linear coefficients, zero-based: element k is the paper's C(k+1), so element
    /// zero is C_1. C_1 and C_2 weight the two tail potentials; C_3 and C_4 the symmetric and
    /// tilt-antisymmetric closure currents; C_5 the ring current; the remaining fourteen are the
    /// eq. (20) expansion. Units are whatever makes the field come out in nT with positions in
    /// Earth radii, which is the paper's convention throughout.
    std::array<T, t89_linear_count> c;
    /// `dx`, R_E — the e-folding length of eq. (20)'s radial envelope.
    T delta_x;
    /// `a_RC`, R_E — the ring current's radial scale length; the radius of its density maximum.
    T a_rc;
    /// `D_0`, R_E — the current sheet half-thickness in the central magnetotail.
    T d_0;
    /// `gamma_RC`, R_E — the ring-current sheet's dayside-to-nightside thickness increment.
    T gamma_rc;
    /// `R_c`, R_E — the "hinging distance": how far out the sheet leaves the dipole equator.
    T r_c;
    /// `G`, R_E — the amplitude of the sheet's transverse (dawn-dusk) bending.
    T g;
    /// `a_T`, R_E — the tail current's radial scale length.
    T a_t;
    /// `D_y`, R_E — the dawn-dusk scale over which the truncation factor `W` falls off.
    T d_y;
    /// `x_0`, R_E — where along the Sun-Earth line `W` makes its transition.
    T x_0;
};

/// The fp64 reference spelling of a coefficient set.
using T89Coefficients = T89Parameters<double>;

/// How many Kp bins the API carries. Seven, because that is what @ref t89_kp_bin returns and what
/// the published `iopt` interface takes — see @ref t89_coefficient_sets for what fills bin 7.
inline constexpr std::size_t t89_bin_count = 7;

/// How many DISTINCT coefficient sets Tsyganenko (1989) Table 1 publishes. Six, not seven.
inline constexpr std::size_t t89_published_set_count = 6;

/**
 * Tsyganenko (1989) Table 1: the fitted parameters, indexed by `bin - 1`.
 *
 * The columns are, in order, `Kp = {0, 0+}`, `{1-, 1, 1+}`, `{2-, 2, 2+}`, `{3-, 3, 3+}`,
 * `{4-, 4, 4+}` and `Kp >= 5-`.
 *
 * **The seventh entry is a stated gap, not a seventh column.** The published table has six columns;
 * the seven-bin split that @ref t89_kp_bin implements — with `{5-, 5, 5+}` and `{>= 6-}` separated —
 * belongs to the post-publication "T89c" revision, whose coefficients appear only in Tsyganenko's
 * GPL-3.0 source and are therefore unavailable to this MIT clean-room implementation. Bins 6 and 7
 * consequently carry the SAME published `Kp >= 5-` set, which is exactly what the 1989 model says
 * about `Kp >= 5-` and is not an extrapolation of anything. A caller evaluating at `Kp >= 6` gets
 * the most disturbed published parameterization, and @ref t89_bin_is_published says so.
 *
 * Every number here is transcribed from the paper and then checked, not trusted: the paper states
 * that `C_16 .. C_19` are not free but are determined from `C_6 .. C_15` and `dx` by
 * `div B = 0` applied to eq. (20), which is four linear identities per bin. All 24 hold to better
 * than 1.5e-3 — see @ref IrbemT89.PublishedCoefficientsAreDivergenceFree, which is a transcription
 * check with no way to pass by accident.
 *
 * @test IrbemT89.PublishedCoefficientsAreDivergenceFree
 * @test IrbemT89.BinSevenRepeatsTheMostDisturbedPublishedSet
 */
inline constexpr std::array<T89Coefficients, t89_bin_count> t89_coefficient_sets{{
    // Kp = 0, 0+
    {{-98.72, -10014.0, 15.03, 76.62, -10237.0, 1.813, 31.10, -0.07464, -0.07764, 0.003303, -1.129,
      0.001663, 0.000988, 18.21, -0.03018, -0.03829, -0.1283, -0.001973, 0.000717},
     24.74, 8.161, 2.08, -0.8799, 9.084, 3.838, 13.55, 26.94, 5.745},
    // Kp = 1-, 1, 1+
    {{-35.64, -12800.0, 14.37, 124.5, -13543.0, 2.316, 35.64, -0.0741, -0.1081, 0.003924, -1.451,
      0.00202, 0.00111, 21.37, -0.04567, -0.05382, -0.1457, -0.002742, 0.001244},
     22.33, 8.119, 1.664, 0.9324, 9.238, 2.426, 13.81, 28.83, 6.052},
    // Kp = 2-, 2, 2+
    {{-77.45, -14588.0, 64.85, 123.9, -16229.0, 2.641, 42.46, -0.07611, -0.1579, 0.004078, -1.391,
      0.00153, 0.000727, 21.86, -0.04199, -0.06523, -0.6412, -0.000948, 0.002276},
     20.90, 6.283, 1.541, 4.183, 9.609, 6.591, 15.08, 30.57, 7.435},
    // Kp = 3-, 3, 3+
    {{-70.12, -16125.0, 90.71, 38.08, -19630.0, 3.181, 47.50, -0.1327, -0.1864, 0.01382, -1.488,
      0.002962, 0.000897, 22.74, -0.04095, -0.09223, -1.059, -0.001766, 0.003034},
     18.64, 6.266, 0.9351, 5.389, 8.573, 5.935, 15.63, 31.47, 8.103},
    // Kp = 4-, 4, 4+
    {{-162.5, -15806.0, 160.6, 5.888, -27534.0, 3.607, 51.10, -0.1006, -0.1927, 0.03353, -1.392,
      0.001594, 0.002439, 22.41, -0.04925, -0.1153, -1.399, 0.000716, 0.002696},
     18.31, 6.196, 0.7677, 5.072, 10.06, 6.668, 16.11, 30.04, 8.260},
    // Kp >= 5- : the last PUBLISHED column, used for bin 6 ...
    {{-128.4, -16184.0, 149.1, 215.5, -36435.0, 4.090, 49.09, -0.0231, -0.1359, 0.01989, -2.298,
      0.004911, 0.003421, 21.79, -0.05447, -0.1149, -0.2214, -0.01355, 0.001185},
     19.48, 5.831, 0.3325, 6.472, 10.47, 9.081, 15.85, 25.27, 7.976},
    // ... and again for bin 7. See the note above: this is the published Kp >= 5- set, not a
    // seventh column, and not an invented one.
    {{-128.4, -16184.0, 149.1, 215.5, -36435.0, 4.090, 49.09, -0.0231, -0.1359, 0.01989, -2.298,
      0.004911, 0.003421, 21.79, -0.05447, -0.1149, -0.2214, -0.01355, 0.001185},
     19.48, 5.831, 0.3325, 6.472, 10.47, 9.081, 15.85, 25.27, 7.976},
}};

/**
 * Whether @p bin has a coefficient set of its own in the published table.
 *
 * False only for bin 7, which repeats bin 6's published `Kp >= 5-` column. Exposed rather than left
 * to a comment because a caller running a storm-time study deserves to know that the two highest
 * bins are not distinguished by this implementation.
 *
 * @param bin the Kp bin, `1..7` as @ref t89_kp_bin returns.
 * @return `true` when the paper publishes a column for that bin; `false` for bin 7 and for any
 *         out-of-range value, which has no published set either.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT89.BinSevenRepeatsTheMostDisturbedPublishedSet
 */
[[nodiscard]] constexpr bool t89_bin_is_published(int bin) {
    return bin >= 1 && bin <= static_cast<int>(t89_published_set_count);
}

/**
 * The parameters for @p bin, in the scalar type the caller's lane evaluates in.
 *
 * The `float` instantiation is not a convenience: it is how the host lane is made to run the same
 * arithmetic as the device kernel, which receives these as `float`. Rounding the coefficients on the
 * way in — rather than evaluating in `double` and rounding the answer — is what makes a
 * host-vs-device disagreement attributable to the DEVICE.
 *
 * @tparam T the scalar type; `double` or `float`.
 * @param bin the Kp bin, `1..7`. Values outside that range are clamped, because this function is
 *        called after @ref t89_field has already decided what status the answer carries and must not
 *        be able to index out of the table.
 * @return the parameter set, by value, converted to @p T.
 * @complexity O(1) — 28 conversions, all of which the compiler folds for a compile-time @p bin.
 * @alloc none; the returned object is inline storage.
 * @test IrbemT89.ParametersRoundTripThroughFloat
 * @test IrbemT89.ParametersClampAnOutOfRangeBin
 */
template <std::floating_point T>
[[nodiscard]] constexpr T89Parameters<T> t89_parameters(int bin) {
    const int index = bin < 1 ? 0
                     : bin > static_cast<int>(t89_bin_count) ? static_cast<int>(t89_bin_count) - 1
                                                             : bin - 1;
    const T89Coefficients& s = t89_coefficient_sets[static_cast<std::size_t>(index)];
    T89Parameters<T> out{};
    for (std::size_t k = 0; k < t89_linear_count; ++k) out.c[k] = static_cast<T>(s.c[k]);
    out.delta_x = static_cast<T>(s.delta_x);
    out.a_rc = static_cast<T>(s.a_rc);
    out.d_0 = static_cast<T>(s.d_0);
    out.gamma_rc = static_cast<T>(s.gamma_rc);
    out.r_c = static_cast<T>(s.r_c);
    out.g = static_cast<T>(s.g);
    out.a_t = static_cast<T>(s.a_t);
    out.d_y = static_cast<T>(s.d_y);
    out.x_0 = static_cast<T>(s.x_0);
    return out;
}

// -------------------------------------------------------------------------------------------
// The evaluator
// -------------------------------------------------------------------------------------------

/**
 * The T89 external field at one GSM point, as three components in nanotesla.
 *
 * This is the whole model, and it is written so that each block is one equation of the paper. The
 * order is the paper's order, and the two coordinate systems are kept apart deliberately (see the
 * file brief): the tail and ring current are evaluated in SM and rotated to GSM; the closure and
 * Chapman-Ferraro terms are evaluated in GSM directly.
 *
 * **On the derivatives.** Every `d/dx` and `d/dy` below is the analytic derivative of the expression
 * above it, and none is optional: the field is the curl of a vector potential whose `x` and `y`
 * dependence runs through `W`, through the sheet surface `z_s` and through the thickness `D`, so
 * dropping any of them does not merely approximate the answer, it breaks `div B = 0` — which is why
 * a finite-difference divergence test catches exactly this class of mistake and is the main test
 * this file carries.
 *
 * With `P = A_phi / rho` the potential-over-radius, the curl of `A = P * (-y, x, 0)` is
 * `B_x = -x dP/dz`, `B_y = -y dP/dz`, `B_z = 2P + x dP/dx + y dP/dy` — which is where the paper's
 * eqs. (14)-(17) and (19) come from, and the form in which they are evaluated here.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device kernel.
 * @param p the fitted parameters for the Kp bin; see @ref t89_parameters.
 * @param sin_tilt `sin(psi)`, the dipole tilt's sine. Taken precomputed because it is a property of
 *        the epoch, not of the point: @ref FieldContext pays for it once per timestamp.
 * @param cos_tilt `cos(psi)`. Must be non-zero — eq. (11) carries `tan(psi)` — which
 *        @ref t89_field checks before it gets here.
 * @param x the GSM x coordinate, R_E.
 * @param y the GSM y coordinate, R_E.
 * @param z the GSM z coordinate, R_E.
 * @return `{B_x, B_y, B_z}` in GSM, nanotesla.
 * @complexity O(1) — about 400 flops, thirteen square roots and one `exp`. No loop, no
 *             branch on data, and no trigonometry: the tilt arrives already resolved into its sine
 *             and cosine.
 * @alloc none.
 * @test IrbemT89.DivergenceVanishesEverywhere
 * @test IrbemT89.ZeroTiltIsMirrorSymmetricAboutTheEquator
 * @test IrbemT89.DawnDuskSymmetryHoldsAtEveryTilt
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, 3> t89_components(const T89Parameters<T>& p, T sin_tilt,
                                                     T cos_tilt, T x, T y, T z) {
    const T half = static_cast<T>(0.5);
    const T one = static_cast<T>(1);
    const T two = static_cast<T>(2);
    const T three = static_cast<T>(3);

    const T c1 = p.c[0];
    const T c2 = p.c[1];
    const T c3 = p.c[2];
    const T c4 = p.c[3];
    const T c5 = p.c[4];

    const T tan_tilt = sin_tilt / cos_tilt;

    // ---- GSM -> SM. A rotation about y by the tilt; SM's z axis IS the dipole axis. ----------
    const T xs = (x * cos_tilt) - (z * sin_tilt);
    const T ys = y;
    const T zs = (x * sin_tilt) + (z * cos_tilt);

    // ---- eq. (11): the warped current sheet surface, and its slopes -------------------------
    // 0.5 tan(psi) (x + R_c - sqrt((x + R_c)^2 + 16)) - G sin(psi) y^4 / (y^4 + L_y^4).
    // The first term hinges: ~0 near the Earth (the sheet lies in the dipole equator, z_SM = 0) and
    // ~tan(psi)(x + R_c) far down-tail (the sheet is parallel to the GSM equator, which in SM is
    // z_SM = x_SM tan psi). The second bends it in the dawn-dusk direction.
    const T hinge = xs + p.r_c;
    const T hinge_root = std::sqrt((hinge * hinge) + static_cast<T>(t89_fixed.hinge_scale2));
    const T y2 = ys * ys;
    const T y4 = y2 * y2;
    const T ly2 = static_cast<T>(t89_fixed.l_y) * static_cast<T>(t89_fixed.l_y);
    const T ly4 = ly2 * ly2;
    const T bend_den = y4 + ly4;
    const T sheet = (half * tan_tilt * (hinge - hinge_root)) - (p.g * sin_tilt * y4 / bend_den);
    const T dsheet_dx = half * tan_tilt * (one - (hinge / hinge_root));
    const T dsheet_dy =
        -p.g * sin_tilt * static_cast<T>(4) * ys * y2 * ly4 / (bend_den * bend_den);

    const T zr = zs - sheet;      // distance from the warped sheet
    const T rho2 = (xs * xs) + (ys * ys);

    // ---- eq. (13): the thickness profiles, and their slopes ---------------------------------
    // h_T and h_RC rise from 0 in the tail to 1 on the dayside (the sheet is thicker sunward);
    // h_1 rises from 0 sunward to 1 beyond x ~ -16 R_E (the far-tail thickening that removes the
    // spurious B_z reversals the paper reports).
    const T lt2 = static_cast<T>(t89_fixed.l_t) * static_cast<T>(t89_fixed.l_t);
    const T lrc2 = static_cast<T>(t89_fixed.l_rc) * static_cast<T>(t89_fixed.l_rc);
    const T rt_root = std::sqrt((xs * xs) + lt2);
    const T rrc_root = std::sqrt((xs * xs) + lrc2);
    const T xh1 = xs + static_cast<T>(t89_fixed.h1_offset);
    const T h1_root = std::sqrt((xh1 * xh1) + static_cast<T>(t89_fixed.h1_scale2));
    const T h_t = half * (one + (xs / rt_root));
    const T h_rc = half * (one + (xs / rrc_root));
    const T h_1 = half * (one - (xh1 / h1_root));
    const T dh_t = half * lt2 / (rt_root * rt_root * rt_root);
    const T dh_rc = half * lrc2 / (rrc_root * rrc_root * rrc_root);
    const T dh_1 = -half * static_cast<T>(t89_fixed.h1_scale2) / (h1_root * h1_root * h1_root);

    const T gamma_t = static_cast<T>(t89_fixed.gamma_t);
    const T gamma_1 = static_cast<T>(t89_fixed.gamma_1);
    const T delta = static_cast<T>(t89_fixed.delta);
    const T d_tail = p.d_0 + (delta * y2) + (gamma_t * h_t) + (gamma_1 * h_1);
    const T d_ring = p.d_0 + (p.gamma_rc * h_rc) + (gamma_1 * h_1);
    const T dd_tail_dx = (gamma_t * dh_t) + (gamma_1 * dh_1);
    const T dd_tail_dy = two * delta * ys;
    const T dd_ring_dx = (p.gamma_rc * dh_rc) + (gamma_1 * dh_1);

    // ---- eq. (12)/(13): the broadened radial coordinates -------------------------------------
    const T xi_t = std::sqrt((zr * zr) + (d_tail * d_tail));
    const T xi_rc = std::sqrt((zr * zr) + (d_ring * d_ring));
    const T u_t = p.a_t + xi_t;
    const T u_rc = p.a_rc + xi_rc;
    const T s_t = std::sqrt(rho2 + (u_t * u_t));
    const T s_rc = std::sqrt(rho2 + (u_rc * u_rc));

    // ---- eq. (13): the truncation factor W, and its slopes -----------------------------------
    const T xw = xs - p.x_0;
    const T dx2 = static_cast<T>(t89_fixed.d_x) * static_cast<T>(t89_fixed.d_x);
    const T w_root = std::sqrt((xw * xw) + dx2);
    const T w_x = half * (one - (xw / w_root));
    const T w_y = one / (one + (y2 / (p.d_y * p.d_y)));
    const T w = w_x * w_y;
    const T dw_dx = w_y * (-half * dx2 / (w_root * w_root * w_root));
    const T dw_dy = w_x * (-two * ys / (p.d_y * p.d_y)) * w_y * w_y;
    const T grad_w = (xs * dw_dx) + (ys * dw_dy);

    // ---- eqs. (14)-(15): the tail current, in SM ---------------------------------------------
    const T q_t = (w / (xi_t * s_t)) * ((c1 / (s_t + u_t)) + (c2 / (s_t * s_t)));
    const T bx_t = q_t * xs * zr;
    const T by_t = q_t * ys * zr;
    const T bz_t = ((w / s_t) * (c1 + (c2 * u_t / (s_t * s_t)))) +
                   (grad_w * (c1 + (c2 / s_t)) / (s_t + u_t)) + (bx_t * dsheet_dx) +
                   (by_t * dsheet_dy) -
                   (q_t * d_tail * ((xs * dd_tail_dx) + (ys * dd_tail_dy)));

    // ---- eqs. (16)-(17): the ring current, in SM ---------------------------------------------
    const T s_rc2 = s_rc * s_rc;
    const T s_rc5 = s_rc2 * s_rc2 * s_rc;
    const T q_rc = three * c5 * u_rc / (xi_rc * s_rc5);
    const T bx_rc = q_rc * xs * zr;
    const T by_rc = q_rc * ys * zr;
    const T bz_rc = (c5 * ((two * u_rc * u_rc) - rho2) / s_rc5) + (bx_rc * dsheet_dx) +
                    (by_rc * dsheet_dy) - (q_rc * d_ring * xs * dd_ring_dx);

    // ---- SM -> GSM. The inverse of the rotation above; a transpose, not a solve. --------------
    const T bx_sm = bx_t + bx_rc;
    const T by_sm = by_t + by_rc;
    const T bz_sm = bz_t + bz_rc;
    T bx = (bx_sm * cos_tilt) + (bz_sm * sin_tilt);
    T by = by_sm;
    T bz = (-bx_sm * sin_tilt) + (bz_sm * cos_tilt);

    // ---- eqs. (18)-(19): the closure currents, in GSM -----------------------------------------
    // Two planar sheets at z = ∓R_T, each an eq. (7) potential with a = D = 0. The paper writes the
    // pair with a ± that resolves |z ± R_T| for a point BETWEEN the sheets; using sgn(z ± R_T)
    // instead is identical there and keeps the function total outside, where the paper's spelling
    // would divide by zero rather than merely being out of the model's domain.
    T f_x[2]{};
    T f_y[2]{};
    T f_z[2]{};
    for (int k = 0; k < 2; ++k) {
        const T side = (k == 0) ? one : -one;
        const T zeta = z + (side * static_cast<T>(t89_fixed.r_t));
        const T s = std::sqrt((zeta * zeta) + (x * x) + (y * y));
        const T sgn = (zeta >= static_cast<T>(0)) ? one : -one;
        const T den = s + (sgn * zeta);
        if (!(den > static_cast<T>(0))) continue;  // only at the sheet's own axis origin
        const T xc = x - static_cast<T>(t89_fixed.x_0c);
        const T sc = std::sqrt((xc * xc) + static_cast<T>(t89_fixed.l_xc2));
        const T wc_x = half * (one - (xc / sc));
        const T wc_y = one / (one + (y * y / static_cast<T>(t89_fixed.d_yc2)));
        const T wc = wc_x * wc_y;
        const T dwc_dx = wc_y * (-half * static_cast<T>(t89_fixed.l_xc2) / (sc * sc * sc));
        const T dwc_dy =
            wc_x * (-two * y / static_cast<T>(t89_fixed.d_yc2)) * wc_y * wc_y;
        f_x[k] = sgn * wc * x / (s * den);
        f_y[k] = sgn * wc * y / (s * den);
        f_z[k] = (wc / s) + (((x * dwc_dx) + (y * dwc_dy)) / den);
    }
    bx += (c3 * (f_x[0] + f_x[1])) + (c4 * (f_x[0] - f_x[1]) * sin_tilt);
    by += (c3 * (f_y[0] + f_y[1])) + (c4 * (f_y[0] - f_y[1]) * sin_tilt);
    bz += (c3 * (f_z[0] + f_z[1])) + (c4 * (f_z[0] - f_z[1]) * sin_tilt);

    // ---- eq. (20): the Chapman-Ferraro and Birkeland systems, in GSM --------------------------
    // A low-order polynomial in y and z under an exp(x / dx) envelope. C_16..C_19 are not free — the
    // paper derives them from C_6..C_15 and dx through div B = 0 — but they are TABULATED, so they
    // are used as published and the four identities are asserted in the suite instead.
    const T env = std::exp(x / p.delta_x);
    const T z2 = z * z;
    bx += env * ((p.c[5] * z * cos_tilt) +
                 (((p.c[6]) + (p.c[7] * y2) + (p.c[8] * z2)) * sin_tilt));
    by += env * ((p.c[9] * y * z * cos_tilt) +
                 (((p.c[10] * y) + (p.c[11] * y * y2) + (p.c[12] * y * z2)) * sin_tilt));
    bz += env * ((((p.c[13]) + (p.c[14] * y2) + (p.c[15] * z2)) * cos_tilt) +
                 (((p.c[16] * z) + (p.c[17] * z * y2) + (p.c[18] * z * z2)) * sin_tilt));

    return {bx, by, bz};
}

/**
 * The T89 external field at one GSM point, in `double` — the reference lane.
 *
 * @param p the position, GSM, in Earth radii.
 * @param sin_tilt `sin(psi)`; @ref HotState::sin_tilt holds it, precomputed per epoch.
 * @param cos_tilt `cos(psi)`; must be non-zero, which @ref t89_field checks.
 * @param bin the Kp bin, `1..7`; out-of-range values are clamped by @ref t89_parameters.
 * @return the external field at @p p, GSM, in nanotesla.
 * @complexity O(1); see @ref t89_components.
 * @alloc none.
 * @test IrbemT89.ReferenceLaneMatchesTheComponentForm
 */
[[nodiscard]] inline FieldVector<Frame::GSM> t89_field_at(Position<Frame::GSM> p, double sin_tilt,
                                                          double cos_tilt, int bin) {
    const std::array<double, 3> b = t89_components<double>(t89_parameters<double>(bin), sin_tilt,
                                                           cos_tilt, p.v[0], p.v[1], p.v[2]);
    return FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
}

/**
 * The T89 external field, with the model's own verdict on whether it should be believed here.
 *
 * The value is **always** returned, including when the status is @ref Status::OutOfValidityRange —
 * that is `status.hpp`'s standing rule and the whole reason a @ref Result exists: extrapolating an
 * empirical fit is a decision only the caller can make. What is refused outright is arithmetic that
 * has no answer: a non-finite input, or a tilt of exactly ±90 degrees, at which eq. (11)'s
 * `tan(psi)` does not exist.
 *
 * Two things are checked against the published envelope, both through `status.hpp` so the rules live
 * in one place: Kp against `0 <= Kp <= 9` (in IRBEM's Kp x 10 scaling), and the position against
 * T89's `r_GEO <= 70 R_E`. The Kp check happens before the binning, so a caller who passes Kp = 12
 * gets the most disturbed published set AND is told the model was never fitted there.
 *
 * There is one more refusal and it is about the OUTPUT, not the input. Eq. (20) carries
 * `exp(x / dx)` with `dx ~ 20 R_E`, so a caller who extrapolates far enough sunward — past about
 * 1e4 R_E, which is nothing a trace can reach but is exactly what a unit-confusion bug produces
 * when kilometres arrive where Earth radii were meant — overflows it. The field is then infinite,
 * and `inf * 0` in the `B_y` assembly makes a NaN out of a component that is EXACTLY zero
 * everywhere the model means anything (@ref IrbemT89.OnTheNoonMidnightMeridianTheFieldStaysInThatPlane).
 * A NaN that escapes here surfaces a hundred RK4 steps later with nothing left pointing at its
 * cause, so a non-finite answer is refused the same way a non-finite input is.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians; positive when the north dipole leans sunward.
 * @param kp_times_ten Kp in IRBEM's `maginput` slot-1 scaling, i.e. Kp x 10, nominally 0..90.
 * @return the field and its caveat. @ref Status::DomainError (with a zero field) for a non-finite
 *         input, a radius inside the Earth, `|psi| >= pi/2`, or a position so far extrapolated that
 *         eq. (20)'s exponential overflows; @ref Status::OutOfValidityRange for a Kp or a radius
 *         outside the published envelope, with the field still computed; otherwise
 *         @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT89.OutOfRangeKpIsReportedButStillEvaluated
 * @test IrbemT89.NonFiniteInputIsADomainError
 * @test IrbemT89.RightAngleTiltIsADomainError
 * @test IrbemT89.AnOverflowingExtrapolationIsADomainErrorNotANaN
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> t89_field(Position<Frame::GSM> p,
                                                               double tilt_rad,
                                                               double kp_times_ten) {
    const FieldVector<Frame::GSM> zero{};
    if (!std::isfinite(p.v[0]) || !std::isfinite(p.v[1]) || !std::isfinite(p.v[2]) ||
        !std::isfinite(tilt_rad) || !std::isfinite(kp_times_ten)) {
        return {Status::DomainError, zero};
    }
    if (!(std::fabs(tilt_rad) < max_tilt_rad)) return {Status::DomainError, zero};

    const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
    const Status where = check_position(ExternalModel::Tsyganenko1989, r, p.v[0]);
    if (where == Status::DomainError) return {Status::DomainError, zero};

    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Kp)] = kp_times_ten;
    const Status drives = check_validity(ExternalModel::Tsyganenko1989, drivers);

    const FieldVector<Frame::GSM> b =
        t89_field_at(p, std::sin(tilt_rad), std::cos(tilt_rad), t89_kp_bin(kp_times_ten));
    if (!std::isfinite(b.v[0]) || !std::isfinite(b.v[1]) || !std::isfinite(b.v[2])) {
        return {Status::DomainError, zero};
    }
    return {first_failure(drives, where), b};
}

/**
 * The T89 external field for a whole epoch's worth of state — the production entry point.
 *
 * Reads the tilt's precomputed sine and cosine and Kp straight out of @ref HotState, so a batch
 * over one timestamp pays for neither the trigonometry nor the driver lookup per point. The tilt
 * ANGLE is still read, because the `|psi| < pi/2` check is on the angle and @ref FieldContext
 * already guarantees `|psi| <= pi/2` but not the strict inequality eq. (11) needs.
 *
 * @param p the position, GSM, in Earth radii.
 * @param ctx the epoch's context; its `hot()` block carries `sin(psi)`, `cos(psi)` and Kp.
 * @return the field and its caveat, exactly as the three-argument overload.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT89.ContextOverloadAgreesWithTheExplicitOne
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> t89_field(Position<Frame::GSM> p,
                                                               const FieldContext& ctx) {
    return t89_field(p, ctx.hot().tilt_rad, ctx.hot().kp);
}

// -------------------------------------------------------------------------------------------
// The batch lanes
// -------------------------------------------------------------------------------------------

/**
 * The T89 field over a whole batch, on the CPU, in `float`.
 *
 * The host twin of `irbem_t89_f32`: the same expressions, in the same order, in the same precision,
 * evaluated from coefficients rounded to `float` FIRST. That is what makes a disagreement between
 * the two lanes attributable to the device — a contraction, a driver's `exp` — rather than to the
 * arithmetic having been written differently on the two sides. It is also the lane a machine with no
 * GPU actually runs.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, GSM, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, nanotesla; overwritten in full.
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`; must be non-zero.
 * @param bin the Kp bin, `1..7`.
 * @return `false` when @p pos is not a whole number of points or @p out is a different length, in
 *         which case nothing is written; `true` otherwise. A bool rather than a throw because this
 *         is the fallback lane of a batch entry point that already reports through @ref Result.
 * @complexity O(N).
 * @alloc none. Not one byte: the loop is over caller-provided spans and one stack parameter set.
 * @test IrbemT89.HostFloatLaneTracksTheReferenceLane
 * @test IrbemT89.HostFloatLaneRejectsMismatchedSpans
 */
[[nodiscard]] inline bool t89_field_host(std::span<const float> pos, std::span<float> out,
                                         float sin_tilt, float cos_tilt, int bin) {
    if (pos.size() % 3 != 0 || out.size() != pos.size()) return false;
    const T89Parameters<float> p = t89_parameters<float>(bin);
    const std::size_t n = pos.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const std::array<float, 3> b =
            t89_components<float>(p, sin_tilt, cos_tilt, pos[(3 * i) + 0], pos[(3 * i) + 1],
                                  pos[(3 * i) + 2]);
        out[(3 * i) + 0] = b[0];
        out[(3 * i) + 1] = b[1];
        out[(3 * i) + 2] = b[2];
    }
    return true;
}

/// How many `float` scalars the device kernel's parameter buffer holds: `sin(psi)`, `cos(psi)`,
/// `C_1..C_19`, then the nine non-linear parameters. Asserted against the kernel registry below.
inline constexpr std::size_t t89_param_count = 2 + t89_linear_count + 9;

/**
 * Pack the epoch's tilt and a Kp bin's parameters into the kernel's parameter buffer.
 *
 * The layout is the kernel's ABI and is stated in exactly two places — here and the comment above
 * `irbem_t89_f32` in `irbem.slang`. A test evaluates both lanes on the same points, which is what
 * actually keeps the two statements in step.
 *
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`.
 * @param bin the Kp bin, `1..7`.
 * @return the parameter block, `t89_param_count` floats, by value.
 * @complexity O(1).
 * @alloc none — the block is the returned object's own inline array.
 * @test IrbemT89.ParameterBlockCarriesTheTiltThenTheCoefficients
 */
[[nodiscard]] inline std::array<float, t89_param_count> t89_param_block(float sin_tilt,
                                                                        float cos_tilt, int bin) {
    const T89Parameters<float> p = t89_parameters<float>(bin);
    std::array<float, t89_param_count> block{};
    block[0] = sin_tilt;
    block[1] = cos_tilt;
    for (std::size_t k = 0; k < t89_linear_count; ++k) block[2 + k] = p.c[k];
    block[21] = p.delta_x;
    block[22] = p.a_rc;
    block[23] = p.d_0;
    block[24] = p.gamma_rc;
    block[25] = p.r_c;
    block[26] = p.g;
    block[27] = p.a_t;
    block[28] = p.d_y;
    block[29] = p.x_0;
    return block;
}

/**
 * The batch's position caveat, accumulated one point at a time.
 *
 * A batch returns ONE @ref Status for N points, so it can only be the worst of them — and it must
 * not be better than the worst, or a caller who checks it learns nothing about the point that was
 * wrong. Computing that honestly must not cost a second pass over the positions: at 2^22 points the
 * input is 100 MB and re-reading it is a measurable fraction of the whole call. So this folds into
 * the loop that is already touching each point, and it does so without a `sqrt`: T89's envelope is
 * a radius band plus a sunward `x` bound, all three monotone, so the smallest radius, the largest
 * radius and the smallest `x` decide the entire batch and the radii can be compared as SQUARES.
 * Two `check_position` calls per CALL then reproduce exactly what N of them would have said.
 *
 * @test IrbemT89.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
struct T89PositionFold {
    /// The smallest `r^2` seen, R_E^2; `+inf` until the first point.
    double r2_lo = std::numeric_limits<double>::infinity();
    /// The largest `r^2` seen, R_E^2; zero until the first point.
    double r2_hi = 0.0;
    /// The smallest GSM `x` seen, R_E; `+inf` until the first point.
    double x_lo = std::numeric_limits<double>::infinity();
    /// False once any point has had a non-finite coordinate. Tracked separately because a NaN
    /// radius compares false against everything and would otherwise slip through both extremes.
    bool finite = true;

    /**
     * Fold one position in.
     *
     * @param p the position, GSM, in Earth radii.
     * @complexity O(1) — one fused radius, three comparisons, no `sqrt` and no branch that a
     *             compiler cannot turn into a select.
     * @alloc none.
     * @test IrbemT89.BatchReportsTheSameEnvelopeTheScalarLaneDoes
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
     *
     * @return @ref Status::DomainError when any point is not finite or is inside the Earth,
     *         @ref Status::OutOfValidityRange when any point is outside T89's published
     *         `r_GEO <= 70 R_E`, otherwise @ref Status::Ok.
     * @complexity O(1) — two square roots and two envelope lookups for the whole batch.
     * @alloc none.
     * @test IrbemT89.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    [[nodiscard]] Status verdict() const {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double lo = finite ? std::sqrt(r2_lo) : nan;
        const double hi = finite ? std::sqrt(r2_hi) : nan;
        return first_failure(check_position(ExternalModel::Tsyganenko1989, lo, x_lo),
                             check_position(ExternalModel::Tsyganenko1989, hi, x_lo));
    }
};

/**
 * The T89 field over a whole batch of GSM points, on the device when that is worth it.
 *
 * **This is the routine to call for more than a handful of points.** @ref t89_field is the reference
 * lane: it is what the batch is verified against, and what runs when there is no device or the batch
 * is too small to pay for one.
 *
 * T89 sits in the regime where the device wins. It is ~400 flops for 24 bytes in and 12 out —
 * ~11 flops/byte, an order of magnitude above the streaming dipole kernel that LOSES 0.69x on this
 * same seam and within a factor of two of IGRF, which wins 8.96x. Nothing is stored per point beyond
 * the answer, there is no branch on data, and the 30-float parameter block is read identically by
 * every lane in a workgroup.
 *
 * **The batch reports the same caveats the scalar entry point does, folded over the whole batch.**
 * One status for N points can only be the worst of them, so that is what it is: if any point is
 * beyond T89's published `r_GEO <= 70 R_E`, the batch says @ref Status::OutOfValidityRange, and if
 * any point is inside the Earth or not finite it says @ref Status::DomainError. An
 * out-of-validity batch is still computed in full — `status.hpp`'s standing rule — and a
 * DOMAIN-ERROR batch is zeroed in full, which is exactly what @ref t89_field does with its one
 * point. Checking positions costs the batch two comparisons per point and two envelope lookups per
 * CALL: the radii are folded as squares, so there is no per-point `sqrt` and the device lane's
 * measured throughput is unchanged.
 *
 * @param points the positions, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians.
 * @param kp_times_ten Kp in IRBEM's slot-1 scaling; binned by @ref t89_kp_bin.
 * @param out receives one field vector per input, GSM, nanotesla; same length as @p points.
 * @return @ref Status::DomainError on a length mismatch, a tilt at which the model has no value, or
 *         a point that is not finite or is inside the Earth, and then every output is zeroed;
 *         @ref Status::OutOfValidityRange when Kp or any point's radius is outside the published
 *         range, with every point still computed;
 *         otherwise @ref Status::Ok. The value is `true` exactly when the device
 *         lane serviced the call — a test asserts this rather than trusting that a GPU was used,
 *         because a silent fallback is what makes a performance claim worthless.
 * @complexity O(N); on the device those N run concurrently over `ceil(N/256)` workgroups.
 * @alloc the device lane stages positions and results into two `std::vector<float>` of `3N`; the
 *        host lane allocates nothing.
 * @test IrbemT89.BatchAgreesWithTheReferenceLane
 * @test IrbemT89.BatchRejectsMismatchedSpans
 * @test IrbemT89.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
[[nodiscard]] inline Result<bool> t89_field_batch(
    std::span<const Position<Frame::GSM>> points, double tilt_rad, double kp_times_ten,
    std::span<FieldVector<Frame::GSM>> out) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(tilt_rad) || !std::isfinite(kp_times_ten)) return {Status::DomainError, false};
    if (!(std::fabs(tilt_rad) < max_tilt_rad)) return {Status::DomainError, false};

    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Kp)] = kp_times_ten;
    const Status drives = check_validity(ExternalModel::Tsyganenko1989, drivers);
    if (n == 0) return {drives, false};

    const int bin = t89_kp_bin(kp_times_ten);
    const double sin_tilt = std::sin(tilt_rad);
    const double cos_tilt = std::cos(tilt_rad);
    // The positions are folded INSIDE the loop that is already reading them — see
    // T89PositionFold. A separate pass would re-read 100 MB at 2^22 points.
    T89PositionFold fold;

#if CHEATAH_SPACE_IRBEM_T89_GPU
    if (gpu::prefer_gpu("irbem_t89_f32", n) &&
        std::filesystem::exists(gpu::shader_path("irbem_t89_f32"))) {
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
        const std::array<float, t89_param_count> block =
            t89_param_block(static_cast<float>(sin_tilt), static_cast<float>(cos_tilt), bin);
        gpu::dispatch_batch("irbem_t89_f32", pos, raw, std::span<const float>(block));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = FieldVector<Frame::GSM>{
                fixarray::vec3d{raw[(3 * i) + 0], raw[(3 * i) + 1], raw[(3 * i) + 2]}};
        }
        return {first_failure(drives, where), true};
    }
#endif

    const T89Parameters<double> p = t89_parameters<double>(bin);
    for (std::size_t i = 0; i < n; ++i) {
        fold.add(points[i]);
        const std::array<double, 3> b = t89_components<double>(
            p, sin_tilt, cos_tilt, points[i].v[0], points[i].v[1], points[i].v[2]);
        out[i] = FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
    }
    const Status where = fold.verdict();
    if (where == Status::DomainError) {
        for (std::size_t i = 0; i < n; ++i) out[i] = FieldVector<Frame::GSM>{};
        return {Status::DomainError, false};
    }
    return {first_failure(drives, where), false};
}

}  // namespace cheatah::space::irbem
