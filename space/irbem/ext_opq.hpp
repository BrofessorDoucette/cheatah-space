#pragma once

/**
 * @file ext_opq.hpp
 * @brief space.irbem — Olson & Pfitzer (1977), the quiet-time EXTERNAL magnetospheric field:
 *        IRBEM's `kext = 5`, and the one model its own `LANDI2LSTAR` fast path hard-wires.
 *
 * ## What the model is
 *
 * Olson, W. P. and Pfitzer, K. A., *Magnetospheric magnetic field modeling*, Annual Scientific
 * Report, AFOSR contract F44620-75-C-0033, McDonnell Douglas Astronautics Co. (January 1977),
 * DTIC AD-A037492. The report's §3 describes the construction and its Appendix publishes the model
 * IN FULL as the subroutine `BXYZMU` (report pp. 64–67): the functional form, the enumeration of
 * terms, the tilt expansion and every one of the 344 fitted coefficients. Nothing here is taken
 * from anywhere else, and nothing is taken from IRBEM's source.
 *
 * The field of the magnetopause, tail and ring currents — no internal field — is a power series in
 * solar-magnetic (SM) Cartesian position with an exponential envelope, fitted to over 600 OGO-3
 * and OGO-5 vector measurements taken under quiet conditions. In the report's own words (p. 64):
 * "a generalized orthonormal least squares program was used to fit the coefficients of a power
 * series (including exponential terms) through fourth order in space and second order in tilt".
 * Written out, with `E = exp(-0.06 r^2)` and `r^2 = x^2 + y^2 + z^2` in Earth radii:
 *
 *     B_x = sum_i (A_i + B_i E) x^p y^{2q} z^s          32 terms, i = 1..32
 *     B_y = y * sum_j (C_j + D_j E) x^p y^{2q} z^s      22 terms, j = 1..22
 *     B_z = sum_i (E_i + F_i E) x^p y^{2q} z^s          32 terms, i = 1..32
 *
 * where the monomials are enumerated by the report's loop (x power outermost, `y^2` power next,
 * `z` power innermost, with the truncation rules that stop each block — see @ref opq_components,
 * which IS that loop). The even-in-`y` factor on `B_x`, `B_z` and the odd one on `B_y` are the
 * dawn-dusk mirror symmetry of a quiet magnetosphere, built in rather than fitted.
 *
 * **Tilt is a polynomial, in degrees.** Every coefficient is a pair from the published tables:
 * `A_i = AA_{2i-1} + AA_{2i} t^2` for a term that is even in the dipole tilt `t`, or
 * `A_i = AA_{2i-1} t + AA_{2i} t^3` for one that is odd, with `t` in DEGREES — the report's
 * `TT(1..4) = {1, t, t^2, t^3}` and its `ITA / ITB / ITC` selector tables (@ref opq_table). Which
 * parity each term has is not arbitrary: at zero tilt the magnetosphere is a mirror about the SM
 * equator, so a term whose `z` power has the wrong parity for its component can only enter
 * multiplied by an odd power of the tilt. @ref IrbemOpq.ParityTablesFollowTheEquatorialMirror
 * derives every one of the 86 selector entries from that rule and compares.
 *
 * **Three radial rules, all published (report p. 66):** the field is set to zero below 2 R_E ("the
 * main field dominates and the variations expressed by this expansion are not sufficiently
 * accurate"), tapered linearly in `r^2` between 2 and 2.5 R_E, and set to zero beyond 15 R_E,
 * where "the field diverges rapidly and a template sets the field to zero". That 15 R_E is the
 * validity limit `status.hpp` carries for `kext = 5`; the template is why the VALUE beyond it is a
 * zero external field and not an extrapolation.
 *
 * **No drivers.** This is a quiet-time model: it reads no Kp, no Dst, no solar wind. The envelope
 * row in `status.hpp` lists no driver, @ref opq_field's context overload reads nothing but the
 * tilt, and @ref IrbemOpq.TheModelReadsNoDrivers plus the oracle sweep in
 * @ref IrbemOpq.HeavyDifferentialOracleIgnoresEveryDriverRegime show that IRBEM's `kext = 5`
 * behaves the same way across the corpus's four activity regimes to the last bit.
 *
 * ## Provenance — MEASURED, NOT ASSUMED: both the form and the coefficients are published
 *
 * The T89 lesson (see `ext_t89.hpp`) is that a model's name is not its equations. So before this
 * file existed the oracle was interrogated as a black box: its external field was isolated as
 * `kext = 5` minus `kext = 0` (the internal IGRF term cancels exactly), and because the published
 * form is LINEAR in every coefficient, a least-squares solve of the published basis against the
 * oracle's field at 400 scattered SM points recovers the folded per-tilt coefficients exactly if —
 * and only if — the oracle evaluates that form. Measured by
 * [`tools/oracle/opq_diff.cpp`](../../../tools/oracle/opq_diff.cpp) against the `-O2` oracle at
 * six dipole tilts from -30.4 to +25.6 degrees: the RMS residual of that solve is **2e-12 nT**
 * (1.5e-13 of the field), i.e. fp64 roundoff. Fitting each recovered coefficient's tilt dependence
 * across the six tilts as a cubic separates the pairs and their parities: the part of every fit
 * the published parity says is zero is **1.1e-9** of the part it says is not, and every one of the
 * 344 recovered values sits within **3.3e-10** (relative) of a six-significant-figure decimal —
 * the decimals printed in the report. IRBEM's `kext = 5` is therefore the published model,
 * coefficient for coefficient, and this header's target is ORACLE PARITY at the internal field's
 * standard (1e-6 relative on B), not a documented gap.
 *
 * And it gets there. The same harness's deviation pass, this header against the oracle at those
 * three tilts, 600 scattered points per region:
 *
 * | region | RMS \|dB\| | max \|dB\| | max relative |
 * |---|---|---|---|
 * | radiation belts, 2.6–10 R_E | 3.5e-14 nT | 2.3e-13 nT | 2.9e-14 |
 * | whole published region, 2.6–15 R_E | 3.5e-14 nT | 3.0e-13 nT | 2.1e-14 |
 * | taper and inner zero, 1.2–2.6 R_E | 1.2e-13 nT | 5.8e-13 nT | 7.1e-12 |
 *
 * against external fields of 40–150 nT: fp64 roundoff in every region, the taper included. The
 * oracle's field is bit-identical across the corpus's four activity regimes (worst
 * \|dB\| = 0.0 nT), and beyond 15 R_E it refuses outright where the published template gives
 * zero — the one difference, by design on both sides, and stated. Carried through the drift-shell
 * chain (`make_lstar` over @ref TotalFieldOpq against IRBEM's `make_lstar1` at `kext = 5`), the
 * worst \|dL*\| over L = 3–6.6 is **0.0062 at matched `options(3,4) = 9`** — the same 0.0066 the
 * internal field alone shows, i.e. the numerics and nothing else — and 0.012 at the matched
 * default resolution, inside IRBEM's own documented 0.010–0.017 discretization error there.
 *
 * The tables below were transcribed from the archive.org scan of the report (pp. 64–66) and
 * every entry — the scan's typeface makes 4 and 9 hard to tell apart — was then checked against
 * that recovery, which is what makes the transcription a measurement rather than a reading. The
 * parity tables recovered from the tilt fit agree with the printed `ITA / ITB / ITC` exactly.
 *
 * ## div B is NOT identically zero here, and the suite says by how much
 *
 * Unlike T89, whose field is the curl of a vector potential, this is a least-squares fit of
 * independent power series to each component, so `div B` is a small residual of the FIT rather
 * than zero. Because the form is a polynomial times an exponential its divergence is available in
 * closed form, and @ref IrbemOpq.StencilDivergenceConvergesToTheAnalyticOne checks that a
 * second-order central difference of the implemented field falls as `h^2` onto that closed form —
 * over three decades of `h` — which is what verifies that every term of the evaluator is the term
 * it claims to be, with no reference model in the loop: measured, the worst gap falls
 * 1.2e-2 → 1.2e-4 → 1.2e-6 → 1.2e-8 nT/R_E as `h` goes 1e-1 → 1e-4. The intrinsic divergence
 * itself is measured and stated by that test rather than hidden, and it is not small: over the
 * belts (r <= 8 R_E) the fit's \|div B\| reaches 12 nT/R_E near the polar axis at 3 R_E (3% of the
 * dipole's own \|B\|/r there, ~1x the EXTERNAL field's \|B\|/r), and past 12 R_E — where the report
 * itself says the series "diverges rapidly" — it exceeds the dipole gradient by an order of
 * magnitude. That is a property of the published fit the oracle shares to the last bit; a caller
 * integrating this model's field lines beyond ~10 R_E should know it.
 *
 * @note Nothing on a hot path allocates, and nothing but the device lane can throw. The evaluator,
 *       the scalar entry points and the fp32 host lane touch nothing but the caller's spans and a
 *       stack parameter block; @ref IrbemOpq.NothingOnTheHeapInTheHotPath counts. The one exception
 *       is @ref opq_field_batch's DEVICE lane, which stages `3N` floats each way and forwards
 *       whatever `gpu::dispatch_batch` throws; its `@alloc` and the dispatch header say so.
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
#include "lstar.hpp"
#include "status.hpp"

// The device lane is opt-in by include path, exactly as ext_t89.hpp's is. Without
// cheatah-gpu-linalg every routine here still compiles and runs on the host.
#if __has_include("cheatah_gpu_linalg/context.hpp")
#  include "gpu/dispatch.hpp"
/// 1 when this translation unit can reach the OPQ device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_OPQ_GPU 1
#else
/// 1 when this translation unit can reach the OPQ device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_OPQ_GPU 0
#endif

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// The shape of the expansion — every number the report states beside its tables
// -------------------------------------------------------------------------------------------

/// How many monomials the `B_x` and `B_z` series carry: the report's `A(32)`, `E(32)`.
inline constexpr std::size_t opq_xz_count = 32;
/// How many monomials the `B_y` series carries (before its overall factor of `y`): `C(22)`.
inline constexpr std::size_t opq_y_count = 22;
/// The rate of the exponential envelope: `E = exp(-0.06 r^2)`, `r` in Earth radii (report p. 67).
inline constexpr double opq_exp_rate = 0.06;
/// Below this `r^2` (2 R_E) the published field is identically zero (report p. 66).
inline constexpr double opq_inner_r2 = 4.0;
/// Below this `r^2` (2.5 R_E) the field is tapered linearly in `r^2` towards the inner zero.
inline constexpr double opq_taper_r2 = 6.25;
/// Above this `r^2` (15 R_E) the published template sets the field to zero (report p. 66).
inline constexpr double opq_template_r2 = 225.0;

/**
 * The published model in full — Olson & Pfitzer (1977), Appendix, subroutine `BXYZMU`, pp. 64–66.
 *
 * Each coefficient of the space expansion is a PAIR of table entries combined with the tilt:
 * entries `2i-1` and `2i` (one-based) of `aa` give `A_i`, and so on. The parity tables say how:
 * `false` (the report's selector value 1) pairs them as `{1, t^2}`, `true` (selector 2) as
 * `{t, t^3}`, `t` being the dipole tilt in degrees. `A`/`B` share `odd_x` (the report's `ITA`),
 * `C`/`D` share `odd_y` (`ITB`), `E`/`F` share `odd_z` (`ITC`).
 *
 * Every value is the six-significant-figure decimal printed in the report, and every value has
 * been checked against a black-box recovery from the IRBEM oracle to better than 1.2e-9 relative
 * (see the file brief) — a transcription check that cannot pass by accident.
 *
 * @test IrbemOpq.TablesHaveThePublishedShape
 * @test IrbemOpq.ParityTablesFollowTheEquatorialMirror
 */
struct OpqTable {
    /// `AA(64)`: the `B_x` non-exponential pairs.
    std::array<double, 2 * opq_xz_count> aa;
    /// `BB(64)`: the `B_x` exponential-envelope pairs.
    std::array<double, 2 * opq_xz_count> bb;
    /// `CC(44)`: the `B_y` non-exponential pairs.
    std::array<double, 2 * opq_y_count> cc;
    /// `DD(44)`: the `B_y` exponential-envelope pairs.
    std::array<double, 2 * opq_y_count> dd;
    /// `EE(64)`: the `B_z` non-exponential pairs.
    std::array<double, 2 * opq_xz_count> ee;
    /// `FF(64)`: the `B_z` exponential-envelope pairs.
    std::array<double, 2 * opq_xz_count> ff;
    /// `ITA`: whether each `B_x` term is odd in the tilt.
    std::array<bool, opq_xz_count> odd_x;
    /// `ITB`: whether each `B_y` term is odd in the tilt.
    std::array<bool, opq_y_count> odd_y;
    /// `ITC`: whether each `B_z` term is odd in the tilt.
    std::array<bool, opq_xz_count> odd_z;
};

/**
 * Olson & Pfitzer (1977) pp. 64–66, as printed — see @ref OpqTable for the layout and the file
 * brief for how each entry was verified.
 *
 * @test IrbemOpq.TablesHaveThePublishedShape
 */
inline constexpr OpqTable opq_table{
    // AA(64), report p. 64
    {-2.26836e-02, -1.01863e-04, 3.42986e+00,  -3.12195e-04, 9.50629e-03,  -2.91512e-06,
     -1.57317e-03, 8.62856e-08,  -4.26478e-05, 1.62924e-08,  -1.27549e-04, 1.90732e-06,
     -1.65983e-02, 8.46680e-09,  -5.55850e-05, 1.37404e-08,  9.91815e-05,  1.59296e-08,
     4.52864e-07,  -7.17669e-09, 4.98627e-05,  3.33662e-10,  -5.97620e-02, 1.60669e-05,
     -2.29457e-01, -1.43777e-04, 1.09403e-03,  -9.15606e-07, 1.60658e-03,  -4.01198e-07,
     -3.15064e-06, 2.03125e-09,  4.92887e-04,  -1.80676e-07, -1.12022e-03, 5.98568e-07,
     -5.90009e-06, 5.16504e-09,  -1.48737e-06, 4.83477e-10,  -7.44379e-04, 3.82472e-06,
     7.41737e-04,  -1.31468e-05, -1.24729e-04, 1.92930e-08,  -1.91764e-04, -5.30371e-08,
     1.38186e-05,  -2.81594e-08, 7.46386e-06,  2.64404e-08,  2.45049e-04,  -1.81802e-07,
     -1.00278e-03, 1.98742e-06,  -1.16425e-05, 1.17556e-08,  -2.46079e-06, -3.45831e-10,
     1.02440e-05,  -1.90716e-08, -4.00855e-05, 1.25818e-07},
    // BB(64), report p. 65
    {9.47753e-02,  1.45981e-04,  -1.82933e+00, 5.54882e-04,  5.03665e-03,  -2.07698e-06,
     1.10959e-01,  -3.45837e-05, -4.40075e-05, 5.06464e-07,  -1.20112e-03, 3.64911e-06,
     1.49849e-01,  -7.44929e-05, 2.46382e-04,  9.65870e-07,  -9.54881e-04, 2.43647e-07,
     3.06520e-04,  3.07836e-07,  6.48301e-03,  1.26251e-06,  -7.09548e-03, -1.55596e-05,
     3.06465e+00,  -7.84893e-05, 4.95145e-03,  3.71921e-06,  -1.52002e-01, 6.81988e-06,
     -8.55686e-05, -9.01230e-08, -3.71458e-04, 1.30476e-07,  -1.82971e-01, 1.51390e-05,
     -1.45912e-04, -2.22778e-07, 6.49278e-05,  -3.72758e-08, -1.59932e-03, 8.04921e-06,
     5.38012e-01,  -1.43182e-04, 1.50000e-04,  5.88020e-07,  -1.59000e-02, 1.60744e-06,
     3.17837e-04,  1.78959e-07,  -8.93794e-03, 6.37549e-06,  1.27887e-03,  -2.45878e-07,
     -1.93210e-01, 6.91233e-06,  -2.80637e-04, -2.57073e-07, 5.78343e-05,  4.52128e-10,
     1.89621e-04,  -4.84911e-08, -1.50058e-02, 6.21772e-06},
    // CC(44), report p. 65
    {-1.88177e-02, -1.92493e-06, -2.89064e-01, -8.49439e-05, -4.76380e-04, -4.52998e-08,
     1.61086e-03,  3.18728e-07,  1.29159e-06,  5.52259e-10,  3.95543e-05,  5.61209e-08,
     1.38287e-03,  5.74237e-07,  1.86489e-06,  7.10175e-10,  1.45243e-07,  -2.97591e-10,
     -2.43029e-03, -6.70000e-07, -2.30624e-02, -6.22193e-06, -2.40815e-05, 2.01689e-08,
     1.76721e-04,  3.78689e-08,  9.88496e-06,  7.33820e-09,  7.32126e-05,  8.43986e-08,
     8.82449e-06,  -6.11708e-08, 1.78881e-04,  8.62589e-07,  3.43724e-06,  2.53783e-09,
     -2.04239e-07, 8.16641e-10,  1.68075e-05,  7.62815e-09,  2.26026e-04,  3.66341e-08,
     3.44637e-07,  2.25531e-10},
    // DD(44), report p. 65
    {2.50143e-03,  1.01200e-06,  3.23821e+00,  1.08589e-05,  -3.39199e-05, -5.27052e-07,
     -9.46161e-02, -1.95413e-09, -4.23614e-06, 1.43153e-08,  -2.62948e-04, 1.05138e-07,
     -2.15784e-01, -2.20717e-07, -2.65687e-05, 1.26370e-08,  5.88917e-07,  -1.13658e-08,
     1.64385e-03,  1.44263e-06,  -1.66045e-01, -1.46096e-05, 1.22811e-04,  3.43922e-08,
     9.66760e-05,  -6.32150e-07, -4.97400e-05, -2.78578e-08, 1.77366e-02,  2.05401e-07,
     -1.91756e-03, -9.49392e-07, -1.99488e-01, -2.07170e-06, -5.40443e-05, 1.59289e-08,
     7.30914e-05,  3.38786e-08,  -1.59537e-04, -1.65504e-07, 1.90940e-02,  2.03238e-06,
     1.01148e-04,  5.20815e-08},
    // EE(64), report p. 65
    {-2.77924e+01, -1.01457e-03, 9.21436e-02,  -8.52177e-06, 5.19106e-01,  8.28881e-05,
     -5.59651e-04, 1.16736e-07,  -2.11206e-03, -5.35469e-07, 4.41990e-01,  -1.33679e-05,
     -7.18642e-04, 6.17358e-08,  -3.51990e-03, -5.29070e-07, 1.88443e-06,  -6.60696e-10,
     -1.34708e-03, 1.02160e-07,  1.58219e-06,  2.05040e-10,  1.18039e+00,  1.58903e-04,
     1.86944e-02,  -4.46477e-06, 5.49869e-02,  4.94690e-06,  -1.18335e-04, 6.95684e-09,
     -2.73839e-04, -9.17883e-08, 2.79126e-02,  -1.02567e-05, -1.25427e-04, 3.07143e-08,
     -5.31826e-04, -2.98476e-08, -4.89899e-05, 4.91480e-08,  3.85563e-01,  4.16966e-05,
     6.74744e-04,  -2.08736e-07, -3.42654e-03, -3.13957e-06, -6.31361e-06, -2.92981e-09,
     -2.63883e-03, -1.32235e-07, -6.19406e-06, 3.54334e-09,  6.65986e-03,  -5.81949e-06,
     -1.88809e-04, 3.62055e-08,  -4.64380e-04, -2.21159e-07, -1.77496e-04, 4.95560e-08,
     -3.18867e-04, -3.17697e-07, -1.05815e-05, 2.22220e-09},
    // FF(64), report pp. 65–66
    {-5.07092e+00, 4.71960e-03,  -3.79851e-03, -3.67309e-06, -6.02439e-01, 1.08490e-04,
     5.09287e-04,  5.62210e-07,  7.05718e-02,  5.13160e-06,  -2.85571e+00, -4.31728e-05,
     1.03185e-03,  1.05332e-07,  1.04106e-02,  1.60749e-05,  4.18031e-05,  3.32759e-08,
     1.20113e-01,  1.40486e-05,  -3.37993e-05, 5.48340e-09,  9.10815e-02,  -4.00608e-04,
     3.75393e-03,  -4.69939e-07, -2.48561e-02, 1.31836e-04,  -2.67755e-04, -7.60285e-08,
     3.04443e-03,  -3.28956e-06, 5.82367e-01,  5.39496e-06,  -6.15261e-04, 4.05316e-08,
     1.13546e-02,  -4.26493e-06, -2.72007e-02, 5.72523e-08,  -2.98576e+00, 3.07325e-05,
     1.51645e-03,  1.25098e-06,  4.07213e-02,  1.05964e-05,  1.04232e-04,  1.77381e-08,
     1.92781e-01,  2.15734e-05,  -1.65741e-05, -1.88683e-09, 2.44803e-01,  1.51316e-05,
     -3.01157e-04, 8.47006e-08,  1.86971e-02,  -6.94074e-06, 9.13198e-03,  -2.38052e-07,
     1.28552e-01,  6.92595e-06,  -8.36792e-05, -6.10021e-08},
    // ITA(32), report p. 64: 2 = odd in the tilt, 1 = even
    {true,  false, true,  false, true,  true,  false, true,  false, true,  false, true,
     false, true,  false, true,  true,  false, true,  true,  true,  false, true,  false,
     true,  false, true,  false, true,  true,  true,  false},
    // ITB(22)
    {true,  false, true,  false, true,  true,  false, true,  true,  true,  false, true,
     false, true,  false, true,  false, true,  true,  true,  false, true},
    // ITC(32)
    {false, true,  false, true,  false, false, true,  false, true,  false, true,  false,
     true,  false, true,  false, false, true,  false, false, false, true,  false, true,
     false, true,  false, true,  false, false, false, true},
};

// -------------------------------------------------------------------------------------------
// The per-tilt fold
// -------------------------------------------------------------------------------------------

/**
 * One tilt's worth of space-expansion coefficients — the report's `A, B, C, D, E, F` arrays
 * after its `DO 1` / `DO 2` loops have folded the tilt in.
 *
 * Templated on the scalar type for the same reason `T89Parameters` is: the device kernel receives
 * these as `float`, so the host lane it is verified against must fold in `double`, round to
 * `float` FIRST and then evaluate, or a disagreement between the lanes cannot be attributed.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` for the device-mirroring one.
 * @test IrbemOpq.ParametersRoundTripThroughFloat
 */
template <std::floating_point T>
struct OpqParameters {
    /// `A_i`: the `B_x` coefficients of the plain monomials.
    std::array<T, opq_xz_count> a;
    /// `B_i`: the `B_x` coefficients of the monomials under the exponential envelope.
    std::array<T, opq_xz_count> b;
    /// `C_j`: the `B_y` coefficients of the plain monomials (times `y`).
    std::array<T, opq_y_count> c;
    /// `D_j`: the `B_y` coefficients of the enveloped monomials (times `y`).
    std::array<T, opq_y_count> d;
    /// `E_i`: the `B_z` coefficients of the plain monomials.
    std::array<T, opq_xz_count> e;
    /// `F_i`: the `B_z` coefficients of the enveloped monomials.
    std::array<T, opq_xz_count> f;
};

/**
 * Fold the dipole tilt into the published pairs — the report's `DO 1` and `DO 2` loops.
 *
 * `TT = {1, t, t^2, t^3}` with `t` the tilt in DEGREES; an even term takes
 * `pair[0] + pair[1] t^2`, an odd term `pair[0] t + pair[1] t^3`. Done in `double` whatever @p T
 * is, then rounded: the rounding on the way in is what makes the `float` lane the device's twin.
 *
 * @tparam T the scalar type of the result; `double` or `float`.
 * @param tilt_deg the dipole tilt `psi` in degrees, positive when the north dipole leans sunward.
 * @return the folded coefficient set, by value.
 * @complexity O(1) — 172 fused pairs, no transcendental.
 * @alloc none; the returned object is inline storage.
 * @test IrbemOpq.FoldedCoefficientsAreEvenOrOddInTheTilt
 * @test IrbemOpq.ParametersRoundTripThroughFloat
 */
template <std::floating_point T>
[[nodiscard]] constexpr OpqParameters<T> opq_parameters(double tilt_deg) {
    const double t = tilt_deg;
    const double t2 = t * t;
    const double t3 = t2 * t;
    const auto fold = [&](const double* pair, bool odd) {
        return odd ? (pair[0] * t) + (pair[1] * t3) : pair[0] + (pair[1] * t2);
    };
    OpqParameters<T> out{};
    for (std::size_t i = 0; i < opq_xz_count; ++i) {
        const std::size_t j = 2 * i;
        out.a[i] = static_cast<T>(fold(&opq_table.aa[j], opq_table.odd_x[i]));
        out.b[i] = static_cast<T>(fold(&opq_table.bb[j], opq_table.odd_x[i]));
        out.e[i] = static_cast<T>(fold(&opq_table.ee[j], opq_table.odd_z[i]));
        out.f[i] = static_cast<T>(fold(&opq_table.ff[j], opq_table.odd_z[i]));
    }
    for (std::size_t i = 0; i < opq_y_count; ++i) {
        const std::size_t j = 2 * i;
        out.c[i] = static_cast<T>(fold(&opq_table.cc[j], opq_table.odd_y[i]));
        out.d[i] = static_cast<T>(fold(&opq_table.dd[j], opq_table.odd_y[i]));
    }
    return out;
}

// -------------------------------------------------------------------------------------------
// The evaluator — the report's loop, in SM
// -------------------------------------------------------------------------------------------

/**
 * The OP-77 external field at one SOLAR-MAGNETIC point, as three SM components in nanotesla.
 *
 * This is the report's subroutine, block for block: the radial rules first (zero inside 2 R_E, a
 * linear-in-`r^2` taper to 2.5 R_E, the zero template beyond 15 R_E), then the envelope
 * `E = exp(-0.06 r^2)`, then the series in the report's enumeration. That enumeration is a fixed
 * pattern that reads as three nested loops with two truncation rules — the trip count depends on
 * nothing but the loop counters, so the work is the same at every point and a device lane does
 * not diverge:
 *
 *  - `x` power `p = 0..4` outermost; `y^2` power `q = 0..2` next, cut off when `p + 2q > 5`;
 *  - `z` power `s` innermost, from 0, cut off when `s > 4` OR the running index `p + 2q + s`
 *    passes 5 — with `B_y` getting one term FEWER than `B_x`/`B_z` in each block, because the
 *    overall factor of `y` costs it one order.
 *
 * Written with the report's own bookkeeping counters (`ijk`, `k`) rather than the closed-form
 * rules above, so that the 32 + 22 + 32 terms come out in exactly the order the tables are laid
 * out in — the counts are asserted, not assumed.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device kernel.
 * @param p the folded coefficients for the epoch's tilt; see @ref opq_parameters.
 * @param x the SM x coordinate, R_E (sunward, in the plane of the dipole axis and the Sun line).
 * @param y the SM y coordinate, R_E (duskward).
 * @param z the SM z coordinate, R_E (along the north dipole axis).
 * @return `{B_x, B_y, B_z}` in SM, nanotesla. Zero — exactly — inside 2 R_E and beyond 15 R_E,
 *         and for any input whose `r^2` is not a number, so that a NaN can never propagate out.
 * @complexity O(1) — 86 fused terms, ~400 flops, one `exp`, no square root and no trigonometry.
 * @alloc none.
 * @test IrbemOpq.MonomialEnumerationHasThePublishedCounts
 * @test IrbemOpq.ZeroTiltIsMirrorSymmetricAboutTheEquator
 * @test IrbemOpq.DawnDuskSymmetryHoldsAtEveryTilt
 * @test IrbemOpq.TheThreeRadialRulesArePublishedOnes
 * @test IrbemOpq.StencilDivergenceConvergesToTheAnalyticOne
 * @test IrbemOpq.TiltIsAContinuousParameter
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, 3> opq_components(const OpqParameters<T>& p, T x, T y, T z) {
    const T zero = static_cast<T>(0);
    const T one = static_cast<T>(1);
    const T y2 = y * y;
    const T r2 = (x * x) + y2 + (z * z);
    // The template beyond 15 R_E and the hard zero inside 2 R_E; the negated comparison also
    // catches a NaN radius, which is the only way a NaN could otherwise reach the series.
    if (!(r2 <= static_cast<T>(opq_template_r2)) || r2 < static_cast<T>(opq_inner_r2)) {
        return {zero, zero, zero};
    }
    // The taper: linear in r^2 from 0 at 2 R_E to 1 at 2.5 R_E, the report's CON.
    const T con = r2 < static_cast<T>(opq_taper_r2)
                      ? (r2 - static_cast<T>(opq_inner_r2)) /
                            (static_cast<T>(opq_taper_r2) - static_cast<T>(opq_inner_r2))
                      : one;
    const T expr = std::exp(-static_cast<T>(opq_exp_rate) * r2);

    T bx = zero;
    T by = zero;
    T bz = zero;
    std::size_t ix = 0;  // the next B_x / B_z term
    std::size_t iy = 0;  // the next B_y term
    T xb = one;          // x^p
    for (int i = 1; i <= 5; ++i) {
        T yexb = xb;     // x^p y^{2q}
        for (int j = 1; j <= 3; ++j) {
            if (i + (2 * j) > 8) break;  // the report's "IF (I+2*J .GT. 8) GO TO 30"
            int ijk = i + (2 * j) + 1;
            int k = 0;
            T zeyexb = yexb;             // x^p y^{2q} z^s
            for (;;) {
                bz += (p.e[ix] + (p.f[ix] * expr)) * zeyexb;
                bx += (p.a[ix] + (p.b[ix] * expr)) * zeyexb;
                ++ix;
                if (ijk > 8) break;      // "IF (IJK .GT. 8) GO TO 20": B_y is one order short
                by += (p.c[iy] + (p.d[iy] * expr)) * zeyexb * y;
                ++iy;
                zeyexb *= z;
                ++ijk;
                ++k;
                if (ijk > 9 || k > 4) break;       // "IF (IJK.LE.9 .AND. K.LE.4) GO TO 10"
            }
            yexb *= y2;
        }
        xb *= x;
    }
    return {bx * con, by * con, bz * con};
}

/**
 * The OP-77 external field at one GSM point, in `double`, from an already-folded parameter set.
 *
 * The rotation in and out is the one every Tsyganenko-family evaluator makes: about `y` by the
 * tilt, GSM to SM, because the report's series is defined in SM (its `z` IS the dipole axis).
 *
 * @param p the position, GSM, in Earth radii.
 * @param par the coefficients folded for this epoch's tilt — @ref opq_parameters. Folding is
 *        172 fused pairs, so a caller evaluating many points at one epoch does it once.
 * @param sin_tilt `sin(psi)`, the dipole tilt's sine; @ref HotState::sin_tilt holds it.
 * @param cos_tilt `cos(psi)`.
 * @return the external field at @p p, GSM, nanotesla.
 * @complexity O(1); see @ref opq_components, plus two rotations.
 * @alloc none.
 * @test IrbemOpq.ReferenceLaneMatchesTheComponentForm
 * @test IrbemOpq.HeavyDifferentialAgreesWithTheIrbemOracle
 * @test IrbemOpq.HeavyDifferentialOracleIgnoresEveryDriverRegime
 */
[[nodiscard]] inline FieldVector<Frame::GSM> opq_field_at(Position<Frame::GSM> p,
                                                          const OpqParameters<double>& par,
                                                          double sin_tilt, double cos_tilt) {
    const double xs = (p.v[0] * cos_tilt) - (p.v[2] * sin_tilt);
    const double zs = (p.v[0] * sin_tilt) + (p.v[2] * cos_tilt);
    const std::array<double, 3> b = opq_components<double>(par, xs, p.v[1], zs);
    return FieldVector<Frame::GSM>{fixarray::vec3d{(b[0] * cos_tilt) + (b[2] * sin_tilt), b[1],
                                                   (-b[0] * sin_tilt) + (b[2] * cos_tilt)}};
}

/**
 * The OP-77 external field at one GSM point, in `double`, folding the tilt on the way.
 *
 * The convenience form of the overload above for a caller with one point and one epoch; a
 * batch folds once and uses the other.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians; the fold converts it to the report's degrees.
 * @return the external field at @p p, GSM, nanotesla.
 * @complexity O(1) — one fold, one evaluation, `sin`, `cos`.
 * @alloc none.
 * @test IrbemOpq.ReferenceLaneMatchesTheComponentForm
 */
[[nodiscard]] inline FieldVector<Frame::GSM> opq_field_at(Position<Frame::GSM> p,
                                                          double tilt_rad) {
    return opq_field_at(p, opq_parameters<double>(tilt_rad * (180.0 / std::numbers::pi)),
                        std::sin(tilt_rad), std::cos(tilt_rad));
}

/**
 * The model's verdict on a point and a tilt — everything @ref opq_field decides apart from the
 * arithmetic, exposed so a composite field can ask without evaluating twice.
 *
 * Three checks, through `status.hpp` so the rules live in one place: finiteness of the inputs,
 * the tilt against `|psi| <= pi/2` (the model has no `tan(psi)` and is total at a right angle, so
 * the bound is the one @ref FieldContext already guarantees and nothing stricter), and the
 * position against the published `r_GEO <= 15 R_E`. There is NO driver check because the model
 * reads no driver; `check_validity` is still consulted so that the envelope table, not this file,
 * is the authority on that.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians.
 * @return @ref Status::DomainError for a non-finite input, a radius inside the Earth or a tilt
 *         beyond a right angle; @ref Status::OutOfValidityRange beyond 15 R_E; otherwise
 *         @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpq.ValidityIsReportedFromBothSides
 * @test IrbemOpq.NonFiniteInputIsADomainError
 * @test IrbemOpq.ATiltBeyondARightAngleIsADomainError
 */
[[nodiscard]] inline Status opq_status(Position<Frame::GSM> p, double tilt_rad) {
    if (!std::isfinite(p.v[0]) || !std::isfinite(p.v[1]) || !std::isfinite(p.v[2]) ||
        !std::isfinite(tilt_rad)) {
        return Status::DomainError;
    }
    if (!(std::fabs(tilt_rad) <= max_tilt_rad)) return Status::DomainError;
    const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
    const Status where = check_position(ExternalModel::OlsonPfitzerQuiet1977, r, p.v[0]);
    const Status drives = check_validity(ExternalModel::OlsonPfitzerQuiet1977, DriverSet{});
    return first_failure(where, drives);
}

/**
 * The OP-77 external field, with the model's own verdict on whether it should be believed here.
 *
 * The value is **always** returned, `status.hpp`'s standing rule — and beyond 15 R_E that value is
 * a ZERO external field, because that is what the published model says there: the report's
 * template "sets the field to zero" where "the power series diverges". A caller who extrapolates
 * this model gets the published extrapolation, which is no external field at all, and is told so.
 * Below 2 R_E the field is likewise the published zero, and the status is @ref Status::Ok because
 * that region is inside the envelope: the report chose zero there on purpose.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians; positive when the north dipole leans sunward.
 * @return the field and its caveat. @ref Status::DomainError (with a zero field) for a non-finite
 *         input, a radius inside the Earth or `|psi| > pi/2`; @ref Status::OutOfValidityRange
 *         beyond 15 R_E, with the published zero; otherwise @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpq.ValidityIsReportedFromBothSides
 * @test IrbemOpq.NonFiniteInputIsADomainError
 * @test IrbemOpq.ATiltBeyondARightAngleIsADomainError
 * @test IrbemOpq.NothingOnTheHeapInTheHotPath
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> opq_field(Position<Frame::GSM> p,
                                                               double tilt_rad) {
    const Status s = opq_status(p, tilt_rad);
    if (s == Status::DomainError) return {s, FieldVector<Frame::GSM>{}};
    return {s, opq_field_at(p, tilt_rad)};
}

/**
 * The OP-77 external field for a whole epoch's worth of state — the production entry point.
 *
 * Reads the tilt out of @ref HotState and nothing else: no driver slot is touched, which is the
 * literal form of "a quiet-time model". Two contexts differing only in their drivers give
 * bit-identical fields through this overload — @ref IrbemOpq.TheModelReadsNoDrivers.
 *
 * @param p the position, GSM, in Earth radii.
 * @param ctx the epoch's context; only `hot().tilt_rad` is read.
 * @return the field and its caveat, exactly as the two-argument overload.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpq.ContextOverloadAgreesWithTheExplicitOne
 * @test IrbemOpq.TheModelReadsNoDrivers
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> opq_field(Position<Frame::GSM> p,
                                                               const FieldContext& ctx) {
    return opq_field(p, ctx.hot().tilt_rad);
}

// -------------------------------------------------------------------------------------------
// The batch lanes
// -------------------------------------------------------------------------------------------

/**
 * The OP-77 field over a whole batch, on the CPU, in `float`.
 *
 * The host twin of `irbem_opq_f32`: the same expressions, in the same order, in the same
 * precision, from coefficients folded in `double` and rounded to `float` FIRST. That is what makes
 * a disagreement between the two lanes attributable to the device — a contraction, a driver's
 * `exp` — rather than to the arithmetic having been written differently on the two sides.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, GSM, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, nanotesla, GSM; overwritten in full.
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`.
 * @param par the coefficients folded for the epoch's tilt, in `float`.
 * @return `false` when @p pos is not a whole number of points or @p out is a different length, in
 *         which case nothing is written; `true` otherwise.
 * @complexity O(N).
 * @alloc none: the loop is over caller-provided spans and the caller's parameter set.
 * @test IrbemOpq.HostFloatLaneTracksTheReferenceLane
 * @test IrbemOpq.HostFloatLaneRejectsMismatchedSpans
 * @test IrbemOpq.DeviceKernelAgreesWithTheHostLane
 */
[[nodiscard]] inline bool opq_field_host(std::span<const float> pos, std::span<float> out,
                                         float sin_tilt, float cos_tilt,
                                         const OpqParameters<float>& par) {
    if (pos.size() % 3 != 0 || out.size() != pos.size()) return false;
    const std::size_t n = pos.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const float x = pos[(3 * i) + 0];
        const float y = pos[(3 * i) + 1];
        const float z = pos[(3 * i) + 2];
        const float xs = (x * cos_tilt) - (z * sin_tilt);
        const float zs = (x * sin_tilt) + (z * cos_tilt);
        const std::array<float, 3> b = opq_components<float>(par, xs, y, zs);
        out[(3 * i) + 0] = (b[0] * cos_tilt) + (b[2] * sin_tilt);
        out[(3 * i) + 1] = b[1];
        out[(3 * i) + 2] = (-b[0] * sin_tilt) + (b[2] * cos_tilt);
    }
    return true;
}

/// How many `float` scalars the device kernel's parameter buffer holds: `sin(psi)`, `cos(psi)`,
/// then `A(32) B(32) C(22) D(22) E(32) F(32)`. Asserted against the kernel registry.
inline constexpr std::size_t opq_param_count = 2 + (4 * opq_xz_count) + (2 * opq_y_count);

/**
 * Pack the epoch's tilt and its folded coefficients into the kernel's parameter buffer.
 *
 * The layout is the kernel's ABI and is stated in exactly two places — here and the comment above
 * `irbem_opq_f32` in `irbem.slang`. A test evaluates both lanes on the same points, which is what
 * actually keeps the two statements in step.
 *
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`.
 * @param tilt_deg the tilt in degrees, for the fold. Carried separately from its sine and cosine
 *        because the fold is a polynomial in the ANGLE and the rotation is not.
 * @return the parameter block, `opq_param_count` floats, by value.
 * @complexity O(1).
 * @alloc none — the block is the returned object's own inline array.
 * @test IrbemOpq.ParameterBlockCarriesTheTiltThenTheCoefficients
 */
[[nodiscard]] inline std::array<float, opq_param_count> opq_param_block(float sin_tilt,
                                                                        float cos_tilt,
                                                                        double tilt_deg) {
    const OpqParameters<float> p = opq_parameters<float>(tilt_deg);
    std::array<float, opq_param_count> block{};
    block[0] = sin_tilt;
    block[1] = cos_tilt;
    std::size_t at = 2;
    for (std::size_t k = 0; k < opq_xz_count; ++k) block[at++] = p.a[k];
    for (std::size_t k = 0; k < opq_xz_count; ++k) block[at++] = p.b[k];
    for (std::size_t k = 0; k < opq_y_count; ++k) block[at++] = p.c[k];
    for (std::size_t k = 0; k < opq_y_count; ++k) block[at++] = p.d[k];
    for (std::size_t k = 0; k < opq_xz_count; ++k) block[at++] = p.e[k];
    for (std::size_t k = 0; k < opq_xz_count; ++k) block[at++] = p.f[k];
    return block;
}

/**
 * The batch's position caveat, accumulated one point at a time — the same fold `ext_t89.hpp`
 * makes, for the same reason: one status for N points can only be the worst of them, and finding
 * it must not cost a second pass over 100 MB of positions. OP-77's envelope is a radius band, so
 * the smallest and largest `r^2` decide the batch and two `check_position` calls per CALL
 * reproduce exactly what N of them would have said.
 *
 * @test IrbemOpq.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
struct OpqPositionFold {
    /// The smallest `r^2` seen, R_E^2; `+inf` until the first point.
    double r2_lo = std::numeric_limits<double>::infinity();
    /// The largest `r^2` seen, R_E^2; zero until the first point.
    double r2_hi = 0.0;
    /// False once any point has had a non-finite coordinate. Tracked separately because a NaN
    /// radius compares false against everything and would otherwise slip through both extremes.
    bool finite = true;

    /**
     * Fold one position in.
     * @param p the position, GSM, in Earth radii.
     * @complexity O(1) — one fused radius, two comparisons, no `sqrt`.
     * @alloc none.
     * @test IrbemOpq.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    constexpr void add(const Position<Frame::GSM>& p) {
        const double r2 = (p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]);
        finite = finite && std::isfinite(r2);
        r2_lo = r2 < r2_lo ? r2 : r2_lo;
        r2_hi = r2 > r2_hi ? r2 : r2_hi;
    }

    /**
     * What the batch's positions say about the model's envelope.
     * @return @ref Status::DomainError when any point is not finite or is inside the Earth,
     *         @ref Status::OutOfValidityRange when any point is beyond 15 R_E, otherwise
     *         @ref Status::Ok.
     * @complexity O(1) — two square roots and two envelope lookups for the whole batch.
     * @alloc none.
     * @test IrbemOpq.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    [[nodiscard]] Status verdict() const {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double lo = finite ? std::sqrt(r2_lo) : nan;
        const double hi = finite ? std::sqrt(r2_hi) : nan;
        // The model publishes no sunward bound, so the x argument cannot change the answer; the
        // radius is what is being asked about.
        return first_failure(check_position(ExternalModel::OlsonPfitzerQuiet1977, lo, 0.0),
                             check_position(ExternalModel::OlsonPfitzerQuiet1977, hi, 0.0));
    }
};

/**
 * The OP-77 field over a whole batch of GSM points, on the device when that is worth it.
 *
 * **This is the routine to call for more than a handful of points.** @ref opq_field is the
 * reference lane: it is what the batch is verified against, and what runs when there is no device
 * or the batch is too small to pay for one.
 *
 * The arithmetic is T89's regime: 86 fused polynomial terms and one `exp` for 24 bytes in and 12
 * out, ~15 flops/byte, no data-dependent branch, and a 174-float parameter block read identically
 * by every lane in a workgroup — see the `irbem_opq_f32` row of `gpu/dispatch.hpp` for the
 * measurement and the crossover derived from it.
 *
 * **The batch reports the same caveats the scalar entry point does, folded over the batch.** If
 * any point is beyond 15 R_E the batch says @ref Status::OutOfValidityRange and is still computed
 * in full (those points come back as the published zero); if any point is inside the Earth or not
 * finite it says @ref Status::DomainError and every output is zeroed.
 *
 * @param points the positions, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians.
 * @param out receives one field vector per input, GSM, nanotesla; same length as @p points.
 * @return @ref Status::DomainError on a length mismatch, a non-finite or right-angle-exceeding
 *         tilt, or a point that is not finite or is inside the Earth, and then every output is
 *         zeroed; @ref Status::OutOfValidityRange when any point's radius is beyond 15 R_E, with
 *         every point still computed; otherwise @ref Status::Ok. The value is `true` exactly when
 *         the device lane serviced the call — asserted by a test rather than trusted, because a
 *         silent fallback is what makes a performance claim worthless.
 * @complexity O(N); on the device those N run concurrently over `ceil(N/256)` workgroups.
 * @alloc the device lane stages positions and results into two `std::vector<float>` of `3N`; the
 *        host lane allocates nothing.
 * @test IrbemOpq.BatchAgreesWithTheReferenceLane
 * @test IrbemOpq.BatchRejectsMismatchedSpans
 * @test IrbemOpq.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 * @test IrbemOpq.BatchFallsBackToTheHostWhenTheShaderWasNeverBuilt
 * @test IrbemOpq.BatchUsesTheDeviceWhenOneIsAvailable
 * @test IrbemOpq.TheDeviceLaneRefusesABadPointBeforeItDispatches
 */
[[nodiscard]] inline Result<bool> opq_field_batch(std::span<const Position<Frame::GSM>> points,
                                                  double tilt_rad,
                                                  std::span<FieldVector<Frame::GSM>> out) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(tilt_rad)) return {Status::DomainError, false};
    if (!(std::fabs(tilt_rad) <= max_tilt_rad)) return {Status::DomainError, false};
    if (n == 0) return {Status::Ok, false};

    const double tilt_deg = tilt_rad * (180.0 / std::numbers::pi);
    const double sin_tilt = std::sin(tilt_rad);
    const double cos_tilt = std::cos(tilt_rad);
    OpqPositionFold fold;

#if CHEATAH_SPACE_IRBEM_OPQ_GPU
    if (gpu::prefer_gpu("irbem_opq_f32", n) &&
        std::filesystem::exists(gpu::shader_path("irbem_opq_f32"))) {
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
        const std::array<float, opq_param_count> block = opq_param_block(
            static_cast<float>(sin_tilt), static_cast<float>(cos_tilt), tilt_deg);
        gpu::dispatch_batch("irbem_opq_f32", pos, raw, std::span<const float>(block));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = FieldVector<Frame::GSM>{
                fixarray::vec3d{raw[(3 * i) + 0], raw[(3 * i) + 1], raw[(3 * i) + 2]}};
        }
        return {where, true};
    }
#endif

    const OpqParameters<double> par = opq_parameters<double>(tilt_deg);
    for (std::size_t i = 0; i < n; ++i) {
        fold.add(points[i]);
        out[i] = opq_field_at(points[i], par, sin_tilt, cos_tilt);
    }
    const Status where = fold.verdict();
    if (where == Status::DomainError) {
        for (std::size_t i = 0; i < n; ++i) out[i] = FieldVector<Frame::GSM>{};
        return {Status::DomainError, false};
    }
    return {where, false};
}

// -------------------------------------------------------------------------------------------
// The total field: IGRF plus OP-77, as one GeoFieldModel
// -------------------------------------------------------------------------------------------

/**
 * IGRF plus Olson & Pfitzer (1977), as a single field — the `TotalFieldT89` shape without the
 * activity level, because this model has none.
 *
 * Satisfies @ref GeoFieldModel, so it drops into `trace_invariant`, `make_lstar` and everything
 * above them. The tilt fold is done ONCE, at construction, and the sine and cosine with it: a
 * trace evaluates the field thousands of times per line, and refolding 172 coefficient pairs at
 * each of them would cost more than the series itself.
 *
 * Beyond 15 R_E the external part is the published zero, so a trace that leaves the model's
 * region simply continues in the internal field — which is the published model's own behaviour
 * and differs from IRBEM's wrapper, which refuses the whole field there. Inside the belts, where
 * every drift shell this library traces lives, the two agree to the internal field's standard.
 *
 * @tparam NMAX the internal field's truncation degree. 10 reproduces IRBEM's own choice, which is
 *         what the differential tests run through; 13 is IGRF-14's full published degree.
 * @test IrbemOpq.TotalFieldSuperposesInternalAndExternal
 * @test IrbemOpq.LstarRunsThroughTheTotalField
 */
template <int NMAX = 10>
    requires(NMAX >= 1)
class TotalFieldOpq {
  public:
    /// The internal part's truncation degree — what generic staging and buffer sizing read.
    static constexpr int degree = NMAX;

    /**
     * @param internal the internal field, already built for the epoch.
     * @param rotations the epoch's frame rotations — built once, reused for every point; their
     *        `dipole_tilt_deg` is the tilt the fold and the rotations use.
     * @complexity O(1) — one fold, `sin`, `cos`.
     * @alloc none.
     * @test IrbemOpq.TotalFieldSuperposesInternalAndExternal
     */
    TotalFieldOpq(const Igrf<NMAX>& internal, const Rotations& rotations)
        : internal_(&internal),
          rotations_(&rotations),
          tilt_rad_(rotations.dipole_tilt_deg * (std::numbers::pi / 180.0)),
          sin_tilt_(std::sin(tilt_rad_)),
          cos_tilt_(std::cos(tilt_rad_)),
          par_(opq_parameters<double>(rotations.dipole_tilt_deg)) {}

    /**
     * The total field at a geographic point.
     *
     * @param p the position, GEO, Earth radii.
     * @return `B_internal + B_external`, in GEO, nT. When the external model refuses the point —
     *         a non-finite input or a radius inside the Earth — the INTERNAL field is returned
     *         alone rather than a zero or a NaN; beyond 15 R_E the external part is the published
     *         zero and the sum is again the internal field. @ref external_status says which.
     * @complexity One IGRF evaluation, one OP-77 evaluation, two 3x3 rotations.
     * @alloc none.
     * @test IrbemOpq.TotalFieldSuperposesInternalAndExternal
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::GEO>& p) const {
        const FieldVector<Frame::GEO> b_int = internal_->evaluate(p);
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        if (opq_status(p_gsm, tilt_rad_) == Status::DomainError) return b_int;
        const FieldVector<Frame::GSM> b_ext = opq_field_at(p_gsm, par_, sin_tilt_, cos_tilt_);
        // A FieldVector transform is a pure rotation — no origin shift, unlike a Position.
        const FieldVector<Frame::GEO> b_ext_geo = transform<Frame::GEO>(b_ext, *rotations_);
        return FieldVector<Frame::GEO>{b_int.v + b_ext_geo.v};
    }

    /**
     * Whether the external model answered at @p p without caveat, and if not, why.
     * @param p the position, GEO, Earth radii.
     * @return the external model's status; @ref Status::Ok when it contributed in full.
     * @complexity One rotation and one envelope check.
     * @alloc none.
     * @test IrbemOpq.TotalFieldReportsWhenTheExternalModelDeclines
     */
    [[nodiscard]] Status external_status(const Position<Frame::GEO>& p) const {
        return opq_status(transform<Frame::GSM>(p, *rotations_), tilt_rad_);
    }

    /// The epoch's frame rotations — what any caller mapping frames needs.
    /// @return the rotations this field was built with.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOpq.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Rotations& rotations() const { return *rotations_; }

    /// The internal part's Gauss coefficient `g(n, m)`, in nT — the INTERNAL field's, for the same
    /// reason `TotalFieldT89::g` forwards: a superposition with a non-potential part has no
    /// expansion of its own, and every caller asking (the dipole moment, the trace step sizing)
    /// means the internal one.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOpq.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double g(int n, int m) const { return internal_->g(n, m); }

    /// The internal part's `h(n, m)`, in nT — see @ref g for why this forwards.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOpq.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double h(int n, int m) const { return internal_->h(n, m); }

    /// The internal field alone — what `dipole_moment` and any device staging need.
    /// @return the internal model.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOpq.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Igrf<NMAX>& internal() const { return *internal_; }

    /// The folded external coefficients this field evaluates with — exposed so a differential
    /// harness can run the bare series through the same numbers the trace sees.
    /// @return the parameter set folded at construction.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOpq.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const OpqParameters<double>& parameters() const { return par_; }

  private:
    const Igrf<NMAX>* internal_;
    const Rotations* rotations_;
    double tilt_rad_;
    double sin_tilt_;
    double cos_tilt_;
    OpqParameters<double> par_;
};


/**
 * Trace a batch of field lines through IGRF plus OP-77 — the entry point `make_lstar` reaches.
 *
 * The shape mirrors the internal field's @ref trace_invariant_batch, and the return value says
 * which lane served the call for the same reason: a silent fallback is what makes a performance
 * claim worthless. **This overload is the host lane only.** The device tracer that composes an
 * external model (`irbem_trace_total_f32`) is written for T89's parameter block; composing OP-77 on
 * the device would need its own tracer kernel, and until one exists the honest answer is `false`
 * here every time, never a quiet substitution of the internal field. Single traces and drift shells
 * still run — in fp64, on the host, at the reference lane's cost.
 *
 * @tparam NMAX the internal field's truncation degree.
 * @param field the superposed model; carries the epoch's rotations and folded coefficients.
 * @param starts the starting positions, GEO, Earth radii.
 * @param pitch_angles_deg the local pitch angle at each start; same length as @p starts.
 * @param out one @ref FieldLine per input. @param statuses one @ref Status per input.
 * @param opt the tracing options.
 * @return @ref Status::Ok when every line closed, @ref Status::OpenFieldLine when any did not,
 *         @ref Status::DomainError on a length mismatch; the value is always `false` (host lane).
 * @complexity O(lines x steps) total-field evaluations, ~900 flops each.
 * @alloc none.
 * @test IrbemOpq.LstarRunsThroughTheTotalField
 * @test IrbemOpq.TotalFieldBatchTraceIsTheHostLaneAndSaysSo
 * @test IrbemOpq.HeavyDifferentialLstarMatchesTheOracleThroughTheTotalField
 */
template <int NMAX>
[[nodiscard]] inline Result<bool> trace_invariant_batch(
    const TotalFieldOpq<NMAX>& field, std::span<const Position<Frame::GEO>> starts,
    std::span<const double> pitch_angles_deg, std::span<FieldLine> out, std::span<Status> statuses,
    const TraceOptions& opt = {}) {
    const std::size_t n = starts.size();
    if (pitch_angles_deg.size() != n || out.size() != n || statuses.size() != n) {
        return {Status::DomainError, false};
    }
    bool all_ok = true;
    for (std::size_t i = 0; i < n; ++i) {
        const Result<FieldLine> r = trace_invariant(field, starts[i], pitch_angles_deg[i], opt);
        out[i] = r.value;
        statuses[i] = r.status;
        all_ok = all_ok && (r.status == Status::Ok);
    }
    return {all_ok ? Status::Ok : Status::OpenFieldLine, false};
}

}  // namespace cheatah::space::irbem
