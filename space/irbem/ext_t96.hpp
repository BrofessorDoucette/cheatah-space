#pragma once

/**
 * @file ext_t96.hpp
 * @brief space.irbem — Tsyganenko (1995/1996), "T96": the first external model driven by the
 *        solar wind and Dst CONTINUOUSLY, and the first whose deployed coefficients were never
 *        published. This header states what that means for a clean room before it computes anything.
 *
 * T89 (@ref ext_t89.hpp) is a step function of Kp: seven coefficient sets, one per bin. T96 is
 * what replaced it — a magnetosphere confined inside a pressure-scaled magnetopause, whose ring
 * current is driven by Dst, whose tail and Birkeland systems are driven by the interplanetary
 * field, and whose every driver enters as a continuous number. It is IRBEM's `kext = 7`, and it is
 * the model most radiation-belt studies of the 1996-2005 literature were done with.
 *
 * ## The provenance question, asked FIRST and answered with measurements
 *
 * The T89 lesson, applied before a line of physics was written: IRBEM's `kext = 4` turned out to
 * be a revision never published as equations, and no refit of the published form could reach it.
 * So the first thing this module did was isolate the oracle's external field —
 * `get_field(kext = 7) - get_field(kext = 0)`, which cancels IGRF exactly — and ask what the two
 * papers actually publish.
 *
 * **What is published.** Tsyganenko, *Modeling the Earth's magnetospheric magnetic field confined
 * within a realistic magnetopause*, JGR **100**(A4):5599 (1995) publishes the mathematical
 * STRUCTURE: an ellipsoid-plus-cylinder magnetopause whose size scales with the solar-wind dynamic
 * pressure; a shielded dipole (Chapman-Ferraro) field represented by "box" harmonics — Cartesian
 * products `exp(x sqrt(p^-2 + r^-2)) cos(y/p) sin(z/r)`, the class the code calls SHLCAR3X3; a
 * ring current and a tail current sheet, each a warped, thickened disc of the family introduced in
 * Tsyganenko (1989), each with its own shielding field; and an interconnection field proportional
 * to the transverse IMF. Tsyganenko & Stern, *Modeling the global magnetic field of the large-scale
 * Birkeland current systems*, JGR **101**:27187 (1996) adds the region-1 and region-2 field-aligned
 * current modules as deformed conical current sheets whose fields are the "conical harmonics" —
 * tangential fields with no radial component, `B_theta ∝ cos(phi)`, `B_phi ∝ sin(phi)`.
 *
 * **What is not.** Neither paper tabulates the fitted coefficients. The 1995 paper's shielding
 * expansions alone carry ~100 numbers per module; the 1996 paper's deformation parameters are
 * likewise stated to exist "in the code". The complete parameter set of T96 was distributed ONLY
 * as Tsyganenko's Fortran (GPL-3.0) and as its LGPL-3.0 re-distribution inside IRBEM, and this
 * MIT clean-room implementation may read neither. **So the verdict for this model is
 * PUBLISHED-FORM-WITH-DOCUMENTED-GAP, decided before implementation, not discovered after.**
 *
 * ## What the oracle's driver structure is — MEASURED, then built in
 *
 * A model whose coefficients cannot be read can still be characterised as a black box, and the
 * driver structure of `kext = 7` turned out to be exact enough to state as theorems. Measured by
 * [`tools/oracle/t96_diff.cpp`](../../../tools/oracle/t96_diff.cpp) (its `probe` pass) against the
 * `-O2` oracle, relative deviations of the stated form from the oracle's own values:
 *
 * - **The field is EXACTLY affine in Dst** at fixed pressure, IMF and tilt: predicting Dst = -63
 *   from Dst = 0 and -100 deviates by 1e-16. There is no Dst x IMF cross term (measured zero to
 *   roundoff) and no Dst x By cross term either.
 * - **The IMF response is EXACTLY homogeneous of degree one in (By, Bz)**: scaling the IMF vector
 *   by 0.25, 0.5, 0.75 scales the response linearly to 1e-16. So it is `B_t · F(theta)` with
 *   `B_t = sqrt(By^2 + Bz^2)` and `theta` the clock angle — which is why the oracle is piecewise
 *   linear in Bz with a kink at Bz = 0 (northward and southward slopes differ by a factor of ~4
 *   in the tail), and even in By to leading order.
 * - **The clock-angle dependence is spanned by `{sin theta, cos theta, sin(theta/2)}`**: on the
 *   noon-midnight meridian exactly (the constant term fits to 1e-5), off it to 1-3% RMS. The
 *   `sin(theta/2)` term is the Birkeland driver: a quantity that is `|Bz|` for pure southward IMF,
 *   zero for pure northward, and `|By|/sqrt(2)` for pure By — see @ref t96_clock_driver.
 * - **The pressure dependence is continuous and, per coefficient, affine in `sqrt(Pdyn)` to
 *   1-3%** over 0.5-10 nPa; the residual is the geometric rescaling of the magnetopause, which is
 *   part of the documented floor below.
 * - **The oracle refuses `r > 40 R_E`** (returns `baddata`), which is the published envelope and
 *   the one `status.hpp` carries.
 *
 * Those five facts fix the form of the driver dependence exactly, so it is built in rather than
 * fitted: with `h = B_t sin(theta/2)` (eq. (D) below),
 *
 *     B(x; psi, Dst, Pdyn, By, Bz) = B_0(x; psi, Pdyn) + Dst B_D(x; psi, Pdyn)
 *                                  + Bz B_Z(x; psi, Pdyn) + h B_A(x; psi, Pdyn)
 *                                  + By B_Y(x; psi, Pdyn)                             (1)
 *
 * and the whole question of "what is T96" reduces to the five SPATIAL fields on the right.
 *
 * ## The spatial form — published modules, one basis, one fit
 *
 * Each of the five fields in (1) is expanded on ONE shared basis of divergence-free fields, all of
 * them the published classes above (`@ref t96_basis` evaluates it, and is the only physics in this
 * file):
 *
 * 1. **Warped-disc current modules** (Tsyganenko 1989 eqs. (7)-(9), (11)-(17); the formalism T96's
 *    tail and ring current are built from). The `A^(1)` and `A^(2)` tail potentials at four scale
 *    lengths and two truncation positions, and the `A^(3)` ring potential at four radii, all on ONE
 *    warped, thickened sheet: `x`-hinged and `y`-bent (eq. 11), thickened sunward and at the flanks
 *    (eq. 13). Evaluated in SM and rotated to GSM, as the paper says. The fixed parameters of
 *    T89 §3 are reused unchanged (@ref t89_fixed). 20 fields.
 * 2. **Toroidal field-aligned-current modules** — `B = ∇T × r` with
 *    `T = t(r) s(theta, r) A_m(theta, phi)`: the Mie/toroidal representation (Backus 1986), of
 *    which Tsyganenko & Stern's conical harmonics are the `t = 1/r`, thin-sheet, undeformed special
 *    case. Here the sheet `s` is a smooth step across a DIPOLE shell `sin^2 theta_0 = r/L` (the
 *    deformation the 1996 paper applies to make the cone follow field lines, done in the potential
 *    where it costs nothing), at six shells and two widths, with radial profiles `r^-1, r^-2,
 *    r^-3`, and azimuthal modes `m = 1` (the dawn-dusk antisymmetric region-1/2 sense) and `m = 2`
 *    (the noon-midnight asymmetry of those systems, measured at half the amplitude of `m = 1` at
 *    2 R_E). `B_r = 0` and `div B = 0` hold identically for ANY `T`. 72 fields.
 * 3. **Box harmonics** — `B = -∇U`, `U = exp(kappa x) {cos, sin}(y/p) {sin, cos}(z/r)` with
 *    `kappa = sqrt(p^-2 + r^-2)` so that `∇²U = 0` exactly: the SHLCAR class of the 1995 paper,
 *    all four `(y, z)` parities because the By-driven field has no parity of its own, three scales
 *    each way. 36 fields.
 *
 * Every basis field is exactly divergence-free — the discs because they are curls of potentials
 * in which every geometric substitution is made (T89's argument), the toroidal terms by
 * construction, the box terms because they are gradients of harmonic functions — so ANY linear
 * combination is, and the finite-difference divergence test has no floor to hide behind
 * (@ref IrbemT96.DivergenceVanishesEverywhere: it falls as `h^2` over three decades).
 *
 * The 128 amplitudes of each of the five fields are made to depend on tilt and pressure through
 * six factors, `{1, sin psi, sin^2 psi} x {1, sqrt(Pdyn)}`, giving 768 numbers per field and 3840
 * in all — @ref t96_tables. **Those numbers are FITTED TO THE ORACLE'S OBSERVABLE OUTPUT** by
 * linear least squares over 5 dipole tilts, 3 pressures and ~360 points inside the magnetopause,
 * with the handful of non-linear geometric parameters (@ref t96_geometry) refined by Nelder-Mead
 * on the same cost. They are derived from a black box's behaviour and from nothing else: no
 * source of Tsyganenko's or IRBEM's was read, and the provenance of every number in this file is
 * the harness that regenerates it. That is the same standing the T89 differential harness has,
 * applied to the coefficients instead of to the residual.
 *
 * ## The documented gap — MEASURED, NOT ASSUMED
 *
 * Measured by the harness's `report` pass over points inside the magnetopause, 1.2-35 R_E, at the
 * five fitting tilts and three fitting pressures, RMS of the vector difference against the RMS of
 * the oracle's own field for that family (this is the free-refit FLOOR: every linear amplitude of
 * the published form is already at its least-squares optimum, so it is what the FORM cannot
 * represent, exactly as T89's 0.44 nT floor was):
 *
 * MEASURED_FAMILY_TABLE
 *
 * Carried through to whole-field differentials over the four corpus regimes (drivers of
 * `tests/irbem_domain_corpus.hpp`, clamped to the published envelope for the two that exceed it):
 *
 * MEASURED_REGIME_TABLE
 *
 * Where the floor sits is itself informative: the interconnection families (`B_Y`, `B_Z`) are
 * potential fields and the box harmonics reproduce them to a few percent; the ring-current and
 * quiet-time families are limited by the disc shapes (T96's ring current is not an `A^(3)` disc
 * and its tail is not an `A^(1)+A^(2)` sheet), and the Birkeland family by the conical-sheet
 * deformation the 1996 paper does with its own stretch functions. A reader who wants the last
 * factor of two must read Tsyganenko's code, which this file may not.
 *
 * @note Nothing on a hot path allocates, and nothing but the device lane can throw. The
 *       amplitude collapse (@ref t96_amplitudes) is a per-epoch, per-driver-set computation of
 *       3840 multiply-adds that @ref t96_field pays per call and a batch pays once; the evaluator
 *       is straight-line arithmetic over one stack basis array. The device lane stages `3N`
 *       floats each way and forwards whatever `gpu::dispatch_batch` throws; its `@alloc` says so.
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
#include "ext_t89.hpp"
#include "frames.hpp"
#include "igrf.hpp"
#include "status.hpp"

// The device lane is opt-in by include path, exactly as ext_t89.hpp's is.
#if __has_include("cheatah_gpu_linalg/context.hpp")
#  include "gpu/dispatch.hpp"
/// 1 when this translation unit can reach the T96 device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_T96_GPU 1
#else
/// 1 when this translation unit can reach the T96 device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_T96_GPU 0
#endif

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// The basis: what it is made of, and how it is indexed
// -------------------------------------------------------------------------------------------

/// How many tail-sheet truncation positions `x_0` the disc modules are evaluated at.
inline constexpr std::size_t t96_tail_truncation_count = 2;
/// How many tail scale lengths `a_T`.
inline constexpr std::size_t t96_tail_scale_count = 4;
/// How many ring-current radii `a_RC`.
inline constexpr std::size_t t96_ring_scale_count = 4;
/// How many dipole shells `L` the field-aligned-current sheets sit on.
inline constexpr std::size_t t96_fac_shell_count = 6;
/// How many angular widths those sheets are given.
inline constexpr std::size_t t96_fac_width_count = 2;
/// How many radial profiles: `r^-1`, `r^-2`, `r^-3`.
inline constexpr std::size_t t96_fac_profile_count = 3;
/// How many azimuthal modes: `m = 1` and `m = 2`.
inline constexpr std::size_t t96_fac_mode_count = 2;
/// How many box-harmonic `(y, z)` parity families: `cos·sin`, `cos·cos`, `sin·sin`, `sin·cos`.
inline constexpr std::size_t t96_box_parity_count = 4;
/// How many box scales in `y` and, separately, in `z`.
inline constexpr std::size_t t96_box_scale_count = 3;

/// The disc modules: two tail potentials per (truncation, scale) plus the ring potentials.
inline constexpr std::size_t t96_disc_count =
    (2 * t96_tail_truncation_count * t96_tail_scale_count) + t96_ring_scale_count;
/// The toroidal field-aligned-current modules.
inline constexpr std::size_t t96_fac_count =
    t96_fac_shell_count * t96_fac_width_count * t96_fac_profile_count * t96_fac_mode_count;
/// The box harmonics.
inline constexpr std::size_t t96_box_count =
    t96_box_parity_count * t96_box_scale_count * t96_box_scale_count;
/// The whole basis: 20 + 72 + 36.
inline constexpr std::size_t t96_basis_count = t96_disc_count + t96_fac_count + t96_box_count;
/// Where each block starts in the basis index.
inline constexpr std::size_t t96_fac_first = t96_disc_count;
/// Where the box harmonics start in the basis index.
inline constexpr std::size_t t96_box_first = t96_disc_count + t96_fac_count;

/// The tilt factors every amplitude carries: `{1, sin psi, sin^2 psi}`.
inline constexpr std::size_t t96_tilt_factor_count = 3;
/// The pressure factors: `{1, sqrt(Pdyn)}`.
inline constexpr std::size_t t96_pressure_factor_count = 2;
/// Tilt x pressure factors per basis function.
inline constexpr std::size_t t96_factor_count = t96_tilt_factor_count * t96_pressure_factor_count;
/// Coefficients per driver family: `t96_basis_count x t96_factor_count`.
inline constexpr std::size_t t96_column_count = t96_basis_count * t96_factor_count;
/// The driver families of eq. (1), in the order the tables and the collapse use them.
inline constexpr std::size_t t96_family_count = 5;

/**
 * Which term of eq. (1) a coefficient table belongs to.
 *
 * The order is the order of @ref t96_tables and of the weights @ref t96_amplitudes applies:
 * `{1, Dst, Bz, h, By}`.
 */
enum class T96Family : std::uint8_t {
    Quiet = 0,     ///< `B_0`: the field at Dst = 0 and zero IMF.
    Dst,           ///< `B_D`: per nT of Dst.
    Bz,            ///< `B_Z`: per nT of IMF Bz (the northward-branch slope).
    Clock,         ///< `B_A`: per nT of `h = B_t sin(theta/2)`, the Birkeland driver.
    By,            ///< `B_Y`: per nT of IMF By.
};

/**
 * The non-linear geometry of the basis — the numbers the linear fit cannot move.
 *
 * Everything the disc modules, the field-aligned-current sheets and the box harmonics need beyond
 * their amplitudes. The sheet geometry (`r_c .. d_y`) is T89's (eqs. 11 and 13) and takes T89's
 * meaning; the scale-length ladders are the disc, shell and box scales the basis is evaluated at.
 * Templated on nothing: the fp32 lane reads these through `static_cast<T>` exactly as it reads
 * @ref t89_fixed.
 *
 * @test IrbemT96.GeometryIsTheFrozenFitGeometry
 */
struct T96Geometry {
    /// `R_c`, R_E — the hinging distance of the warped sheet, T89 eq. (11).
    double r_c;
    /// `G`, R_E — the amplitude of the sheet's dawn-dusk bending, T89 eq. (11).
    double g;
    /// `D_0`, R_E — the sheet half-thickness in the central tail, T89 eq. (13).
    double d_0;
    /// `gamma_RC`, R_E — the ring sheet's dayside thickness increment, T89 eq. (13).
    double gamma_rc;
    /// `D_y`, R_E — the dawn-dusk scale of the tail truncation factor `W`, T89 eq. (13).
    double d_y;
    /// The tail truncation positions `x_0`, R_E, T89 eq. (13) — one `W` per entry.
    std::array<double, t96_tail_truncation_count> x_0;
    /// The tail scale lengths `a_T`, R_E — one `A^(1)` and one `A^(2)` per entry per truncation.
    std::array<double, t96_tail_scale_count> a_t;
    /// The ring radii `a_RC`, R_E — one `A^(3)` per entry.
    std::array<double, t96_ring_scale_count> a_rc;
    /// The dipole shells `L`, R_E, the field-aligned-current sheets follow: `sin^2 theta_0 = r/L`.
    std::array<double, t96_fac_shell_count> fac_shell;
    /// The angular half-widths of those sheets, radians.
    std::array<double, t96_fac_width_count> fac_width;
    /// The box `y` scales `p_i`, R_E.
    std::array<double, t96_box_scale_count> box_p;
    /// The box `z` scales `r_k`, R_E.
    std::array<double, t96_box_scale_count> box_r;
};

// ---- BEGIN GENERATED: written by tools/oracle/t96_diff.cpp `fit`; do not edit by hand ----------

/**
 * The frozen geometry — the starting ladders refined by the harness's Nelder-Mead pass.
 *
 * @test IrbemT96.GeometryIsTheFrozenFitGeometry
 */
inline constexpr T96Geometry t96_geometry{
    /* r_c */ 10.8228765, /* g */ 5.40080041, /* d_0 */ 1.53053164, /* gamma_rc */ -0.105105896, /* d_y */ 27.3054127,
    /* x_0 */ {5.64362667, -4.70755261},
    /* a_t */ {2.11466257, 6.72687226, 19.8766187, 39.7464194},
    /* a_rc */ {2.1440289, 3.7687063, 6.05167693, 17.0039816},
    /* fac_shell */ {3.96624426, 6.55293792, 9.15563147, 11.2477354, 21.156181, 46.3161443},
    /* fac_width */ {0.0923025449, 0.485290252},
    /* box_p */ {3.14070776, 12.0837768, 19.3949163},
    /* box_r */ {4.54044129, 6.32833314, 20.7723996}};

/**
 * The fitted amplitude tables: `t96_tables[family][basis * 6 + tilt_factor * 2 + pressure_factor]`.
 *
 * Fitted to the IRBEM oracle's black-box output by the harness; see the file brief for what that
 * means and what it does not. Families in @ref T96Family order.
 *
 * @test IrbemT96.TablesAreNotEmpty
 */
inline constexpr std::array<std::array<double, t96_column_count>, t96_family_count> t96_tables{{
    // B_0 (quiet)
    {{
     3358.94997, -3676.54386, 137.286899, -656.043263, 887.585874, -786.009852, 
     -1606.276, 3134.80957, -105.39714, 1060.93066, -801.317995, 1877.08766, 
     -7405.40058, 8583.96587, -475.236712, 1686.18404, -2255.63649, 1394.49649, 
     -1011.21551, -2749.4887, -737.375079, -1761.02347, -193.059096, -1184.72541, 
     846.073484, -5042.01828, 828.942975, -247.498934, 1522.076, 618.192863, 
     120.941364, 1270.93122, 105.497479, 296.370654, -20.0422864, 541.241911, 
     4874.87273, -1735.94357, -723.720764, -1967.88814, 1774.43134, -3953.91831, 
     -416.793487, -861.882311, 77.4607948, 74.8943079, -38.4475348, -1.00758434, 
     1612.84801, -423.501434, 341.925563, 505.607156, 841.582523, -1393.76821, 
     -8096.53067, 5851.66107, -849.27458, -522.02218, -547.881936, 1672.0263, 
     1585.75536, -3060.56634, -93.7202997, -1053.66289, -1273.63872, 1343.92541, 
     -546.850506, -3628.43083, -136.812119, -33.7525672, -202.572203, -1725.66435, 
     -1694.80956, 2582.72312, -219.034861, -1530.40713, 112.15356, 3689.28277, 
     721.926708, 1009.85367, -24.4507965, -357.171365, 128.732993, 601.295144, 
     -2941.9969, 3059.55469, -80.762601, 3743.34255, -897.114761, -2715.4154, 
     -335.97708, -638.353664, 38.5928515, -55.287915, -22.5495327, 77.30434, 
     -1762.84366, 1004.55092, -154.235415, -100.994974, -52.4162286, 203.287238, 
     10222.9838, -2754.70757, 442.188217, 997.598027, -305.990924, -2420.93264, 
     -15178.5296, -3501.64024, 885.65529, -3291.14325, -647.777429, 6516.70564, 
     -406.606866, 11267.9399, 1659.93544, 4098.30521, 717.735906, 4114.94289, 
     2.67211504, -2.02799758, 5.5408399, -3.82300611, -29.5634115, -1.30410053, 
     -0.41007159, -1.28169274, 7.29584849, -5.72239677, -8.04320649, 8.58435684, 
     -13.2528425, 14.5178334, -29.4839639, 18.3715509, 128.014767, 6.57043067, 
     6.73906472, 2.36027099, -29.4697881, 28.1177594, 20.5111826, -27.6502534, 
     14.4143891, -19.0257921, 29.2213024, -19.6842376, -116.705448, -11.0451322, 
     -11.4945866, 0.947390003, 26.3236261, -29.963678, 17.5459002, -5.08214574, 
     -3.19244702, -4.72604904, 38.4958067, -12.9450337, -13.3438958, 11.3831917, 
     -45.2186492, 30.7530028, 3.43360886, 7.00249034, 18.972694, -13.5970759, 
     63.7003761, -99.3924749, -228.692752, 33.3493197, 159.836957, -2.62537543, 
     260.670924, -131.230954, -70.08812, -42.0917353, -21.4069149, 18.1153375, 
     -268.04862, 337.100481, 261.440254, 96.2787722, -348.472936, -76.5307542, 
     -414.756974, 93.7969911, 110.059088, 108.214903, -125.612342, 107.303319, 
     5.67588684, -7.03656241, 9.09560583, -3.61714678, 1.75254135, -0.963451268, 
     1.03664322, -0.281940192, 9.6939961, -4.86808242, -21.346708, 11.9075356, 
     -28.5654189, 39.6173304, -56.2389987, 34.2633223, 19.878621, 25.191751, 
     14.7160535, -4.99714376, -37.1568393, 20.3774017, 97.1661302, -62.5118264, 
     36.6782974, -56.9300835, 74.4340362, -58.0162564, -31.8840151, -43.8956861, 
     -32.4476081, 7.76562269, 20.4024764, -11.3580942, -41.1124987, 47.0623723, 
     -63.9506144, 51.8461252, 0.774146757, -90.5160875, -65.2783867, 32.7969922, 
     -3.00703016, 70.3636472, -71.9653953, 35.0539985, -23.8555359, 53.3922276, 
     162.969008, -453.724553, -44.3155223, 852.684161, 276.421651, -337.39628, 
     -124.631509, -743.787138, 379.353803, -18.2816878, -90.4837335, -383.329301, 
     1893.9501, -454.888158, -1026.68814, -849.261994, 990.504554, 1145.25024, 
     1093.12718, 1679.62873, -451.953991, -705.907172, 1507.97135, 430.488995, 
     -7.92162967, 8.44995013, -12.1496, 15.373907, 15.083988, 2.21311516, 
     -3.69163248, 8.06843812, -1.31864214, 7.59316875, 10.0992162, -23.0970782, 
     43.556728, -46.0614562, 91.3996252, -83.0692422, -89.9405063, 20.0388569, 
     20.8621672, -59.559146, 9.26528648, -28.0151404, -92.5757341, 207.235083, 
     -59.7631249, 61.9684029, -146.001628, 96.798094, 143.515302, -88.9734968, 
     45.6742518, 73.8971978, 2.48561736, 0.382755077, 58.2398711, -360.215839, 
     228.82514, -460.515731, -479.877946, 384.34081, 149.748021, -436.562184, 
     -30.5778261, -107.159964, -134.352683, -31.8725049, 384.002072, -257.772906, 
     -1793.9233, 4515.34901, 3972.28404, -3130.85166, 669.276944, 3580.53372, 
     527.385475, 1228.2744, 1698.71634, -103.294986, -1019.43388, 296.824406, 
     -2327.23055, -4781.92347, 1351.72825, -198.172305, 69.876895, -354.708087, 
     -809.192934, -3722.83478, 619.296317, 1006.53204, -272.875414, -2397.02461, 
     5.83907641, -2.27518625, 11.8333331, -0.258213926, -10.5268166, 12.4316992, 
     2.43074731, -4.87530149, -8.60785107, 14.8434828, -48.1789145, 59.3201469, 
     -68.3698627, 60.3521157, -110.554407, 30.3278847, 30.0897017, -101.515838, 
     -5.12103541, 41.7770009, 40.1846035, -68.4277319, 234.586205, -402.708378, 
     54.5710951, -96.0968848, 180.277027, -69.171177, -8.21305374, 139.997069, 
     -72.0435033, -66.0429209, -84.8358992, 113.636072, -12.1400552, 508.968312, 
     -45.1056342, 90.0740881, 412.067262, -278.291958, -64.8711912, 1037.80773, 
     184.678601, -133.215187, 195.585746, 31.1408311, -311.66173, 105.247269, 
     1213.86798, -1780.65713, -4715.46534, 2847.54363, -3221.06348, -7743.36923, 
     -1610.46159, 1012.78699, -3071.50061, 160.615979, -274.847379, 2164.56582, 
     262.016337, 2575.7723, -576.668065, 635.950026, -621.735225, -1913.54551, 
     -54.5355346, -99.3613868, -809.933251, -70.9315881, -51.4318151, -861.713208, 
     5.80536912, -12.5020779, -49.289788, 31.8605445, -129.791957, 74.0075651, 
     -11.6038803, 6.60366688, -54.2414292, 32.9821192, -0.865883354, -19.1706022, 
     -44.2933701, 75.3181319, 259.937346, -140.166982, 555.696524, -335.177202, 
     -6.33482384, -6.88078807, 179.41698, -96.1433216, -170.465287, 204.725, 
     46.9993964, -81.3534869, -314.206675, 129.097588, -580.516051, 282.94327, 
     16.1310444, -5.14072904, -85.0988726, 10.8084708, -211.441816, -138.058531, 
     -247.1838, 487.022525, 352.511872, -181.059333, 89.0132114, -1013.52022, 
     -100.656151, 189.786448, 172.685184, -171.306641, -47.2625407, 0.399501372, 
     1035.48547, -4864.9316, 458.472893, -483.185184, 2768.36438, 8322.52507, 
     1031.40839, -2234.44238, 1357.59038, 228.770614, 1101.21845, -1347.85523, 
     318.134678, 6238.38253, 366.315388, 1289.75794, 511.486909, 1506.98361, 
     866.410352, 3787.14463, 187.407171, -57.7046496, 161.14756, 1855.04442, 
     47.447137, -64.829269, 134.130706, -65.6201196, -206.120623, 144.00331, 
     66.6215753, -59.2678878, 56.8094218, -8.43295869, -298.359534, 206.03435, 
     -206.038684, 278.907471, -623.107445, 288.951338, 631.808003, -464.436801, 
     -203.415164, 171.46564, -331.347278, 118.311111, 752.501335, -387.183817, 
     184.679019, -278.361882, 616.531595, -284.222196, -346.150925, 305.968586, 
     202.473173, -109.607083, 66.6164543, 22.9519877, 271.775894, -394.841238, 
     56.784888, -60.5225744, -422.215961, 203.986292, 276.701288, 144.817435, 
     -57.9155536, -1.38997292, -178.465955, 95.8890386, 355.520119, -139.411884, 
     -327.584155, 2205.28503, 1023.60841, -266.299007, -1977.00894, -2995.17014, 
     85.5709336, 728.403198, -113.611725, -208.029879, -566.866671, -251.570867, 
     -182.485522, -3408.61363, -814.770738, -771.518931, 379.731214, -872.644693, 
     -829.091825, -1642.2144, 390.335361, -379.88992, -1276.69342, 1184.32918, 
     82.9353105, -68.1145016, 0.668238534, -10.2871031, -67.7717029, 46.5712307, 
     -225.542076, 185.476782, 1.0055006, 25.8111636, 203.837212, -123.699466, 
     356.70328, -290.458306, 1.57224705, -39.0407131, -357.746835, 192.916711, 
     -339.836572, 280.192375, 127.272095, 27.7340942, 263.789065, -906.944572, 
     644.893578, -534.954315, -499.51428, 51.4779077, -910.503557, 1911.62769, 
     -214.1207, 221.925372, 470.437863, -16.8008796, 1623.03378, -2252.72897, 
     323.682496, -267.325945, -136.488766, -3.92200731, -227.547445, 909.084823, 
     -583.522506, 492.218823, 515.395453, -103.99807, 764.33607, -1888.56948, 
     -214.526283, -284.658474, -495.597239, 63.2301368, -1147.22806, 2217.25782, 
     7.1252019, -8.4610834, 6.51702087, -13.0939622, -25.0676665, 8.96108652, 
     -14.2341258, 15.1591646, -2.38652876, 22.1330194, 55.589085, -16.9902717, 
     8.22409782, -7.7941457, -7.17818045, -6.63503871, -32.4290607, 9.09095605, 
     69.2644351, -18.5396151, -298.005329, 353.23937, 92.3952304, -87.0573223, 
     -181.339807, 97.6638376, 303.883624, -570.515768, -241.108641, 162.114693, 
     49.6684576, -27.8320206, 333.610046, 89.0795178, 289.935804, -132.966973, 
     -46.9506862, 5.40214651, 319.266047, -370.76756, -86.141094, 88.6765349, 
     141.94361, -74.3113161, -365.205932, 641.55601, 224.453672, -163.941412, 
     -27.5644981, 15.4753838, -605.335503, -366.53913, -355.651685, 148.51231, 
     -5.30333301, 10.4393068, 31.7883345, -35.1801019, 17.857224, -1.45209352, 
     -0.414250986, -20.0522426, -73.719891, 71.281777, -65.7708978, 15.6365559, 
     36.7901429, 7.69192602, 149.370926, -118.541983, 131.987995, -37.4425971, 
     143.639257, 55.988359, -885.055805, 1304.33258, 11.0999111, -375.35518, 
     -790.726031, 472.143151, -292.667268, -493.404871, 179.14291, 398.17674, 
     660.171102, -608.207546, 109.324531, 198.106631, 6.73620363, -287.679185, 
     -14.9436101, -262.45341, 1490.55962, -2088.56689, 43.3217586, 539.060823, 
     770.722922, -339.738159, 272.024923, 847.802233, -382.921618, -543.432127, 
     -497.316846, 475.073736, -185.166471, -241.29317, 13.6061565, 414.328628, 
     -0.991303426, 3.4734524, 0.183941266, -1.08374809, 18.129226, -12.3070762, 
     1.43501037, -6.63902432, -1.26200596, 3.4477297, -38.3583464, 28.2562417, 
     3.15235425, 0.118585137, -2.13476879, 0.571077168, 18.8118958, -15.1642883, 
     -587.292378, 519.524922, 332.341271, -62.0118312, 161.374041, -251.917759, 
     255.966739, -209.017696, -230.213606, -100.607172, -11.9823791, 128.88139, 
     -44.3001335, 13.2109693, 56.3067272, -11.4476233, -12.5815256, 21.7746965, 
     958.612374, -856.709802, -630.418959, 177.490655, -305.01148, 406.538799, 
     -429.541922, 373.811739, 528.201947, 26.2320902, 119.220525, -217.912733, 
     48.4115034, -13.4823612, -121.922918, 50.5594627, -29.1049707, -32.7875685}},
    // B_D (per nT Dst)
    {{
     23.1030096, 53.2412512, -9.41875417, 0.790335975, -5.82011928, 9.24072496, 
     -140.042528, -29.9340437, 3.841324, -9.81233017, -12.2267388, -7.90275795, 
     68.4701038, -132.387852, 23.1332249, 5.31751056, 18.3310898, -17.2227826, 
     48.4950371, -56.4761859, 2.00969226, 8.2653301, 1.22486665, -15.1094212, 
     -87.1065014, 102.496351, -13.3581288, -11.5562769, -0.326665297, 20.1117112, 
     1.46063576, 14.0142753, 0.288280736, -1.19062718, 1.80897961, 4.43727865, 
     -90.7751237, 13.3294874, -2.2356293, 1.20549708, -19.3395394, -13.9615764, 
     -0.34932274, -1.41592874, -0.673565919, -0.242356659, 0.355883316, 0.0121185759, 
     -39.522599, -32.3263532, -5.42715067, -5.89587644, -7.82846899, -2.56956918, 
     176.259182, -10.265491, 28.9234215, 19.2505706, 38.6654607, 0.193201805, 
     -51.3674071, 117.616802, -3.32251034, 0.789395298, -8.47646326, 11.9844594, 
     -24.9200311, -26.738768, -2.48831108, -2.09758102, -6.40002793, -5.19652289, 
     87.2080986, -111.772631, 2.31303456, 10.9939876, 12.432596, -20.7183675, 
     -1.49364067, 11.5102627, -0.457246402, 0.400198442, 0.326502247, 1.26338782, 
     90.129845, -11.1137664, 11.4055459, -1.84171838, 9.63868177, 10.8164112, 
     0.460181207, -9.38721304, -0.328053203, 0.546782178, 0.0452713065, -2.10034996, 
     42.4913536, -17.4754998, 1.85303221, -1.73800392, 9.84545247, -5.5779597, 
     -371.875668, 65.8970602, -6.50862639, 11.9034285, -51.0983722, 33.8390291, 
     828.170941, 11.2804871, -0.11307122, -17.0720338, 109.97543, -73.7810333, 
     -49.6989533, -436.232891, -5.18787468, -26.5425217, -12.9422854, -82.5097223, 
     0.0714433945, -0.0364819726, 0.0152737595, -0.00471038865, 0.161786108, 0.0114162886, 
     0.0277367728, -0.0244117033, 0.0415196041, 0.00869837148, 0.256875934, -0.0203875269, 
     -0.389011121, 0.183364755, -0.0172627531, -0.0122312632, -0.565333415, -0.0645794464, 
     -0.122917067, 0.105983109, -0.286043131, -0.0274263929, -1.33239429, 0.12876329, 
     0.449259195, -0.198167245, 0.0467995882, 0.00940771743, 0.24537367, 0.13918747, 
     0.121210904, -0.100412925, 0.390268841, 0.00464692076, 1.50090491, -0.0955847904, 
     0.862313454, -0.407706112, -0.0884122516, 0.075462282, -0.333823705, 0.113680195, 
     0.00857910589, -0.157274807, 0.203685013, -0.139684333, 0.113177504, 0.0431248587, 
     -4.87697545, 2.41374673, 0.771561232, -0.365121161, 0.705791934, -0.727334752, 
     0.397893228, 1.11770543, -0.75588936, 0.896737281, -1.50988772, 0.0447934299, 
     6.97257175, -3.5264834, -0.175825104, -0.586853756, -1.11584876, 2.23797231, 
     -1.32624948, -2.382877, -0.436677294, -1.158009, 1.50430381, -0.532319163, 
     -0.00728854169, 0.00467232301, -0.0686539139, 0.0132244748, 0.40023942, -0.0890764098, 
     0.0220546795, -0.0107240102, -0.0158336806, 0.027602299, -0.209230443, 0.025045902, 
     -0.138458312, 0.0495980747, 0.229833224, -0.102192981, -1.6593538, 0.250862328, 
     -0.238474277, 0.124356649, -0.110174951, -0.130807634, 1.07304279, -0.198721169, 
     0.228103887, -0.0429715559, -0.0944007444, 0.136898337, 1.15701947, -0.159918708, 
     0.370255291, -0.126210686, 0.323822622, 0.0972368775, -1.80224574, 0.375492462, 
     -0.358213363, -0.174356317, 0.791657588, 0.152882383, -0.767452792, 0.00501712527, 
     1.4733282, -0.0265273725, -0.617189874, 0.220019848, 0.387970063, -0.356965447, 
     2.71560429, 1.74522137, -4.7484288, -3.0965286, 6.00678621, 0.328514904, 
     -13.7726374, -1.37739477, 4.03579603, -2.60093972, -1.34369537, -0.287854563, 
     -10.7658569, 1.04773259, 8.36939335, 6.11152925, -10.8591024, 0.279851979, 
     25.8306947, 11.7343971, 4.66845249, 6.76589622, -2.46793451, 7.76777527, 
     0.0487514519, -0.040430686, -0.0524046407, -0.0463157074, 0.14651431, -0.0663345794, 
     -0.0728531267, -0.00478812661, -0.0320118923, -0.0695297762, 0.544161283, -0.20034375, 
     -0.410692556, 0.279108876, 0.378672022, 0.194757744, -1.13316992, 0.119599191, 
     0.322246404, 0.156444353, 0.373237889, 0.187250116, -3.28979313, 0.642291334, 
     0.614832152, -0.386282038, -0.519333265, -0.17678367, 1.56009469, 0.145446341, 
     -0.554989213, -0.271821043, -0.858426785, 0.0677176904, 5.59271912, -0.419762832, 
     0.0177167086, 1.42320121, 2.0093199, -1.72223832, -3.73679144, 4.35832645, 
     -0.392824696, -0.0225141225, 2.90242587, 0.995081367, -4.27398051, 2.99686817, 
     -0.308129282, -14.1913866, -16.9351876, 15.8908547, 2.07168971, -23.8120225, 
     7.9817041, -0.283256727, -25.4701126, -9.2492045, 15.0212799, -11.6203591, 
     -4.98624722, 15.4241053, -12.3434233, -5.41558973, -5.6194075, -2.88210804, 
     -21.2432699, -12.8650926, -6.98021298, -4.27779843, -7.87640767, -4.27737354, 
     -0.0482306972, -0.0114105807, 0.041206672, -0.00477779475, 0.0603834522, -0.0445532718, 
     0.0144177977, 0.00348495146, 0.0858692963, -0.11595165, 0.191497529, -0.277980749, 
     0.483037507, -0.0523937664, -0.222194656, 0.013138319, 0.303257122, 0.045628945, 
     -0.25229914, -0.0135207234, -0.97146334, 0.744245448, -0.223764736, 1.45148602, 
     -0.811222414, 0.23662355, 0.321026621, 0.0275748433, -1.09046727, -0.152009957, 
     0.409498473, 0.0930091938, 2.08748786, -1.26599503, -1.40559212, -1.75677403, 
     -2.38986247, 0.151412278, -3.17824773, 1.32382727, 3.12823339, -4.8987026, 
     -1.53355427, -0.0127636922, -2.95671674, -1.57528348, 2.13990595, -2.38243296, 
     18.0436366, 2.58635163, 31.8895545, -12.6930907, 14.9338698, 26.6224834, 
     12.3908986, 4.83325901, 31.131157, 16.4932409, 10.2235659, 8.72976466, 
     2.23915424, -2.09084974, 0.114819845, -4.37959201, 1.98215147, 5.71046776, 
     -23.3594143, -10.3750393, 3.93801019, -2.55144319, -7.30807227, -1.11032766, 
     0.0157951678, 0.0149803878, -0.122027218, 0.0434083755, 0.773033368, -0.283916663, 
     0.199726629, -0.0830682729, 0.180656186, -0.0156743184, 0.238045471, 0.0746719962, 
     -0.107373328, -0.0176839946, 0.624770899, -0.169727452, -3.89191675, 1.04244686, 
     -1.15374579, 0.420000437, 0.0248252741, -0.44984049, 1.02105621, -1.21988231, 
     0.301570385, -0.0181454516, -0.556656461, 0.210222489, 3.96584852, -0.774500839, 
     1.31561475, -0.417351776, -1.50626037, 1.35694149, -1.77680066, 1.71939185, 
     2.4595161, -1.29237926, -1.17263755, 0.89570887, -3.25531138, 1.66144852, 
     -0.125071621, 0.530789762, -0.810896828, 0.964393702, -2.05221297, 0.993323268, 
     -21.6739305, 11.7380398, -9.5326682, -2.63533899, -14.4367267, -5.15873062, 
     0.39469908, -10.8664246, -7.3130796, -5.57138274, -16.3482692, -0.54181782, 
     7.34964226, -19.9890384, 2.73281086, 5.36441213, 6.60330658, -4.52672857, 
     16.4127969, 29.4510748, -2.66446559, -6.44227937, 4.643211, 2.83802949, 
     -0.12150547, 0.115628921, -1.08041068, 0.371859919, -0.220818674, -0.036854496, 
     -0.342057734, 0.203979022, -0.704260221, 0.144005231, -0.365012213, -0.194319273, 
     0.815703199, -0.578153808, 4.75817569, -1.63215027, 4.22953455, -0.998415808, 
     2.21403114, -1.03636245, 3.29058225, -0.165628215, 3.70396245, -0.163383428, 
     -0.917740869, 0.651559447, -3.99760236, 1.3676293, -6.99899689, 1.76831396, 
     -3.19368938, 1.20239234, -0.943706399, -1.3067496, -3.78270607, 1.02178751, 
     -0.709416136, 0.263397667, 2.99306571, -1.12582256, 3.61509607, -0.713742822, 
     0.692061552, -0.390437408, 1.7065565, -0.450042661, 2.93256325, -0.669532514, 
     6.31797326, -4.19357465, -7.48107146, 4.69735307, -6.18059922, 2.21260914, 
     -8.1117275, 6.81483127, -3.87659255, -0.121787073, -6.71721841, 2.93567606, 
     -0.993511215, 8.91217125, 6.33053387, -2.73816467, 9.71565507, -1.62987895, 
     5.20021193, -15.9475289, 1.91380998, 8.71498706, 12.9600316, -5.44637407, 
     -0.713455525, 0.368336424, 0.322503625, -0.0303928987, 1.45885443, -0.908690494, 
     1.85126768, -0.970069764, -1.14621616, 0.215473998, -4.39839505, 2.61460979, 
     -2.86038848, 1.52651248, 2.23240663, -0.575856224, 7.64065891, -4.52694331, 
     6.38155831, -3.0865823, 1.07694715, -0.99714312, 5.58257776, 3.39365416, 
     -11.4586462, 4.7388969, 0.15376482, 1.63601709, 1.60994116, -6.29270613, 
     6.05334594, -1.43331308, -2.77082572, -1.25608016, -17.6943741, 7.38729988, 
     -6.51537233, 3.20951704, -1.02732139, 0.737007114, -6.79868029, -3.52848462, 
     11.9499607, -4.99802346, -0.138905482, -1.02877772, 1.48504321, 6.17031378, 
     -5.60942057, 1.70438542, 3.17565702, 0.374176807, 12.7106498, -6.84501085, 
     -0.0872151572, 0.0637322219, -0.0732506763, 0.113837731, 0.175749447, -0.108233049, 
     0.183187421, -0.123409899, 0.0172142537, -0.214620792, -0.375288757, 0.225921402, 
     -0.108328954, 0.0681985803, 0.109092554, 0.0915136165, 0.225490792, -0.136255738, 
     0.665077773, -0.460586337, 6.20385068, -3.16997724, -1.1831282, 1.08642518, 
     -0.6083151, 0.37587419, -11.2796573, 5.12397666, 0.752617439, -1.04760773, 
     0.443767019, -0.232306863, 1.38270163, -1.38723481, -1.41851232, 0.599372435, 
     -0.772098733, 0.501072487, -6.71159558, 3.35196143, 1.16802809, -1.11414671, 
     0.796738914, -0.450950332, 13.126647, -5.55416109, -0.669202931, 1.12378233, 
     -0.609028712, 0.309736741, -2.30611343, 1.8296022, 2.17489452, -0.831752142, 
     0.0220278194, -0.0175107416, -0.464687579, 0.246715131, -0.444978361, 0.139082818, 
     0.00351299993, 0.0381600782, 1.40076842, -0.650100895, 1.20687423, -0.415564033, 
     -0.165160703, 0.00633850333, -2.52765663, 1.08097574, -2.26406136, 0.76905767, 
     -0.211033815, -1.75535865, -10.3589878, 1.90925362, 10.5005614, 0.729206721, 
     2.91128673, 0.29952203, 11.898422, -4.46822145, -0.710721297, -2.94069131, 
     -3.67458134, 1.80452237, 4.67879067, -1.66232868, -2.87538784, -0.723882488, 
     -0.563170303, 3.29108447, 17.4618135, -3.27955354, -17.7657378, -1.07825499, 
     -2.32391331, -1.85355494, -21.6312511, 7.88800131, 2.82754595, 5.22961892, 
     2.92547117, -1.01232272, -2.49586513, 0.635088064, 4.51334129, -0.832714393, 
     0.00644378147, -0.0155653738, 0.0617054481, -0.019713149, -0.433870023, 0.205151238, 
     -0.0787818144, 0.0747559329, -0.0667926258, 0.0344520855, 0.975571514, -0.592938556, 
     0.0586446376, -0.047374691, 0.0372930891, -0.0287519184, -0.503681043, 0.352072809, 
     8.13820829, -5.1782997, -9.01177885, 1.81666137, -0.0582298703, 11.4295696, 
     -7.83170385, 4.29764157, 5.27620598, -0.978682817, -3.03964191, -7.87391154, 
     1.66564381, -0.675676638, 1.73817125, -0.108185623, -0.591467501, 0.965374917, 
     -12.2232399, 7.97451074, 15.7731438, -3.61738135, 0.167609392, -17.7737727, 
     11.245419, -6.42610401, -10.8461657, 2.64942834, 3.73409028, 12.4026753, 
     -2.10001274, 0.870647079, -1.94260647, -0.09754426, 1.75285637, -1.6657141}},
    // B_Z (per nT Bz)
    {{
     -2.59788009, 8.14982255, -0.59029987, 2.39858458, -2.99416339, 7.7115394, 
     0.452626405, 4.72608627, -0.272537363, -1.98269735, -3.22582051, -4.24317345, 
     11.895373, -34.2431464, -2.0008476, -2.62073191, 2.51839887, -9.98529919, 
     1.19696614, -4.08382889, -0.568326252, 2.78991206, 0.0686039563, 1.73966945, 
     -15.2061032, 40.845242, 8.56337207, -3.7488696, 4.31549761, 5.79603762, 
     1.60841654, 3.12223875, 0.121944636, -0.950933182, 0.156450225, 0.270407561, 
     4.70102646, -2.81603398, -2.82470192, 2.55127757, 2.38879844, -9.98414943, 
     -2.37349835, -0.0104236451, 0.561149311, 0.108452624, -0.0384468881, 0.617769284, 
     -5.02760639, 11.500624, 0.936230775, -1.74769103, 3.35684727, -2.75204268, 
     11.9460893, -44.0140713, 3.84319753, -3.35131105, 6.57861519, -11.0216341, 
     1.43471187, 4.22031201, -0.745738186, 3.08647393, -4.29768754, 1.17962424, 
     0.410315709, 13.1699526, -0.810955666, 1.99344976, 0.380528119, 1.15854519, 
     -7.27055283, -22.2134132, 2.38112475, -3.44529524, -2.84703395, 9.74608362, 
     -0.0250294288, 0.527809509, 0.537353211, -0.815581312, 0.131117952, 0.250100797, 
     24.1221854, -8.71390236, -10.5054496, 7.12573598, -5.40308506, -5.0373569, 
     -2.23839218, -1.69838982, 0.476258047, -0.195507388, 0.0182535795, 0.659994651, 
     1.56895311, -1.13037841, -0.232846246, 0.159671106, -2.00051004, 0.417927369, 
     -7.27791803, 4.46807823, 1.19119933, -0.0951139846, 10.4171149, 4.43672539, 
     1.83466084, -4.92696745, -2.30449514, -1.48745061, -2.54102731, -32.1163353, 
     51.7983957, 116.010026, -0.396448345, 1.91595294, 7.28956732, 13.0585229, 
     0.126042533, -0.0880935393, 0.130390021, -0.100579295, -0.15854529, 0.176663576, 
     0.0302825865, -0.0149154868, 0.0510006623, -0.0487931185, -0.0572228992, 0.0629207201, 
     -0.557791678, 0.390284317, -0.547196598, 0.413581167, 0.731278872, -0.806247416, 
     -0.0994290121, 0.0428849953, -0.205130607, 0.203218789, 0.157694288, -0.247011709, 
     0.537116683, -0.375772289, 0.543322336, -0.392735646, -0.708797714, 0.787744698, 
     0.0439620239, -0.00131760037, 0.216162307, -0.215904127, 0.0279072956, 0.168819797, 
     0.218396834, -0.0873025365, -0.195604038, 0.116013805, -0.212699058, 0.18431742, 
     0.146490151, -0.145610464, -0.132664102, 0.0796726266, 0.00200147884, 0.0756118302, 
     -2.08953206, 0.542513954, 0.50474558, -0.00253741241, 1.73460996, -1.63905651, 
     -1.17800971, 0.982865046, 0.718185969, -0.251533028, -0.0954339693, -0.605419682, 
     4.32658584, 0.317591578, 0.0209082597, -1.06741943, -1.56562757, 2.11472558, 
     2.54923606, -1.82853319, -0.96305957, 0.0258396442, 0.526431993, 1.29262162, 
     0.0109691163, -0.0156956437, 0.0415118684, -0.0304266938, -0.087457786, 0.112557288, 
     0.012623279, 0.00170101986, 0.0137352144, -0.0140229679, -0.0852101168, 0.0200321917, 
     -0.249450033, 0.166342778, -0.235876637, 0.157112243, 0.717607024, -0.864397216, 
     0.0110703592, -0.0780132751, -0.199740493, 0.167246342, 0.0289934475, 0.149881659, 
     0.306644311, -0.142930709, 0.238307523, -0.105991323, -0.835817435, 1.05991897, 
     -0.212638903, 0.298577298, 0.324144361, -0.289160595, 0.781078437, -0.822211414, 
     0.434147512, -1.72825954, -0.014237169, 0.121585691, 0.228914773, 0.108182727, 
     0.392284506, -0.715794046, 0.214214316, 0.0627842602, -0.140308316, 0.31080258, 
     -0.948220526, 18.3762997, 0.361670241, -1.62745496, -0.481798217, -0.93721108, 
     -1.44099307, 5.76200273, -1.65604634, -1.43625167, 1.47845623, -1.64726869, 
     -18.0523496, -38.1779249, -0.526378848, 3.70099968, -2.77178803, -6.00172432, 
     -11.1618204, -3.12505318, 4.92563363, 4.29460324, 2.02097416, -10.0364978, 
     0.0204627904, -0.0237628787, 0.102067134, -0.114154616, 0.10588191, -0.153584319, 
     0.113043492, -0.107473404, -0.0677153927, 0.0489137782, -0.0761937457, 0.0207693375, 
     -0.165890554, 0.125680373, -0.703091157, 0.733261607, -0.610909882, 0.86033104, 
     -0.785867875, 0.792615718, 0.198162496, -0.219945278, 1.51322243, -1.2942555, 
     0.178442767, -0.0789343247, 0.980633633, -0.931538049, 0.75169183, -1.10356296, 
     1.19607365, -1.25445254, -0.187715231, 0.331077511, -3.42186775, 3.32490532, 
     -2.05600965, 4.67450376, 0.634085764, -0.661259689, 0.551471054, -0.945356453, 
     -2.23978779, 2.94249303, 0.834535337, -1.01777905, 0.747857875, -1.42073253, 
     10.4113591, -45.2562898, -7.66860599, 6.7227594, -0.660258356, 2.90466006, 
     15.5910491, -26.5853041, -9.40201173, 10.7792421, -3.32259051, 9.84843753, 
     21.9481887, 74.6492086, -1.06225619, 3.52661992, 4.7649748, 13.9444221, 
     9.25105796, 23.6942367, -3.1732695, -5.59544523, 3.58179687, 5.60609572, 
     0.051237799, -0.0848103477, -0.158254376, 0.128000429, -0.289447333, 0.290825986, 
     -0.029249354, 0.033778014, -0.15384722, 0.0657822098, -0.0956615512, -0.0202273449, 
     -0.282243334, 0.510600117, 0.746987802, -0.598456903, 1.80717389, -1.91324864, 
     0.563557124, -0.539392966, 0.581186149, -0.094014427, -0.359266258, 1.09159274, 
     0.354459974, -0.515844672, -0.871923703, 0.720592442, -2.04833005, 2.42193957, 
     -1.11984291, 1.15474964, -0.384170455, -0.412351792, 2.05004494, -2.93823396, 
     -0.0293394862, -0.85682361, -0.880442784, 1.02644453, -1.06936797, 0.477054098, 
     0.88394811, -1.32877949, -1.01633306, 1.30478852, -0.239207189, 0.506651075, 
     2.02961903, 7.49422631, 10.5759919, -9.98011197, -0.720317162, 1.81369921, 
     -7.45657728, 14.1486026, 13.2003851, -13.169922, -6.1544856, 1.63196804, 
     -4.56056977, 12.1971361, 2.53494836, -2.90855795, 0.975276479, -1.59697206, 
     -5.26065305, -4.70663046, 1.28864939, -7.2492554, 0.943565684, 0.98660033, 
     -0.176398973, 0.144178931, -0.248579736, 0.178857154, 0.0554215901, -0.0505307774, 
     -0.177233304, 0.150966792, -0.476884588, 0.462189659, 0.0868599582, -0.0505159548, 
     0.989373916, -0.639602996, 1.31982541, -0.931408821, -0.430858238, 0.467306144, 
     1.35588483, -1.12499357, 2.53118333, -2.93550204, -2.43411697, 2.2604711, 
     -1.00570667, 0.73648994, -1.49298006, 1.15671935, 0.161839001, -0.192737425, 
     -1.53731064, 1.3324496, -2.97946581, 3.95922445, 3.54168655, -3.18470726, 
     0.139696587, -1.09675017, 0.207174166, -0.640441597, 1.68154609, -0.442659688, 
     -0.105975204, -0.141984731, 0.910758627, -1.51923777, 1.36710934, -0.0490670962, 
     -4.1756962, 22.1447862, -1.67310791, 5.73353759, -4.70575082, 0.457724877, 
     -1.74487269, 4.33214338, -5.326222, 11.936158, 0.232512632, -11.0263251, 
     -21.0752386, -78.3001843, -4.01326241, -9.22961105, -0.19398552, -20.4879217, 
     -1.17709113, -20.1293269, -2.02489239, 6.15596072, -1.57045873, -0.4669383, 
     -0.0163083692, 0.0251679528, 0.79757201, -0.63389125, -1.42480191, 1.26988794, 
     0.246261527, -0.128315139, 0.892151184, -0.767306908, -1.71821574, 1.19635202, 
     0.411453447, -0.223898545, -3.6538116, 2.96500047, 6.76485313, -6.25682486, 
     0.102714533, -0.292955174, -4.27892112, 4.42607388, 5.3872129, -4.16839585, 
     -0.557623105, 0.529351526, 4.16767964, -3.46616869, -7.1476508, 6.70479951, 
     -1.29264421, 1.40280223, 5.14405696, -6.39193125, -3.02961515, 2.69031571, 
     1.20307615, -0.801367071, -0.459060906, 0.639802071, 0.748123386, -1.12769971, 
     0.71978319, -0.541106796, -1.07475845, 1.35164101, 0.2139667, -0.65148554, 
     -5.1621995, -3.80879481, 1.10563251, -3.68330304, -4.50053152, 6.20297423, 
     -4.89371778, 2.54354821, 3.86060261, -9.43531604, 3.54656429, 3.99677165, 
     17.4507511, 29.2891559, -0.620784475, 9.07536115, 8.86330846, 2.14807491, 
     8.70156133, 3.17453395, -2.20351232, 5.4116161, -5.42634145, 3.38499004, 
     -0.0539775887, -0.0322254631, -0.0536074695, 0.0539610405, -0.150712573, 0.151096189, 
     0.131517123, 0.0793817417, 0.126746441, -0.1412108, 0.43549327, -0.47758666, 
     -0.250587163, -0.0915646675, -0.182307358, 0.218263101, -0.821730146, 0.948167315, 
     1.14975657, 0.79863832, 1.48812401, -1.03463991, 3.14025371, -0.434961448, 
     -1.67614298, -1.93281562, -2.59684251, 1.9744734, -6.35143503, 0.990828967, 
     -0.590840704, 0.579429498, 1.99764325, -1.84263166, 8.5388883, -3.32546282, 
     -1.55701807, -0.594344679, -1.44394205, 0.984825431, -3.25545043, 0.366775097, 
     3.14996812, 1.45434655, 2.44948141, -1.80315398, 6.74135093, -0.728592435, 
     -9.84571131, -0.445379387, -1.50814506, 1.37093181, -9.9281068, 2.86537628, 
     -0.0997357429, 0.0530775273, 0.020695809, 0.0231790955, -0.0854421502, 0.110630828, 
     0.202672126, -0.110477448, 0.000917940332, -0.0902163446, 0.144012478, -0.208426509, 
     -0.104009643, 0.05838935, -0.012022858, 0.0606201619, -0.0550480732, 0.0950020641, 
     2.08496268, -0.986124775, -2.19298474, 0.547095989, 2.88294607, -2.52986402, 
     -4.27836991, 2.1569349, 3.40577865, 0.174360159, -5.00914875, 4.52842702, 
     2.33923851, -1.21735358, -2.88289781, -0.262696478, 2.23817282, -2.14965026, 
     -1.93924101, 0.876958646, 2.37126318, -0.617927138, -2.88909511, 2.50654643, 
     3.99467977, -1.93226607, -4.08582664, 0.0370423891, 5.15534862, -4.59709952, 
     -2.50818488, 1.23394511, 5.57128522, -0.24461186, -2.76597788, 2.63722593, 
     0.0374425745, -0.0438264144, 0.0339932918, -0.0383990241, 0.0721860494, -0.0800091435, 
     -0.111137257, 0.147408539, -0.122743853, 0.143397549, -0.292345388, 0.274063776, 
     0.190698514, -0.257554937, 0.266664112, -0.315525521, 0.691421989, -0.591254941, 
     -1.52983878, -0.0194252359, -1.26946009, 0.828546885, -0.606294895, 1.69089276, 
     1.67772742, -0.567426729, 2.12067269, -1.38633087, -0.688652946, -1.49590081, 
     1.0066895, -0.684861823, -2.71309973, 2.05124432, -1.35402403, 2.44707305, 
     2.78274127, -0.226327993, 1.96743575, -1.33079971, 1.01314396, -2.50360383, 
     -3.70883022, 1.5973476, -3.0219974, 2.00640072, 1.36219846, 1.71712946, 
     0.392697829, -0.286607058, 3.32456242, -2.44744853, 0.748751529, -2.32741195, 
     -0.0477481791, 0.0332239885, 0.0302848331, -0.0243310268, -0.0192314284, 0.0316600709, 
     0.0704320759, -0.0400735364, -0.0532278347, 0.050608849, -0.00479469014, -0.0334992587, 
     -0.0216435467, 0.00216528134, 0.0226124886, -0.025578609, 0.0168863984, 0.00308323081, 
     0.846435319, -1.40619401, -0.0460705305, -0.788680575, 0.764580555, -0.848733487, 
     -0.917946029, 2.12976011, 1.00402268, 0.0959327158, 0.292100643, 1.27416288, 
     -0.467869229, 0.0296428735, -0.475143071, 0.249275706, -0.401158213, -0.210433859, 
     -1.40585469, 2.09342682, -0.135460992, 1.39095719, -1.0251945, 1.06032868, 
     1.8116707, -3.3335241, -1.29961789, -0.351363984, -0.575434381, -1.70808978, 
     0.457313669, 0.0541902727, 0.665915773, -0.334581139, 0.523110004, 0.357437065}},
    // B_A (per nT h)
    {{
     -243.895211, -630.778687, 43.1486538, 26.4901177, -31.9532056, -110.61527, 
     467.8273, 1027.8308, 41.0551008, 43.8333497, 43.9647467, 168.907988, 
     275.809664, 798.10517, -100.507697, -88.1185957, 38.720895, 204.314957, 
     372.662724, 1289.44013, -85.223861, -152.536978, 71.1885261, 218.056374, 
     1.60173653, -245.58241, 105.402503, 123.908994, -6.26295824, -65.1137445, 
     -54.3809961, -207.246419, 17.8837263, 29.8327172, -8.98931485, -24.9898608, 
     -12.1113223, -0.291087467, -42.0998439, -38.7454956, -14.6172345, -128.170198, 
     14.9348708, 54.898934, 0.355056277, -2.97491495, 2.47827238, 10.6754691, 
     100.802011, 342.846026, 26.9007162, 103.198404, 35.7158089, 10.9705383, 
     -218.763224, -822.116168, -209.931819, -319.475182, -2.00606623, -107.538033, 
     -199.946382, -662.733705, 16.585194, -52.164092, -87.0560009, -104.712468, 
     315.370583, 1011.58437, -16.9693053, -9.51904775, 49.9345155, 147.199177, 
     128.262415, 461.438095, -47.1586121, -160.078321, 53.7277765, 176.178676, 
     -44.6314565, -165.504525, 8.63809284, 0.0309746408, -8.9821533, -21.5956032, 
     -57.6196331, -65.4136484, -22.1017306, 106.024779, 17.9072176, -20.4748554, 
     19.0752206, 76.0347748, -4.98182775, -14.4261823, 3.67848445, 15.8028217, 
     -8.25485883, -11.0177725, 22.1612495, -6.53721125, -61.1789254, 24.9188828, 
     -57.3920243, 31.0190861, -82.2599457, 21.8859459, 103.043557, 1.57426202, 
     245.405087, -61.105811, 4.61327858, -0.735797052, 59.1920294, -188.337012, 
     -376.000904, -569.322921, 52.9355991, 167.701941, -64.5405467, -123.475234, 
     2.82890285, -1.50166586, 2.03542229, -1.68836662, -6.17431086, 2.39681679, 
     1.14036894, -0.369043905, 1.91572034, -1.0176111, 4.24711622, -1.47911777, 
     -15.6501552, 8.04033265, -8.29022607, 7.13740765, 30.9981793, -11.0983373, 
     -4.44925332, 1.52175168, -9.42852289, 5.2674712, -26.6687371, 9.29080796, 
     18.1746295, -7.92486974, 7.60249354, -7.42191591, -35.817533, 9.87388541, 
     5.12748766, -2.86120008, 12.1588301, -7.41699294, 39.473623, -13.9432516, 
     12.778126, -10.780765, -0.899538093, 4.07559341, 7.57411136, 2.09167599, 
     -11.0692341, 7.0415278, -2.74966409, 4.04373544, -2.59724837, 2.6599612, 
     -78.0577803, 67.8501824, -3.80290446, -23.4877988, -42.2002298, -27.2252637, 
     78.1289937, -54.8368792, 4.68148479, -24.6750847, 25.4504918, -15.9409408, 
     82.5398307, -95.1032044, 32.2886856, 35.8443843, 60.5861293, 82.6816078, 
     -143.890349, 112.594683, 20.2394881, 36.4476752, -11.5761524, -7.7902621, 
     0.817543398, -0.587720811, -0.360675345, 0.635011139, -1.48555453, 3.24340541, 
     -1.29886346, 0.575493992, -2.82622091, 2.71330918, -7.76320715, 3.79245829, 
     -6.00430237, 1.51487486, 1.16452465, -3.19125068, 13.1188353, -16.810466, 
     14.7998501, -5.86254565, 20.3853464, -19.0603808, 47.6790152, -29.6917439, 
     16.6608717, -4.44489, 3.7003044, 0.355189312, -23.6995943, 11.3864545, 
     -27.1731505, 13.6329317, -31.3337548, 25.9520148, -45.7739946, 39.9441627, 
     -18.145821, -0.325357851, -3.03826212, -4.93638481, 5.31355058, -15.3946446, 
     16.799, -19.4734007, -1.6177016, -4.47225383, 26.9177146, -11.5606484, 
     140.60594, 18.4436936, 44.1841365, 33.3445531, 33.4556289, 125.609972, 
     -142.280823, 220.312796, 31.7391923, 40.0392299, -187.344962, 78.8895284, 
     -379.862932, 36.8432839, -113.143775, -38.840678, -176.234962, -372.242947, 
     248.745576, -659.192167, -48.1811389, -160.002409, 109.866698, -59.6936277, 
     -1.07878221, 0.346417586, -3.75804301, 1.40760744, 13.7474967, -5.03932336, 
     2.27874201, -0.630806666, 0.542463409, 0.0715349602, 3.50772113, -3.60089298, 
     4.6933686, -0.977673829, 31.7696044, -10.1612745, -81.4837838, 18.8338427, 
     -21.6498911, 4.62367311, -3.35635197, -0.388845083, -25.9712485, 32.7561392, 
     -0.350490168, -6.47043018, -54.2426202, 9.38476281, 107.168642, -13.8530134, 
     86.4202511, -15.2605844, 4.90071445, 1.1409609, -8.28015316, -59.2241347, 
     -21.3999904, 33.5530995, 1.13838661, -2.63806504, -7.33062293, 14.8363092, 
     -11.3915754, 17.9087891, -3.97688642, -5.19389543, -2.92154195, 15.235514, 
     119.832664, -321.994312, 7.24865346, 43.3485519, -6.93212478, -93.5141136, 
     42.8239125, -220.027809, 50.1062424, 29.2798658, -15.3746335, -151.964632, 
     406.811652, 453.820365, -62.8165554, -136.232856, 110.709174, 40.3650187, 
     343.721231, 664.265349, -37.5495354, 128.556215, 101.378683, 102.293229, 
     0.425227191, -0.541248646, 3.6682955, -0.62773237, 9.23606845, -2.82666774, 
     -4.11065413, 2.64124786, 1.02967503, -0.862047859, 6.50154886, -0.37197365, 
     -2.29550575, 3.78261107, -36.320588, 8.5864954, -52.5073856, 14.319866, 
     51.6163465, -24.5406998, -4.18605308, 7.98804116, -55.7072837, -14.5247767, 
     -21.6075951, 4.60581298, 58.7848916, -13.817679, 69.0641642, -35.2927944, 
     -98.3006629, 7.41560923, 8.91734833, -24.1590047, 189.889538, 47.3184319, 
     5.62553181, -10.1954979, 8.63117451, 3.6124987, -13.2838973, 0.0295058563, 
     3.67105674, -2.6306809, 8.45961107, 8.78609107, -19.5541978, -16.0761626, 
     0.0540232531, 59.7740487, -92.2698087, -29.09016, 15.3316093, 12.3366552, 
     -5.34808005, -17.6936546, -81.1050455, -89.769086, 156.549765, 189.183354, 
     120.797756, -7.95699356, 24.2049995, -38.7790965, 128.083248, 146.191577, 
     -26.3284947, 342.427086, -96.9023771, 91.9495303, 15.814283, 82.0690355, 
     0.885517573, -0.59726338, -5.22517995, 4.86223307, 4.78527381, -1.11003988, 
     0.104593664, 0.940506922, -2.30030882, 4.25435811, 14.5396838, -6.48154376, 
     -13.9942714, 13.1222595, 47.3048135, -36.2257808, 31.341875, -35.2538656, 
     -17.1020049, -4.64382292, 7.91708256, -30.7716989, -91.3563089, 56.8753933, 
     25.3158827, -15.4286523, -82.5347232, 51.6794711, -92.3718286, 77.6682863, 
     2.37417425, 28.4730742, -13.4504461, 61.0544301, -114.061007, -74.1593984, 
     19.1947286, -5.78532537, -5.49539584, -1.00966399, -24.5712101, -16.2296082, 
     9.35683683, -6.43291743, -0.119991681, -8.06778797, -24.2656892, 6.52429558, 
     -362.328674, 156.027105, 78.6455082, -44.8042549, 56.3809617, 139.184295, 
     -126.164866, 143.721058, 28.3272744, 89.5584171, 149.767459, -106.343946, 
     -542.297892, -717.533418, 291.42514, 403.751281, -13.52809, 191.348517, 
     -756.43018, -937.122582, -37.8025134, -69.9482494, -272.102898, -226.387811, 
     -1.72438362, -0.224528318, 9.35649888, -2.00610514, -9.43723303, -4.6077739, 
     -4.85201299, -0.433576839, 5.10423739, -1.14362494, -12.205105, -3.87567667, 
     23.1417021, -11.3935059, -50.4723586, 13.5862339, 4.99048418, 56.3799801, 
     49.7903013, -0.183881109, -15.200104, 8.42783198, 52.7332299, 47.2599263, 
     -29.155819, 22.035689, 54.608221, -25.6884708, 82.0884105, -118.496259, 
     11.4849042, -43.6062963, -74.1833101, -7.47140091, 178.437927, -29.3277204, 
     -0.32851458, -1.8305314, -5.33655273, -1.58433925, 23.6870266, 22.0220235, 
     -1.13096784, 1.03664572, -3.3592274, 0.820347907, 15.2435919, 14.2475809, 
     190.05288, 2.51935415, -21.2720633, 40.5180593, -10.5350177, -179.870188, 
     81.5841368, -42.8614709, -29.9889979, -15.8505349, -35.4576897, -92.7789526, 
     302.996177, 339.628745, -158.421846, -239.666415, -208.868912, -22.4570359, 
     352.780672, 489.586794, 293.175888, -76.0515834, -178.284942, 196.636614, 
     0.18308098, 0.160753187, -0.308976148, 0.621664529, -1.33477898, 0.412347056, 
     -0.3713289, -0.841317537, 1.50897566, -2.20313269, 5.35468323, -1.29969924, 
     0.431813107, 2.12021432, -3.48942451, 4.41419183, -12.3853855, 2.76828396, 
     -9.16365314, 9.959007, -5.76622905, 2.38865689, -25.7145861, 1.73927542, 
     21.5997499, -19.9813655, 7.8744416, -1.71542263, 47.6951749, -2.78728874, 
     -35.3670624, 19.0653413, 1.44074848, -7.08748372, -19.9738215, -15.1833783, 
     7.84178375, -9.55363197, 5.57143247, -2.39773052, 24.4744204, 0.303740064, 
     -18.346271, 19.6429836, -8.43701386, 2.41777254, -45.6771543, -2.70664503, 
     38.9968896, -27.2328444, 1.24375902, 5.3475874, 14.5367981, 30.635117, 
     -0.447904336, 0.306638708, -0.188932225, 0.0946958706, -0.216017298, -0.04481478, 
     0.725366645, -0.545705476, 0.567248403, -0.334653949, 0.42652999, 0.3233904, 
     -0.296904003, 0.261363323, -0.416209599, 0.301205125, -0.236064791, -0.299963349, 
     11.9834448, -5.95391521, 4.23001547, 2.06905364, 0.525856498, -8.04975314, 
     -18.1442729, 9.02831095, -15.7845032, 0.175029189, 2.88042792, 12.6431511, 
     6.86167462, -3.61096677, 12.1140896, -0.743453187, -3.04377412, -5.88239662, 
     -11.3808964, 5.64177903, -4.19911778, -2.84373439, -1.63001437, 8.21555797, 
     16.9314353, -8.41853637, 15.5857552, 2.98702661, -0.398066033, -13.6470992, 
     -6.73880145, 3.65679656, -8.30064791, -12.4479216, 1.4183731, 8.77127086, 
     1.04036256, -0.505091693, 1.07919173, -0.695334366, -0.920841149, -0.191118148, 
     -3.09890081, 1.6088565, -3.11817854, 2.0883611, 2.96845628, 0.258056615, 
     5.47655032, -2.99287021, 6.31631325, -4.43196112, -5.95346206, 0.134244388, 
     -5.7169951, -1.44063475, -11.38144, 13.5047562, 14.3099461, 0.448571392, 
     8.41524344, 3.83061319, -9.82783199, -1.09390713, 6.9591866, -11.2330983, 
     -15.8364337, 2.69772481, 0.137763669, 0.84095439, -16.9675894, -4.53555599, 
     8.60949014, 1.41778277, 20.8371296, -24.2450121, -27.8197159, 1.49960911, 
     -10.8581682, -5.15821796, 10.8164955, 7.4057093, -0.946291537, 15.3502744, 
     19.7455749, -3.93147851, -2.72653249, -2.89492634, 18.4835952, 5.99398731, 
     0.163011745, -0.0172298704, 0.325113341, -0.217077001, -0.843377049, 0.388854662, 
     -0.281529069, 0.0528962478, -0.472031457, 0.198529066, 1.79817695, -0.711970056, 
     0.102250219, -0.0187085281, 0.243936872, -0.0619026332, -0.903068587, 0.245685405, 
     -2.67191953, -5.79232433, -11.2845264, 20.8322335, -0.383934858, -4.33495949, 
     9.55902112, 4.11032754, -4.41408954, -11.5043856, -25.4058039, 20.0308449, 
     -8.29192324, 2.29197691, 3.06265836, 0.764447769, 26.2204056, -12.0985837, 
     -0.589138857, 11.2435945, 17.6168634, -32.8918098, 11.7140476, 1.41918877, 
     -5.40325064, -10.9492895, 6.25025566, 19.0641967, 17.7633449, -21.3372067, 
     8.72118074, -1.66845702, -3.74427766, -2.01726756, -32.0058644, 15.350751}},
    // B_Y (per nT By)
    {{
     3.82271982, -3.46408523, 0.274120901, -4.03009236, 2.88359869, 1.16834643, 
     -6.95384499, 0.613201643, 1.25264998, 2.21654439, -5.19150745, -4.07818136, 
     -2.82709821, 15.8388917, -5.76411796, 12.0992769, -3.95416081, -3.50527693, 
     -1.48666422, -6.64814571, -1.51769589, -2.81640436, 1.10099776, -1.45076907, 
     -5.97293403, -47.3093629, 4.9524205, -4.25768214, -0.449557322, 5.07354852, 
     -0.0919067869, -0.835519967, -0.391188621, 1.22030314, -0.350172065, 0.069287078, 
     1.99158824, 50.4281322, 10.2025151, -13.7750233, -0.0695778545, 1.39510119, 
     -0.630743393, -5.64442537, -0.136331205, 0.143744999, 0.0976875526, 0.0936183642, 
     -0.861346115, -6.34850952, 3.35422574, -1.99909923, 1.25890769, 0.608448182, 
     3.16576248, 17.5866459, -4.79967355, 9.48760346, -0.262306449, 4.60690464, 
     -3.86762806, 7.66830533, 1.88465042, -4.92533731, -3.63025486, -2.26039132, 
     -7.0280758, -15.7881123, 0.633981509, -2.77415922, -0.896775235, -4.21897419, 
     22.4939232, 11.5098151, -7.62848392, 6.19505156, 7.28227993, 11.3891142, 
     2.20337531, 4.55397229, 0.305374337, -0.236005596, 0.177554892, 1.63077866, 
     -20.5012969, -20.814075, -6.53391389, 9.44637288, -2.97403049, -19.7685497, 
     1.15194424, -0.89664251, -0.176835162, -0.0727233089, 0.457186368, 0.955486936, 
     0.233947624, -0.259626415, -0.25960213, 0.604129269, 0.780031936, -0.410671346, 
     -0.499380608, 1.59864675, 2.85219014, -5.52586743, -2.83458873, 0.969205717, 
     -0.443937356, -2.74638656, -7.92464822, 14.1796337, 0.321390244, 4.63410957, 
     1.94509744, -0.839491703, -4.45945721, -0.191324104, -0.576551465, -2.4747094, 
     -0.0176114831, 0.00488328852, 0.00246119972, 0.00630255149, 0.0191098167, -0.000807965326, 
     -0.00244112432, 0.000802759699, 0.00191406104, 0.00361071556, 0.0368136566, -0.0329603442, 
     0.0724394253, -0.0240997402, -0.00640825146, -0.0296698306, -0.0596209521, -0.00563292149, 
     0.00346540928, 0.00299998816, -0.0226156711, -0.00646964736, -0.202710183, 0.160966213, 
     -0.0660181133, 0.028874035, -0.00108605562, 0.0307779788, 0.0695217043, -0.00917653188, 
     0.008966412, -0.0165698035, 0.0417389804, -0.00764313634, 0.260347861, -0.16981399, 
     0.040099872, -0.0884808731, -0.0299977894, -0.0107591314, -0.0428705576, 0.0531439966, 
     -0.0555901685, 0.0500886384, -0.0558927142, 0.017119652, -0.0493198231, 0.00525772654, 
     -0.105301127, 0.660326272, 0.124849099, 0.0537620814, 0.213116885, -0.623977043, 
     0.546441302, -0.420714871, 0.393153988, -0.0923205852, 0.0269764512, 0.179884662, 
     -0.28439432, -1.06574846, 0.089116546, -0.314833897, 0.0344784161, 1.32312493, 
     -1.19149906, 0.742474791, -0.73263061, 0.312605562, 0.0770929829, 0.155238188, 
     0.00537190864, -0.00489473309, 1.03915246e-05, 0.00254144979, -0.0190195286, -0.020517061, 
     -0.0223519153, 0.0261649254, 0.0302678197, -0.0214762461, 0.0724585087, -0.114479119, 
     0.0186773336, -0.0106034227, -0.0183007939, -0.0138851725, -0.0188519198, 0.185891612, 
     0.0947967243, -0.138635534, -0.181733996, 0.101375822, -0.446306057, 0.754552789, 
     0.00622254774, 0.0124872704, -0.0394260707, 0.0415403823, 0.0156315658, -0.230579815, 
     -0.0263289649, 0.121665617, 0.225003786, -0.119468962, 0.419016355, -0.967027052, 
     0.0874386024, 0.149731507, 0.0515452139, -0.0960466531, -0.0260614184, 0.224839189, 
     0.193860499, -0.0786731878, -0.0865941279, 0.032791066, -0.236411665, 0.323851245, 
     -0.949426422, -1.66305956, -0.689335556, 0.426073307, -0.366091261, -0.919948092, 
     -1.75747242, 0.320558221, 0.296185715, 0.112326196, 0.794281555, -0.738340082, 
     -0.439396063, 8.43924995, 1.20223133, 3.7853435, 0.995041218, -3.6465043, 
     4.55218345, -0.239310189, 1.99539008, -0.381747647, 0.678914949, -4.49036659, 
     -0.0176072415, -0.00474060473, 0.0307643757, -0.0184982962, 0.115539258, -0.00251172217, 
     0.0136192366, -0.0623758515, -0.0237093061, 0.00956554055, -0.16693175, 0.451406937, 
     0.276425882, -0.0896819206, -0.387238968, 0.240488474, -0.947051913, 0.0645363, 
     0.0302656281, 0.415328803, 0.103980991, -0.0633173552, 0.515194723, -3.02716152, 
     -0.499141846, 0.191803413, 0.647079908, -0.432581509, 1.48505166, -0.0764828414, 
     -0.222565338, -0.589799607, -0.104079518, 0.0260934312, -0.039895359, 4.42587083, 
     -2.65027874, 2.1438076, 0.0396524155, 1.03749435, 1.71412518, -1.35670993, 
     -0.789217906, 0.297203379, 0.103886345, 1.05990919, 1.83698404, -0.833936098, 
     20.3497548, -15.3638482, 0.268988683, -9.25489673, -2.30489184, -0.289251804, 
     4.36735409, 0.484010796, -0.0451426089, -10.3988144, -6.0013113, -3.77070219, 
     3.76484094, -11.227907, -2.53130116, -6.7426552, -0.466375426, -5.39180868, 
     -0.430225968, -3.54219445, 1.061453, -2.53643702, -0.283468267, 0.524397403, 
     0.0983140082, -0.0698194608, -0.211257019, 0.133572371, -0.360932078, 0.114251348, 
     -0.00381068872, 0.0908223533, -0.134577615, 0.0828556071, -0.208985716, -0.368576551, 
     -0.558593272, 0.306754548, 1.03100019, -0.71792446, 1.79367581, -0.108017013, 
     -0.11556437, -0.636598823, 0.328850803, -0.286970422, 0.487811173, 3.70625311, 
     0.790241187, -0.355232608, -1.27777364, 0.919711249, -2.36040508, -0.101391149, 
     0.425608597, 0.84341318, -0.172049549, 0.358515523, -0.549402661, -5.99129035, 
     3.78547908, -3.10569067, -0.459919402, -1.08403607, -2.96778657, 1.06290606, 
     1.21591866, -0.578920363, 0.108467834, -1.52087759, -2.65996309, 0.765246057, 
     -29.7765569, 24.0463766, 1.03424329, 14.2028919, -0.498260254, 16.9460172, 
     -6.81645208, 0.379735144, -4.72315013, 18.215348, 4.36166466, 14.4587169, 
     -3.86982107, 1.63567296, -2.24878438, -0.174810898, 1.3990424, 3.0884603, 
     -3.19782027, 0.35126812, -1.40959009, 2.12256781, 2.39063444, 7.24186455, 
     -0.0332374476, -0.0265886029, -0.195207647, 0.14076606, -0.369237875, 0.579941216, 
     -0.074213869, 0.0188980371, -0.448643968, 0.23283228, -0.733193886, 0.763861788, 
     0.294811413, -0.00226610389, 0.288491457, -0.218940212, 0.804611584, -1.98760118, 
     0.448805624, -0.0732985552, 2.24211939, -1.00348529, 3.61288524, -4.52845829, 
     -0.439425349, 0.0752503314, 0.0764366223, -0.093605447, -0.279089916, 1.85457061, 
     -0.649675024, -0.0193376197, -2.49301771, 0.88669691, -4.03180176, 6.2966039, 
     -1.54932954, 1.26770767, 1.01733739, -0.397494338, 3.91382713, -2.08063151, 
     -0.80946176, 0.166596129, 0.985809299, -0.240133741, 3.11679827, -0.928959808, 
     13.0108116, -9.27380159, 0.409202203, -7.87718921, 0.931449788, -20.8792353, 
     5.46897377, 2.1979887, 3.53663605, -10.4553921, -0.625875227, -22.6836874, 
     3.4613656, 1.54291911, 4.20841295, 7.53408665, 1.99868385, 3.96573971, 
     -0.4201876, 1.52174532, 0.350126735, 1.86483127, -0.203896987, -0.665408732, 
     0.0248452673, -0.00639511502, -0.333611995, 0.265297011, -0.234268155, -0.03760866, 
     -0.116340896, 0.0361146719, 0.0613831973, -0.109246445, 0.303115323, -0.448597436, 
     -0.242289945, 0.0512645811, 1.75279293, -1.38998522, 0.162608494, 1.3568327, 
     0.152391761, 0.409474219, -1.33833019, 0.767418208, -2.41654604, 2.24978413, 
     0.307620789, -0.0504378624, -2.0784722, 1.67729803, 0.679021206, -2.37538663, 
     0.224947003, -0.836563989, 1.95041784, -0.539569712, 4.36695734, -3.8904757, 
     0.230210679, -0.255546312, 0.0618122621, 0.048893424, -1.75676181, 1.44665208, 
     0.45882522, 0.0227075449, -0.531376198, 0.448507056, -1.33462639, 0.454099311, 
     -2.3956427, 1.34092054, -3.74416972, 4.50894255, 0.318649102, 6.32894891, 
     -2.4517628, -2.90616788, -0.622068452, 3.12302589, -0.0363096345, 13.1248788, 
     -2.73147583, 0.7929943, 1.91005574, -6.18271351, -3.58271972, 1.53912699, 
     0.948228912, 1.63689939, -0.689471834, -1.99534029, -3.13349192, -2.38473353, 
     -0.183132685, 0.116289616, 0.0475569864, -0.0465456942, -0.106159717, 0.114269385, 
     0.445858045, -0.244285672, -0.144980729, 0.136415011, 0.172362877, -0.24800588, 
     -0.668487786, 0.290400541, 0.277837594, -0.248730352, 0.00518513732, 0.218605657, 
     2.66122728, -2.74232094, -0.558493984, 0.603403468, 1.04740455, -0.548866331, 
     -4.54212919, 4.82013898, 1.24889232, -1.35213819, -2.7923894, 1.93219445, 
     3.74816706, -3.89095017, -1.89955248, 1.70601511, 2.15145543, -1.95549281, 
     -2.67674433, 2.76262212, 0.577927215, -0.621589902, -0.813768807, 0.295128135, 
     4.73223461, -5.0449904, -1.30082996, 1.4031503, 2.13464012, -1.15008352, 
     -4.86787326, 5.16070625, 2.17305229, -1.8565269, -1.01257143, 0.449032395, 
     -0.00802922209, 0.0164796251, -0.0107731223, 0.0147062545, -0.0575934148, 0.0499924308, 
     -0.0520468352, 0.0080024174, 0.038155737, -0.0431608603, 0.129917249, -0.125648801, 
     0.0558765505, -0.0212916674, -0.0315343061, 0.0323314644, -0.0782406663, 0.0821318319, 
     1.98502892, -1.39045245, -0.254267331, 0.00640824504, 0.659416773, 0.170599001, 
     -2.14517487, 1.68155952, 0.0968652175, 0.429441658, -0.995714242, -0.40870363, 
     0.00952448372, -0.161811339, 0.765758294, -1.12316722, 0.358644593, 0.146498567, 
     -1.96627825, 1.36123151, 0.217190061, 0.0312895436, -0.643702829, -0.202554384, 
     2.14227795, -1.64082999, 0.0569871133, -0.616446064, 0.830499048, 0.638028859, 
     0.105113684, 0.0686219271, -1.69167879, 2.20950224, 0.120471608, -0.669565846, 
     0.247223399, -0.238625772, -0.0465166869, 0.0373261346, -0.148882966, 0.119152895, 
     -0.680996852, 0.621896326, 0.0909994893, -0.0648589428, 0.447008539, -0.332472792, 
     1.12991719, -0.971306627, -0.0756692707, 0.0322601545, -0.800725869, 0.547959409, 
     -2.26451535, 3.08462007, 1.4696812, -1.2362992, -0.615585637, -0.0896663339, 
     3.69252952, -5.53005117, -2.88359101, 2.62800992, 0.446326716, 1.02627485, 
     -8.32656351, 8.75372037, 1.94496195, -2.32038142, -0.428269937, -1.41386952, 
     2.40290676, -3.50537679, -1.89437091, 1.49639708, 1.04418621, -0.206033614, 
     -2.58474565, 5.23362989, 3.61454203, -3.10981408, -0.811942952, -0.809729005, 
     7.96246911, -8.66597911, -2.05971599, 2.45813447, 0.892487127, 1.08138293, 
     -0.0347834227, -0.0351512882, -0.00753530526, 0.00678293213, -0.0221204728, 0.00610591388, 
     0.117100707, 0.0780372824, 0.00687872452, -0.0118325405, 0.0536912563, -0.0203283201, 
     -0.104543877, -0.0227579945, -0.00299315642, 0.00802425022, -0.017149884, 0.0048587618, 
     -1.33980147, 1.47423181, 0.693360507, -0.217842326, 0.549733221, -0.646977384, 
     -3.00622456, -1.39378575, -0.445477948, 0.076003909, -2.21100925, 1.80497817, 
     10.5452248, -0.815643842, 0.120981963, -0.127124242, 0.432155154, -0.441397337, 
     2.30514141, -1.82112551, -1.12276466, 0.382839086, -0.499620717, 0.804462761, 
     3.93882701, 0.56886126, 0.752162162, -0.174293179, 2.75457865, -2.37856023, 
     -24.4783378, 2.28143058, -0.2064243, 0.212276127, -0.354216287, 0.463449153}},
}};

// ---- END GENERATED ---------------------------------------------------------------------------

/// The basis evaluated at one point: one field vector per basis function.
template <std::floating_point T>
using T96Basis = std::array<std::array<T, 3>, t96_basis_count>;

/**
 * Every basis field at one GSM point.
 *
 * This is the whole of the physics in this file, and the fitting harness calls THIS function, so
 * the shipped evaluator is the thing that was fitted rather than a copy of it. Three blocks, each
 * one published class; see the file brief for the citations.
 *
 * **Disc modules, in SM.** The GSM point is rotated about `y` by the tilt, the warped sheet surface
 * and its slopes are formed (T89 eq. 11), the thickness profiles and their slopes (T89 eq. 13),
 * then for each truncation position a `W` and its gradient, and for each tail scale the `A^(1)`
 * and `A^(2)` fields (T89 eqs. 14-15, separated by coefficient), and for each ring radius the
 * `A^(3)` field (eqs. 16-17). With `P = A_phi / rho`, `B_x = -x dP/dz`, `B_y = -y dP/dz`,
 * `B_z = 2P + x dP/dx + y dP/dy`; every `d/dx`, `d/dy` below is the analytic derivative of the
 * line it follows and none is optional — dropping one breaks `div B = 0`.
 *
 * **Toroidal modules, in SM.** `B = ∇T × r` has, in spherical components, `B_r = 0`,
 * `B_theta = (1/sin theta) dT/dphi`, `B_phi = -dT/dtheta`. With `T = t(r) s(theta, r) A_m`,
 * `A_1 = sin theta sin phi = y/r`, `A_2 = sin^2 theta sin 2phi = 2xy/r^2` (both smooth on the
 * axis), `s = S(theta - theta_0) + S(pi - theta - theta_0)`, `S(u) = (1 + tanh(u/w)) / 2`, and
 * `theta_0(r) = asin(sqrt(min(1, r/L)))` — the colatitude at which the dipole shell `L` is crossed
 * at radius `r`. No radial derivative of `T` is ever needed, which is why the shell may depend on
 * `r` at no cost.
 *
 * **Box harmonics, in GSM.** `B = -∇U` for the four parity families; `kappa^2 = p^-2 + r^-2`
 * makes each `U` harmonic.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device.
 * @param g the geometry — @ref t96_geometry in production, a candidate in the fitting harness.
 * @param sin_tilt `sin(psi)`, precomputed per epoch.
 * @param cos_tilt `cos(psi)`; must be non-zero — the sheet warp carries `tan(psi)`.
 * @param x the GSM x coordinate, R_E.
 * @param y the GSM y coordinate, R_E.
 * @param z the GSM z coordinate, R_E.
 * @param out receives one `{B_x, B_y, B_z}` in GSM per basis function, in the index order the
 *            constants above lay out: discs, then toroidal, then box.
 * @complexity O(1): about 2 200 flops, ~20 square roots, 9 exponentials, 12 sine/cosine pairs,
 *             one `atan2`, six `asin` and 24 `tanh`. No loop with a data-dependent trip count and
 *             no branch on data except the axis guards, which are selects.
 * @alloc none.
 * @test IrbemT96.DivergenceVanishesEverywhere
 * @test IrbemT96.EveryBasisFieldIsDivergenceFreeOnItsOwn
 * @test IrbemT96.DawnDuskSymmetryOfTheDiscAndBoxBlocks
 */
template <std::floating_point T>
inline void t96_basis(const T96Geometry& g, T sin_tilt, T cos_tilt, T x, T y, T z,
                      T96Basis<T>& out) {
    const T half = static_cast<T>(0.5);
    const T one = static_cast<T>(1);
    const T two = static_cast<T>(2);
    const T three = static_cast<T>(3);
    const T zero = static_cast<T>(0);

    const T tan_tilt = sin_tilt / cos_tilt;

    // ---- GSM -> SM ---------------------------------------------------------------------------
    const T xs = (x * cos_tilt) - (z * sin_tilt);
    const T ys = y;
    const T zs = (x * sin_tilt) + (z * cos_tilt);

    // ---- T89 eq. (11): the warped sheet, its slopes -------------------------------------------
    const T hinge = xs + static_cast<T>(g.r_c);
    const T hinge_root = std::sqrt((hinge * hinge) + static_cast<T>(t89_fixed.hinge_scale2));
    const T y2 = ys * ys;
    const T y4 = y2 * y2;
    const T ly2 = static_cast<T>(t89_fixed.l_y) * static_cast<T>(t89_fixed.l_y);
    const T ly4 = ly2 * ly2;
    const T bend_den = y4 + ly4;
    const T sheet = (half * tan_tilt * (hinge - hinge_root)) -
                    (static_cast<T>(g.g) * sin_tilt * y4 / bend_den);
    const T dsheet_dx = half * tan_tilt * (one - (hinge / hinge_root));
    const T dsheet_dy = -static_cast<T>(g.g) * sin_tilt * static_cast<T>(4) * ys * y2 * ly4 /
                        (bend_den * bend_den);
    const T zr = zs - sheet;
    const T rho2 = (xs * xs) + (ys * ys);

    // ---- T89 eq. (13): the thickness profiles, their slopes ----------------------------------
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
    const T d_tail = static_cast<T>(g.d_0) + (delta * y2) + (gamma_t * h_t) + (gamma_1 * h_1);
    const T d_ring = static_cast<T>(g.d_0) + (static_cast<T>(g.gamma_rc) * h_rc) + (gamma_1 * h_1);
    const T dd_tail_dx = (gamma_t * dh_t) + (gamma_1 * dh_1);
    const T dd_tail_dy = two * delta * ys;
    const T dd_ring_dx = (static_cast<T>(g.gamma_rc) * dh_rc) + (gamma_1 * dh_1);
    const T xi_t = std::sqrt((zr * zr) + (d_tail * d_tail));
    const T xi_rc = std::sqrt((zr * zr) + (d_ring * d_ring));
    const T sigma = (xs * dd_tail_dx) + (ys * dd_tail_dy);

    // ---- the disc modules, in SM, then rotated -----------------------------------------------
    // A local (bx, by, bz) in SM is committed to `out` by this lambda so the rotation is written
    // once. The index is advanced by the caller's loop; the lambda captures by reference.
    std::size_t j = 0;
    const auto commit = [&](T bx_sm, T by_sm, T bz_sm) {
        out[j][0] = (bx_sm * cos_tilt) + (bz_sm * sin_tilt);
        out[j][1] = by_sm;
        out[j][2] = (-bx_sm * sin_tilt) + (bz_sm * cos_tilt);
        ++j;
    };

    const T dx2 = static_cast<T>(t89_fixed.d_x) * static_cast<T>(t89_fixed.d_x);
    const T d_y2 = static_cast<T>(g.d_y) * static_cast<T>(g.d_y);
    for (std::size_t a = 0; a < t96_tail_truncation_count; ++a) {
        // T89 eq. (13): the truncation factor W and its slopes, for this x_0.
        const T xw = xs - static_cast<T>(g.x_0[a]);
        const T w_root = std::sqrt((xw * xw) + dx2);
        const T w_x = half * (one - (xw / w_root));
        const T w_y = one / (one + (y2 / d_y2));
        const T w = w_x * w_y;
        const T dw_dx = w_y * (-half * dx2 / (w_root * w_root * w_root));
        const T dw_dy = w_x * (-two * ys / d_y2) * w_y * w_y;
        const T grad_w = (xs * dw_dx) + (ys * dw_dy);
        for (std::size_t b = 0; b < t96_tail_scale_count; ++b) {
            const T u_t = static_cast<T>(g.a_t[b]) + xi_t;
            const T s_t = std::sqrt(rho2 + (u_t * u_t));
            // T89 eqs. (14)-(15), the C_1 (A^(1)) coefficient ...
            const T q1 = w / (xi_t * s_t * (s_t + u_t));
            const T bx1 = q1 * xs * zr;
            const T by1 = q1 * ys * zr;
            commit(bx1, by1,
                   (w / s_t) + (grad_w / (s_t + u_t)) + (bx1 * dsheet_dx) + (by1 * dsheet_dy) -
                       (q1 * d_tail * sigma));
            // ... and the C_2 (A^(2)) coefficient.
            const T q2 = w / (xi_t * s_t * s_t * s_t);
            const T bx2 = q2 * xs * zr;
            const T by2 = q2 * ys * zr;
            commit(bx2, by2,
                   (w * u_t / (s_t * s_t * s_t)) + (grad_w / (s_t * (s_t + u_t))) +
                       (bx2 * dsheet_dx) + (by2 * dsheet_dy) - (q2 * d_tail * sigma));
        }
    }
    for (std::size_t b = 0; b < t96_ring_scale_count; ++b) {
        // T89 eqs. (16)-(17): the A^(3) disc, the only one with a finite moment.
        const T u_rc = static_cast<T>(g.a_rc[b]) + xi_rc;
        const T s_rc = std::sqrt(rho2 + (u_rc * u_rc));
        const T s_rc2 = s_rc * s_rc;
        const T s_rc5 = s_rc2 * s_rc2 * s_rc;
        const T q3 = three * u_rc / (xi_rc * s_rc5);
        const T bx3 = q3 * xs * zr;
        const T by3 = q3 * ys * zr;
        commit(bx3, by3,
               (((two * u_rc * u_rc) - rho2) / s_rc5) + (bx3 * dsheet_dx) + (by3 * dsheet_dy) -
                   (q3 * d_ring * xs * dd_ring_dx));
    }

    // ---- the toroidal field-aligned-current modules, in SM, then rotated --------------------
    {
        const T r = std::sqrt(rho2 + (zs * zs));
        const T rho = std::sqrt(rho2);
        const T cos_th = zs / r;
        const T sin_th = rho / r;
        const T theta = std::atan2(rho, zs);
        // On the axis phi is undefined and A_1, A_2 vanish with sin theta; any direction gives the
        // same (continuous) field there, so noon is chosen.
        const bool off_axis = rho > zero;
        const T cos_ph = off_axis ? xs / rho : one;
        const T sin_ph = off_axis ? ys / rho : zero;
        const T sin2ph = two * sin_ph * cos_ph;
        const T cos2ph = (cos_ph * cos_ph) - (sin_ph * sin_ph);
        const T inv_r = one / r;
        const T pi = static_cast<T>(std::numbers::pi);
        for (std::size_t li = 0; li < t96_fac_shell_count; ++li) {
            const T ratio = r / static_cast<T>(g.fac_shell[li]);
            const T sin_th0 = std::sqrt(ratio < one ? ratio : one);
            const T th0 = std::asin(sin_th0);
            for (std::size_t wi = 0; wi < t96_fac_width_count; ++wi) {
                const T wdt = static_cast<T>(g.fac_width[wi]);
                const T a1 = std::tanh((theta - th0) / wdt);
                const T a2 = std::tanh(((pi - theta) - th0) / wdt);
                const T s = (half * (one + a1)) + (half * (one + a2));
                const T ds = (half * (one - (a1 * a1)) / wdt) - (half * (one - (a2 * a2)) / wdt);
                T t = inv_r;
                for (std::size_t pi_ = 0; pi_ < t96_fac_profile_count; ++pi_) {
                    // t = r^-1, r^-2, r^-3 in turn.
                    for (std::size_t mi = 0; mi < t96_fac_mode_count; ++mi) {
                        T b_theta = zero;
                        T b_phi = zero;
                        if (mi == 0) {
                            // A_1 = sin theta sin phi: dA/dphi / sin theta = cos phi,
                            // dA/dtheta = cos theta sin phi.
                            b_theta = t * s * cos_ph;
                            b_phi = -t * ((ds * sin_th * sin_ph) + (s * cos_th * sin_ph));
                        } else {
                            // A_2 = sin^2 theta sin 2phi: dA/dphi / sin theta = 2 sin theta cos 2phi,
                            // dA/dtheta = 2 sin theta cos theta sin 2phi.
                            b_theta = t * s * two * sin_th * cos2ph;
                            b_phi = -t * ((ds * sin_th * sin_th * sin2ph) +
                                          (s * two * sin_th * cos_th * sin2ph));
                        }
                        // theta-hat = (cos th cos ph, cos th sin ph, -sin th); phi-hat = (-sin ph, cos ph, 0)
                        commit((b_theta * cos_th * cos_ph) - (b_phi * sin_ph),
                               (b_theta * cos_th * sin_ph) + (b_phi * cos_ph), -b_theta * sin_th);
                    }
                    t *= inv_r;
                }
            }
        }
    }

    // ---- the box harmonics, in GSM -------------------------------------------------------------
    for (std::size_t fam = 0; fam < t96_box_parity_count; ++fam) {
        for (std::size_t pi_ = 0; pi_ < t96_box_scale_count; ++pi_) {
            const T p = static_cast<T>(g.box_p[pi_]);
            const T cy = std::cos(y / p);
            const T sy = std::sin(y / p);
            for (std::size_t rk = 0; rk < t96_box_scale_count; ++rk) {
                const T rr = static_cast<T>(g.box_r[rk]);
                const T cz = std::cos(z / rr);
                const T sz = std::sin(z / rr);
                const T kappa = std::sqrt((one / (p * p)) + (one / (rr * rr)));
                const T e = std::exp(kappa * x);
                switch (fam) {
                    case 0:  // U = e cos(y/p) sin(z/r)
                        out[j] = {-kappa * e * cy * sz, (e / p) * sy * sz, -(e / rr) * cy * cz};
                        break;
                    case 1:  // U = e cos(y/p) cos(z/r)
                        out[j] = {-kappa * e * cy * cz, (e / p) * sy * cz, (e / rr) * cy * sz};
                        break;
                    case 2:  // U = e sin(y/p) sin(z/r)
                        out[j] = {-kappa * e * sy * sz, -(e / p) * cy * sz, -(e / rr) * sy * cz};
                        break;
                    default:  // U = e sin(y/p) cos(z/r)
                        out[j] = {-kappa * e * sy * cz, -(e / p) * cy * cz, (e / rr) * sy * sz};
                        break;
                }
                ++j;
            }
        }
    }
}

// -------------------------------------------------------------------------------------------
// The drivers, collapsed to one amplitude per basis function
// -------------------------------------------------------------------------------------------

/**
 * The Birkeland driver `h = B_t sin(theta/2)`, eq. (D): the measured third clock-angle function.
 *
 * Written without the angle: `sin^2(theta/2) = (1 - cos theta) / 2` and `cos theta = Bz / B_t`, so
 * `h = sqrt(B_t (B_t - Bz) / 2)`. It is `|Bz|` for pure southward IMF, `0` for pure northward,
 * `|By| / sqrt(2)` for pure By, and homogeneous of degree one in `(By, Bz)` — which is what makes
 * the oracle exactly linear along every ray from the origin of the IMF plane and kinked only at
 * the northward ray. Even in By by construction.
 *
 * @param by_nt IMF By, GSM, nT.
 * @param bz_nt IMF Bz, GSM, nT.
 * @return `h` in nT; zero for a zero IMF.
 * @complexity O(1): two square roots.
 * @alloc none.
 * @test IrbemT96.ClockDriverIsTheMeasuredThirdFunction
 */
[[nodiscard]] inline double t96_clock_driver(double by_nt, double bz_nt) {
    const double bt = std::sqrt((by_nt * by_nt) + (bz_nt * bz_nt));
    const double arg = 0.5 * bt * (bt - bz_nt);
    return std::sqrt(arg > 0.0 ? arg : 0.0);  // roundoff can make bt - bz slightly negative
}

/**
 * The 128 amplitudes for one epoch and one driver set, plus the tilt's sine and cosine.
 *
 * Everything @ref t96_components needs beyond the point. The `float` instantiation is how the host
 * lane is made to run the device kernel's arithmetic: the collapse is done in `double` and rounded
 * ONCE, so both lanes multiply the same `float` amplitudes against the same `float` basis.
 *
 * @tparam T the scalar type; `double` or `float`.
 * @test IrbemT96.AmplitudesCollapseTheTablesLinearly
 */
template <std::floating_point T>
struct T96Amplitudes {
    /// One amplitude per basis function, in @ref t96_basis's index order.
    std::array<T, t96_basis_count> c;
    /// `sin(psi)`.
    T sin_tilt;
    /// `cos(psi)`.
    T cos_tilt;
};

/**
 * Collapse the tables, the tilt, the pressure and the drivers into one amplitude per basis field.
 *
 * Eq. (1) evaluated for its coefficients: with the family weights `{1, Dst, Bz, h, By}`, the tilt
 * factors `{1, sin psi, sin^2 psi}` and the pressure factors `{1, sqrt(Pdyn)}`,
 *
 *     c_j = sum_f w_f sum_{t,p} tables[f][6 j + 2 t + p] T_t(psi) P_p(Pdyn).
 *
 * This is the per-epoch cost of the model: 3 840 multiply-adds, paid once per batch and once per
 * @ref t96_field call. Pdyn is taken as given: a negative pressure has no square root and no
 * meaning, and @ref t96_field refuses it before calling here.
 *
 * @tparam T the scalar type the amplitudes are wanted in.
 * @param tilt_rad the dipole tilt, radians.
 * @param dst_nt Dst, nT.
 * @param pdyn_npa solar-wind dynamic pressure, nPa; must be `>= 0`.
 * @param by_nt IMF By, nT.
 * @param bz_nt IMF Bz, nT.
 * @return the amplitudes, computed in `double` and converted to @p T once.
 * @complexity O(`t96_family_count x t96_column_count`) — 3 840 multiply-adds.
 * @alloc none; the returned object is inline storage.
 * @test IrbemT96.AmplitudesCollapseTheTablesLinearly
 * @test IrbemT96.AmplitudesRoundTripThroughFloat
 */
template <std::floating_point T>
[[nodiscard]] inline T96Amplitudes<T> t96_amplitudes(double tilt_rad, double dst_nt,
                                                    double pdyn_npa, double by_nt,
                                                    double bz_nt) {
    const double sp = std::sin(tilt_rad);
    const double cp = std::cos(tilt_rad);
    const std::array<double, t96_tilt_factor_count> tf{1.0, sp, sp * sp};
    const std::array<double, t96_pressure_factor_count> pf{1.0, std::sqrt(pdyn_npa)};
    const std::array<double, t96_family_count> weight{1.0, dst_nt, bz_nt,
                                                      t96_clock_driver(by_nt, bz_nt), by_nt};
    std::array<double, t96_factor_count> factor{};
    for (std::size_t t = 0; t < t96_tilt_factor_count; ++t) {
        for (std::size_t p = 0; p < t96_pressure_factor_count; ++p) {
            factor[(t * t96_pressure_factor_count) + p] = tf[t] * pf[p];
        }
    }
    T96Amplitudes<T> out{};
    for (std::size_t j = 0; j < t96_basis_count; ++j) {
        double acc = 0.0;
        for (std::size_t f = 0; f < t96_family_count; ++f) {
            double inner = 0.0;
            for (std::size_t k = 0; k < t96_factor_count; ++k) {
                inner += t96_tables[f][(j * t96_factor_count) + k] * factor[k];
            }
            acc += weight[f] * inner;
        }
        out.c[j] = static_cast<T>(acc);
    }
    out.sin_tilt = static_cast<T>(sp);
    out.cos_tilt = static_cast<T>(cp);
    return out;
}

// -------------------------------------------------------------------------------------------
// The evaluator
// -------------------------------------------------------------------------------------------

/**
 * The T96 external field at one GSM point: the basis dotted with the amplitudes.
 *
 * The sum runs in basis order, `j = 0 .. 127`, with each product added as it is formed — the
 * same order and the same association the device kernel uses, so the two fp32 lanes differ only by
 * a driver's contraction choices.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device.
 * @param a the collapsed amplitudes and the tilt, from @ref t96_amplitudes.
 * @param x the GSM x coordinate, R_E.
 * @param y the GSM y coordinate, R_E.
 * @param z the GSM z coordinate, R_E.
 * @return `{B_x, B_y, B_z}` in GSM, nanotesla.
 * @complexity O(1); one @ref t96_basis plus 384 multiply-adds.
 * @alloc none — one stack basis array of `128 x 3` scalars.
 * @test IrbemT96.ReferenceLaneMatchesTheComponentForm
 * @test IrbemT96.HostFloatLaneTracksTheReferenceLane
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, 3> t96_components(const T96Amplitudes<T>& a, T x, T y, T z) {
    T96Basis<T> basis{};
    t96_basis<T>(t96_geometry, a.sin_tilt, a.cos_tilt, x, y, z, basis);
    T bx = static_cast<T>(0);
    T by = static_cast<T>(0);
    T bz = static_cast<T>(0);
    for (std::size_t j = 0; j < t96_basis_count; ++j) {
        bx += a.c[j] * basis[j][0];
        by += a.c[j] * basis[j][1];
        bz += a.c[j] * basis[j][2];
    }
    return {bx, by, bz};
}

/**
 * The T96 external field at one GSM point, in `double`, for amplitudes already collapsed.
 *
 * The entry a trace uses: the amplitudes are a property of the epoch and the drivers, so a loop
 * over points collapses once and calls this.
 *
 * @param p the position, GSM, in Earth radii.
 * @param a the amplitudes for the epoch and driver set.
 * @return the external field at @p p, GSM, in nanotesla.
 * @complexity O(1); see @ref t96_components.
 * @alloc none.
 * @test IrbemT96.ReferenceLaneMatchesTheComponentForm
 */
[[nodiscard]] inline FieldVector<Frame::GSM> t96_field_at(Position<Frame::GSM> p,
                                                          const T96Amplitudes<double>& a) {
    const std::array<double, 3> b = t96_components<double>(a, p.v[0], p.v[1], p.v[2]);
    return FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
}

/**
 * The T96 external field, with the model's own verdict on whether it should be believed here.
 *
 * The value is **always** returned, including when the status is @ref Status::OutOfValidityRange —
 * `status.hpp`'s standing rule. Refused outright, with a zero field, is only arithmetic that has no
 * answer: a non-finite input, a negative pressure (no square root), a tilt of exactly ±90 degrees
 * (the sheet warp carries `tan(psi)`), a point inside the Earth, or an extrapolation so far sunward
 * that a box harmonic's `exp(kappa x)` overflows and the assembly makes a NaN — which is refused
 * for T89's reason: a NaN that escapes here surfaces a hundred RK4 steps later with nothing left
 * pointing at its cause.
 *
 * Everything else is checked against the published envelope through `status.hpp`, so the rules
 * live in one place: the four drivers against `-100 <= Dst <= 20`, `0.5 <= Pdyn <= 10`,
 * `|By| <= 10`, `|Bz| <= 10`, and the position against `r_GEO <= 40 R_E`. An out-of-envelope
 * driver is still evaluated — the form is affine in Dst and homogeneous in the IMF, so what a
 * caller at Dst = -300 gets is exactly the published extrapolation, and is told so.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians; positive when the north dipole leans sunward.
 * @param dst_nt Dst, nT.
 * @param pdyn_npa solar-wind dynamic pressure, nPa.
 * @param by_nt IMF By, GSM, nT.
 * @param bz_nt IMF Bz, GSM, nT.
 * @return the field and its caveat. @ref Status::DomainError (with a zero field) for the refusals
 *         above; @ref Status::OutOfValidityRange for a driver or a radius outside the published
 *         envelope, with the field still computed; otherwise @ref Status::Ok.
 * @complexity O(1): one amplitude collapse (3 840 multiply-adds) and one evaluation.
 * @alloc none.
 * @test IrbemT96.OutOfRangeDriversAreReportedButStillEvaluated
 * @test IrbemT96.ValidityIsCheckedFromBothSidesOfEveryBound
 * @test IrbemT96.NonFiniteInputIsADomainError
 * @test IrbemT96.RightAngleTiltAndNegativePressureAreDomainErrors
 * @test IrbemT96.AnOverflowingExtrapolationIsADomainErrorNotANaN
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> t96_field(Position<Frame::GSM> p,
                                                               double tilt_rad, double dst_nt,
                                                               double pdyn_npa, double by_nt,
                                                               double bz_nt) {
    const FieldVector<Frame::GSM> zero{};
    if (!std::isfinite(p.v[0]) || !std::isfinite(p.v[1]) || !std::isfinite(p.v[2]) ||
        !std::isfinite(tilt_rad) || !std::isfinite(dst_nt) || !std::isfinite(pdyn_npa) ||
        !std::isfinite(by_nt) || !std::isfinite(bz_nt)) {
        return {Status::DomainError, zero};
    }
    if (!(std::fabs(tilt_rad) < max_tilt_rad)) return {Status::DomainError, zero};
    if (pdyn_npa < 0.0) return {Status::DomainError, zero};

    const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
    const Status where = check_position(ExternalModel::Tsyganenko1996, r, p.v[0]);
    if (where == Status::DomainError) return {Status::DomainError, zero};

    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Dst)] = dst_nt;
    drivers[static_cast<std::size_t>(Driver::Pdyn)] = pdyn_npa;
    drivers[static_cast<std::size_t>(Driver::ByIMF)] = by_nt;
    drivers[static_cast<std::size_t>(Driver::BzIMF)] = bz_nt;
    const Status drives = check_validity(ExternalModel::Tsyganenko1996, drivers);

    const FieldVector<Frame::GSM> b =
        t96_field_at(p, t96_amplitudes<double>(tilt_rad, dst_nt, pdyn_npa, by_nt, bz_nt));
    if (!std::isfinite(b.v[0]) || !std::isfinite(b.v[1]) || !std::isfinite(b.v[2])) {
        return {Status::DomainError, zero};
    }
    return {first_failure(drives, where), b};
}

/**
 * The T96 external field for a whole epoch's worth of state — the production entry point.
 *
 * Reads the tilt and the four drivers straight out of @ref HotState. The tilt ANGLE is read rather
 * than its sine and cosine because the `|psi| < pi/2` check is on the angle.
 *
 * @param p the position, GSM, in Earth radii.
 * @param ctx the epoch's context; its `hot()` block carries the tilt, Dst, Pdyn, By and Bz.
 * @return the field and its caveat, exactly as the six-argument overload.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT96.ContextOverloadAgreesWithTheExplicitOne
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> t96_field(Position<Frame::GSM> p,
                                                               const FieldContext& ctx) {
    const HotState& h = ctx.hot();
    return t96_field(p, h.tilt_rad, h.dst, h.pdyn, h.by_imf, h.bz_imf);
}

// -------------------------------------------------------------------------------------------
// The batch lanes
// -------------------------------------------------------------------------------------------

/**
 * The T96 field over a whole batch, on the CPU, in `float`.
 *
 * The host twin of `irbem_t96_f32`: the same basis, the same amplitudes rounded to `float` once,
 * the same accumulation order. What runs on a machine with no device, and what the device is held
 * to.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, GSM, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, nanotesla; overwritten in full.
 * @param a the amplitudes, in `float`, from `t96_amplitudes<float>`.
 * @return `false` when @p pos is not a whole number of points or @p out is a different length, in
 *         which case nothing is written; `true` otherwise.
 * @complexity O(N).
 * @alloc none.
 * @test IrbemT96.HostFloatLaneTracksTheReferenceLane
 * @test IrbemT96.HostFloatLaneRejectsMismatchedSpans
 */
[[nodiscard]] inline bool t96_field_host(std::span<const float> pos, std::span<float> out,
                                         const T96Amplitudes<float>& a) {
    if (pos.size() % 3 != 0 || out.size() != pos.size()) return false;
    const std::size_t n = pos.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const std::array<float, 3> b =
            t96_components<float>(a, pos[(3 * i) + 0], pos[(3 * i) + 1], pos[(3 * i) + 2]);
        out[(3 * i) + 0] = b[0];
        out[(3 * i) + 1] = b[1];
        out[(3 * i) + 2] = b[2];
    }
    return true;
}

/// How many `float` scalars the device kernel's parameter buffer holds: `sin(psi)`, `cos(psi)`,
/// then the 128 collapsed amplitudes. The geometry is compiled into the kernel, as T89's fixed
/// parameters are.
inline constexpr std::size_t t96_param_count = 2 + t96_basis_count;

/**
 * Pack the amplitudes into the kernel's parameter buffer.
 *
 * The layout is the kernel's ABI and is stated in exactly two places — here and above
 * `irbem_t96_f32` in `irbem.slang`; a test evaluates both lanes on the same points.
 *
 * @param a the amplitudes, in `float`.
 * @return the parameter block, `t96_param_count` floats, by value.
 * @complexity O(`t96_basis_count`).
 * @alloc none — the block is the returned object's own inline array.
 * @test IrbemT96.ParameterBlockCarriesTheTiltThenTheAmplitudes
 */
[[nodiscard]] inline std::array<float, t96_param_count> t96_param_block(
    const T96Amplitudes<float>& a) {
    std::array<float, t96_param_count> block{};
    block[0] = a.sin_tilt;
    block[1] = a.cos_tilt;
    for (std::size_t j = 0; j < t96_basis_count; ++j) block[2 + j] = a.c[j];
    return block;
}

/**
 * The batch's position caveat, folded one point at a time — @ref T89PositionFold's device, for
 * T96's envelope: a radius band with no sunward bound, so the smallest and largest radii decide
 * the whole batch, compared as squares.
 *
 * @test IrbemT96.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
struct T96PositionFold {
    /// The smallest `r^2` seen, R_E^2; `+inf` until the first point.
    double r2_lo = std::numeric_limits<double>::infinity();
    /// The largest `r^2` seen, R_E^2; zero until the first point.
    double r2_hi = 0.0;
    /// The smallest GSM `x` seen, R_E; `+inf` until the first point.
    double x_lo = std::numeric_limits<double>::infinity();
    /// False once any point has had a non-finite coordinate.
    bool finite = true;

    /**
     * Fold one position in.
     * @param p the position, GSM, in Earth radii.
     * @complexity O(1) — no `sqrt`, no branch that is not a select.
     * @alloc none.
     * @test IrbemT96.BatchReportsTheSameEnvelopeTheScalarLaneDoes
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
     *         @ref Status::OutOfValidityRange when any point is beyond `r_GEO = 40 R_E`,
     *         otherwise @ref Status::Ok.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemT96.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    [[nodiscard]] Status verdict() const {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double lo = finite ? std::sqrt(r2_lo) : nan;
        const double hi = finite ? std::sqrt(r2_hi) : nan;
        return first_failure(check_position(ExternalModel::Tsyganenko1996, lo, x_lo),
                             check_position(ExternalModel::Tsyganenko1996, hi, x_lo));
    }
};

/**
 * The T96 field over a whole batch of GSM points, on the device when that is worth it.
 *
 * The same contract as @ref t89_field_batch, driver for driver: one status for N points is the
 * worst of them; an out-of-validity batch is computed in full; a domain-error batch — a bad tilt or
 * pressure, a point inside the Earth or not finite — is zeroed in full; the value says whether the
 * DEVICE serviced the call. The amplitudes are collapsed ONCE per batch on the host: a batch
 * shares one epoch and one driver set, and 3 840 multiply-adds per thread would be 10^5 redundant
 * copies of a calculation done once.
 *
 * @param points the positions, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians.
 * @param dst_nt Dst, nT.
 * @param pdyn_npa solar-wind dynamic pressure, nPa.
 * @param by_nt IMF By, nT.
 * @param bz_nt IMF Bz, nT.
 * @param out receives one field vector per input, GSM, nanotesla; same length as @p points.
 * @return @ref Status::DomainError on a length mismatch, a non-finite or right-angle tilt, a
 *         non-finite or negative pressure, a non-finite driver, or a point that is not finite or
 *         is inside the Earth, and then every output is zeroed; @ref Status::OutOfValidityRange
 *         when a driver or any point's radius is outside the published envelope, with every point
 *         still computed; otherwise @ref Status::Ok. The value is `true` exactly when the device
 *         lane serviced the call.
 * @complexity O(N); on the device those N run concurrently over `ceil(N/256)` workgroups.
 * @alloc the device lane stages positions and results into two `std::vector<float>` of `3N`; the
 *        host lane allocates nothing.
 * @test IrbemT96.BatchAgreesWithTheReferenceLane
 * @test IrbemT96.BatchRejectsMismatchedSpans
 * @test IrbemT96.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
[[nodiscard]] inline Result<bool> t96_field_batch(std::span<const Position<Frame::GSM>> points,
                                                  double tilt_rad, double dst_nt, double pdyn_npa,
                                                  double by_nt, double bz_nt,
                                                  std::span<FieldVector<Frame::GSM>> out) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(tilt_rad) || !std::isfinite(dst_nt) || !std::isfinite(pdyn_npa) ||
        !std::isfinite(by_nt) || !std::isfinite(bz_nt)) {
        return {Status::DomainError, false};
    }
    if (!(std::fabs(tilt_rad) < max_tilt_rad) || pdyn_npa < 0.0) return {Status::DomainError, false};

    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Dst)] = dst_nt;
    drivers[static_cast<std::size_t>(Driver::Pdyn)] = pdyn_npa;
    drivers[static_cast<std::size_t>(Driver::ByIMF)] = by_nt;
    drivers[static_cast<std::size_t>(Driver::BzIMF)] = bz_nt;
    const Status drives = check_validity(ExternalModel::Tsyganenko1996, drivers);
    if (n == 0) return {drives, false};

    T96PositionFold fold;

#if CHEATAH_SPACE_IRBEM_T96_GPU
    if (gpu::prefer_gpu("irbem_t96_f32", n) &&
        std::filesystem::exists(gpu::shader_path("irbem_t96_f32"))) {
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
        const std::array<float, t96_param_count> block = t96_param_block(
            t96_amplitudes<float>(tilt_rad, dst_nt, pdyn_npa, by_nt, bz_nt));
        gpu::dispatch_batch("irbem_t96_f32", pos, raw, std::span<const float>(block));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = FieldVector<Frame::GSM>{
                fixarray::vec3d{raw[(3 * i) + 0], raw[(3 * i) + 1], raw[(3 * i) + 2]}};
        }
        return {first_failure(drives, where), true};
    }
#endif

    const T96Amplitudes<double> a = t96_amplitudes<double>(tilt_rad, dst_nt, pdyn_npa, by_nt, bz_nt);
    for (std::size_t i = 0; i < n; ++i) {
        fold.add(points[i]);
        out[i] = t96_field_at(points[i], a);
    }
    const Status where = fold.verdict();
    if (where == Status::DomainError) {
        for (std::size_t i = 0; i < n; ++i) out[i] = FieldVector<Frame::GSM>{};
        return {Status::DomainError, false};
    }
    return {first_failure(drives, where), false};
}

// -------------------------------------------------------------------------------------------
// The total field — IGRF plus T96, as one GeoFieldModel
// -------------------------------------------------------------------------------------------

/**
 * IGRF plus Tsyganenko 1996, as a single field a tracer can follow.
 *
 * The same shape as `TotalFieldT89`, for the same reason: `trace_invariant`, `make_lstar` and
 * everything above them take one @ref GeoFieldModel and follow the TOTAL field without knowing
 * there are two. The amplitudes are collapsed ONCE, at construction — they are a property of the
 * epoch and the drivers, and a trace evaluates the field ~10^5 times — so a field evaluation is
 * one basis evaluation and one dot product. The external model's refusals fall back to the
 * internal field alone (the best available answer there, and never a NaN into an integrator), and
 * @ref external_status says when that happened.
 *
 * Traces through this field run the host lane: the device tracer is T89's, and a T96 tracer is a
 * second kernel this file does not own.
 *
 * @tparam NMAX the internal field's truncation degree; 10 reproduces IRBEM's own choice.
 * @test IrbemT96.TotalFieldSuperposesInternalAndExternal
 * @test IrbemT96.TotalFieldTracesAndReportsWhenTheExternalModelDeclines
 */
template <int NMAX = 10>
class TotalFieldT96 {
  public:
    /// The internal part's truncation degree — what generic staging and buffer sizing read.
    static constexpr int degree = NMAX;

    /**
     * @param internal the internal field, already built for the epoch.
     * @param rotations the epoch's frame rotations — built once, reused for every point.
     * @param dst_nt Dst, nT.
     * @param pdyn_npa solar-wind dynamic pressure, nPa.
     * @param by_nt IMF By, GSM, nT.
     * @param bz_nt IMF Bz, GSM, nT.
     */
    TotalFieldT96(const Igrf<NMAX>& internal, const Rotations& rotations, double dst_nt,
                  double pdyn_npa, double by_nt, double bz_nt)
        : internal_(&internal),
          rotations_(&rotations),
          tilt_rad_(rotations.dipole_tilt_deg * (std::numbers::pi / 180.0)),
          dst_nt_(dst_nt),
          pdyn_npa_(pdyn_npa),
          by_nt_(by_nt),
          bz_nt_(bz_nt),
          amplitudes_(t96_amplitudes<double>(tilt_rad_, dst_nt, pdyn_npa >= 0.0 ? pdyn_npa : 0.0,
                                             by_nt, bz_nt)) {}

    /**
     * The total field at a geographic point.
     * @param p the position, GEO, Earth radii.
     * @return `B_internal + B_external` in GEO, nT; the internal field alone when the external
     *         model refuses the point (see @ref external_status).
     * @complexity One IGRF evaluation, one T96 evaluation, two 3x3 rotations.
     * @alloc none.
     * @test IrbemT96.TotalFieldSuperposesInternalAndExternal
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::GEO>& p) const {
        const FieldVector<Frame::GEO> b_int = internal_->evaluate(p);
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        if (external_status(p) == Status::DomainError) return b_int;
        const FieldVector<Frame::GSM> b_ext = t96_field_at(p_gsm, amplitudes_);
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
     * @complexity One radius, two envelope checks.
     * @alloc none.
     * @test IrbemT96.TotalFieldTracesAndReportsWhenTheExternalModelDeclines
     */
    [[nodiscard]] Status external_status(const Position<Frame::GEO>& p) const {
        const Position<Frame::GSM> g = transform<Frame::GSM>(p, *rotations_);
        if (!std::isfinite(g.v[0]) || !std::isfinite(g.v[1]) || !std::isfinite(g.v[2])) {
            return Status::DomainError;
        }
        if (!(std::fabs(tilt_rad_) < max_tilt_rad) || pdyn_npa_ < 0.0) return Status::DomainError;
        const double r = std::sqrt((g.v[0] * g.v[0]) + (g.v[1] * g.v[1]) + (g.v[2] * g.v[2]));
        const Status where = check_position(ExternalModel::Tsyganenko1996, r, g.v[0]);
        if (where == Status::DomainError) return where;
        DriverSet drivers{};
        drivers[static_cast<std::size_t>(Driver::Dst)] = dst_nt_;
        drivers[static_cast<std::size_t>(Driver::Pdyn)] = pdyn_npa_;
        drivers[static_cast<std::size_t>(Driver::ByIMF)] = by_nt_;
        drivers[static_cast<std::size_t>(Driver::BzIMF)] = bz_nt_;
        return first_failure(check_validity(ExternalModel::Tsyganenko1996, drivers), where);
    }

    /// The epoch's frame rotations.
    /// @return the rotations this field was built with.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT96.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Rotations& rotations() const { return *rotations_; }

    /// The collapsed amplitudes this field evaluates with.
    /// @return the amplitudes, computed at construction.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT96.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const T96Amplitudes<double>& amplitudes() const { return amplitudes_; }

    /// The internal part's Gauss coefficient `g(n, m)`, in nT — forwarded, for `TotalFieldT89`'s
    /// reason: the superposition has no expansion of its own and every caller means the internal one.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT96.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double g(int n, int m) const { return internal_->g(n, m); }

    /// The internal part's `h(n, m)`, in nT — see @ref g.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT96.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double h(int n, int m) const { return internal_->h(n, m); }

    /// The internal field alone.
    /// @return the internal model.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT96.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Igrf<NMAX>& internal() const { return *internal_; }

  private:
    const Igrf<NMAX>* internal_;
    const Rotations* rotations_;
    double tilt_rad_;
    double dst_nt_;
    double pdyn_npa_;
    double by_nt_;
    double bz_nt_;
    T96Amplitudes<double> amplitudes_;
};

}  // namespace cheatah::space::irbem
