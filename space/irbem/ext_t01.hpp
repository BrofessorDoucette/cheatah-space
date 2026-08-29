#pragma once

/**
 * @file ext_t01.hpp
 * @brief space.irbem — Tsyganenko (2002), "T01": the near-magnetosphere model with a dawn-dusk
 *        asymmetry, driven CONTINUOUSLY by the solar wind — as published, with the unpublished
 *        parts named and measured rather than invented.
 *
 * T89 (`ext_t89.hpp`) is a closed-form fit sorted into Kp bins. T01 is a different kind of object:
 * every current system is built by taking a simple analytic field, DEFORMING space around it so its
 * geometry matches observations, and letting the solar wind drive its amplitude and shape
 * continuously — the dynamic pressure `P_dyn`, the corrected `Dst*`, and two history integrals of the
 * interplanetary field, `G1` and `G2`, which are hour-long averages of the coupling functions
 * defined in `tests/irbem_domain_corpus.hpp` and IRBEM's `maginput` table. Nothing here is binned:
 * nearby drivers give nearby fields, and the suite asserts that smoothness rather than the bin
 * identity T89's suite asserts.
 *
 * ## The papers, and what each equation number below refers to
 *
 *  - **Paper 1**: Tsyganenko, *A model of the near magnetosphere with a dawn-dusk asymmetry.
 *    1. Mathematical structure*, J. Geophys. Res. **107**(A8), 1179, doi:10.1029/2001JA000219
 *    (2002). Equation numbers written "P1 eq. (n)".
 *  - **Paper 2**: Tsyganenko, *... 2. Parameterization and fitting to observations*, J. Geophys.
 *    Res. **107**(A8), 1176, doi:10.1029/2001JA000220 (2002). Written "P2 eq. (n)".
 *
 * P1 eq. (1) is the whole plan: `B_E = B_CF + B_RC + B_T + B_FAC + B_INT` — the Chapman-Ferraro
 * field, the ring current, the cross-tail current, the field-aligned (Birkeland) currents and the
 * IMF interconnection field. P2 eq. (11) is the same sum with every amplitude written out as a
 * function of the drivers, and P2 Table 1 gives the 24 linear coefficients and 18 non-linear
 * parameters of that expression.
 *
 * ## What is implemented here, and how each piece is built
 *
 * **The cross-tail current, `B_T`** (P1 §2.2). One vector potential, `A = A_phi e_phi` coaxial
 * with `z_GSM`, P1 eqs. (2)-(5): a thick current disc whose profile is a five-term sum over
 * Table 1's `f_i, b_i, c_i`, with the half-thickness `D(X, Y)` of eq. (6). Two MODULES are cut
 * from it by the affine substitution of P1 ¶16 — `{X, Y, Z} -> {eta X - (eta - 1) X_m, eta Y,
 * eta Z}`, `D -> eta D`, plus a sunward shift — with different scale factors `eta`, so that one
 * ("short", `eta = 1.1`) carries the near-tail current that falls off quickly and the other
 * ("long", `eta = 0.25`) the current that persists down-tail. Their amplitudes `t_1, t_2` are P2
 * eq. (2). The dipole tilt enters by two DEFORMATIONS of the untilted field, P1 §2.2.2: a warping of
 * the sheet across the tail (eqs. 7-10, the `phi -> F` map in the `y z` plane) and a radially
 * dependent rotation about `y` (eqs. 11-14), which hinges the sheet from the dipole equator near the
 * Earth to the GSM equator far down-tail.
 *
 * **The Birkeland currents, `B_FAC`** (P1 §2.3). Each of Region 1 and Region 2 is a conical current
 * sheet of P1 eqs. (16)-(17), Tsyganenko (1991)'s analytic field, with two longitudinal Fourier modes
 * `m = 1, 2`. That cone is DEFORMED in `(r, theta)` by eqs. (18)-(19) with the 30 coefficients per
 * mode of Table 2, which is what turns radial cones into currents that follow dipolar field lines;
 * the southern system is the northern one rotated by 180 degrees about `x` with its polarity
 * reversed, eq. (20); the noon-midnight asymmetry and the tilt come from the azimuthal deformation
 * of eq. (21) (its constants `Delta phi`, `b`, `rho_0`, `beta`, `R_H`, `epsilon` are given in
 * ¶47-48); and the whole system is expanded or contracted by the scaling of eq. (25) with the
 * `G2`-dependent factors of P2 eq. (9). Amplitudes are P2 eq. (8). Each family is multiplied by ONE
 * unit constant the papers omit, measured against the oracle and documented at
 * @ref T01UnitConstants.
 *
 * **The interconnection field, `B_INT`** (P2 §3.5): a uniform field along the transverse IMF,
 * `epsilon(theta) B_perp`, with the penetration efficiency of P2 eq. (10) rising from 0.07 for
 * northward IMF to 0.62 for southward. This is the ONE term whose whole structure can be read
 * straight off the oracle, because `B_y` and `B_z` drive nothing else in P2 eq. (11): the
 * differential's pass 3 shows the oracle's response to them is a single UNIFORM vector to 8e-14 nT
 * over scattered points, and that it is exactly linear in `sin^2(theta / 2)` — the published form,
 * confirmed to eleven digits. Its two coefficients are a documented disagreement; see
 * @ref T01Coefficients.
 *
 * **Every deformation is ONE operation**, @ref detail::t01::push_forward: if `T` maps a point `r`
 * to `r*` with Jacobian `J = d r* / d r`, the field `B'(r) = adj(J) B*(T r)` — the adjugate of the
 * Jacobian applied to the undeformed field at the image point — is divergence-free whenever `B*`
 * is, because it is the pull-back of the flux two-form. P1 eqs. (8)-(10), (14) and (22)-(24) are
 * that identity written out in cylindrical and Cartesian coordinates; here it is written once and
 * every Jacobian is an analytic derivative of its map. That is also why `div B = 0` is the test
 * that carries this file: a wrong partial derivative anywhere in six Jacobians breaks it, and
 * nothing else does.
 *
 * ## What is NOT implemented, and why — the provenance verdict, MEASURED
 *
 * Three pieces of P1 eq. (1) are absent, for three different reasons, and the differential harness
 * `tools/oracle/t01_diff.cpp` measures what each absence costs against IRBEM's `kext = 9`:
 *
 *  1. **Every shielding field, including all of `B_CF`.** P1 ¶13 and ¶29 say that each module
 *     "is provided with its own shielding field" — a curl-free expansion whose coefficients were
 *     obtained by least squares against the model magnetopause — and `B_CF` is nothing but the
 *     shielding field of the Earth's dipole. Those coefficient arrays (hundreds of numbers) appear
 *     in neither paper; they exist only in Tsyganenko's distributed source, which is GPL-3.0 and
 *     which this MIT clean room does not read. The unshielded modules are the published physics;
 *     the shielding is what confines it, and its absence is largest near the magnetopause and in
 *     the dayside `B_z`.
 *  2. **The ring current, `B_RC` (symmetric and partial).** P1 §2.1 defers its analytic form to
 *     Tsyganenko, J. Geophys. Res. **105**, 27739 (2000), doi:10.1029/2000JA000138 — a paper that is
 *     published but closed-access with no open copy (Unpaywall, Semantic Scholar and NTRS were all
 *     queried and returned none), so its equations could not be read in this room. P2 eqs. (5)-(7)
 *     give the amplitudes, scaling factors and rotation angle that would drive it; they multiply a
 *     basis this file does not have and are therefore not carried.
 *  3. **P1 §2.4-2.5** (the magnetopause and the shielding representation) are the pages after
 *     the author's own PDF ends; the interconnection term's final form is nevertheless fully stated
 *     in P2 §3.5 and is implemented from there.
 *
 * Seven numeric facts the published text uses without stating were resolved from the vendored
 * IRBEM `source/Tsyganenko01.f` DATA statements, as rule 1 permits for a single number, and are
 * cited at their definitions in @ref t01_tail_fixed and @ref t01_tail_modules: the two modules'
 * `eta`, `X_m`, sunward shift and `Delta D_x`, the eq. (6) scales `Delta X = 7`, `Delta Y = 20`,
 * the eq. (7) warp scale `L = 20`, and that the fitted inner-edge shift `Delta X_0 + Delta X_1 G2`
 * attaches to the short module alone.
 *
 * Two kinds of thing the papers leave open were settled by a BLACK-BOX comparison with the oracle
 * (never by reading its logic), and both are stated as measurements at
 * @ref t01_tail_modules and @ref T01UnitConstants: the placement of the sunward shift relative to the
 * `eta` scaling, and the two unit constants — one per current family — that turn the printed
 * potentials into nanotesla, which neither paper states (@ref T01UnitConstants).
 *
 * The numbers those measurements produced, and the residual the published form leaves against the
 * oracle with and without a free re-fit of its amplitudes, are in the file's differential section
 * below the API and in `tools/oracle/t01_diff.cpp`. The verdict is stated once, here:
 * **PUBLISHED-FORM-WITH-DOCUMENTED-GAP** — the form and the amplitude coefficients are published
 * and implemented; the shielding coefficients are not published and the ring-current form is
 * unobtainable, and the differential says how much field that is.
 *
 * @note Nothing on a hot path allocates. The evaluator touches its stack and the caller's spans;
 *       the device lane of @ref t01_field_batch stages `3N` floats each way and forwards whatever
 *       `gpu::dispatch_batch` throws — its `@alloc` says so, as `ext_t89.hpp`'s does.
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

// The device lane is opt-in by include path, exactly as ext_t89.hpp's is.
#if __has_include("cheatah_gpu_linalg/context.hpp")
#  include <filesystem>
#  include "gpu/dispatch.hpp"
/// 1 when this translation unit can reach the T01 device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_T01_GPU 1
#else
/// 1 when this translation unit can reach the T01 device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_T01_GPU 0
#endif

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// The published structure: the tail disc, its two modules, the Birkeland cones
// -------------------------------------------------------------------------------------------

/// How many terms the cross-tail disc potential of P1 eq. (2) sums: `N = 5` (P1 ¶14).
inline constexpr std::size_t t01_tail_term_count = 5;

/**
 * P1 Table 1 — the cross-tail current disc, `f_i`, `b_i`, `c_i` of eqs. (2)-(5).
 *
 * Fitted once by Tsyganenko to give the current profile of P1 Fig. 1 (a steep inner edge at
 * `~10 R_E` and a slow tailward decay) and FIXED thereafter; the two tail modules are cut from this
 * one disc by coordinate scaling, not by re-fitting it.
 *
 * @test IrbemT01.PublishedTablesAreTranscribed
 */
struct T01TailDisc {
    /// `f_i`, the term amplitudes.
    std::array<double, t01_tail_term_count> f;
    /// `b_i`, R_E — the radial scale of each term; the disc's inner edge sits near the smallest.
    std::array<double, t01_tail_term_count> b;
    /// `c_i`, R_E — each term's vertical offset, added to the broadened `zeta`.
    std::array<double, t01_tail_term_count> c;
};

/// P1 Table 1, transcribed.
inline constexpr T01TailDisc t01_tail_disc{
    {-71.093466, -1014.3086, -1272.9394, -3224.9359, -44546.862},
    {10.901012, 12.683939, 13.517920, 14.867750, 15.123064},
    {0.79540700, 0.67166018, 1.17486632, 2.56524992, 10.0198679}};

/**
 * The tail parameters P1 states or uses as fixed constants — the eq. (6) thickness scales, the
 * eq. (7) warp scale and the eq. (13) hinge coefficient.
 *
 * `delta_x`, `delta_y` and `l` are used by P1 eqs. (6)-(7) without a stated value ("`Delta Y`
 * enters as a factor by `Delta D_y` and hence can be kept fixed", ¶15; "`G` and `L` were considered
 * as constant parameters", ¶20, but P2 Table 1 fits only `G`). They were resolved as numeric facts
 * from the vendored IRBEM `source/Tsyganenko01.f`, lines 802 (`XL = 20`), 984 (`exp(X / 7)`) and
 * 985 (`(Y / 20)^2`), per rule 1. `rh2` is stated in P1 ¶23.
 *
 * @test IrbemT01.FixedParametersAreThePublishedOnes
 */
struct T01TailFixed {
    /// `Delta X`, R_E — the sunward e-folding scale of the eq. (6) dayside flaring.
    double delta_x;
    /// `Delta Y`, R_E — the dawn-dusk scale of the eq. (6) flank thickening.
    double delta_y;
    /// `L`, R_E — the scale of the eq. (7) warping across the tail.
    double l;
    /// `R_H2`, R_E — the eq. (13) coefficient making the hinging distance depend on `Z^2 / r^2`;
    /// "-5.2 R_E was specified to provide the observed amplitude of the tilt-related shift of the
    /// magnetotail boundary at `X ~ -30 R_E`" (P1 ¶23).
    double rh2;
};

/// The fixed tail parameters; see @ref T01TailFixed for where each comes from.
inline constexpr T01TailFixed t01_tail_fixed{7.0, 20.0, 20.0, -5.2};

/**
 * One tail module's coordinate substitution, P1 ¶16: `X' = eta (X - shift) - (eta - 1) X_m`,
 * `Y' = eta Y`, `Z' = eta Z`, `D -> eta D`, plus its own dayside-flaring amplitude in eq. (6).
 *
 * P1 ¶17 and its Fig. 1 caption state `eta = 1.1, X_m = -12, X_s = 6` for the short module and
 * `eta = 0.25, X_m = -12, X_s = 4` for the long one, presenting them as the plotted values; that
 * they are also the model's, and that `Delta D_x` is 1 for the short module and 0 for the long one,
 * was resolved from `source/Tsyganenko01.f` lines 875-876 and 912 (DATA statements) per rule 1.
 *
 * **Where the shift sits is AMBIGUOUS, and the measurement does not settle it.** The paper says
 * "one can introduce a variable shift `X_s` of the current sheet along the Sun-Earth line" without
 * saying whether it is applied before or after the `eta` scaling (the two differ by a factor `eta`
 * in the effective displacement: 0.6 R_E for the short module, 3 R_E for the long one).
 * `tools/oracle/t01_diff.cpp`'s pass 4 fits both readings to the oracle's response to `G1` — which
 * isolates the tail modules, because `G1` drives nothing else in P2 eq. (11) — and they come back
 * at residuals of **0.5226** and **0.5170**, six parts in a thousand apart on a residual of one
 * half. That is not a discrimination, and this file does not pretend it is one: the reading
 * implemented is the one the sentence reads most naturally in, the shift is in PHYSICAL
 * coordinates and the scaling is applied to the shifted point, and the alternative is left
 * recorded here rather than silently chosen. It matters little: the difference between the two
 * fields is a fraction of the gap the missing shielding already leaves (see @ref
 * T01UnitConstants).
 *
 * @test IrbemT01.FixedParametersAreThePublishedOnes
 */
struct T01TailModule {
    /// `eta`, the scale factor about the fixed point `X_m`; > 1 compresses the disc sunward.
    double eta;
    /// `X_m`, R_E — the fixed point of the scaling: the current-density peak stays there.
    double x_m;
    /// `X_s`, R_E — the sunward shift of the whole module, in PHYSICAL coordinates before scaling.
    double x_s;
    /// `Delta D_x`, R_E — the eq. (6) dayside-flaring amplitude for this module.
    double delta_d_x;
};

/**
 * The two tail modules, P1 ¶17: index 0 the "short" near-tail module, 1 the "long" one.
 *
 * The shift is applied in physical coordinates, before the scaling; the harness's `G1` probe does
 * NOT distinguish that from the other reading, and @ref T01TailModule records both residuals.
 *
 * @test IrbemT01.FixedParametersAreThePublishedOnes
 */
inline constexpr std::array<T01TailModule, 2> t01_tail_modules{{{1.1, -12.0, 6.0, 1.0},
                                                                 {0.25, -12.0, 4.0, 0.0}}};

/**
 * One Birkeland mode's published parameters — P1 Table 2, one column: the conical amplitude
 * `B_m`, the cone's apex colatitude `theta_0` (eq. 15, Fig. 3), and the 29 coefficients of the
 * `(r, theta)` deformation, eqs. (18)-(19).
 *
 * @test IrbemT01.PublishedTablesAreTranscribed
 */
struct T01Cone {
    /// Which Fourier mode, `m = 1` (dawn-dusk antisymmetric) or `m = 2`; the exponent in eq. (17).
    int m;
    /// `B_m`, the mode's overall magnitude coefficient, P1 Table 2.
    double b_m;
    /// `theta_0`, radians — the colatitude of the cone's sheet at the ionosphere.
    double theta_0;
    /// `Delta theta`, radians — half the sheet's angular thickness (P1 ¶36: 3.5 degrees for
    /// Region 1 and 5.2 for Region 2 in total).
    double delta_theta;
    /// `a_1 .. a_9` of eq. (18).
    std::array<double, 9> a;
    /// `b_1 .. b_6` of eq. (18), R_E.
    std::array<double, 6> b;
    /// `c_1 .. c_10` of eq. (19).
    std::array<double, 10> c;
    /// `d_1 .. d_4` of eq. (19), R_E.
    std::array<double, 4> d;
};

/// The Region 1 and Region 2 cones' half-thickness, radians: half of P1 ¶36's 3.5 and 5.2 degrees.
inline constexpr double t01_r1_half_thickness = 1.75 * std::numbers::pi / 180.0;
/// See @ref t01_r1_half_thickness.
inline constexpr double t01_r2_half_thickness = 2.6 * std::numbers::pi / 180.0;

/**
 * P1 Table 2, transcribed: index 0 Region 1 `m = 1`, 1 Region 1 `m = 2`, 2 Region 2 `m = 1`,
 * 3 Region 2 `m = 2`.
 *
 * @test IrbemT01.PublishedTablesAreTranscribed
 */
inline constexpr std::array<T01Cone, 4> t01_cones{{
    // Region 1, m = 1
    {1, 0.1618068, 0.7113545, t01_r1_half_thickness,
     {-0.1797958, 2.999643, -0.9322709, -0.6811060, 0.2099057, -8.358816, -14.86034, 0.3838363,
      -16.30946},
     {4.537023, 2.685836, 27.97833, 6.330871, 1.876532, 18.95619},
     {0.9651528, 0.4217195, -0.0895777, -1.823556, 0.7457045, -0.5785917, -1.010201, 0.0111239,
      0.0957293, -0.3599292},
     {8.713701, 0.9763933, 3.834603, 2.492118}},
    // Region 1, m = 2
    {2, 0.7058027, 0.5567714, t01_r1_half_thickness,
     {-0.2845939, 5.715471, -2.472821, -0.7738802, 0.3478294, -11.37654, -38.64769, 0.6932928,
      -212.4017},
     {4.944205, 3.071270, 33.05882, 7.387534, 2.366769, 79.22573},
     {0.6154290, 0.5592051, -0.1796585, -1.654932, 0.7309109, -0.4926293, -1.130266, -0.00961398,
      0.1484586, -0.2215347},
     {7.883593, 0.02768252, 2.950281, 1.212635}},
    // Region 2, m = 1
    {1, 0.1278764, 0.8867880, t01_r2_half_thickness,
     {-0.2320034, 1.805623, -32.37241, -0.9931491, 0.3175086, -2.492466, -16.21600, 0.2695393,
      -6.752691},
     {3.971795, 14.54478, 41.10158, 7.912890, 1.258297, 9.583548},
     {1.014142, 0.5104135, -0.1790431, -1.756358, 0.7561987, -0.6775248, -0.04014016, 0.0144680,
      0.1200522, -0.2203585},
     {4.508964, 0.8221624, 1.779934, 1.102650}},
    // Region 2, m = 2
    {2, 0.4036015, 0.7247997, t01_r2_half_thickness,
     {-0.3302974, 2.827731, -45.44406, -1.611104, 0.4927112, -0.00325846, -49.59015, 0.3796217,
      -233.7884},
     {4.312667, 18.05052, 28.95320, 11.09948, 0.7471650, 67.10246},
     {0.5667097, 0.6468520, -0.1560666, -1.460805, 0.7719654, -0.6658989, 0.251518e-5, 0.0242602,
      0.1195003, -0.2625739},
     {4.377173, 0.2421191, 2.503483, 1.071587}},
}};

/**
 * The azimuthal deformation of the Birkeland field, P1 eq. (21), and its constants (¶47-48).
 *
 * `F = phi - (Delta phi + b (rho^2 - 1) / (rho^2 + rho_0^2)) sin phi - beta Psi / [1 + (r / R_H)^
 * epsilon]^(1 / epsilon)`, in cylindrical coordinates about `y_GSM` with `phi = -arctan(z / x)`.
 * The first term shifts the ovals to the nightside, the second drags them with the dipole tilt.
 *
 * @test IrbemT01.FixedParametersAreThePublishedOnes
 */
struct T01FacDeformation {
    /// `Delta phi`, radians, per region: "0.055 and 0.030 for the Region 1 and 2 currents".
    std::array<double, 2> delta_phi;
    /// `b`, the amplitude of the distance-dependent nightside shift; 0.5.
    double b;
    /// `rho_0`, R_E — where that shift saturates; 7.0.
    double rho_0;
    /// `beta`, the "slippage" between the Earth and the oval under tilt; 0.9.
    double beta;
    /// `R_H`, R_E — the hinging distance of the tilt dragging; 10.
    double r_h;
    /// `epsilon`, the radial sharpness of the dragging profile; 3.
    double epsilon;
};

/// P1 ¶47-48, as stated.
inline constexpr T01FacDeformation t01_fac_deformation{{0.055, 0.030}, 0.5, 7.0, 0.9, 10.0, 3.0};

/**
 * The two unit constants the papers omit — MEASURED against the oracle, and stated as such.
 *
 * Both papers print their building blocks without the constant that turns them into nanotesla.
 * P1 ¶35 normalizes the Birkeland currents to "10^6 A per half-wave of each sinusoidal mode" and
 * Table 2 gives `B_m` of order 0.1-0.7, but eqs. (16)-(17) as printed have no unit factor, and with
 * `B_m` alone the Region 1 field would be ~0.1 nT where Table 2's own `<B>` column says ~13 nT. The
 * tail disc of eqs. (2)-(5) with Table 1's `f_i` gives ~-30 nT per unit amplitude at the Earth for
 * the short module, while P2 ¶50's worked example puts the WHOLE shielded tail at -22.8 nT there
 * with `t_1 ~ 4.7`. So each family carries one scalar that P2's fitted amplitudes were expressed
 * against and that the text does not state.
 *
 * Each is measured ONCE, as a least-squares scale of this file's basis against the oracle's
 * response to the driver that isolates that family in P2 eq. (11) — `G1` drives only the two tail
 * amplitudes (and the oracle's response to it is linear to 2e-14, which is pass 2's check that the
 * isolation is exact), `G2` drives the four Birkeland amplitudes, their two scaling factors and the
 * short module's inner-edge shift — over the 999 points of `tools/oracle/t01_diff.cpp`'s pass 4,
 * scattered at `3 <= r <= 15 R_E`, `x >= -14 R_E`, over four epochs. Its printed output IS these
 * numbers:
 *
 *  - **tail**: fitting `a_1 M_1 + a_2 M_2` to `dB/dG1` gives `a_1 = 0.023797` against the published
 *    `t_1^(2) = 0.319`, hence `0.07460` per unit amplitude, with a **33.7%** residual. The short
 *    module alone leaves 34.2%, the long module alone 59.3%: the near-tail module carries the
 *    response and the long one adds almost nothing to it. Worse, the long module's fitted weight
 *    (`a_2 = +0.000466`) is ten times smaller than the published ratio `t_2^(2) / t_1^(2) = -0.191`
 *    would make it, AND of the opposite sign — its UNSHIELDED field, a current ring stretched to
 *    four times the disc, is simply not a proxy for the confined field the oracle carries. The one
 *    constant is applied to both modules, and the differential section says what that costs.
 *  - **Birkeland**: fitting the two regions with the published amplitude ratios of P2 eq. (8) to
 *    `dB/dG2` gives a common scale of `629.0` with a **35.6%** residual; freeing the two regions
 *    independently gives 629.9 and 637.0 and lowers the residual by 0.0001. That is the strongest
 *    statement this file can make about P2 eq. (8): the published amplitude RATIOS are what the
 *    oracle carries, and only the unit is unpublished.
 *
 * A single scalar per family cannot repair a wrong spatial pattern: the pattern agreement in those
 * fits is what validates the implementation, and the scalar is the unit. Two CONTROLS in the same
 * pass say the fit is measuring the model rather than flattering it — switching the tilt
 * deformations off raises the tail residual from 33.7% to 52.3%, so P1 eqs. (7)-(14) are doing
 * real work; and the two readings of P1's sunward shift differ by 0.6% of a 52% residual, so that
 * ambiguity is NOT resolved by this measurement (see @ref T01TailModule).
 *
 * @test IrbemT01.UnitConstantsAreTheMeasuredOnes
 */
struct T01UnitConstants {
    /// Nanotesla per unit tail-module amplitude, applied to `t_1 B_T1 + t_2 B_T2`.
    double tail;
    /// Nanotesla per unit conical amplitude, applied to every `b_l^(m) B_Rl^(m)`.
    double fac;
};

/// The measured unit constants; see @ref T01UnitConstants for the measurement.
inline constexpr T01UnitConstants t01_units{0.07460, 629.0};

// -------------------------------------------------------------------------------------------
// The published parameterization — P2 Table 1
// -------------------------------------------------------------------------------------------

/**
 * P2 Table 1: the 24 linear coefficients and the non-linear parameters this file reads.
 *
 * The ring-current entries (`s`, `p`, `xi_S0`, `beta_S`, `xi_P0`, `beta_P`, `Dst_0*`) and the
 * magnetopause compression index `kappa` are published in the same table and NOT carried here,
 * because they multiply bases this file does not have (see the file brief). Carrying numbers that
 * nothing reads would be a claim of completeness the code cannot back.
 *
 * @test IrbemT01.PublishedTablesAreTranscribed
 */
struct T01Coefficients {
    /// `t_1^(0..3)`: the short tail module's amplitude, P2 eq. (2).
    std::array<double, 4> t1;
    /// `t_2^(0..3)`: the long tail module's, P2 eq. (2).
    std::array<double, 4> t2;
    /// `b_1^(10), b_1^(11), b_1^(20), b_1^(21)`: Region 1, modes 1 and 2, P2 eq. (8).
    std::array<double, 4> b1;
    /// `b_2^(10), b_2^(11), b_2^(20), b_2^(21)`: Region 2, likewise.
    std::array<double, 4> b2;
    /// `epsilon_0`, `epsilon_1` of P2 eq. (10) — printed as `epsilon_1`, `epsilon_2` in the table.
    std::array<double, 2> eps;
    /// `alpha_1`, `alpha_2`: the pressure exponents of P2 eq. (2).
    double alpha1;
    double alpha2;
    /// `Delta X_0`, `Delta X_1`: the inner-edge shift of P2 eq. (3), R_E and R_E per unit `G2`.
    double dx0;
    double dx1;
    /// `D_0`, R_E — the central half-thickness of the tail sheet, P1 eq. (6).
    double d0;
    /// `R_H`, R_E — the tail hinging distance `R_H0` of P1 eq. (13).
    double rh0;
    /// `G`, R_E — the tail warping amplitude of P1 eq. (7).
    double g;
    /// `Delta D_y`, R_E — the flank thickening of P1 eq. (6).
    double delta_d_y;
    /// `varsigma_1^(0), varsigma_1^(1)`: the Region 1 scaling of P2 eq. (9).
    std::array<double, 2> sigma1;
    /// `varsigma_2^(0), varsigma_2^(1)`: Region 2, likewise.
    std::array<double, 2> sigma2;
    /// `P_d0`, nPa — the reference pressure of P2 eq. (2), "assumed equal to 2 nPa" (P2 ¶30).
    double pdyn0;
};

/// P2 Table 1, transcribed, with `P_d0 = 2` from P2 ¶30.
inline constexpr T01Coefficients t01_coefficients{
    {2.483, 0.583, 0.319, -0.088},   {-1.173, 3.575, -0.061, -0.011},
    {0.281, 0.166, -0.029, 0.026},   {-0.236, -0.077, 0.091, -0.025},
    {0.068, 0.554},
    1.244, 0.380, 0.689, -0.046, 2.36, 8.94, 28.3, 3.900,
    {1.13, 0.014}, {1.03, 0.030}, 2.0};

// -------------------------------------------------------------------------------------------
// Drivers, and the per-epoch state they resolve to
// -------------------------------------------------------------------------------------------

/**
 * The six drivers T01 reads, in IRBEM's `maginput` units.
 *
 * `G1` and `G2` are the hour averages defined in IRBEM's `maginput` table and P2 eqs. (1) and (4):
 * `G1 = <V h(B_perp) sin^3(theta / 2)>` with `h = (B_perp / 40)^2 / (1 + B_perp / 40)`, and
 * `G2 = 0.005 <V B_s>` with `B_s` the southward IMF magnitude. They are integrals over history, so
 * they are inputs here, never derived from the instantaneous `B_y, B_z`.
 *
 * @test IrbemT01.NearbyDriversGiveNearbyFields
 */
struct T01Drivers {
    /// `Dst`, nT. The corrected `Dst* = 0.8 Dst - 13 sqrt(P_dyn)` (P2 ¶30) is formed internally.
    double dst;
    /// `P_dyn`, nPa; must be positive, since it enters as `(P_dyn / P_d0)^alpha` and `sqrt(P_dyn)`.
    double pdyn;
    /// IMF `B_y`, GSM, nT.
    double by_imf;
    /// IMF `B_z`, GSM, nT.
    double bz_imf;
    /// `G1`.
    double g1;
    /// `G2`.
    double g2;
};

/**
 * Everything the evaluator needs that is a property of the EPOCH and the DRIVERS rather than of the
 * point: the tilt, the module amplitudes of P2 eq. (11), the interconnection vector, the short
 * module's inner-edge shift and the Birkeland scaling factors.
 *
 * Resolved once by @ref t01_state and read by every point; this is also, float for float, the
 * device kernel's parameter block (@ref t01_param_block). The `float` instantiation is how the
 * host fp32 lane is made to run the same arithmetic as the device: the amplitudes are formed in
 * `double` and ROUNDED, so a host-device disagreement is attributable to the device.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the kernel.
 * @test IrbemT01.StateCarriesThePublishedAmplitudes
 */
template <std::floating_point T>
struct T01State {
    /// `sin(Psi)` — the only trigonometric function of the tilt the model reads.
    T sin_tilt;
    /// `Psi`, radians — the tilt ANGLE, which P1 eq. (21) uses directly as `beta Psi`.
    T tilt;
    /// `t_1`, the short tail module's amplitude.
    T t1;
    /// `t_2`, the long tail module's amplitude.
    T t2;
    /// Birkeland amplitudes `b_l^(m)`: Region 1 `m = 1`, Region 1 `m = 2`, Region 2 `m = 1`,
    /// Region 2 `m = 2` — P2's numbers; the unit constant is applied by the evaluator.
    std::array<T, 4> b_fac;
    /// The interconnection field, `epsilon B_perp`: its `y` and `z` components, nT.
    T int_y;
    T int_z;
    /// The short module's fitted inner-edge shift `Delta X`, R_E, P2 eq. (3).
    T shift1;
    /// `varsigma_1`, `varsigma_2` — the Region 1 and 2 scaling factors, P2 eq. (9).
    T zeta1;
    T zeta2;
};

/**
 * The corrected `Dst*` of P2 ¶30: `0.8 Dst - 13 sqrt(P_dyn)` — the part of `Dst` attributable to
 * magnetospheric sources once the induction and magnetopause contributions are removed.
 *
 * @param dst `Dst`, nT.
 * @param pdyn `P_dyn`, nPa; a non-positive value makes the square root NaN, which the caller checks.
 * @return `Dst*`, nT.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01.StateCarriesThePublishedAmplitudes
 */
[[nodiscard]] inline double t01_dst_star(double dst, double pdyn) {
    return (0.8 * dst) - (13.0 * std::sqrt(pdyn));
}

/**
 * Resolve the drivers and the tilt into the evaluator's state — P2 eqs. (2), (3), (8), (9), (10).
 *
 * @tparam T the scalar type of the state; the arithmetic is `double` and rounded on the way out.
 * @param tilt_rad the dipole tilt `Psi`, radians.
 * @param d the drivers. `pdyn` must be positive; the caller (@ref t01_field) guarantees it.
 * @return the state.
 * @complexity O(1): two `pow`, three `sqrt`, one `sin`.
 * @alloc none.
 * @test IrbemT01.StateCarriesThePublishedAmplitudes
 */
template <std::floating_point T>
[[nodiscard]] inline T01State<T> t01_state(double tilt_rad, const T01Drivers& d) {
    const T01Coefficients& k = t01_coefficients;
    const double dst_star = t01_dst_star(d.dst, d.pdyn);
    const double pr = d.pdyn / k.pdyn0;
    // P2 eq. (2): the two tail amplitudes.
    const double t1 = k.t1[0] + (k.t1[1] * std::pow(pr, k.alpha1)) + (k.t1[2] * d.g1) +
                      (k.t1[3] * dst_star);
    const double t2 = k.t2[0] + (k.t2[1] * std::pow(pr, k.alpha2)) + (k.t2[2] * d.g1) +
                      (k.t2[3] * dst_star);
    // P2 eq. (10): the penetration efficiency, from the clock angle theta = atan2(By, Bz) —
    // sin^2(theta / 2) = (1 - cos theta) / 2 = (1 - Bz / B_perp) / 2. A vanishing transverse IMF
    // has no clock angle and no interconnection field, so epsilon is then irrelevant.
    const double bperp = std::sqrt((d.by_imf * d.by_imf) + (d.bz_imf * d.bz_imf));
    const double sin2_half = bperp > 0.0 ? 0.5 * (1.0 - (d.bz_imf / bperp)) : 0.0;
    const double eps = k.eps[0] + (k.eps[1] * sin2_half);
    T01State<T> s{};
    s.sin_tilt = static_cast<T>(std::sin(tilt_rad));
    s.tilt = static_cast<T>(tilt_rad);
    s.t1 = static_cast<T>(t1);
    s.t2 = static_cast<T>(t2);
    // P2 eq. (8): b_l^(m) = b_l^(m0) + b_l^(m1) G2.
    s.b_fac[0] = static_cast<T>(k.b1[0] + (k.b1[1] * d.g2));
    s.b_fac[1] = static_cast<T>(k.b1[2] + (k.b1[3] * d.g2));
    s.b_fac[2] = static_cast<T>(k.b2[0] + (k.b2[1] * d.g2));
    s.b_fac[3] = static_cast<T>(k.b2[2] + (k.b2[3] * d.g2));
    s.int_y = static_cast<T>(eps * d.by_imf);
    s.int_z = static_cast<T>(eps * d.bz_imf);
    // P2 eq. (3), attached to the short module (see T01TailModule).
    s.shift1 = static_cast<T>(k.dx0 + (k.dx1 * d.g2));
    // P2 eq. (9).
    s.zeta1 = static_cast<T>(k.sigma1[0] + (k.sigma1[1] * d.g2));
    s.zeta2 = static_cast<T>(k.sigma2[0] + (k.sigma2[1] * d.g2));
    return s;
}

// -------------------------------------------------------------------------------------------
// The evaluator
// -------------------------------------------------------------------------------------------

namespace detail::t01 {

/// A three-vector in the scalar type of the lane; an aggregate so it lives in registers.
template <std::floating_point T>
struct Vec3 {
    T x;
    T y;
    T z;
};

/// A Jacobian `d r* / d r`, row-major: `m[3 i + j] = d r*_i / d r_j`.
template <std::floating_point T>
struct Jac {
    std::array<T, 9> m;
};

/**
 * The divergence-free push-forward of a field through a deformation: `B'(r) = adj(J) B*(T r)`.
 *
 * With `J = d r* / d r`, the adjugate `adj(J) = det(J) J^-1` applied to `B*` at the image point is
 * the pull-back of the flux two-form `B . dS`, so `div B' = det(J) (div B*)(T r) = 0` whenever `B*`
 * is divergence-free. P1 eqs. (8)-(10), (14) and (22)-(24) are this identity in cylindrical and
 * Cartesian coordinates; writing it once means every deformation in the model is exact by
 * construction and the only thing that can go wrong is a Jacobian entry — which the suite's
 * `div B` test catches.
 *
 * @tparam T the scalar type.
 * @param j the Jacobian of the map at the ORIGINAL point.
 * @param b the undeformed field at the IMAGE point.
 * @return `adj(j) b`.
 * @complexity O(1): nine 2x2 determinants and one 3x3 product.
 * @alloc none.
 * @test IrbemT01.PushForwardIsTheAdjugate
 */
template <std::floating_point T>
[[nodiscard]] constexpr Vec3<T> push_forward(const Jac<T>& j, const Vec3<T>& b) {
    const std::array<T, 9>& a = j.m;
    // Cofactors C_ki; (adj J b)_i = sum_k C_ki b_k.
    const T c00 = (a[4] * a[8]) - (a[5] * a[7]);
    const T c01 = -((a[3] * a[8]) - (a[5] * a[6]));
    const T c02 = (a[3] * a[7]) - (a[4] * a[6]);
    const T c10 = -((a[1] * a[8]) - (a[2] * a[7]));
    const T c11 = (a[0] * a[8]) - (a[2] * a[6]);
    const T c12 = -((a[0] * a[7]) - (a[1] * a[6]));
    const T c20 = (a[1] * a[5]) - (a[2] * a[4]);
    const T c21 = -((a[0] * a[5]) - (a[2] * a[3]));
    const T c22 = (a[0] * a[4]) - (a[1] * a[3]);
    return {(c00 * b.x) + (c10 * b.y) + (c20 * b.z), (c01 * b.x) + (c11 * b.y) + (c21 * b.z),
            (c02 * b.x) + (c12 * b.y) + (c22 * b.z)};
}

/**
 * The cross-tail disc field, P1 eqs. (2)-(5), in the module's SCALED coordinates.
 *
 * The potential `A = A_phi e_phi = P (X, Y, Z) (-Y, X, 0)` with `P = A_phi / rho`, so
 * `B_x = -X dP/dZ`, `B_y = -Y dP/dZ`, `B_z = 2P + X dP/dX + Y dP/dY` — the same curl T89 uses.
 * Each term of eq. (2) collapses to `P_i = 2 f_i b_i sqrt(v_i - 4 b_i^2) / v_i^(3/2)` with
 * `v_i = S_i^(1)^2 S_i^(2)^2 = (rho^2 + b_i^2 + w_i^2)^2 - 4 b_i^2 rho^2`, `w_i = zeta + c_i`,
 * `zeta = sqrt(Z^2 + D^2)` — and `v_i - 4 b_i^2 >= 4 b_i^2 (w_i^2 - 1) > 0` because every
 * `w_i >= eta D_0 + c_i > 1`, so the square root is always real. `D` depends on `X` and `Y`
 * (eq. 6) and that dependence is differentiated, not dropped: dropping it breaks `div B = 0`.
 *
 * @tparam T the scalar type.
 * @param x the scaled coordinates, R_E.
 * @param y see @p x.
 * @param z see @p x.
 * @param d the sheet half-thickness at `(x, y)`, already scaled by `eta`, R_E.
 * @param d_x `dd / dx`, and @param d_y `dd / dy`.
 * @return the disc's field, in the scaled frame's axes (which are the GSM axes).
 * @complexity O(1): five terms, each two square roots.
 * @alloc none.
 * @test IrbemT01.DivergenceVanishesEverywhere
 */
template <std::floating_point T>
[[nodiscard]] inline Vec3<T> tail_disc(T x, T y, T z, T d, T d_x, T d_y) {
    const T two = static_cast<T>(2);
    const T four = static_cast<T>(4);
    const T rho2 = (x * x) + (y * y);
    const T zeta = std::sqrt((z * z) + (d * d));
    const T dzeta_dx = d * d_x / zeta;
    const T dzeta_dy = d * d_y / zeta;
    const T dzeta_dz = z / zeta;
    T p = static_cast<T>(0);
    T p_x = static_cast<T>(0);
    T p_y = static_cast<T>(0);
    T p_z = static_cast<T>(0);
    for (std::size_t i = 0; i < t01_tail_term_count; ++i) {
        const T f = static_cast<T>(t01_tail_disc.f[i]);
        const T b = static_cast<T>(t01_tail_disc.b[i]);
        const T c = static_cast<T>(t01_tail_disc.c[i]);
        const T w = zeta + c;
        const T b2 = b * b;
        const T u = rho2 + b2 + (w * w);
        const T v = (u * u) - (four * b2 * rho2);
        const T q = std::sqrt(v - (four * b2));
        const T sv = std::sqrt(v);
        p += two * f * b * q / (v * sv);
        const T dp_dv = f * b * ((static_cast<T>(12) * b2) - (two * v)) / (v * v * sv * q);
        const T dv_drho2 = (two * u) - (four * b2);
        const T dv_dw = four * u * w;
        p_x += dp_dv * ((dv_drho2 * two * x) + (dv_dw * dzeta_dx));
        p_y += dp_dv * ((dv_drho2 * two * y) + (dv_dw * dzeta_dy));
        p_z += dp_dv * dv_dw * dzeta_dz;
    }
    return {-x * p_z, -y * p_z, (two * p) + (x * p_x) + (y * p_y)};
}

/**
 * One tail module's untilted field at a GSM point: the disc of @ref tail_disc under P1 ¶16's
 * substitution, `X' = eta (X - shift) - (eta - 1) X_m`, `Y' = eta Y`, `Z' = eta Z`, `D -> eta D`,
 * with eq. (6)'s `D` evaluated in the scaled coordinates.
 *
 * The substitution is affine, so `B(r) = B_disc(T r)` is divergence-free with no Jacobian factor;
 * the paper writes the substituted formulas and the amplitudes `t_1, t_2` absorb any constant.
 *
 * @tparam T the scalar type.
 * @param mod the module's `eta`, `X_m`, `Delta D_x` (its own `X_s` is in @p shift).
 * @param shift the total sunward shift, R_E: the module's `X_s` plus, for the short module, the
 *        fitted `Delta X` of P2 eq. (3).
 * @param d0 `D_0`, and @param ddy `Delta D_y`, of eq. (6), R_E.
 * @param x the GSM point, R_E; @param y and @param z likewise.
 * @return the module's field, GSM axes, per unit amplitude.
 * @complexity O(1): one `exp` plus @ref tail_disc.
 * @alloc none.
 * @test IrbemT01.DivergenceVanishesEverywhere
 */
template <std::floating_point T>
[[nodiscard]] inline Vec3<T> tail_module(const T01TailModule& mod, T shift, T d0, T ddy, T x, T y,
                                         T z) {
    const T eta = static_cast<T>(mod.eta);
    const T xs = (eta * (x - shift)) - ((eta - static_cast<T>(1)) * static_cast<T>(mod.x_m));
    const T ys = eta * y;
    const T zs = eta * z;
    const T dx = static_cast<T>(t01_tail_fixed.delta_x);
    const T dy = static_cast<T>(t01_tail_fixed.delta_y);
    const T ddx = static_cast<T>(mod.delta_d_x);
    const T flare = ddx * std::exp(xs / dx);
    const T yn = ys / dy;
    // eq. (6), then the eta scaling of D, with both partials carried.
    const T d = eta * (d0 + (ddy * yn * yn) + flare);
    const T d_x = eta * flare / dx;
    const T d_y = eta * static_cast<T>(2) * ddy * yn / dy;
    return tail_disc<T>(xs, ys, zs, d, d_x, d_y);
}

/**
 * The tilted cross-tail field at a GSM point: both modules, warped (P1 eqs. 7-10) and then bent
 * (eqs. 11-14), each deformation through @ref push_forward.
 *
 * The bending map is `X* = X cos Psi* - Z sin Psi*`, `Z* = X sin Psi* + Z cos Psi*` with
 * `sin Psi* = R_H sin Psi / (R_H^3 + r^3)^(1/3)` and `R_H = R_H0 + R_H2 Z^2 / r^2` — the rotation
 * that takes the dipole equator `z = -x tan Psi` onto `z* = 0` near the Earth and fades to nothing
 * far out. The warping map is a rotation of `(Y, Z)` by the position-dependent angle
 * `g = G sin Psi Y rho^2 / (rho^4 + L^4)` (eq. 7 with `rho cos phi = Y`, so no `atan2` is needed
 * and the axis `rho = 0` is an ordinary point).
 *
 * Order: the warp is P1's "first deformation" (¶19) and the bending its "second" (¶22), applied to
 * the warped field, so `B = Bend[Warp[B_0]]`: the point is first bent, then warped, then the
 * undeformed field is read; the Jacobians are applied in the reverse order.
 *
 * @tparam T the scalar type.
 * @param s the epoch's state (amplitudes, tilt, the short module's shift).
 * @param x the GSM point, R_E; @param y and @param z likewise.
 * @return `t_1 B_T1 + t_2 B_T2`, tilted, GSM, nT.
 * @complexity O(1): one `pow`, one `sin`/`cos` pair, two @ref tail_module evaluations.
 * @alloc none.
 * @test IrbemT01.DivergenceVanishesEverywhere
 * @test IrbemT01.TheWarpedSheetFollowsTheDipoleEquatorNearEarth
 */
template <std::floating_point T>
[[nodiscard]] inline Vec3<T> tail_field(const T01State<T>& s, T x, T y, T z) {
    const T zero = static_cast<T>(0);
    const T one = static_cast<T>(1);
    const T two = static_cast<T>(2);
    const T third = static_cast<T>(1.0 / 3.0);

    // ---- eqs. (11)-(13): the bending, and its Jacobian --------------------------------------
    const T r2 = (x * x) + (y * y) + (z * z);
    const T r = std::sqrt(r2);
    const T zr2 = r2 > zero ? (z * z) / r2 : zero;
    const T rh0 = static_cast<T>(t01_coefficients.rh0);
    const T rh2 = static_cast<T>(t01_tail_fixed.rh2);
    const T rh = rh0 + (rh2 * zr2);
    const T w = std::pow((rh * rh * rh) + (r2 * r), -third);
    const T sp = rh * s.sin_tilt * w;
    const T cp = std::sqrt(one - (sp * sp));
    const T w4 = w * w * w * w;
    // dR_H/dx_i = R_H2 * 2 Z (delta_iz r^2 - Z x_i) / r^4; dw/dx_i = -w^4 (R_H^2 dR_H/dx_i + r x_i).
    const T inv_r4 = r2 > zero ? one / (r2 * r2) : zero;
    const std::array<T, 3> pos{x, y, z};
    std::array<T, 3> dsp{};
    std::array<T, 3> dcp{};
    for (int i = 0; i < 3; ++i) {
        const T delta_iz = i == 2 ? one : zero;
        const T drh = rh2 * two * z * ((delta_iz * r2) - (z * pos[static_cast<std::size_t>(i)])) * inv_r4;
        const T dw = -w4 * ((rh * rh * drh) + (r * pos[static_cast<std::size_t>(i)]));
        dsp[static_cast<std::size_t>(i)] = s.sin_tilt * ((drh * w) + (rh * dw));
        dcp[static_cast<std::size_t>(i)] = -(sp / cp) * dsp[static_cast<std::size_t>(i)];
    }
    Jac<T> jb{};
    for (int i = 0; i < 3; ++i) {
        const T dix = i == 0 ? one : zero;
        const T diy = i == 1 ? one : zero;
        const T diz = i == 2 ? one : zero;
        const std::size_t k = static_cast<std::size_t>(i);
        jb.m[k] = (dix * cp) + (x * dcp[k]) - (diz * sp) - (z * dsp[k]);
        jb.m[3 + k] = diy;
        jb.m[6 + k] = (dix * sp) + (x * dsp[k]) + (diz * cp) + (z * dcp[k]);
    }
    const T x1 = (x * cp) - (z * sp);
    const T y1 = y;
    const T z1 = (x * sp) + (z * cp);

    // ---- eqs. (7)-(10): the warping at the bent point, and its Jacobian -----------------------
    const T l = static_cast<T>(t01_tail_fixed.l);
    const T l4 = l * l * l * l;
    const T gg = static_cast<T>(t01_coefficients.g) * s.sin_tilt;
    const T rho2 = (y1 * y1) + (z1 * z1);
    const T den = (rho2 * rho2) + l4;
    const T g = gg * y1 * rho2 / den;
    const T dcore = two * (l4 - (rho2 * rho2)) / (den * den);
    const T g_y = gg * ((rho2 / den) + (y1 * y1 * dcore));
    const T g_z = gg * y1 * z1 * dcore;
    const T cg = std::cos(g);
    const T sg = std::sin(g);
    const T x2 = x1;
    const T y2 = (y1 * cg) - (z1 * sg);
    const T z2 = (y1 * sg) + (z1 * cg);
    Jac<T> jw{};
    jw.m = {one, zero, zero, zero, cg - (z2 * g_y), -sg - (z2 * g_z), zero, sg + (y2 * g_y),
            cg + (y2 * g_z)};

    // ---- the untilted modules at the doubly deformed point -----------------------------------
    const T d0 = static_cast<T>(t01_coefficients.d0);
    const T ddy = static_cast<T>(t01_coefficients.delta_d_y);
    const Vec3<T> m1 = tail_module<T>(t01_tail_modules[0],
                                      static_cast<T>(t01_tail_modules[0].x_s) + s.shift1, d0, ddy,
                                      x2, y2, z2);
    const Vec3<T> m2 = tail_module<T>(t01_tail_modules[1], static_cast<T>(t01_tail_modules[1].x_s),
                                      d0, ddy, x2, y2, z2);
    const T unit = static_cast<T>(t01_units.tail);
    const Vec3<T> b0{unit * ((s.t1 * m1.x) + (s.t2 * m2.x)), unit * ((s.t1 * m1.y) + (s.t2 * m2.y)),
                     unit * ((s.t1 * m1.z) + (s.t2 * m2.z))};
    return push_forward<T>(jb, push_forward<T>(jw, b0));
}

/**
 * The floor under a cylindrical radius that sits on an axis. The conical harmonics and the
 * azimuthal deformations are regular ON their axes — the physical field there is finite — but their
 * spherical/cylindrical spellings have `0 / 0` there; evaluating an infinitesimal step off the axis
 * gives the limit to roundoff. `1e-9 R_E` is 6 mm.
 */
template <std::floating_point T>
inline constexpr T axis_floor = static_cast<T>(1e-9);

/**
 * One deformed conical Birkeland mode, northern current system only: P1 eqs. (16)-(19).
 *
 * The undeformed field is the curl of `A_r = B_m T^(m)(theta) sin(m phi)` — eq. (16) — so
 * `B_r = 0`, `B_theta = m B_m T^(m) cos(m phi) / (r sin theta)`, `B_phi = -B_m dT^(m)/dtheta
 * sin(m phi) / r`, with the piecewise `T^(m)` of eq. (17): `tan^m(theta / 2)` inside the cone,
 * `cot^m(theta / 2)` (times a matching constant) outside, and the current-carrying blend between
 * `theta_0 -+ Delta theta`. Both outer branches are curl-free — a check the suite makes — so all the
 * current flows in the sheet. The blend is evaluated in `t = tan(theta / 2)`; the outer branch in
 * `u = cot(theta / 2) = sin theta / (1 - cos theta)`, which is regular where `t` overflows at
 * `theta = pi`.
 *
 * That cone is then read at `(r', theta', phi)` of eqs. (18)-(19) and pushed forward through the
 * Cartesian Jacobian of that map, which is assembled by the chain rule from `dr'/dr, dr'/dtheta,
 * dtheta'/dr, dtheta'/dtheta` and the spherical-coordinate gradients.
 *
 * @tparam T the scalar type.
 * @param c the mode's Table 2 column.
 * @param x the GSM point, R_E; @param y and @param z likewise.
 * @return the northern system's field, per unit amplitude.
 * @complexity O(1): two `atan2`, one `sin`/`cos` pair, ~a dozen square roots.
 * @alloc none.
 * @test IrbemT01.DivergenceVanishesEverywhere
 * @test IrbemT01.ConeBranchesAreCurlFreeOutsideTheSheet
 */
template <std::floating_point T>
[[nodiscard]] inline Vec3<T> cone_mode(const T01Cone& c, T x, T y, T z) {
    const T zero = static_cast<T>(0);
    const T one = static_cast<T>(1);
    const T two = static_cast<T>(2);
    const T half = static_cast<T>(0.5);
    const T three = static_cast<T>(3);

    // ---- spherical coordinates of the point, with the axis floored ----------------------------
    const T rho_raw = std::sqrt((x * x) + (y * y));
    const T rho = rho_raw > axis_floor<T> ? rho_raw : axis_floor<T>;
    const T r = std::sqrt((rho * rho) + (z * z));
    const T th = std::atan2(rho, z);
    const T ph = std::atan2(y, x);
    const T st = std::sin(th);
    const T ct = std::cos(th);
    const T s2t = two * st * ct;
    const T c2t = (ct * ct) - (st * st);
    const T s3t = (s2t * ct) + (c2t * st);
    const T c3t = (c2t * ct) - (s2t * st);

    // ---- eq. (18): r' and its partials ---------------------------------------------------------
    const auto rad = [&](std::size_t k) {
        const T bk = static_cast<T>(c.b[k]);
        return std::sqrt((r * r) + (bk * bk));
    };
    const auto ak = [&](std::size_t k) { return static_cast<T>(c.a[k]); };
    const auto bk = [&](std::size_t k) { return static_cast<T>(c.b[k]); };
    // Terms of the form r / sqrt(r^2 + b^2): derivative b^2 / (r^2 + b^2)^(3/2);
    // r / (r^2 + b^2): derivative (b^2 - r^2) / (r^2 + b^2)^2.
    const auto s_term = [&](std::size_t k, T& val, T& der) {
        const T q = rad(k);
        val = r / q;
        der = bk(k) * bk(k) / (q * q * q);
    };
    const auto p_term = [&](std::size_t k, T& val, T& der) {
        const T q2 = (r * r) + (bk(k) * bk(k));
        val = r / q2;
        der = ((bk(k) * bk(k)) - (r * r)) / (q2 * q2);
    };
    T v1{};
    T d1{};
    T v2{};
    T d2{};
    T v3{};
    T d3{};
    T v4{};
    T d4{};
    T v5{};
    T d5{};
    T v6{};
    T d6{};
    s_term(0, v1, d1);
    p_term(1, v2, d2);
    s_term(2, v3, d3);
    p_term(3, v4, d4);
    s_term(4, v5, d5);
    p_term(5, v6, d6);
    const T inv_r = one / r;
    const T g0 = (ak(0) * inv_r) + (ak(1) * v1) + (ak(2) * v2);
    const T g0_r = (-ak(0) * inv_r * inv_r) + (ak(1) * d1) + (ak(2) * d2);
    const T g1 = ak(3) + (ak(4) * inv_r) + (ak(5) * v3) + (ak(6) * v4);
    const T g1_r = (-ak(4) * inv_r * inv_r) + (ak(5) * d3) + (ak(6) * d4);
    const T g2 = (ak(7) * v5) + (ak(8) * v6);
    const T g2_r = (ak(7) * d5) + (ak(8) * d6);
    const T r1 = r + g0 + (g1 * ct) + (g2 * c2t);
    const T r1_r = one + g0_r + (g1_r * ct) + (g2_r * c2t);
    const T r1_th = (-g1 * st) - (g2 * two * s2t);

    // ---- eq. (19): theta' and its partials -----------------------------------------------------
    const auto ck = [&](std::size_t k) { return static_cast<T>(c.c[k]); };
    const auto dk = [&](std::size_t k) { return static_cast<T>(c.d[k]); };
    const auto sd_term = [&](std::size_t k, T& val, T& der) {
        const T q = std::sqrt((r * r) + (dk(k) * dk(k)));
        val = r / q;
        der = dk(k) * dk(k) / (q * q * q);
    };
    const auto pd_term = [&](std::size_t k, T& val, T& der) {
        const T q2 = (r * r) + (dk(k) * dk(k));
        val = r / q2;
        der = ((dk(k) * dk(k)) - (r * r)) / (q2 * q2);
    };
    T e1{};
    T f1{};
    T e2{};
    T f2{};
    T e3{};
    T f3{};
    T e4{};
    T f4{};
    sd_term(0, e1, f1);
    sd_term(1, e2, f2);
    pd_term(2, e3, f3);
    pd_term(3, e4, f4);
    const T h1 = ck(0) + (ck(1) * inv_r) + (ck(2) * inv_r * inv_r) + (ck(3) * e1);
    const T h1_r = (-ck(1) * inv_r * inv_r) - (two * ck(2) * inv_r * inv_r * inv_r) + (ck(3) * f1);
    const T h2 = ck(4) + (ck(5) * e2) + (ck(6) * e3);
    const T h2_r = (ck(5) * f2) + (ck(6) * f3);
    const T h3 = ck(7) + (ck(8) * inv_r) + (ck(9) * e4);
    const T h3_r = (-ck(8) * inv_r * inv_r) + (ck(9) * f4);
    const T th1 = th + (h1 * st) + (h2 * s2t) + (h3 * s3t);
    const T th1_r = (h1_r * st) + (h2_r * s2t) + (h3_r * s3t);
    const T th1_th = one + (h1 * ct) + (h2 * two * c2t) + (h3 * three * c3t);

    // ---- eqs. (16)-(17): the cone at (r', theta', phi) ----------------------------------------
    const T st1 = std::sin(th1);
    const T ct1 = std::cos(th1);
    const T tminus = std::tan(half * static_cast<T>(c.theta_0 - c.delta_theta));
    const T tplus = std::tan(half * static_cast<T>(c.theta_0 + c.delta_theta));
    const int m = c.m;
    const T mm = static_cast<T>(m);
    const T twom1 = (two * mm) + one;
    // t^m and t^(m-1) for m in {1, 2}, spelled without pow.
    const auto pow_m = [m](T t) { return m == 1 ? t : t * t; };
    const auto pow_m1 = [m](T t) { return m == 1 ? static_cast<T>(1) : t; };
    const auto pow_2m1 = [m](T t) {
        const T t2 = t * t;
        return m == 1 ? t2 * t : t2 * t2 * t;
    };
    T t_over_sin{};   // T^(m)(theta') / sin theta'
    T t_prime{};      // dT^(m)/dtheta at theta'
    const T t = st1 / (one + ct1);   // tan(theta'/2); finite away from theta' = pi
    const T dt = half * (one + (t * t));
    if (th1 < static_cast<T>(c.theta_0 - c.delta_theta)) {
        // Inside the cone: T = t^m.
        t_over_sin = pow_m1(t) * (one + (t * t)) * half;
        t_prime = mm * pow_m1(t) * dt;
    } else if (th1 <= static_cast<T>(c.theta_0 + c.delta_theta)) {
        // In the sheet: the eq. (17) blend, in t.
        const T dtt = tplus - tminus;
        const T tm = pow_m(t);
        const T tm1 = pow_m1(t);
        const T tail = (pow_2m1(t) - pow_2m1(tminus)) / (twom1 * tm);
        const T tt = ((tm * (tplus - t)) + tail) / dtt;
        const T dtail = (((mm + one) * tm) + (mm * pow_2m1(tminus) / (tm * t))) / twom1;
        const T dtt_dt = ((mm * tm1 * (tplus - t)) - tm + dtail) / dtt;
        t_over_sin = tt * (one + (t * t)) / (two * t);
        t_prime = dtt_dt * dt;
    } else {
        // Outside the cone: T = K cot^m(theta/2), in u = cot(theta'/2), regular at theta' = pi.
        const T kk = (pow_2m1(tplus) - pow_2m1(tminus)) / (twom1 * (tplus - tminus));
        const T u = st1 / (one - ct1);
        t_over_sin = kk * pow_m1(u) * (one + (u * u)) * half;
        t_prime = -kk * mm * pow_m1(u) * (one + (u * u)) * half;
    }
    const T cph = std::cos(ph);
    const T sph = std::sin(ph);
    const T cmph = m == 1 ? cph : (cph * cph) - (sph * sph);
    const T smph = m == 1 ? sph : two * sph * cph;
    const T bm = static_cast<T>(c.b_m);
    const T b_th = mm * bm * t_over_sin * cmph / r1;
    const T b_ph = -bm * t_prime * smph / r1;
    const Vec3<T> that{ct1 * cph, ct1 * sph, -st1};
    const Vec3<T> phat{-sph, cph, zero};
    const Vec3<T> shat{st1 * cph, st1 * sph, ct1};
    const Vec3<T> bstar{(b_th * that.x) + (b_ph * phat.x), (b_th * that.y) + (b_ph * phat.y),
                        (b_th * that.z) + (b_ph * phat.z)};

    // ---- the Cartesian Jacobian of (x, y, z) -> r' s-hat(theta', phi) --------------------------
    const std::array<T, 3> dq_dr{(r1_r * shat.x) + (r1 * th1_r * that.x),
                                 (r1_r * shat.y) + (r1 * th1_r * that.y),
                                 (r1_r * shat.z) + (r1 * th1_r * that.z)};
    const std::array<T, 3> dq_dth{(r1_th * shat.x) + (r1 * th1_th * that.x),
                                  (r1_th * shat.y) + (r1 * th1_th * that.y),
                                  (r1_th * shat.z) + (r1 * th1_th * that.z)};
    const std::array<T, 3> dq_dph{r1 * st1 * phat.x, r1 * st1 * phat.y, zero};
    const T xr = rho * cph;
    const T yr = rho * sph;
    const T inv_r2 = inv_r * inv_r;
    const std::array<T, 3> dr_dq{xr * inv_r, yr * inv_r, z * inv_r};
    const std::array<T, 3> dth_dq{xr * z * inv_r2 / rho, yr * z * inv_r2 / rho, -rho * inv_r2};
    const std::array<T, 3> dph_dq{-yr / (rho * rho), xr / (rho * rho), zero};
    Jac<T> j{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t k = 0; k < 3; ++k) {
            j.m[(3 * i) + k] =
                (dq_dr[i] * dr_dq[k]) + (dq_dth[i] * dth_dq[k]) + (dq_dph[i] * dph_dq[k]);
        }
    }
    return push_forward<T>(j, bstar);
}

/**
 * One region's Birkeland field at a GSM point: both modes, both hemispheres (P1 eq. 20), the
 * azimuthal and tilt deformation of eq. (21) through @ref push_forward, and the scaling of eq. (25).
 *
 * Eq. (21) is a rotation of `(X, Z)` about `y` by the position-dependent angle `h = a Z / rho -
 * tau` with `a = Delta phi + b (rho^2 - 1) / (rho^2 + rho_0^2)`, `tau = beta Psi / [1 + (r /
 * R_H)^3]^(1/3)`: `X* = X cos h + Z sin h`, `Z* = Z cos h - X sin h` (the paper's `phi = -arctan(Z
 * / X)` and `X* = rho cos F`, `Z* = -rho sin F` say exactly this). Eq. (25)'s scaling reads the
 * whole thing at `zeta r` and multiplies by `zeta`.
 *
 * @tparam T the scalar type.
 * @param region 0 for Region 1, 1 for Region 2 — selects the cones, `Delta phi` and `zeta`.
 * @param s the epoch's state.
 * @param x the GSM point, R_E; @param y and @param z likewise.
 * @return the region's field, GSM, nT.
 * @complexity O(1): four @ref cone_mode evaluations and one `pow`.
 * @alloc none.
 * @test IrbemT01.DivergenceVanishesEverywhere
 * @test IrbemT01.BirkelandCurrentsAreAntisymmetricAcrossNoonMidnightAtZeroTilt
 */
template <std::floating_point T>
[[nodiscard]] inline Vec3<T> fac_region(int region, const T01State<T>& s, T x, T y, T z) {
    const T zero = static_cast<T>(0);
    const T one = static_cast<T>(1);
    const T two = static_cast<T>(2);
    const T third = static_cast<T>(1.0 / 3.0);
    const std::size_t reg = static_cast<std::size_t>(region);
    const T zeta = region == 0 ? s.zeta1 : s.zeta2;
    const T unit = static_cast<T>(t01_units.fac);
    const T amp1 = unit * s.b_fac[2 * reg];
    const T amp2 = unit * s.b_fac[(2 * reg) + 1];
    const T01Cone& c1 = t01_cones[2 * reg];
    const T01Cone& c2 = t01_cones[(2 * reg) + 1];

    // ---- eq. (25): the point the scaled system is read at ------------------------------------
    const T xs = zeta * x;
    const T ys = zeta * y;
    const T zs = zeta * z;

    // ---- eq. (21): the azimuthal + tilt deformation, and its Jacobian --------------------------
    const T01FacDeformation& fd = t01_fac_deformation;
    const T rho_raw = std::sqrt((xs * xs) + (zs * zs));
    const T rho = rho_raw > axis_floor<T> ? rho_raw : axis_floor<T>;
    const T rho2 = rho * rho;
    const T r = std::sqrt(rho2 + (ys * ys));
    const T rho0 = static_cast<T>(fd.rho_0);
    const T bb = static_cast<T>(fd.b);
    const T qden = rho2 + (rho0 * rho0);
    const T a = static_cast<T>(fd.delta_phi[reg]) + (bb * (rho2 - one) / qden);
    const T a_coef = two * bb * ((rho0 * rho0) + one) / (qden * qden);   // da/dx_i = a_coef x_i
    const T rh = static_cast<T>(fd.r_h);
    const T rr = r / rh;
    const T tau_base = one + (rr * rr * rr);
    const T tau = static_cast<T>(fd.beta) * s.tilt * std::pow(tau_base, -third);
    // dtau/dx_i = -beta Psi (1 + (r/R_H)^3)^(-4/3) r x_i / R_H^3.
    const T tau_coef = -static_cast<T>(fd.beta) * s.tilt * std::pow(tau_base, -static_cast<T>(4) * third) *
                       r / (rh * rh * rh);
    const T zr = zs / rho;
    const T h = (a * zr) - tau;
    const T h_x = (a_coef * xs * zr) + (a * (-xs * zs / (rho2 * rho))) - (tau_coef * xs);
    const T h_y = -(tau_coef * ys);
    const T h_z = (a_coef * zs * zr) + (a * xs * xs / (rho2 * rho)) - (tau_coef * zs);
    const T ch = std::cos(h);
    const T sh = std::sin(h);
    const T xd = (xs * ch) + (zs * sh);
    const T zd = (zs * ch) - (xs * sh);
    Jac<T> j{};
    j.m = {ch + (zd * h_x), zd * h_y, sh + (zd * h_z), zero, one, zero, -sh - (xd * h_x), -(xd * h_y),
           ch - (xd * h_z)};

    // ---- eq. (20): northern system plus its mirrored, polarity-reversed southern twin ----------
    const Vec3<T> n1 = cone_mode<T>(c1, xd, ys, zd);
    const Vec3<T> s1 = cone_mode<T>(c1, xd, -ys, -zd);
    const Vec3<T> n2 = cone_mode<T>(c2, xd, ys, zd);
    const Vec3<T> s2 = cone_mode<T>(c2, xd, -ys, -zd);
    const Vec3<T> sum{(amp1 * (n1.x - s1.x)) + (amp2 * (n2.x - s2.x)),
                      (amp1 * (n1.y + s1.y)) + (amp2 * (n2.y + s2.y)),
                      (amp1 * (n1.z + s1.z)) + (amp2 * (n2.z + s2.z))};
    const Vec3<T> b = push_forward<T>(j, sum);
    return {zeta * b.x, zeta * b.y, zeta * b.z};
}

}  // namespace detail::t01

/**
 * The T01 external field at one GSM point, as three components in nanotesla — P1 eq. (1) with the
 * terms this file carries: `B_T + B_FAC(R1) + B_FAC(R2) + B_INT`.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device kernel.
 * @param s the epoch's state, from @ref t01_state.
 * @param x the GSM x coordinate, R_E; @param y and @param z likewise.
 * @return `{B_x, B_y, B_z}` in GSM, nanotesla.
 * @complexity O(1): roughly 3 000 flops, ~60 square roots, ~20 transcendentals; no data-dependent
 *             loop and one three-way branch per cone evaluation (which one of eq. (17)'s pieces the
 *             deformed colatitude lands in).
 * @alloc none.
 * @test IrbemT01.DivergenceVanishesEverywhere
 * @test IrbemT01.ZeroTiltIsMirrorSymmetricAboutTheEquator
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, 3> t01_components(const T01State<T>& s, T x, T y, T z) {
    const detail::t01::Vec3<T> tail = detail::t01::tail_field<T>(s, x, y, z);
    const detail::t01::Vec3<T> r1 = detail::t01::fac_region<T>(0, s, x, y, z);
    const detail::t01::Vec3<T> r2 = detail::t01::fac_region<T>(1, s, x, y, z);
    return {tail.x + r1.x + r2.x, tail.y + r1.y + r2.y + s.int_y, tail.z + r1.z + r2.z + s.int_z};
}

/**
 * The T01 external field at one GSM point, in `double` — the reference lane, from a resolved state.
 *
 * @param p the position, GSM, in Earth radii.
 * @param s the epoch's state; see @ref t01_state.
 * @return the external field at @p p, GSM, in nanotesla.
 * @complexity O(1); see @ref t01_components.
 * @alloc none.
 * @test IrbemT01.ReferenceLaneMatchesTheComponentForm
 */
[[nodiscard]] inline FieldVector<Frame::GSM> t01_field_at(Position<Frame::GSM> p,
                                                          const T01State<double>& s) {
    const std::array<double, 3> b = t01_components<double>(s, p.v[0], p.v[1], p.v[2]);
    return FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
}

/**
 * Whether the drivers can be evaluated at all, as opposed to being outside the fit.
 *
 * `P_dyn` enters as `(P_dyn / P_d0)^alpha` and `sqrt(P_dyn)`, neither of which has a real value
 * for a non-positive pressure; that is a domain error, not an extrapolation.
 *
 * @param d the drivers.
 * @return `true` when every driver is finite and `P_dyn > 0`.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01.NonFiniteInputIsADomainError
 */
[[nodiscard]] inline bool t01_drivers_evaluable(const T01Drivers& d) {
    return std::isfinite(d.dst) && std::isfinite(d.pdyn) && std::isfinite(d.by_imf) &&
           std::isfinite(d.bz_imf) && std::isfinite(d.g1) && std::isfinite(d.g2) && d.pdyn > 0.0;
}

/**
 * The published envelope's verdict on a driver set: `status.hpp`'s T01 row, `-50 <= Dst <= 20`,
 * `0.5 <= P_dyn <= 5`, `|B_y|, |B_z| <= 5`, `0 <= G1, G2 <= 10`.
 *
 * @param d the drivers.
 * @return @ref Status::OutOfValidityRange when any is outside its published closed interval,
 *         @ref Status::DomainError when any is not finite, otherwise @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01.EveryDriverIsCheckedFromBothSides
 */
[[nodiscard]] inline Status t01_check_drivers(const T01Drivers& d) {
    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Dst)] = d.dst;
    drivers[static_cast<std::size_t>(Driver::Pdyn)] = d.pdyn;
    drivers[static_cast<std::size_t>(Driver::ByIMF)] = d.by_imf;
    drivers[static_cast<std::size_t>(Driver::BzIMF)] = d.bz_imf;
    drivers[static_cast<std::size_t>(Driver::G1)] = d.g1;
    drivers[static_cast<std::size_t>(Driver::G2)] = d.g2;
    return check_validity(ExternalModel::Tsyganenko2001, drivers);
}

/**
 * The T01 external field, with the model's own verdict on whether it should be believed here.
 *
 * The value is **always** returned when it exists, including under @ref Status::OutOfValidityRange
 * — `status.hpp`'s standing rule. What is refused is arithmetic with no answer: a non-finite input,
 * a non-positive `P_dyn`, `|Psi| >= pi / 2` (the tilt is a rotation angle with a definite sign and
 * the deformations were fitted for `|Psi| <= ~35 degrees`; the mathematics survives up to but not
 * including a right angle, where `sin Psi* = sin Psi` on the axis and the bending Jacobian's
 * `1 / cos Psi*` does not exist), a point inside the Earth, or a non-finite answer.
 *
 * Two envelopes are checked through `status.hpp`, so the rules live in one place: the six drivers
 * against the T01 row, and the position against `x_GSM >= -15 R_E` — the model was fitted to data
 * sunward of that plane (P2 ¶46) and has no tail current beyond it to speak of; the tail modules
 * still evaluate there and the status says so.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `Psi`, radians; positive when the north dipole leans sunward.
 * @param d the drivers.
 * @return the field and its caveat.
 * @complexity O(1): one @ref t01_state and one @ref t01_components.
 * @alloc none.
 * @test IrbemT01.OutOfRangeDriversAreReportedButStillEvaluated
 * @test IrbemT01.TheTailwardBoundaryIsCheckedFromBothSides
 * @test IrbemT01.NonFiniteInputIsADomainError
 * @test IrbemT01.RightAngleTiltIsADomainError
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> t01_field(Position<Frame::GSM> p,
                                                               double tilt_rad,
                                                               const T01Drivers& d) {
    const FieldVector<Frame::GSM> zero{};
    if (!std::isfinite(p.v[0]) || !std::isfinite(p.v[1]) || !std::isfinite(p.v[2]) ||
        !std::isfinite(tilt_rad) || !t01_drivers_evaluable(d)) {
        return {Status::DomainError, zero};
    }
    if (!(std::fabs(tilt_rad) < max_tilt_rad)) return {Status::DomainError, zero};
    const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
    const Status where = check_position(ExternalModel::Tsyganenko2001, r, p.v[0]);
    if (where == Status::DomainError) return {Status::DomainError, zero};
    const Status drives = t01_check_drivers(d);
    const FieldVector<Frame::GSM> b = t01_field_at(p, t01_state<double>(tilt_rad, d));
    if (!std::isfinite(b.v[0]) || !std::isfinite(b.v[1]) || !std::isfinite(b.v[2])) {
        return {Status::DomainError, zero};
    }
    return {first_failure(drives, where), b};
}

/**
 * The six T01 drivers out of a context's `maginput` vector.
 *
 * `Dst`, `P_dyn`, `B_y`, `B_z` are on the hot block; `G1` and `G2` are read from the cold driver
 * vector, two loads that no per-point loop should notice next to the evaluator's ~3 000 flops.
 *
 * @param ctx the epoch's context.
 * @return the drivers.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01.ContextOverloadAgreesWithTheExplicitOne
 */
[[nodiscard]] inline T01Drivers t01_drivers_of(const FieldContext& ctx) {
    const HotState& h = ctx.hot();
    return {h.dst, h.pdyn, h.by_imf, h.bz_imf, ctx.drivers()[static_cast<std::size_t>(Driver::G1)],
            ctx.drivers()[static_cast<std::size_t>(Driver::G2)]};
}

/**
 * The T01 external field for a whole epoch's worth of state — the context entry point.
 *
 * @param p the position, GSM, in Earth radii.
 * @param ctx the epoch's context.
 * @return the field and its caveat, exactly as the three-argument overload.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01.ContextOverloadAgreesWithTheExplicitOne
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> t01_field(Position<Frame::GSM> p,
                                                               const FieldContext& ctx) {
    return t01_field(p, ctx.hot().tilt_rad, t01_drivers_of(ctx));
}

// -------------------------------------------------------------------------------------------
// The batch lanes
// -------------------------------------------------------------------------------------------

/**
 * The T01 field over a whole batch, on the CPU, in `float` — the host twin of `irbem_t01_f32`.
 *
 * The same expressions, in the same order, in the same precision, from a state rounded to `float`
 * FIRST; what makes a device disagreement attributable to the device.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, GSM, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, nanotesla; overwritten in full.
 * @param s the epoch's state, in `float`.
 * @return `false` when @p pos is not a whole number of points or @p out is a different length, in
 *         which case nothing is written; `true` otherwise.
 * @complexity O(N).
 * @alloc none.
 * @test IrbemT01.HostFloatLaneTracksTheReferenceLane
 * @test IrbemT01.HostFloatLaneRejectsMismatchedSpans
 */
[[nodiscard]] inline bool t01_field_host(std::span<const float> pos, std::span<float> out,
                                         const T01State<float>& s) {
    if (pos.size() % 3 != 0 || out.size() != pos.size()) return false;
    const std::size_t n = pos.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const std::array<float, 3> b =
            t01_components<float>(s, pos[(3 * i) + 0], pos[(3 * i) + 1], pos[(3 * i) + 2]);
        out[(3 * i) + 0] = b[0];
        out[(3 * i) + 1] = b[1];
        out[(3 * i) + 2] = b[2];
    }
    return true;
}

/// How many `float` scalars the device kernel's parameter buffer holds — the @ref T01State fields
/// in declaration order. Asserted against the kernel registry.
inline constexpr std::size_t t01_param_count = 13;

/**
 * Pack a state into the kernel's parameter buffer.
 *
 * The layout is the kernel's ABI and is stated in exactly two places — here and above
 * `irbem_t01_f32` in `irbem.slang`: `sin Psi, Psi, t_1, t_2, b_R1^(1), b_R1^(2), b_R2^(1),
 * b_R2^(2), int_y, int_z, shift_1, zeta_1, zeta_2`.
 *
 * @param s the state, in `float`.
 * @return the parameter block, by value.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemT01.ParameterBlockCarriesTheStateInOrder
 */
[[nodiscard]] inline std::array<float, t01_param_count> t01_param_block(const T01State<float>& s) {
    return {s.sin_tilt, s.tilt,  s.t1,     s.t2,     s.b_fac[0], s.b_fac[1], s.b_fac[2],
            s.b_fac[3], s.int_y, s.int_z, s.shift1, s.zeta1,    s.zeta2};
}

/**
 * The batch's position caveat, accumulated one point at a time — @ref T89PositionFold's shape for
 * T01's envelope, which is a tailward plane rather than a radius: the smallest `x_GSM` and the
 * smallest `r` decide the batch, and the radius is compared as a square, so there is no per-point
 * `sqrt`.
 *
 * @test IrbemT01.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
struct T01PositionFold {
    /// The smallest `r^2` seen, R_E^2; `+inf` until the first point.
    double r2_lo = std::numeric_limits<double>::infinity();
    /// The smallest GSM `x` seen, R_E; `+inf` until the first point.
    double x_lo = std::numeric_limits<double>::infinity();
    /// False once any point has had a non-finite coordinate.
    bool finite = true;

    /**
     * Fold one position in.
     * @param p the position, GSM, in Earth radii.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemT01.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    constexpr void add(const Position<Frame::GSM>& p) {
        const double r2 = (p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]);
        finite = finite && std::isfinite(r2);
        r2_lo = r2 < r2_lo ? r2 : r2_lo;
        x_lo = p.v[0] < x_lo ? p.v[0] : x_lo;
    }

    /**
     * What the batch's positions say about the model's envelope.
     * @return @ref Status::DomainError when any point is not finite or is inside the Earth,
     *         @ref Status::OutOfValidityRange when any point is tailward of `x_GSM = -15 R_E`,
     *         otherwise @ref Status::Ok.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemT01.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    [[nodiscard]] Status verdict() const {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return check_position(ExternalModel::Tsyganenko2001, finite ? std::sqrt(r2_lo) : nan, x_lo);
    }
};

/**
 * The T01 field over a whole batch of GSM points, on the device when that is worth it.
 *
 * The routine to call for more than a handful of points; @ref t01_field is what it is verified
 * against. T01 is ~3 000 flops for 24 bytes in and 12 out — an order of magnitude denser than T89,
 * which already wins 15x at 2^22 points on this seam — so the device lane is the expected home of
 * any real batch; the registry row in `gpu/dispatch.hpp` records the measured crossover.
 *
 * The batch reports the worst caveat of its points, folded in the staging loop, exactly as
 * @ref t89_field_batch does: a domain error anywhere zeroes every output; an out-of-validity point
 * or driver caveats the whole batch and every point is still computed.
 *
 * @param points the positions, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `Psi`, radians.
 * @param d the drivers.
 * @param out receives one field vector per input; same length as @p points.
 * @return the status, and `true` exactly when the device lane serviced the call.
 * @complexity O(N); on the device those N run concurrently over `ceil(N / 256)` workgroups.
 * @alloc the device lane stages two `std::vector<float>` of `3N`; the host lane allocates nothing.
 * @test IrbemT01.BatchAgreesWithTheReferenceLane
 * @test IrbemT01.BatchRejectsMismatchedSpans
 * @test IrbemT01.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
[[nodiscard]] inline Result<bool> t01_field_batch(std::span<const Position<Frame::GSM>> points,
                                                  double tilt_rad, const T01Drivers& d,
                                                  std::span<FieldVector<Frame::GSM>> out) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(tilt_rad) || !t01_drivers_evaluable(d)) return {Status::DomainError, false};
    if (!(std::fabs(tilt_rad) < max_tilt_rad)) return {Status::DomainError, false};
    const Status drives = t01_check_drivers(d);
    if (n == 0) return {drives, false};
    T01PositionFold fold;

#if CHEATAH_SPACE_IRBEM_T01_GPU
    if (gpu::prefer_gpu("irbem_t01_f32", n) &&
        std::filesystem::exists(gpu::shader_path("irbem_t01_f32"))) {
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
        const std::array<float, t01_param_count> block =
            t01_param_block(t01_state<float>(tilt_rad, d));
        gpu::dispatch_batch("irbem_t01_f32", pos, raw, std::span<const float>(block));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = FieldVector<Frame::GSM>{
                fixarray::vec3d{raw[(3 * i) + 0], raw[(3 * i) + 1], raw[(3 * i) + 2]}};
        }
        return {first_failure(drives, where), true};
    }
#endif

    const T01State<double> s = t01_state<double>(tilt_rad, d);
    for (std::size_t i = 0; i < n; ++i) {
        fold.add(points[i]);
        const std::array<double, 3> b =
            t01_components<double>(s, points[i].v[0], points[i].v[1], points[i].v[2]);
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
// The total field: IGRF plus T01, as one GeoFieldModel
// -------------------------------------------------------------------------------------------

/**
 * IGRF plus Tsyganenko 2001, as a single field — `TotalFieldT89`'s shape for a continuously driven
 * model.
 *
 * Satisfies @ref GeoFieldModel, so it drops into `trace_invariant`, `make_lstar` and everything
 * above them. The state is resolved ONCE at construction — the amplitudes are functions of the
 * drivers, not of the point — so a trace pays for the driver arithmetic never and for the
 * evaluator once per field call.
 *
 * @tparam NMAX the internal field's truncation degree; 10 reproduces IRBEM's own choice.
 * @test IrbemT01.TotalFieldSuperposesInternalAndExternal
 */
template <int NMAX = 10>
class TotalFieldT01 {
  public:
    /// The internal part's truncation degree, so `M::degree` means the same thing as on @ref Igrf.
    static constexpr int degree = NMAX;

    /**
     * @param internal the internal field, already built for the epoch.
     * @param rotations the epoch's frame rotations — built once, reused for every point.
     * @param drivers the six drivers; resolved to the model's state here, once.
     */
    TotalFieldT01(const Igrf<NMAX>& internal, const Rotations& rotations, const T01Drivers& drivers)
        : internal_(&internal),
          rotations_(&rotations),
          drivers_(drivers),
          state_(t01_state<double>(rotations.dipole_tilt_deg * (std::numbers::pi / 180.0), drivers)) {}

    /**
     * The total field at a geographic point.
     *
     * @param p the position, GEO, Earth radii.
     * @return `B_internal + B_external`, GEO, nT. When the external model refuses the point
     *         (@ref Status::DomainError) the INTERNAL field is returned alone, for the reason
     *         `TotalFieldT89::evaluate` gives: it is still the best available answer, and a NaN
     *         would fail hundreds of RK4 steps later with nothing pointing at its cause.
     * @complexity One IGRF evaluation, one T01 evaluation, two 3x3 rotations.
     * @alloc none.
     * @test IrbemT01.TotalFieldSuperposesInternalAndExternal
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::GEO>& p) const {
        const FieldVector<Frame::GEO> b_int = internal_->evaluate(p);
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        const Result<FieldVector<Frame::GSM>> b_ext =
            t01_field(p_gsm, rotations_->dipole_tilt_deg * (std::numbers::pi / 180.0), drivers_);
        if (b_ext.status == Status::DomainError) return b_int;
        const FieldVector<Frame::GEO> b_ext_geo = transform<Frame::GEO>(b_ext.value, *rotations_);
        return FieldVector<Frame::GEO>{b_int.v + b_ext_geo.v};
    }

    /**
     * Whether the external model answered at @p p, and if not, why.
     * @param p the position, GEO, Earth radii.
     * @return the external model's status; @ref Status::Ok when it contributed without caveat.
     * @complexity One T01 evaluation and one rotation.
     * @alloc none.
     * @test IrbemT01.TotalFieldReportsWhenTheExternalModelDeclines
     */
    [[nodiscard]] Status external_status(const Position<Frame::GEO>& p) const {
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        return t01_field(p_gsm, rotations_->dipole_tilt_deg * (std::numbers::pi / 180.0), drivers_)
            .status;
    }

    /// The drivers this field was built for.
    /// @return the value passed to the constructor.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const T01Drivers& drivers() const { return drivers_; }

    /// The resolved state — what the device staging of a total-field trace would upload.
    /// @return the state, in `double`.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const T01State<double>& state() const { return state_; }

    /// The epoch's frame rotations.
    /// @return the rotations this field was built with.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Rotations& rotations() const { return *rotations_; }

    /// The internal part's Gauss coefficient `g(n, m)`, nT — forwarded for the reason
    /// `TotalFieldT89::g` gives: a superposition has no expansion of its own, and every caller
    /// asking means the internal field's.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double g(int n, int m) const { return internal_->g(n, m); }

    /// The internal part's `h(n, m)`, nT — see @ref g.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double h(int n, int m) const { return internal_->h(n, m); }

    /// The internal field alone.
    /// @return the internal model.
    /// @complexity O(1). @alloc none.
    /// @test IrbemT01.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Igrf<NMAX>& internal() const { return *internal_; }

  private:
    const Igrf<NMAX>* internal_;
    const Rotations* rotations_;
    T01Drivers drivers_;
    T01State<double> state_;
};

}  // namespace cheatah::space::irbem
