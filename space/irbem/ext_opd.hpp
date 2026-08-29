#pragma once

/**
 * @file ext_opd.hpp
 * @brief space.irbem — the Olson-Pfitzer DYNAMIC external field (IRBEM `kext = 6`): a quiet
 *        magnetosphere compressed by the solar-wind dynamic pressure and deepened by the ring
 *        current the pressure-corrected Dst measures.
 *
 * T89 (`ext_t89.hpp`) is a model of the magnetosphere's activity BINS: Kp sorts the data, one
 * coefficient set per bin, and the field is a step function of the driver. This model is the other
 * kind. It is driven by three CONTINUOUS quantities — the solar-wind proton density `n`, the
 * solar-wind speed `V`, and the Dst index — and it answers with a field that is a smooth function
 * of all three. That is the whole reason it exists beside T89: a storm's main phase is a continuous
 * ramp in `n V^2` and Dst, not a staircase of Kp, and a drift shell followed through a ramp must
 * not jump when the driver crosses a bin edge. @ref IrbemOpd.TheFieldIsSmoothInEveryDriver pins
 * the smoothness, and @ref IrbemOpd.SolarWindEntersOnlyThroughDynamicPressure pins the one
 * combination of `n` and `V` the model is allowed to see.
 *
 * ## What is published, and what this file therefore is — MEASURED, NOT ASSUMED
 *
 * The reference IRBEM cites is Pfitzer, Olson & Mogstad, *A time dependent source driven
 * magnetospheric magnetic field model*, EOS Trans. AGU **69**, 426 (1988). **That is a conference
 * abstract.** It states the model's idea and no equation. No later paper published the functional
 * form or a single coefficient; the model exists as distributed code, and IRBEM's `kext = 6`
 * evaluates that code. Its 1977 quiet-time parent (Olson & Pfitzer, *Magnetospheric magnetic
 * field modeling*, AFOSR annual report, McDonnell Douglas, 1977) is likewise a report, not a
 * paper, and this MIT clean room reads neither it nor the LGPL Fortran that carries its
 * coefficients. The provenance protocol's third branch — "form published, coefficients not" — does
 * not even apply; here NEITHER is, and the honest deliverable is different in kind from T89's:
 *
 * 1. **The documented STRUCTURE is implemented, from the published descriptions of the model.**
 *    SPENVIS, *Background: magnetic field models* (spenvis.oma.be/help/background/magfield),
 *    describes the 1988 model in three sentences: the solar-wind pressure sets the scale and
 *    strength of the magnetopause currents, the tail currents are scaled the same way, and the
 *    ring current is driven by "a modified Dst index" with the magnetopause contribution removed.
 *    Jordan, *Empirical models of the magnetospheric magnetic field*, Rev. Geophys. **32**, 139
 *    (1994) §4 says the same. So this file evaluates, in one line,
 *
 *        B(r; n, V, Dst) = s^3 B_q(s r)  +  ΔDst* · R(r),        s = (P / P_0)^(1/6),
 *
 *    a quiet field `B_q` compressed self-similarly by the pressure ratio, plus a symmetric ring
 *    current whose amplitude is the pressure-corrected Dst's departure from its quiet value.
 * 2. **Every constant comes from a published source, and each is cited where it is defined.**
 *    The dynamic pressure `P = m_p n V^2` (@ref opd_pressure_npa); the standoff law
 *    `R_mp ∝ P^(-1/6)` from Chapman-Ferraro pressure balance and the self-similar `s^3 B(s r)`
 *    scaling of a boundary-current field (Mead, JGR **69**, 1181, 1964) (@ref opd_scale); the
 *    pressure correction `Dst* = Dst - 7.26 sqrt(P) + 11` (O'Brien & McPherron, JGR **105**, 7707,
 *    2000) (@ref opd_dst_star); the Dessler-Parker-Sckopke normalisation that makes the ring's
 *    field at the Earth's centre equal to its Dst* increment (@ref opd_ring_coefficient); and the
 *    quiet field itself, which is **Tsyganenko (1989) at its quietest published parameter set**
 *    (Table 1, `Kp = 0, 0+`) — the closest thing to the unpublished 1977 quiet model that a clean
 *    room can hold, already implemented, divergence-free and tilt-aware in `ext_t89.hpp`, and used
 *    here exactly as published with its ring geometry (`a_RC = 8.161`, `D_0 = 2.08 R_E`) reused
 *    for the Dst-driven ring.
 * 3. **The gap to IRBEM's `kext = 6` is measured, not hidden**, by
 *    [`tools/oracle/opd_diff.cpp`](../../../tools/oracle/opd_diff.cpp) against the `-O2` oracle,
 *    with the external field isolated as `kext = 6` minus `kext = 0` so IGRF cancels exactly and
 *    the tilt read from the oracle itself. The same harness establishes, as black-box facts about
 *    the oracle that this implementation reproduces BY CONSTRUCTION, that: `n` and `V` enter it
 *    only through `n V^2` (bit-identical fields for equal products); it is exactly affine in Dst;
 *    its `∂B/∂Dst` does not depend on the pressure; and it refuses (returns `baddata`) anything
 *    strictly outside `5 <= n <= 50`, `300 <= V <= 500`, `-100 <= Dst <= 20`, `r <= 60 R_E` —
 *    which is the envelope `status.hpp` reports, with the difference that this library still
 *    returns the value there. The same harness then asks whether the oracle IS the structure
 *    above, and the answer is: in the tail, yes — on the nightside out to 15 R_E the best fit of
 *    `a B(s r; P_ref)` to the oracle at `n = 10, 20, 40, 50` finds `s = 1.110, 1.260, 1.420,
 *    1.475` against `(P/P_ref)^(1/6) = 1.122, 1.260, 1.414, 1.468`, with `a` within 4% of `s^3`
 *    — and inside 6.5 R_E, no: there the family is nearly a pure amplitude scaling
 *    (`s ≈ 1.00`, `a ≈ (P/P_ref)^0.65`) that no published relation predicts, and over `n V^2` the
 *    family has numerical rank ~4-5. So the structure is right where the boundary and tail
 *    currents live and is not the distributed arithmetic in the inner magnetosphere, and the
 *    difference in the quiet base (T89's bin 1 against the unpublished 1977 model) is on top of
 *    that. Measured, inside the Shue (1997) magnetopause, at three tilts, drivers swept
 *    continuously around each corpus regime clipped to the envelope:
 *
 *    | regime (n, V, Dst) | belts 3-10 R_E: RMS abs | RMS rel | box 3-35 R_E: RMS abs | RMS rel |
 *    |---|---|---|---|---|
 *    | quiet (5, 380, -8) | 15.2 nT | 52% | 25.0 nT | 78% |
 *    | moderate (8, 450, -42) | 26.4 nT | 47% | 34.4 nT | 66% |
 *    | storm, clipped (20, 500, -100) | 49.7 nT | 47% | 64.1 nT | 62% |
 *    | extreme, clipped (45, 500, -100) | 67.7 nT | 51% | 117 nT | 76% |
 *
 *    and along the Dst axis at geosynchronous (`n = 10, V = 450`): RMS 14.6-24.4 nT, 33-43% of
 *    the field, with `∂B_z/∂Dst` at noon 0.33 here against the oracle's 0.30. The floor
 *    experiment — one shared set of the structure's 20 amplitudes (the 19 quiet coefficients
 *    under the compression plus the ring) let float over every regime, tilt and the whole box —
 *    reaches 67.8% of the field against the shipped constants' 72.7%: no re-weighting of this form
 *    gets it close to the distributed model, because the quiet base's SHAPE differs. Read those
 *    numbers before trusting either model at a specific point; the differential test in the suite
 *    asserts them as a REGRESSION envelope, never as agreement.
 *
 * ## Two frames, as in T89
 *
 * The compressed quiet field is T89's and is evaluated in T89's own two frames. The ring is a flat
 * disc in the DIPOLE equator — a symmetric ring current sits on the dipole's equatorial plane, not
 * the GSM one — so it is evaluated in SM and rotated back, the same rotation T89 does for its own
 * ring and tail. GSM in, GSM out.
 *
 * ## Divergence
 *
 * Both pieces are curls: `B_q` of Tsyganenko's vector potentials, `R` of the `A^(3)` disc potential
 * (Tsyganenko 1989, eq. 9). A similarity transform preserves `div B = 0` — for `B'(r) = s^3 B(s r)`,
 * `div B'(r) = s^4 (div B)(s r)` — and the sum of two solenoidal fields is solenoidal, so the
 * second-order stencil residual falls as `h^2` everywhere. @ref IrbemOpd.DivergenceVanishesEverywhere
 * measures it across all four activity regimes and five tilts; it cannot be fooled by a wrong
 * reference because it uses none.
 *
 * @note Nothing on a hot path allocates and nothing but the device lane can throw: the evaluator,
 *       the scalar entry points and the fp32 host lane touch only the caller's spans and one stack
 *       parameter block. The device lane stages `3N` floats each way and forwards whatever
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
#include "ext_t89.hpp"
#include "frames.hpp"
#include "igrf.hpp"
#include "status.hpp"

// The device lane is opt-in by include path, exactly as ext_t89.hpp's is.
#if __has_include("cheatah_gpu_linalg/context.hpp")
#  include "gpu/dispatch.hpp"
/// 1 when this translation unit can reach the OPD device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_OPD_GPU 1
#else
/// 1 when this translation unit can reach the OPD device kernel; 0 when it is host-only.
#  define CHEATAH_SPACE_IRBEM_OPD_GPU 0
#endif

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// The published constants
// -------------------------------------------------------------------------------------------

/**
 * The constants of the dynamic model that are not the quiet field's own, each with its source.
 *
 * Gathered into one `constexpr` object so that every number in the evaluator has a name, a unit
 * and a citation, and so the two arithmetic lanes read one definition.
 *
 * @test IrbemOpd.TheConstantsAreThePublishedOnes
 */
struct OpdConstants {
    /// `m_p` in the units that make `n [cm^-3] V^2 [km^2 s^-2]` come out in nPa: the CODATA 2018
    /// proton mass `1.67262192369e-27 kg` times `1e6 m^-3 / cm^-3` times `1e6 m^2 s^-2 / km^2 s^-2`
    /// times `1e9 nPa / Pa`. The solar wind is treated as protons, which is what the IRBEM driver
    /// table calls "solar wind density".
    double pressure_npa_per_cc_km2s2;
    /// The reference density, cm^-3: the quiet end of the model's documented envelope
    /// (`5 <= Dsw <= 50`, IRBEM `general_information.rst`), at which the quiet field is
    /// uncompressed.
    double reference_density;
    /// The reference speed, km/s: the canonical nominal solar wind, in the middle of the documented
    /// `300 <= Vsw <= 500`.
    double reference_velocity;
    /// `b` in `Dst* = Dst - b sqrt(P) + c`, nT per sqrt(nPa) — O'Brien & McPherron, JGR 105,
    /// 7707 (2000), the pressure correction that removes the magnetopause-current contribution
    /// from the index, which is exactly what the model's "modified Dst" is described as.
    double dst_pressure_coefficient;
    /// `c` in the same relation, nT — the quiet-time offset of the same paper.
    double dst_offset;
};

/**
 * The constants, as their sources state them.
 * @test IrbemOpd.TheConstantsAreThePublishedOnes
 */
inline constexpr OpdConstants opd_constants{
    /* pressure_npa_per_cc_km2s2 */ 1.67262192369e-6,
    /* reference_density */ 5.0,
    /* reference_velocity */ 400.0,
    /* dst_pressure_coefficient */ 7.26,
    /* dst_offset */ 11.0};

/// The T89 Kp bin whose published parameter set is the quiet field: `Kp = 0, 0+`, Table 1's first
/// column. One bin, always, so nothing in this model is ever a step function of anything.
inline constexpr int opd_quiet_bin = 1;

/**
 * The solar-wind dynamic pressure, in nPa.
 *
 * `P = m_p n V^2`, the ram pressure of a proton wind. This is the ONLY combination of the two
 * solar-wind drivers the model reads — a fact about the published idea (Chapman-Ferraro pressure
 * balance knows nothing but `P`) that the oracle turns out to share bit for bit.
 *
 * @param density_cc the proton density, cm^-3.
 * @param velocity_kms the bulk speed, km/s.
 * @return the dynamic pressure, nPa; 1.338 nPa at the reference `(5, 400)`.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpd.DynamicPressureIsTheProtonRamPressure
 */
[[nodiscard]] constexpr double opd_pressure_npa(double density_cc, double velocity_kms) {
    return opd_constants.pressure_npa_per_cc_km2s2 * density_cc * velocity_kms * velocity_kms;
}

/**
 * The reference pressure `P_0`, at which the quiet field is used uncompressed.
 * @return `opd_pressure_npa(5, 400)`, ~1.338 nPa.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpd.ReferenceConditionsReduceToTheQuietField
 */
[[nodiscard]] constexpr double opd_reference_pressure_npa() {
    return opd_pressure_npa(opd_constants.reference_density, opd_constants.reference_velocity);
}

/**
 * The compression ratio `s = R_0 / R_mp = (P / P_0)^(1/6)`.
 *
 * Pressure balance at the subsolar magnetopause, `P = B_mp^2 / 2 mu_0` with `B_mp ∝ R_mp^-3`,
 * puts the boundary at `R_mp ∝ P^(-1/6)` (Chapman & Ferraro 1931; the textbook statement is
 * Kivelson & Russell, *Introduction to Space Physics*, 1995, §6.2). A self-similar compression of
 * the whole boundary-current system by the factor `s` then scales its field as `s^3 B(s r)` (Mead,
 * JGR 69, 1181, 1964, whose boundary-field coefficients all carry `R_mp^-3`). This function is the
 * `s`; the evaluator applies the transform.
 *
 * @param density_cc the proton density, cm^-3; must be positive.
 * @param velocity_kms the bulk speed, km/s; must be positive.
 * @return `s`, exactly 1 at the reference conditions, 1.58 at the envelope's most compressed
 *         corner `(50, 500)`, 0.86 at its least `(5, 300)`.
 * @complexity O(1) — one `pow`.
 * @alloc none.
 * @test IrbemOpd.CompressionFollowsTheSixthRootOfPressure
 * @test IrbemOpd.CompressionScalesTheCentralFieldAsTheCubeOfTheStandoff
 */
[[nodiscard]] inline double opd_scale(double density_cc, double velocity_kms) {
    return std::pow(opd_pressure_npa(density_cc, velocity_kms) / opd_reference_pressure_npa(),
                    1.0 / 6.0);
}

/**
 * The pressure-corrected Dst, `Dst* = Dst - b sqrt(P) + c`.
 *
 * O'Brien & McPherron (2000), with `b = 7.26 nT/sqrt(nPa)` and `c = 11 nT`: the index with the
 * magnetopause currents' contribution removed, leaving the part the ring current is responsible
 * for. This is the "modified Dst" the model's description names as its ring-current driver.
 *
 * @param dst_nt the Dst index, nT.
 * @param pressure_npa the dynamic pressure, nPa; must be non-negative.
 * @return `Dst*`, nT.
 * @complexity O(1) — one `sqrt`.
 * @alloc none.
 * @test IrbemOpd.DstIsPressureCorrectedBeforeItDrivesTheRing
 */
[[nodiscard]] inline double opd_dst_star(double dst_nt, double pressure_npa) {
    return dst_nt - (opd_constants.dst_pressure_coefficient * std::sqrt(pressure_npa)) +
           opd_constants.dst_offset;
}

/**
 * The ring current's Dst* increment over the quiet field's own, `ΔDst* = Dst*(P, Dst) - Dst*(P_0, 0)`.
 *
 * The quiet field already contains a quiet ring current (T89's `C_5` term), so what the Dst driver
 * adds is the DIFFERENCE from quiet — which is what makes the model reduce to the quiet field
 * exactly at the reference conditions rather than double-counting a ring. The offset `c` cancels
 * in the difference; only `Dst` and the pressure correction survive.
 *
 * @param density_cc the proton density, cm^-3.
 * @param velocity_kms the bulk speed, km/s.
 * @param dst_nt the Dst index, nT.
 * @return `Dst - 7.26 (sqrt(P) - sqrt(P_0))`, nT; zero at `(5, 400, 0)`.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpd.DstIsPressureCorrectedBeforeItDrivesTheRing
 */
[[nodiscard]] inline double opd_ring_increment(double density_cc, double velocity_kms,
                                               double dst_nt) {
    return opd_dst_star(dst_nt, opd_pressure_npa(density_cc, velocity_kms)) -
           opd_dst_star(0.0, opd_reference_pressure_npa());
}

/**
 * The amplitude of the `A^(3)` ring disc whose field at the Earth's centre is `delta_dst_star`.
 *
 * The ring is the finite-moment disc potential of Tsyganenko (1989) eq. (9), with radius `a` and
 * half-thickness `D`. On its axis at the origin its field is `B_z(0) = 2 C / (a + D)^3` (set
 * `rho = z = 0` in eq. 16-17's `B_z = C (2 u^2 - rho^2) / S^5` with `u = a + D`, `S = u`) — `z`
 * being the DIPOLE axis, since the ring lies in the dipole equator; in GSM the central field
 * points along `(sin psi, 0, cos psi)`. The Dessler-Parker-Sckopke relation (Dessler & Parker,
 * JGR 64, 2239, 1959; Sckopke, JGR 71, 3125, 1966) says the symmetric ring's field at the centre
 * IS the pressure-corrected Dst, so the amplitude that makes `B_z(0) = ΔDst*` is
 * `C = ΔDst* (a + D)^3 / 2`.
 *
 * @param delta_dst_star the ring's Dst* increment, nT; see @ref opd_ring_increment.
 * @param a the disc's radial scale, R_E — T89's quiet `a_RC = 8.161`.
 * @param d the disc's half-thickness, R_E — T89's quiet `D_0 = 2.08`.
 * @return `C`, in the units that make eq. (16)-(17) come out in nT with positions in R_E.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpd.TheRingIsNormalisedToDstAtTheCentre
 */
[[nodiscard]] constexpr double opd_ring_coefficient(double delta_dst_star, double a, double d) {
    const double u = a + d;
    return 0.5 * delta_dst_star * u * u * u;
}

// -------------------------------------------------------------------------------------------
// The parameter block
// -------------------------------------------------------------------------------------------

/**
 * Everything the evaluator needs for one epoch and one driver triple, in the lane's scalar type.
 *
 * As with @ref T89Parameters, the `float` instantiation is not a convenience: the device kernel
 * receives these as `float`, so the host lane verified against it must round them FIRST — the
 * scale, the ring coefficient and the quiet set alike — and then evaluate, or a disagreement between
 * the two lanes cannot be attributed to the device. @ref opd_parameters builds both.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` for the device-mirroring one.
 * @test IrbemOpd.ParametersRoundTripThroughFloat
 */
template <std::floating_point T>
struct OpdParameters {
    /// The quiet field's published set — T89 Table 1, `Kp = 0, 0+` — whose `a_rc` and `d_0` also
    /// fix the Dst-driven ring's geometry.
    T89Parameters<T> quiet;
    /// The compression ratio `s`; see @ref opd_scale.
    T scale;
    /// The ring disc's amplitude `C`; see @ref opd_ring_coefficient.
    T ring;
};

/**
 * The parameters for one driver triple, computed in `double` and rounded once to @p T.
 *
 * Computing in `double` and rounding the RESULT is what keeps the `float` block the nearest float
 * to the true parameters rather than the product of float intermediates: `s` is a sixth root and
 * `C` a cube, and both amplify a rounding in their argument.
 *
 * @tparam T the scalar type; `double` or `float`.
 * @param density_cc the proton density, cm^-3; positive, which @ref opd_field has already checked.
 * @param velocity_kms the bulk speed, km/s; positive, likewise.
 * @param dst_nt the Dst index, nT.
 * @return the block, by value.
 * @complexity O(1) — one `pow`, two `sqrt`, 30 conversions.
 * @alloc none; the returned object is inline storage.
 * @test IrbemOpd.ParametersRoundTripThroughFloat
 */
template <std::floating_point T>
[[nodiscard]] inline OpdParameters<T> opd_parameters(double density_cc, double velocity_kms,
                                                     double dst_nt) {
    const T89Coefficients& quiet = t89_coefficient_sets[static_cast<std::size_t>(opd_quiet_bin - 1)];
    OpdParameters<T> out{};
    out.quiet = t89_parameters<T>(opd_quiet_bin);
    out.scale = static_cast<T>(opd_scale(density_cc, velocity_kms));
    out.ring = static_cast<T>(opd_ring_coefficient(
        opd_ring_increment(density_cc, velocity_kms, dst_nt), quiet.a_rc, quiet.d_0));
    return out;
}

// -------------------------------------------------------------------------------------------
// The evaluator
// -------------------------------------------------------------------------------------------

/**
 * The Dst-driven ring's field at one GSM point, per unit amplitude, in nanotesla.
 *
 * A flat, symmetric `A^(3)` disc in the dipole equator: Tsyganenko (1989) eqs. (16)-(17) with the
 * sheet unwarped (`z_s = 0`) and its thickness constant, so every `dz_s/dx`, `dz_s/dy` and `dD/dx`
 * term of the published form is identically zero and what is left is
 *
 *     q = 3 u / (xi S^5),  B_x = q x z,  B_y = q y z,  B_z = (2 u^2 - rho^2) / S^5,
 *
 * with `xi = sqrt(z^2 + D^2)`, `u = a + xi`, `S = sqrt(rho^2 + u^2)`, all in SM, and the
 * result rotated back to GSM. Split out of @ref opd_components so that the unit shape — the
 * model's `∂B/∂Dst` — can be tested and normalised on its own.
 *
 * @tparam T the scalar type.
 * @param a the disc's radial scale, R_E.
 * @param d the disc's half-thickness, R_E.
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`.
 * @param x the GSM x coordinate, R_E.
 * @param y the GSM y coordinate, R_E.
 * @param z the GSM z coordinate, R_E.
 * @return `{B_x, B_y, B_z}` in GSM per unit `C`; multiply by @ref opd_ring_coefficient's `C`.
 * @complexity O(1) — about 40 flops and three square roots.
 * @alloc none.
 * @test IrbemOpd.TheRingIsNormalisedToDstAtTheCentre
 * @test IrbemOpd.TheDstGradientIsTheRingShapeAtEveryPressure
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, 3> opd_ring_unit(T a, T d, T sin_tilt, T cos_tilt, T x, T y,
                                                    T z) {
    const T two = static_cast<T>(2);
    const T three = static_cast<T>(3);
    // GSM -> SM: the ring lies in the dipole equator.
    const T xs = (x * cos_tilt) - (z * sin_tilt);
    const T ys = y;
    const T zs = (x * sin_tilt) + (z * cos_tilt);
    const T rho2 = (xs * xs) + (ys * ys);
    const T xi = std::sqrt((zs * zs) + (d * d));
    const T u = a + xi;
    const T s = std::sqrt(rho2 + (u * u));
    const T s2 = s * s;
    const T s5 = s2 * s2 * s;
    const T q = three * u / (xi * s5);
    const T bx_sm = q * xs * zs;
    const T by_sm = q * ys * zs;
    const T bz_sm = ((two * u * u) - rho2) / s5;
    // SM -> GSM: the transpose.
    return {(bx_sm * cos_tilt) + (bz_sm * sin_tilt), by_sm,
            (-bx_sm * sin_tilt) + (bz_sm * cos_tilt)};
}

/**
 * The Olson-Pfitzer dynamic external field at one GSM point, as three components in nanotesla.
 *
 * The one line of the file brief: `s^3 B_q(s r) + C R(r)`. The quiet field is T89's evaluator on
 * the SCALED position — a point at `r` in a magnetosphere compressed by `s` sees what the quiet
 * magnetosphere has at `s r` — with the `s^3` that keeps the transform a similarity of the
 * currents; the ring is @ref opd_ring_unit at the UNSCALED position, because the Dst-driven
 * increment is a property of the inner magnetosphere and not of the boundary. Both are
 * solenoidal, so their sum is.
 *
 * @tparam T the scalar type; `double` for the reference lane, `float` to mirror the device kernel.
 * @param p the parameters for the epoch's drivers; see @ref opd_parameters.
 * @param sin_tilt `sin(psi)`, precomputed per epoch.
 * @param cos_tilt `cos(psi)`; must be non-zero, which @ref opd_field checks before it gets here.
 * @param x the GSM x coordinate, R_E.
 * @param y the GSM y coordinate, R_E.
 * @param z the GSM z coordinate, R_E.
 * @return `{B_x, B_y, B_z}` in GSM, nanotesla.
 * @complexity O(1) — T89's ~400 flops plus ~50; no loop, no branch on data.
 * @alloc none.
 * @test IrbemOpd.DivergenceVanishesEverywhere
 * @test IrbemOpd.ReferenceConditionsReduceToTheQuietField
 * @test IrbemOpd.CompressionScalesTheCentralFieldAsTheCubeOfTheStandoff
 * @test IrbemOpd.DawnDuskSymmetryHoldsAtEveryTilt
 */
template <std::floating_point T>
[[nodiscard]] inline std::array<T, 3> opd_components(const OpdParameters<T>& p, T sin_tilt,
                                                     T cos_tilt, T x, T y, T z) {
    const T s = p.scale;
    const T s3 = s * s * s;
    const std::array<T, 3> quiet =
        t89_components<T>(p.quiet, sin_tilt, cos_tilt, s * x, s * y, s * z);
    const std::array<T, 3> ring =
        opd_ring_unit<T>(p.quiet.a_rc, p.quiet.d_0, sin_tilt, cos_tilt, x, y, z);
    return {(s3 * quiet[0]) + (p.ring * ring[0]), (s3 * quiet[1]) + (p.ring * ring[1]),
            (s3 * quiet[2]) + (p.ring * ring[2])};
}

/**
 * The dynamic external field at one GSM point, in `double` — the reference lane.
 *
 * @param p the position, GSM, in Earth radii.
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`; must be non-zero, which @ref opd_field checks.
 * @param density_cc the solar-wind proton density, cm^-3; positive.
 * @param velocity_kms the solar-wind speed, km/s; positive.
 * @param dst_nt the Dst index, nT.
 * @return the external field at @p p, GSM, in nanotesla.
 * @complexity O(1); see @ref opd_components, plus @ref opd_parameters once.
 * @alloc none.
 * @test IrbemOpd.ReferenceLaneMatchesTheComponentForm
 */
[[nodiscard]] inline FieldVector<Frame::GSM> opd_field_at(Position<Frame::GSM> p, double sin_tilt,
                                                          double cos_tilt, double density_cc,
                                                          double velocity_kms, double dst_nt) {
    const std::array<double, 3> b =
        opd_components<double>(opd_parameters<double>(density_cc, velocity_kms, dst_nt), sin_tilt,
                               cos_tilt, p.v[0], p.v[1], p.v[2]);
    return FieldVector<Frame::GSM>{fixarray::vec3d{b[0], b[1], b[2]}};
}

/**
 * The dynamic external field, with the model's own verdict on whether it should be believed here.
 *
 * The value is **always** returned, including under @ref Status::OutOfValidityRange — `status.hpp`'s
 * standing rule. The oracle, by contrast, returns `baddata` strictly outside the same envelope
 * (measured: `n = 50.0001`, `V = 500.001`, `Dst = 20.001` and `-100.001`, `r = 60.05` all refuse;
 * the bounds themselves are accepted), and a caller porting from it should expect a number with a
 * caveat where it got a sentinel.
 *
 * Refused outright, as @ref Status::DomainError with a zero field, is arithmetic with no answer:
 * a non-finite input; a tilt of `|psi| >= pi/2`, at which the quiet field's `tan(psi)` does not
 * exist; a radius inside the Earth; a density or speed that is not positive, for which there is no
 * pressure and no magnetopause to compress (`s` would be zero or complex); and a non-finite answer,
 * which the quiet field's `exp(x / dx)` produces far enough sunward that only a unit-confusion bug
 * gets there, and which must not be handed to a tracer as a NaN.
 *
 * @param p the position, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians.
 * @param density_cc the solar-wind proton density, cm^-3 — IRBEM's `maginput(3)`.
 * @param velocity_kms the solar-wind speed, km/s — `maginput(4)`.
 * @param dst_nt the Dst index, nT — `maginput(2)`.
 * @return the field and its caveat: @ref Status::DomainError (zero field) for the refusals above;
 *         @ref Status::OutOfValidityRange for a driver or a radius outside the documented envelope,
 *         with the field still computed; otherwise @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpd.ValidityIsReportedFromBothSidesOfEveryBound
 * @test IrbemOpd.NonFiniteInputIsADomainError
 * @test IrbemOpd.AVanishingSolarWindIsADomainError
 * @test IrbemOpd.RightAngleTiltIsADomainError
 * @test IrbemOpd.AnOverflowingExtrapolationIsADomainErrorNotANaN
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> opd_field(Position<Frame::GSM> p,
                                                               double tilt_rad, double density_cc,
                                                               double velocity_kms, double dst_nt) {
    const FieldVector<Frame::GSM> zero{};
    if (!std::isfinite(p.v[0]) || !std::isfinite(p.v[1]) || !std::isfinite(p.v[2]) ||
        !std::isfinite(tilt_rad) || !std::isfinite(density_cc) || !std::isfinite(velocity_kms) ||
        !std::isfinite(dst_nt)) {
        return {Status::DomainError, zero};
    }
    if (!(std::fabs(tilt_rad) < max_tilt_rad)) return {Status::DomainError, zero};
    if (!(density_cc > 0.0) || !(velocity_kms > 0.0)) return {Status::DomainError, zero};

    const double r = std::sqrt((p.v[0] * p.v[0]) + (p.v[1] * p.v[1]) + (p.v[2] * p.v[2]));
    const Status where = check_position(ExternalModel::OlsonPfitzerDynamic1988, r, p.v[0]);
    if (where == Status::DomainError) return {Status::DomainError, zero};

    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Dsw)] = density_cc;
    drivers[static_cast<std::size_t>(Driver::Vsw)] = velocity_kms;
    drivers[static_cast<std::size_t>(Driver::Dst)] = dst_nt;
    const Status drives = check_validity(ExternalModel::OlsonPfitzerDynamic1988, drivers);

    const FieldVector<Frame::GSM> b = opd_field_at(p, std::sin(tilt_rad), std::cos(tilt_rad),
                                                   density_cc, velocity_kms, dst_nt);
    if (!std::isfinite(b.v[0]) || !std::isfinite(b.v[1]) || !std::isfinite(b.v[2])) {
        return {Status::DomainError, zero};
    }
    return {first_failure(drives, where), b};
}

/**
 * The dynamic external field for a whole epoch's worth of state — the production entry point.
 *
 * The tilt and Dst come out of @ref HotState; the density and speed are the two `maginput` slots
 * the hot block does not carry, so they are read through @ref FieldContext::driver. The tilt ANGLE is
 * what is checked, as in T89, because the `|psi| < pi/2` requirement is strict and the context
 * only guarantees `<=`.
 *
 * @param p the position, GSM, in Earth radii.
 * @param ctx the epoch's context.
 * @return the field and its caveat, exactly as the five-argument overload.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpd.ContextOverloadAgreesWithTheExplicitOne
 */
[[nodiscard]] inline Result<FieldVector<Frame::GSM>> opd_field(Position<Frame::GSM> p,
                                                               const FieldContext& ctx) {
    return opd_field(p, ctx.hot().tilt_rad, ctx.driver(Driver::Dsw), ctx.driver(Driver::Vsw),
                     ctx.hot().dst);
}

// -------------------------------------------------------------------------------------------
// The batch lanes
// -------------------------------------------------------------------------------------------

/**
 * The dynamic field over a whole batch, on the CPU, in `float`.
 *
 * The host twin of `irbem_opd_f32`: the same expressions, in the same order, in the same precision,
 * from a parameter block rounded to `float` FIRST. What makes a host-vs-device disagreement
 * attributable to the device.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, GSM, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, nanotesla; overwritten in full.
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`; must be non-zero.
 * @param p the parameter block, already rounded; see @ref opd_parameters.
 * @return `false` when @p pos is not a whole number of points or @p out is a different length, in
 *         which case nothing is written; `true` otherwise.
 * @complexity O(N).
 * @alloc none.
 * @test IrbemOpd.HostFloatLaneTracksTheReferenceLane
 * @test IrbemOpd.HostFloatLaneRejectsMismatchedSpans
 */
[[nodiscard]] inline bool opd_field_host(std::span<const float> pos, std::span<float> out,
                                         float sin_tilt, float cos_tilt,
                                         const OpdParameters<float>& p) {
    if (pos.size() % 3 != 0 || out.size() != pos.size()) return false;
    const std::size_t n = pos.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const std::array<float, 3> b = opd_components<float>(
            p, sin_tilt, cos_tilt, pos[(3 * i) + 0], pos[(3 * i) + 1], pos[(3 * i) + 2]);
        out[(3 * i) + 0] = b[0];
        out[(3 * i) + 1] = b[1];
        out[(3 * i) + 2] = b[2];
    }
    return true;
}

/// How many `float` scalars the device kernel's parameter buffer holds: T89's thirty (tilt sine
/// and cosine, then the quiet set) followed by `s` and `C`. Asserted against the kernel registry.
inline constexpr std::size_t opd_param_count = t89_param_count + 2;

/**
 * Pack the epoch's tilt and the drivers' parameters into the kernel's parameter buffer.
 *
 * The first thirty floats are EXACTLY @ref t89_param_block for the quiet bin, because the kernel
 * hands them to the shared `t89_eval` untouched; the last two are the scale and the ring
 * amplitude. Stated here and above `irbem_opd_f32` in `irbem.slang`, and kept in step by a test
 * that runs both lanes on the same points.
 *
 * @param sin_tilt `sin(psi)`.
 * @param cos_tilt `cos(psi)`.
 * @param p the parameter block, rounded to `float`.
 * @return the buffer, `opd_param_count` floats, by value.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemOpd.ParameterBlockIsTheQuietBlockThenTheScaleAndTheRing
 */
[[nodiscard]] inline std::array<float, opd_param_count> opd_param_block(
    float sin_tilt, float cos_tilt, const OpdParameters<float>& p) {
    const std::array<float, t89_param_count> quiet = t89_param_block(sin_tilt, cos_tilt, opd_quiet_bin);
    std::array<float, opd_param_count> block{};
    for (std::size_t k = 0; k < t89_param_count; ++k) block[k] = quiet[k];
    block[t89_param_count] = p.scale;
    block[t89_param_count + 1] = p.ring;
    return block;
}

/**
 * The batch's position caveat, accumulated one point at a time — @ref T89PositionFold's method
 * against this model's envelope.
 *
 * The envelope is `r <= 60 R_E` and nothing else, so the largest and smallest squared radii decide
 * the whole batch with no per-point `sqrt`; the finiteness flag is kept apart because a NaN
 * radius compares false against both extremes.
 *
 * @test IrbemOpd.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
struct OpdPositionFold {
    /// The smallest `r^2` seen, R_E^2; `+inf` until the first point.
    double r2_lo = std::numeric_limits<double>::infinity();
    /// The largest `r^2` seen, R_E^2; zero until the first point.
    double r2_hi = 0.0;
    /// The smallest GSM `x` seen, R_E — carried so the envelope check is the same call the scalar
    /// lane makes, though this model publishes no `x` limit.
    double x_lo = std::numeric_limits<double>::infinity();
    /// False once any point has had a non-finite coordinate.
    bool finite = true;

    /**
     * Fold one position in.
     * @param p the position, GSM, in Earth radii.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemOpd.BatchReportsTheSameEnvelopeTheScalarLaneDoes
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
     *         @ref Status::OutOfValidityRange when any point is beyond `60 R_E`, otherwise
     *         @ref Status::Ok.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemOpd.BatchReportsTheSameEnvelopeTheScalarLaneDoes
     */
    [[nodiscard]] Status verdict() const {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double lo = finite ? std::sqrt(r2_lo) : nan;
        const double hi = finite ? std::sqrt(r2_hi) : nan;
        return first_failure(check_position(ExternalModel::OlsonPfitzerDynamic1988, lo, x_lo),
                             check_position(ExternalModel::OlsonPfitzerDynamic1988, hi, x_lo));
    }
};

/**
 * The dynamic field over a whole batch of GSM points, on the device when that is worth it.
 *
 * The shape of @ref t89_field_batch exactly: one status for N points, the worst of them; an
 * out-of-validity batch computed in full and a domain-error batch zeroed in full; the returned
 * value `true` exactly when the device serviced the call. The kernel is T89's plus a ring, at the
 * same ~11 flops/byte, and runs the SAME shared `t89_eval` — so the device wins by the same margin
 * and the crossover is the T89 row's (see `gpu/dispatch.hpp`).
 *
 * @param points the positions, GSM, in Earth radii.
 * @param tilt_rad the dipole tilt `psi`, radians.
 * @param density_cc the solar-wind proton density, cm^-3.
 * @param velocity_kms the solar-wind speed, km/s.
 * @param dst_nt the Dst index, nT.
 * @param out receives one field vector per input, GSM, nanotesla; same length as @p points.
 * @return @ref Status::DomainError on a length mismatch, a non-finite or non-positive driver, a
 *         tilt at which the model has no value, or a point that is not finite or is inside the
 *         Earth — every output then zeroed; @ref Status::OutOfValidityRange when a driver or any
 *         point's radius is outside the documented envelope, every point still computed; otherwise
 *         @ref Status::Ok. The value is `true` exactly when the device lane serviced the call.
 * @complexity O(N); on the device those N run concurrently over `ceil(N/256)` workgroups.
 * @alloc the device lane stages positions and results into two `std::vector<float>` of `3N`; the
 *        host lane allocates nothing.
 * @test IrbemOpd.BatchAgreesWithTheReferenceLane
 * @test IrbemOpd.BatchRejectsMismatchedSpans
 * @test IrbemOpd.BatchReportsTheSameEnvelopeTheScalarLaneDoes
 */
[[nodiscard]] inline Result<bool> opd_field_batch(std::span<const Position<Frame::GSM>> points,
                                                  double tilt_rad, double density_cc,
                                                  double velocity_kms, double dst_nt,
                                                  std::span<FieldVector<Frame::GSM>> out) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(tilt_rad) || !std::isfinite(density_cc) || !std::isfinite(velocity_kms) ||
        !std::isfinite(dst_nt)) {
        return {Status::DomainError, false};
    }
    if (!(std::fabs(tilt_rad) < max_tilt_rad)) return {Status::DomainError, false};
    if (!(density_cc > 0.0) || !(velocity_kms > 0.0)) return {Status::DomainError, false};

    DriverSet drivers{};
    drivers[static_cast<std::size_t>(Driver::Dsw)] = density_cc;
    drivers[static_cast<std::size_t>(Driver::Vsw)] = velocity_kms;
    drivers[static_cast<std::size_t>(Driver::Dst)] = dst_nt;
    const Status drives = check_validity(ExternalModel::OlsonPfitzerDynamic1988, drivers);
    if (n == 0) return {drives, false};

    const double sin_tilt = std::sin(tilt_rad);
    const double cos_tilt = std::cos(tilt_rad);
    OpdPositionFold fold;

#if CHEATAH_SPACE_IRBEM_OPD_GPU
    if (gpu::prefer_gpu("irbem_opd_f32", n) &&
        std::filesystem::exists(gpu::shader_path("irbem_opd_f32"))) {
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
        const std::array<float, opd_param_count> block =
            opd_param_block(static_cast<float>(sin_tilt), static_cast<float>(cos_tilt),
                            opd_parameters<float>(density_cc, velocity_kms, dst_nt));
        gpu::dispatch_batch("irbem_opd_f32", pos, raw, std::span<const float>(block));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = FieldVector<Frame::GSM>{
                fixarray::vec3d{raw[(3 * i) + 0], raw[(3 * i) + 1], raw[(3 * i) + 2]}};
        }
        return {first_failure(drives, where), true};
    }
#endif

    const OpdParameters<double> p = opd_parameters<double>(density_cc, velocity_kms, dst_nt);
    for (std::size_t i = 0; i < n; ++i) {
        fold.add(points[i]);
        const std::array<double, 3> b = opd_components<double>(
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

// -------------------------------------------------------------------------------------------
// The total field — IGRF plus this model, as one GeoFieldModel
// -------------------------------------------------------------------------------------------

/**
 * IGRF plus the Olson-Pfitzer dynamic field, as a single field a tracer can follow.
 *
 * The same shape as `TotalFieldT89`, for the same reason: `trace_invariant`, `make_lstar` and
 * everything above them take one @ref GeoFieldModel and follow the TOTAL field without knowing
 * there are two. The frames are handled once, through the epoch's @ref Rotations; the external
 * model's refusals fall back to the internal field alone (the best available answer there, and
 * never a NaN into an integrator), and @ref external_status says when that happened.
 *
 * @tparam NMAX the internal field's truncation degree; 10 reproduces IRBEM's own choice.
 * @test IrbemOpd.TotalFieldSuperposesInternalAndExternal
 * @test IrbemOpd.TotalFieldTracesAndReportsWhenTheExternalModelDeclines
 */
template <int NMAX = 10>
class TotalFieldOpd {
  public:
    /// The internal part's truncation degree — what generic staging and buffer sizing read.
    static constexpr int degree = NMAX;

    /**
     * @param internal the internal field, already built for the epoch.
     * @param rotations the epoch's frame rotations — built once, reused for every point.
     * @param density_cc the solar-wind proton density, cm^-3.
     * @param velocity_kms the solar-wind speed, km/s.
     * @param dst_nt the Dst index, nT.
     * @test IrbemOpd.TotalFieldSuperposesInternalAndExternal
     */
    constexpr TotalFieldOpd(const Igrf<NMAX>& internal, const Rotations& rotations,
                            double density_cc, double velocity_kms, double dst_nt)
        : internal_(&internal),
          rotations_(&rotations),
          density_cc_(density_cc),
          velocity_kms_(velocity_kms),
          dst_nt_(dst_nt) {}

    /**
     * The total field at a geographic point.
     * @param p the position, GEO, Earth radii.
     * @return `B_internal + B_external` in GEO, nT; the internal field alone when the external
     *         model refuses the point (see @ref external_status).
     * @complexity One IGRF evaluation, one OPD evaluation, two 3x3 rotations.
     * @alloc none.
     * @test IrbemOpd.TotalFieldSuperposesInternalAndExternal
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::GEO>& p) const {
        const FieldVector<Frame::GEO> b_int = internal_->evaluate(p);
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        const double tilt_rad = rotations_->dipole_tilt_deg * (std::numbers::pi / 180.0);
        const Result<FieldVector<Frame::GSM>> b_ext =
            opd_field(p_gsm, tilt_rad, density_cc_, velocity_kms_, dst_nt_);
        if (b_ext.status == Status::DomainError) return b_int;
        const FieldVector<Frame::GEO> b_ext_geo = transform<Frame::GEO>(b_ext.value, *rotations_);
        return FieldVector<Frame::GEO>{b_int.v + b_ext_geo.v};
    }

    /**
     * Whether the external model answered at @p p, and if not, why.
     * @param p the position, GEO, Earth radii.
     * @return the external model's status; @ref Status::Ok when it contributed without caveat.
     * @complexity One OPD evaluation and one rotation.
     * @alloc none.
     * @test IrbemOpd.TotalFieldTracesAndReportsWhenTheExternalModelDeclines
     */
    [[nodiscard]] Status external_status(const Position<Frame::GEO>& p) const {
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        const double tilt_rad = rotations_->dipole_tilt_deg * (std::numbers::pi / 180.0);
        return opd_field(p_gsm, tilt_rad, density_cc_, velocity_kms_, dst_nt_).status;
    }

    /// The epoch's frame rotations.
    /// @return the rotations this field was built with.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOpd.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Rotations& rotations() const { return *rotations_; }

    /// The internal part's Gauss coefficient `g(n, m)`, nT — the internal field's, for the same
    /// reason `TotalFieldT89::g` forwards: a superposition has no harmonic expansion of its own.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOpd.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double g(int n, int m) const { return internal_->g(n, m); }

    /// The internal part's `h(n, m)`, nT.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOpd.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr double h(int n, int m) const { return internal_->h(n, m); }

    /// The internal field alone.
    /// @return the internal model.
    /// @complexity O(1). @alloc none.
    /// @test IrbemOpd.TotalFieldSuperposesInternalAndExternal
    [[nodiscard]] constexpr const Igrf<NMAX>& internal() const { return *internal_; }

  private:
    const Igrf<NMAX>* internal_;
    const Rotations* rotations_;
    double density_cc_;
    double velocity_kms_;
    double dst_nt_;
};

}  // namespace cheatah::space::irbem
