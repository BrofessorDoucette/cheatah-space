#pragma once

/**
 * @file field.hpp
 * @brief space.irbem — the magnetic field at a point, its first derivatives, and everything
 *        algebraic that follows from them.
 *
 * @ref lstar.hpp is the compute-bound half of this module: one point in, ~10⁵ field evaluations,
 * six scalars out. This file is the other half — the **streaming** half, where the answer is as
 * large as the question and the work per point is a few hundred flops. Four IRBEM routines live
 * here, and they are ordered by how much arithmetic they do per byte moved, because that ordering
 * is the whole design:
 *
 * | routine | field evaluations / point | flops / byte | lane |
 * |---|---|---|---|
 * | `GET_FIELD_MULTI` (@ref field_batch) | 1 | ~20 | device above 128 points, **measured 8.96×** |
 * | `GET_HEMI_MULTI` (@ref hemisphere_batch) | 3 | ~60 | device above 43 points |
 * | `GET_BDERIVS` (@ref bderivs_batch) | 4 | ~80 | device above 32 points |
 * | `COMPUTE_GRAD_CURV_CURL` (@ref grad_curv_curl_batch) | **0** | ~0.4 | **host, always** |
 *
 * The last row is the one worth stating out loud. `COMPUTE_GRAD_CURV_CURL` touches no field model
 * at all — it is pure algebra over `GET_BDERIVS`'s outputs — and it is therefore the one routine
 * in this module a GPU can only make slower. Measured: 136 bytes per point cross the bus (16 fp32
 * in, 18 out) for **11.6 ns of host work**; at this machine's ~12 GB/s of effective PCIe bandwidth
 * the transfer alone is ~11.3 ns/point, so the copy costs what the computation costs *before* the
 * ~30 µs submit floor is paid and before the kernel runs. There is no batch size that recovers
 * that: the ratio is fixed, not amortizable. @ref grad_curv_curl_batch is a host loop, deliberately
 * and permanently, and no kernel for it exists in `gpu/irbem.slang`.
 *
 * ## The finite-difference step is the whole design of GET_BDERIVS
 *
 * `GET_BDERIVS` has no closed form here. IGRF's Jacobian *is* analytically available — the Legendre
 * recursion carries `A'ⁿₘ` already — but IRBEM's routine is defined as a **finite difference with a
 * caller-supplied step `dX`**, and reproducing it means differencing. That puts the accuracy of
 * every derivative in this file on one number, pulled in two directions:
 *
 *  - **too large** and truncation dominates. A one-sided difference has error `(h/2)·|∂²B|`, and
 *    for a field falling as `r⁻³` that is a relative error of about `2h/r` on `∇|B|`. Measured on
 *    IGRF-14 against a Richardson-extrapolated reference, over 60 points spread across `r = 1.05
 *    … 11.5 Re`: `h = 10⁻³ Re` costs `1.2 × 10⁻³`, `h = 10⁻⁴` costs `1.2 × 10⁻⁴`. Dead linear, as
 *    the theory says.
 *  - **too small** and cancellation destroys it. The difference of two nearby field magnitudes
 *    loses every digit the two share, so the roundoff term grows as `ε·r/h` — it *diverges* as the
 *    step shrinks. In fp64 that floor sits at `h ≈ 3 × 10⁻⁸ … 10⁻⁷ Re`; in fp32 it sits four
 *    orders of magnitude higher.
 *
 * The two curves cross at `h ≈ r·√ε`, which is why the automatic step here is **proportional to
 * the radius** rather than absolute: both terms scale with `r`, so one ratio serves LEO and the
 * outer belt alike. Measured optima, same 60 points, error relative to `|∇|B||`:
 *
 * | step | fp64 host lane, max / median | fp32 device lane, max / median |
 * |---|---|---|
 * | `h = r · 5 × 10⁻⁸` (@ref host_step_ratio) | **1.1e-07 / 7.3e-08** | 7.5e-01 / 3.6e-01 |
 * | `h = r · 10⁻⁴` (@ref device_step_ratio) | 2.2e-04 / 1.5e-04 | **4.1e-04 / 2.2e-04** |
 * | `h = 10⁻³ Re` (absolute, IRBEM-style) | 1.2e-03 / 2.2e-04 | 1.2e-03 / 3.0e-04 |
 *
 * Read the middle row across: at the *device's* step the fp64 lane is already at `2.2e-04`, and
 * fp32 only doubles it. **The step, not the precision, is what limits the device lane** — an fp32
 * kernel is not the problem, being unable to take a `10⁻⁷` step in fp32 is. And read the first
 * column down: a step chosen for fp64 is catastrophic in fp32, a 75% error. There is therefore no
 * single default that serves both lanes, and this file does not pretend otherwise — @ref auto_step
 * takes the lane as an argument and the two constants are named after their lanes.
 *
 * A caller who supplies `dX` explicitly gets exactly that step on both lanes, which is what a
 * differential comparison against IRBEM needs: at `dX = 10⁻³` the two implementations are being
 * compared at matched resolution, the same discipline `docs/ERROR_BUDGET.md` §2(a) imposes on L\*.
 *
 * ## What the reference actually computes, established as a black box
 *
 * IRBEM is LGPL-3.0 and this repository is MIT, so its Fortran is never read. Its *behaviour* is
 * fair game, and the following were measured by running the compiled `-O2` oracle against inputs
 * chosen to separate the candidate definitions. Each is reproducible from
 * `tools/oracle/` and each is asserted by a test here:
 *
 *  - `GET_BDERIVS` is a **forward** (one-sided) difference, never central. Feeding the oracle's own
 *    `GET_FIELD_MULTI` outputs at `x` and `x + dX·ê_j` through `(B(x+dX·ê_j) − B(x))/dX` reproduces
 *    its `gradBmag` and `diffB` **bit for bit** — relative difference exactly `0.0` at
 *    `dX = 10⁻¹, 10⁻², 10⁻³` over 300 points. A central difference disagrees at the 10⁻³ level, so
 *    this is not a coincidence of tolerance.
 *  - `diffB(i,j) = ∂Bᵢ/∂xⱼ`, Fortran column-major, so component `i` is the fast index. @ref
 *    BDerivatives::diff_b uses `(row, col) = (i, j)`, the same convention in mathematical notation.
 *  - `COMPUTE_GRAD_CURV_CURL`'s `curvature` is `Â − (Â·B̂)B̂` with `Â = (B̂·∇)B/|B|` — it projects
 *    using `Â`'s **own** parallel component and never touches `gradBmag`. For a consistent Jacobian
 *    that equals the documented `(B̂·∇)B̂` identically; for the *inconsistent* inputs a black-box
 *    probe can feed it, the two differ by a multiple of `B̂`, and the oracle follows the first.
 *    Reproducing it to `3.7 × 10⁻¹⁵` while the alternative is off by 10³ settles which.
 *  - `GET_HEMI_MULTI` returns the sign of `d|B|/ds` along `+B̂` — the northern hemisphere is the
 *    side of the magnetic equator where the field is *rising* in the direction B points. Agreed
 *    with the oracle on 400 / 400 random points spanning `r = 1.2 … 8 Re`.
 *  - `options(2) = 0` means "initialize IGRF once per year (year.5)" — the published options table
 *    says so, and it is why every oracle comparison here uses epoch **2015.5** exactly rather than
 *    the day-of-year fraction. With that epoch and degree **10** — IRBEM's internal truncation,
 *    identified by sweeping degrees 8…13 and finding the residual collapse — @ref field_batch
 *    reproduces `GET_FIELD_MULTI` to `7.3 × 10⁻¹⁴` relative, four hundred million times inside the
 *    `1 × 10⁻⁶` `Bgeo` budget. At degree 13 (IGRF-14 as IAGA published it) the same comparison is
 *    `1.6 × 10⁻⁴`, which is the *model* difference, not an error in either implementation.
 *
 * ## Sources
 *
 *  - Alken et al., *International Geomagnetic Reference Field: the fourteenth generation*, Earth
 *    Planets Space (2025) — the field itself; evaluated by @ref Igrf.
 *  - Northrop, *The Adiabatic Motion of Charged Particles*, Interscience (1963), ch. 1 — the
 *    guiding-centre drifts these derivative products are the geometric factors of: the gradient
 *    drift `(B̂ × ∇⊥B)/B`, the curvature drift `B̂ × (B̂·∇)B̂`, and the radius of curvature.
 *  - Roederer, *Dynamics of Geomagnetically Trapped Radiation*, Springer (1970), §1.3.
 */

#include <cmath>
#include <cstddef>
#include <cstdint>

#include <array>
#include <concepts>
#include <limits>
#include <span>
#include <vector>

#include "frames.hpp"
#include "igrf.hpp"
#include "policy.hpp"
#include "status.hpp"

// The device lane is opt-in by include path, exactly as @ref lstar.hpp's is. Without
// cheatah-gpu-linalg every routine below still compiles and runs its host loop.
#if __has_include("cheatah_gpu_linalg/context.hpp")
#  include "gpu/dispatch.hpp"
#  define CHEATAH_SPACE_IRBEM_FIELD_GPU 1
#endif

namespace cheatah::space::irbem {

/**
 * A magnetic field model that can be asked for `B` at a geographic Cartesian point.
 *
 * @ref Igrf satisfies it, and so will an internal-plus-external superposition when the external
 * models land — which is the point of naming the requirement rather than hard-coding @ref Igrf into
 * the single-point routines. The batch routines are *not* written against this concept, because the
 * device lane has to upload a specific model's coefficients and cannot dispatch on a callable.
 *
 * @tparam M the candidate model type.
 */
template <class M>
concept GeoFieldModel = requires(const M& model, const Position<Frame::GEO>& p) {
    { model.evaluate(p) } -> std::same_as<FieldVector<Frame::GEO>>;
};

// ---------------------------------------------------------------------------------------------
// The finite-difference step
// ---------------------------------------------------------------------------------------------

/**
 * Which arithmetic a finite difference will be taken in — the argument @ref auto_step needs and
 * the one thing a caller cannot infer from the step alone.
 *
 * Not a lane *selector*: @ref bderivs_batch still decides where to run from the batch size and the
 * measured crossover. This says which precision the differenced values will have arrived in, which
 * is what moves the optimal step by four orders of magnitude (see the file brief).
 */
enum class DifferenceLane : std::uint8_t {
    Fp64Host,   ///< The `double` host lane — @ref host_step_ratio.
    Fp32Device, ///< The `float` device lane — @ref device_step_ratio.
};

/**
 * The measured optimal step ratio for a difference of `double` field values, as a fraction of the
 * geocentric radius.
 *
 * `h = r · 5 × 10⁻⁸` sits in the trough between truncation (`~2h/r`) and fp64 cancellation
 * (`~ε₆₄·r/h`). Measured max / median relative error on `∇|B|`, 60 points over `r = 1.05 … 11.5`,
 * against a Richardson-extrapolated reference: **1.1 × 10⁻⁷ / 7.3 × 10⁻⁸**.
 */
inline constexpr double host_step_ratio = 5.0e-8;

/**
 * The measured optimal step ratio for a difference of `float` field values.
 *
 * Three thousand times larger than @ref host_step_ratio, and that factor is `√(ε₃₂/ε₆₄)` and
 * nothing else. Measured max / median on the same 60 points: **4.1 × 10⁻⁴ / 2.2 × 10⁻⁴** — of which
 * `2.2 × 10⁻⁴ / 1.5 × 10⁻⁴` is truncation the fp64 lane would also pay at this step. The device is
 * limited by the step it can afford, not by the arithmetic it does.
 */
inline constexpr double device_step_ratio = 1.0e-4;

/**
 * The finite-difference step for a point at radius @p radius_re, in Earth radii.
 *
 * Proportional to the radius on purpose: truncation error scales as `h/r` and cancellation as
 * `r/h`, so their crossing point moves with `r` and a single absolute step is optimal at exactly
 * one altitude. IRBEM takes an absolute `dX` for the whole batch, which is why this is a *default*
 * rather than the only option — @ref bderivs and @ref bderivs_batch accept an explicit step, and a
 * differential comparison must use one so both implementations are differenced identically.
 *
 * @param radius_re the geocentric distance of the point, Earth radii.
 * @param lane which precision the differenced field values will be in.
 * @return the step in Earth radii; floored at the ratio itself, so a point at or inside `r = 1`
 *         still gets a positive step rather than a zero one that would divide by zero.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemField.AutoStepTracksTheRadiusAndTheLane
 */
[[nodiscard]] inline double auto_step(double radius_re, DifferenceLane lane) {
    const double ratio = (lane == DifferenceLane::Fp32Device) ? device_step_ratio : host_step_ratio;
    const double r = (radius_re > 1.0 && std::isfinite(radius_re)) ? radius_re : 1.0;
    return ratio * r;
}

// ---------------------------------------------------------------------------------------------
// GET_BDERIVS
// ---------------------------------------------------------------------------------------------

/**
 * The field and its first derivatives at one point — IRBEM's `GET_BDERIVS` outputs, typed.
 *
 * Fixed size, trivially copyable, no allocation: this is what a million-point ephemeris is stored
 * as, so it is 128 bytes of plain data and not a handle to anything.
 */
struct BDerivatives {
    /// `B` at the point, GEO, nT.
    FieldVector<Frame::GEO> b{};
    /// `|B|`, nT — IRBEM's `Bmag`, and the `Blocal` of the tracing routines.
    double b_mag = 0.0;
    /// `∇|B|` in GEO, nT/Re — IRBEM's `gradBmag`.
    fixarray::vec3d grad_b_mag{};
    /// `∂Bᵢ/∂xⱼ` at `(row, col) = (i, j)`, nT/Re — IRBEM's `diffB(i,j)`, whose Fortran
    /// column-major storage makes `i` the fast index and therefore matches this indexing.
    fixarray::mat3d diff_b{};
};

/**
 * The field and its first derivatives at one point, by forward differences — the fp64 reference
 * lane.
 *
 * Four field evaluations: the base point, then `x + h·ê_x`, `x + h·ê_y`, `x + h·ê_z`. **Forward,
 * not central**, because that is what the reference does — established black-box, and exactly (the
 * measured relative difference is `0.0`, not "small"). A central difference would cost seven
 * evaluations and buy two orders of magnitude; it is offered by nobody here because agreeing with
 * IRBEM is the contract, and the accuracy that matters is quantified in the file brief instead of
 * silently improved.
 *
 * The gradient of the magnitude and the Jacobian of the vector are differenced from the SAME four
 * evaluations. That is not merely an economy: it is what makes `grad_b_mag` and `diff_b` mutually
 * consistent to the same order, which is what @ref grad_curv_curl's `grad_par` and `curvature`
 * silently assume when they compute the same parallel derivative two different ways.
 *
 * @tparam M the field model type.
 * @param model the field model.
 * @param p the point, GEO, Earth radii.
 * @param step_re the difference step `dX` in Earth radii; `0.0` (the default) selects
 *        @ref auto_step for @ref DifferenceLane::Fp64Host. Pass an explicit step to compare
 *        against another implementation at matched resolution.
 * @return the derivatives, or @ref Status::DomainError for a non-finite point, a point at the
 *         origin, or a non-finite step. The value is zero-filled in the failure case rather than
 *         left indeterminate.
 * @complexity Four field evaluations — ~7 600 flops at IGRF degree 13 — plus ~40 for the
 *             differencing.
 * @alloc none.
 * @test IrbemField.BderivsMatchesTheAnalyticDipoleJacobian
 * @test IrbemField.BderivsIsAForwardDifferenceNotACentralOne
 */
template <GeoFieldModel M>
[[nodiscard]] inline Result<BDerivatives> bderivs(const M& model, const Position<Frame::GEO>& p,
                                                  double step_re = 0.0) {
    BDerivatives out{};
    const double r = fixarray::norm(p.v);
    if (!(r > 0.0) || !std::isfinite(r) || !std::isfinite(step_re) || step_re < 0.0) {
        return {Status::DomainError, out};
    }
    const double h = (step_re > 0.0) ? step_re : auto_step(r, DifferenceLane::Fp64Host);

    out.b = model.evaluate(p);
    out.b_mag = out.b.magnitude();

    for (std::size_t j = 0; j < 3; ++j) {
        Position<Frame::GEO> q = p;
        q.v[j] += h;
        const FieldVector<Frame::GEO> bj = model.evaluate(q);
        out.grad_b_mag[j] = (bj.magnitude() - out.b_mag) / h;
        for (std::size_t i = 0; i < 3; ++i) { out.diff_b(i, j) = (bj.v[i] - out.b.v[i]) / h; }
    }
    return {Status::Ok, out};
}

// ---------------------------------------------------------------------------------------------
// COMPUTE_GRAD_CURV_CURL
// ---------------------------------------------------------------------------------------------

/**
 * The guiding-centre geometry that follows from @ref BDerivatives — IRBEM's
 * `COMPUTE_GRAD_CURV_CURL` outputs, typed.
 *
 * Every member is a *geometric factor*, not a velocity: the gradient drift of a particle is
 * `(m v⊥²/2q)·grad_drift` and the curvature drift is `(m v∥²/q)·curv_drift` (Northrop 1963, §1.3).
 * Separating the geometry from the particle is what lets one field evaluation serve every energy
 * and pitch angle at that point, which is exactly the shape a radiation-belt calculation has.
 */
struct GradCurvCurl {
    /// `∇|B| · B̂`, nT/Re — the field's rate of change ALONG itself. Its sign is the magnetic
    /// hemisphere (see @ref hemisphere) and its zero is the magnetic equator.
    double grad_par = 0.0;
    /// `∇|B| − grad_par·B̂`, nT/Re — the part of the gradient that drives drift.
    fixarray::vec3d grad_perp{};
    /// `(B̂ × grad_perp)/|B|`, 1/Re — the geometric factor of the gradient drift.
    fixarray::vec3d grad_drift{};
    /// `(B̂·∇)B̂`, 1/Re — the field-line curvature vector, pointing toward the centre of curvature.
    fixarray::vec3d curvature{};
    /// `1/|curvature|`, Re — the radius of curvature; infinite for a straight field line.
    double r_curv = 0.0;
    /// `B̂ × curvature`, 1/Re — the geometric factor of the curvature drift.
    fixarray::vec3d curv_drift{};
    /// `∇ × B`, nT/Re — proportional to the local current density through Ampère's law, and zero
    /// wherever the field is a pure potential field.
    fixarray::vec3d curl_b{};
    /// `∇ · B`, nT/Re — zero by Maxwell, so a nonzero value measures the differencing error and
    /// nothing else. That makes it the free residual diagnostic this struct carries.
    double div_b = 0.0;
};

/**
 * The guiding-centre geometry at one point, from its field derivatives — pure algebra.
 *
 * No field model, no evaluations, ~50 flops. This is the routine the file brief measures at
 * 11.6 ns/point on the host and rules off the device permanently.
 *
 * The one subtlety is `curvature`. The identity `(B̂·∇)B̂ = [Â − (Â·B̂)B̂]` with `Â = (B̂·∇)B/|B|`
 * holds because `B̂·B̂ = 1` forces the derivative of `B̂` to be perpendicular to it. Writing it that
 * way — projecting out `Â`'s own parallel part rather than subtracting `grad_par·B̂/|B|` — is what
 * the reference does, and the difference is not cosmetic: the two agree only when `grad_b_mag` and
 * `diff_b` are mutually consistent, which finite differences make them only to first order. The
 * form used here needs `grad_b_mag` not at all, so `curvature` and `r_curv` are exactly
 * perpendicular to `B̂` by construction rather than approximately.
 *
 * @param d the field and its derivatives at the point, as @ref bderivs produces them.
 * @return the geometry, or @ref Status::DomainError when `|B|` is zero or non-finite — the point
 *         where `B̂` does not exist and every output below is meaningless rather than merely large.
 * @complexity O(1) — ~50 flops, one square root for `|curvature|`, no transcendentals.
 * @alloc none.
 * @test IrbemField.GradCurvCurlReproducesTheOracleAlgebra
 * @test IrbemField.CurvatureIsPerpendicularToTheField
 */
[[nodiscard]] inline Result<GradCurvCurl> grad_curv_curl(const BDerivatives& d) {
    GradCurvCurl g{};
    if (!(d.b_mag > 0.0) || !std::isfinite(d.b_mag)) return {Status::DomainError, g};

    const fixarray::vec3d bhat = d.b.v / d.b_mag;

    g.grad_par = fixarray::dot(bhat, d.grad_b_mag);
    g.grad_perp = d.grad_b_mag - (bhat * g.grad_par);
    g.grad_drift = fixarray::cross(bhat, g.grad_perp) / d.b_mag;

    // Â = (B̂·∇)B / |B|, then the perpendicular projection that IS (B̂·∇)B̂.
    fixarray::vec3d a_hat{};
    for (std::size_t i = 0; i < 3; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 3; ++j) { s += bhat[j] * d.diff_b(i, j); }
        a_hat[i] = s / d.b_mag;
    }
    g.curvature = a_hat - (bhat * fixarray::dot(a_hat, bhat));

    const double kappa = fixarray::norm(g.curvature);
    // A straight field line has infinite radius of curvature. That is the physics, not an overflow:
    // returning a huge finite number instead would make a uniform-field test look like a tight one.
    g.r_curv = kappa > 0.0 ? 1.0 / kappa : std::numeric_limits<double>::infinity();
    g.curv_drift = fixarray::cross(bhat, g.curvature);

    g.curl_b = fixarray::vec3d{d.diff_b(2, 1) - d.diff_b(1, 2), d.diff_b(0, 2) - d.diff_b(2, 0),
                               d.diff_b(1, 0) - d.diff_b(0, 1)};
    g.div_b = d.diff_b(0, 0) + d.diff_b(1, 1) + d.diff_b(2, 2);
    return {Status::Ok, g};
}

// ---------------------------------------------------------------------------------------------
// GET_HEMI_MULTI
// ---------------------------------------------------------------------------------------------

/**
 * Which magnetic hemisphere a point lies in — IRBEM's `xHEMI`, with its integer values.
 *
 * "Magnetic" and not geographic: the boundary is the magnetic equator of the field line through the
 * point, which at the surface is displaced by up to ~15° from the geographic equator and, in the
 * South Atlantic Anomaly, is somewhere a geographic test would not put it at all.
 */
enum class Hemisphere : std::int8_t {
    South = -1,   ///< South of the magnetic equator along this field line.
    Invalid = 0,  ///< No magnetic field to speak of — IRBEM's `0`.
    North = 1,    ///< North of the magnetic equator along this field line.
};

/**
 * Which magnetic hemisphere @p p lies in.
 *
 * The criterion is the sign of `d|B|/ds` along `+B̂`: the field falls to a minimum at the magnetic
 * equator and rises toward both feet, and `B` points from the southern foot over the equator to the
 * northern one, so a *rising* field in the direction B points means the equator is behind you and
 * you are north of it. That is one signed scalar and it is exactly @ref GradCurvCurl::grad_par,
 * which is why this routine and the derivative routines are in the same file.
 *
 * Computed with a **central** difference along `B̂` rather than the forward differences of
 * @ref bderivs — three evaluations instead of four, and symmetric, which matters here specifically:
 * near the equator the true derivative passes through zero, and a one-sided step can overshoot the
 * minimum and report the wrong side of a boundary the answer is a discrete function of. A
 * derivative that is merely inaccurate is fine; a *sign* that is wrong is a different hemisphere.
 *
 * @tparam M the field model type.
 * @param model the field model.
 * @param p the point, GEO, Earth radii.
 * @param step_re the step along `B̂` in Earth radii; `0.0` selects @ref auto_step for the host lane.
 *        Unlike @ref bderivs this is not accuracy-critical — only the sign survives — so the
 *        default is generous rather than optimal.
 * @return the hemisphere, with @ref Status::DomainError for a non-finite or origin point and
 *         @ref Hemisphere::Invalid where the field vanishes or the derivative is exactly zero.
 * @complexity Three field evaluations.
 * @alloc none.
 * @test IrbemField.HemisphereAgreesWithTheOracleGoldens
 * @test IrbemField.HemisphereFlipsAcrossTheDipoleEquator
 */
template <GeoFieldModel M>
[[nodiscard]] inline Result<Hemisphere> hemisphere(const M& model, const Position<Frame::GEO>& p,
                                                   double step_re = 0.0) {
    const double r = fixarray::norm(p.v);
    if (!(r > 0.0) || !std::isfinite(r) || !std::isfinite(step_re) || step_re < 0.0) {
        return {Status::DomainError, Hemisphere::Invalid};
    }
    const FieldVector<Frame::GEO> b = model.evaluate(p);
    const double bmag = b.magnitude();
    if (!(bmag > 0.0)) return {Status::Ok, Hemisphere::Invalid};

    // The step is measured ALONG the field, so it is scaled by the radius exactly as a coordinate
    // step is: the same 1e-4-ish fraction of the local scale height.
    const double h = (step_re > 0.0) ? step_re : auto_step(r, DifferenceLane::Fp32Device);
    const fixarray::vec3d bhat = b.v / bmag;
    Position<Frame::GEO> fwd{p.v + (bhat * h)};
    Position<Frame::GEO> bwd{p.v - (bhat * h)};
    const double slope = model.evaluate(fwd).magnitude() - model.evaluate(bwd).magnitude();
    if (slope > 0.0) return {Status::Ok, Hemisphere::North};
    if (slope < 0.0) return {Status::Ok, Hemisphere::South};
    return {Status::Ok, Hemisphere::Invalid};
}

// ---------------------------------------------------------------------------------------------
// The device lane
// ---------------------------------------------------------------------------------------------

#ifdef CHEATAH_SPACE_IRBEM_FIELD_GPU
namespace detail {

/**
 * Pack a model's coefficients and the Legendre normalisation into the two `float` buffers
 * `irbem_igrf_f32` binds.
 *
 * Done ONCE per batch on the host, never per thread: the 26-epoch IGRF table is 5 460 doubles and
 * interpolating it on the device would be one redundant copy of that work per point. The
 * normalisation is `constexpr` here and therefore already computed at compile time.
 *
 * @tparam NMAX the truncation degree.
 * @param model the model, already built for the epoch.
 * @param coef receives `g[slots]` then `h[slots]`.
 * @param norm receives `e[slots]`, `f[slots]`, `diagonal[NMAX+1]`.
 * @complexity O(NMAX²) — 105 slots at degree 13.
 * @alloc the two vectors are sized by the caller; this fills them.
 * @test IrbemField.FieldBatchAgreesWithTheReferenceLane
 */
template <int NMAX>
inline void pack_igrf_tables(const Igrf<NMAX>& model, std::vector<float>& coef,
                             std::vector<float>& norm) {
    constexpr int kSlots = ((NMAX + 1) * (NMAX + 2)) / 2;
    coef.assign(2 * static_cast<std::size_t>(kSlots), 0.0F);
    norm.assign(2 * static_cast<std::size_t>(kSlots) + NMAX + 1, 0.0F);
    for (int deg = 1; deg <= NMAX; ++deg) {
        for (int m = 0; m <= deg; ++m) {
            const std::size_t k = (static_cast<std::size_t>(deg) * (deg + 1)) / 2 + m;
            coef[k] = static_cast<float>(model.g(deg, m));
            coef[static_cast<std::size_t>(kSlots) + k] = static_cast<float>(model.h(deg, m));
        }
    }
    constexpr auto kNorm = ::cheatah::space::irbem::detail::make_legendre_normalisation<NMAX, double>();
    for (int k = 0; k < kSlots; ++k) {
        norm[static_cast<std::size_t>(k)] = static_cast<float>(kNorm.e[static_cast<std::size_t>(k)]);
        norm[static_cast<std::size_t>(kSlots + k)] =
            static_cast<float>(kNorm.f[static_cast<std::size_t>(k)]);
    }
    for (int deg = 0; deg <= NMAX; ++deg) {
        norm[static_cast<std::size_t>((2 * kSlots) + deg)] =
            static_cast<float>(kNorm.diagonal[static_cast<std::size_t>(deg)]);
    }
}

/**
 * Evaluate IGRF at `pos.size()/3` points on the device, through the existing `irbem_igrf_f32`
 * entry point.
 *
 * The launcher for the five-binding IGRF shape, which `gpu::dispatch_batch` cannot express (it
 * binds pos/out/params/dims and nothing else) and which no shared launcher provides. It lives here
 * rather than in `gpu/dispatch.hpp` so that adding a *consumer* of an existing kernel does not
 * touch the file every kernel author shares.
 *
 * @tparam NMAX the truncation degree.
 * @param model the model, already built for the epoch.
 * @param pos the points, xyz-interleaved, `3N` floats, GEO, Earth radii.
 * @param out receives the field, xyz-interleaved, `3N` floats, nT.
 * @return `false` when there is no device or the kernel was never compiled — the caller's cue to
 *         run the host loop, not an error.
 * @complexity One dispatch over `ceil(N/256)` workgroups; O(N·NMAX²) concurrent flops.
 * @alloc five device buffers, returned to the context's pool on scope exit, plus the two staged
 *        coefficient vectors.
 * @test IrbemField.FieldBatchUsesTheDeviceWhenOneIsAvailable
 */
template <int NMAX>
[[nodiscard]] inline bool igrf_on_device(const Igrf<NMAX>& model, std::span<const float> pos,
                                         std::span<float> out) {
    const std::size_t n = pos.size() / 3;
    if (n == 0) return true;
    if (!gpu::available() || !std::filesystem::exists(gpu::shader_path("irbem_igrf_f32"))) {
        return false;
    }

    std::vector<float> coef;
    std::vector<float> norm;
    pack_igrf_tables(model, coef, norm);
    const std::array<std::uint32_t, 2> dims{static_cast<std::uint32_t>(n),
                                            static_cast<std::uint32_t>(NMAX)};

    namespace gl = gpu::detail::gl;
    gl::detail::Context& c = gl::detail::ctx();
    gpu::detail::Leases lease;
    // Positions and results are large and touched only by the kernel, so they are device-local;
    // the two tables and dims are small and host-written once, so they take the mapped path.
    gl::detail::Buffer* b_pos = lease.add(c.new_data_buffer(pos.size() * sizeof(float)));
    gl::detail::Buffer* b_out = lease.add(c.new_data_buffer(out.size() * sizeof(float)));
    gl::detail::Buffer* b_cf = lease.add(c.new_buffer(coef.size() * sizeof(float)));
    gl::detail::Buffer* b_nr = lease.add(c.new_buffer(norm.size() * sizeof(float)));
    gl::detail::Buffer* b_dm = lease.add(c.new_buffer(dims.size() * sizeof(std::uint32_t)));
    c.upload(b_pos, pos.data(), pos.size() * sizeof(float));
    c.upload(b_cf, coef.data(), coef.size() * sizeof(float));
    c.upload(b_nr, norm.data(), norm.size() * sizeof(float));
    c.upload(b_dm, dims.data(), dims.size() * sizeof(std::uint32_t));
    {
        const gpu::detail::SpvDirScope scope(gpu::shader_dir().string());
        c.dispatch_1d("irbem_igrf_f32", lease.data(), 5, n);
    }
    c.download(b_out, out.data(), out.size() * sizeof(float));
    return true;
}

}  // namespace detail
#endif  // CHEATAH_SPACE_IRBEM_FIELD_GPU

namespace detail {

/**
 * Fill @p pos with the `3N` floats `irbem_igrf_f32` wants, from typed positions.
 * @param points the points, GEO, Earth radii.
 * @param pos receives `3N` floats, xyz-interleaved; sized by the caller.
 * @complexity O(N).
 * @alloc none — the span is the caller's.
 * @test IrbemField.FieldBatchAgreesWithTheReferenceLane
 */
inline void interleave(std::span<const Position<Frame::GEO>> points, std::span<float> pos) {
    for (std::size_t i = 0; i < points.size(); ++i) {
        pos[(3 * i) + 0] = static_cast<float>(points[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(points[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(points[i].v[2]);
    }
}

}  // namespace detail

// ---------------------------------------------------------------------------------------------
// GET_FIELD_MULTI
// ---------------------------------------------------------------------------------------------

/**
 * The field at every point of a batch — IRBEM's `GET_FIELD_MULTI`.
 *
 * The reference is a bare loop over points, and it is the most trivially parallel routine in the
 * whole library: one point in, one vector out, no state carried between iterations. It is also the
 * *least* arithmetically intense thing worth offloading — ~1 900 flops for 24 bytes in and 24 out,
 * about 20 flops/byte — which is exactly why the crossover is consulted rather than assumed. The
 * measurement that sets it, on an RTX 3070 Ti against `igrf.hpp`'s fp64 host lane built
 * `-O3 -march=native -ffp-contract=off`, 2²⁰ points:
 *
 *     device 36.6 ns/eval (27.3 Mpts/s) | host 306.5 ns/eval (3.3 Mpts/s) | 8.96x
 *
 * and the crossover that follows from those two throughputs and the measured ~30 µs synchronous
 * dispatch floor is ~113 points, rounded to 128 in `gpu/dispatch.hpp`'s registry. Below it the
 * submit cost alone exceeds the whole computation and the host wins; a per-point loop calling
 * @ref Igrf::evaluate can never reach the device at all, which is the reason this entry point
 * takes the whole batch.
 *
 * The device lane returns fp32. Measured maximum relative deviation against the fp64 host lane over
 * 2²⁰ points: **8.8 × 10⁻⁷**, inside the `1 × 10⁻⁶` `Bgeo` budget of `docs/ERROR_BUDGET.md` §4 —
 * and it is inside it because nothing here accumulates. One point, one thread, no reduction; the
 * budget's accumulation concern bites the integrals in @ref lstar.hpp, not this.
 *
 * @tparam NMAX the IGRF truncation degree. Degree 10 is IRBEM's internal truncation and the one a
 *         differential comparison must use; degree 13 is IGRF-14 as IAGA published it.
 * @param model the internal field model, already built for the epoch.
 * @param points the points, GEO, Earth radii.
 * @param b receives one field vector per point, GEO, nT; same length as @p points.
 * @param b_mag receives `|B|` per point, nT; same length as @p points. Computed from the returned
 *        vector, so the two are consistent to the last bit on both lanes.
 * @return @ref Status::DomainError on a length mismatch, @ref Status::Ok otherwise. The value is
 *         `true` when the device serviced the call — asserted by a test rather than assumed,
 *         because a silent fallback to the host is what makes a speed claim worthless.
 * @complexity O(N·NMAX²); on the device those run concurrently.
 * @alloc none on the host lane. The device lane stages `3N` floats in and `3N` out.
 * @test IrbemField.FieldBatchAgreesWithTheReferenceLane
 * @test IrbemField.FieldBatchMatchesTheOracleGoldens
 * @test IrbemField.FieldBatchUsesTheDeviceWhenOneIsAvailable
 */
template <int NMAX>
[[nodiscard]] inline Result<bool> field_batch(const Igrf<NMAX>& model,
                                              std::span<const Position<Frame::GEO>> points,
                                              std::span<FieldVector<Frame::GEO>> b,
                                              std::span<double> b_mag) {
    const std::size_t n = points.size();
    if (b.size() != n || b_mag.size() != n) return {Status::DomainError, false};
    if (n == 0) return {Status::Ok, false};

#ifdef CHEATAH_SPACE_IRBEM_FIELD_GPU
    if (gpu::prefer_gpu("irbem_igrf_f32", n)) {
        std::vector<float> pos(3 * n);
        std::vector<float> out(3 * n);
        detail::interleave(points, pos);
        if (detail::igrf_on_device(model, pos, out)) {
            for (std::size_t i = 0; i < n; ++i) {
                b[i] = FieldVector<Frame::GEO>{fixarray::vec3d{
                    out[(3 * i) + 0], out[(3 * i) + 1], out[(3 * i) + 2]}};
                b_mag[i] = b[i].magnitude();
            }
            return {Status::Ok, true};
        }
        // No device or no compiled SPIR-V after all: a deployment problem, not a reason to refuse
        // to compute. Fall through to the host loop.
    }
#endif

    for (std::size_t i = 0; i < n; ++i) {
        b[i] = model.evaluate(points[i]);
        b_mag[i] = b[i].magnitude();
    }
    return {Status::Ok, false};
}

/**
 * The field and its first derivatives at every point of a batch — IRBEM's `GET_BDERIVS`.
 *
 * **One dispatch, not four.** Each point needs four field evaluations — the base and three
 * one-sided neighbours — and the obvious implementation issues four batched dispatches of `N`
 * points each. That pays the ~30 µs submit floor four times for work that has no dependency between
 * the four groups whatsoever. Building the `4N` points up front and dispatching once pays it once,
 * and at the same time quadruples the occupancy of the launch, which is what pulls the crossover
 * down: this routine reaches the device at **32 points**, because 32 points are 128 evaluations and
 * 128 evaluations is `irbem_igrf_f32`'s measured crossover.
 *
 * The layout is `[all N base points][all N +x][all N +y][all N +z]` rather than four consecutive
 * points per input point. Same dispatch either way, but this way each of the four groups is a
 * contiguous run whose lanes read neighbouring coefficients in the same order, and the host-side
 * differencing walks four unit-stride streams instead of one stride-4 one.
 *
 * On the device the differenced values are fp32, and the file brief's table is the consequence:
 * the achievable accuracy is ~4 × 10⁻⁴ relative against ~1 × 10⁻⁷ on the host, and the step that
 * achieves it is three thousand times larger. That is why @p step_re defaults to *the lane's* step
 * rather than to a constant, and why a caller comparing the two lanes must pass an explicit step to
 * both. It is a real limitation, quantified here rather than hidden behind an average.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param points the points, GEO, Earth radii.
 * @param out receives one @ref BDerivatives per point; same length as @p points.
 * @param step_re the difference step `dX` in Earth radii; `0.0` (the default) selects
 *        @ref auto_step for whichever lane actually runs.
 * @return @ref Status::DomainError on a length mismatch or a non-finite step; @ref Status::Ok
 *         otherwise. The value is `true` when the device serviced the call.
 * @complexity O(4·N·NMAX²) field evaluations plus ~40 flops per point of differencing.
 * @alloc none on the host lane. The device lane stages `12N` floats in and `12N` out.
 * @test IrbemField.BderivsBatchAgreesWithTheReferenceLane
 * @test IrbemField.BderivsBatchMatchesTheOracleGoldens
 * @test IrbemField.BderivsBatchUsesTheDeviceWhenOneIsAvailable
 */
template <int NMAX>
[[nodiscard]] inline Result<bool> bderivs_batch(const Igrf<NMAX>& model,
                                                std::span<const Position<Frame::GEO>> points,
                                                std::span<BDerivatives> out,
                                                double step_re = 0.0) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(step_re) || step_re < 0.0) return {Status::DomainError, false};
    if (n == 0) return {Status::Ok, false};

#ifdef CHEATAH_SPACE_IRBEM_FIELD_GPU
    // The dispatch is 4N points wide, so the crossover is asked about 4N — the size the kernel
    // actually sees. In points of THIS routine that is a crossover of 32.
    if (gpu::prefer_gpu("irbem_igrf_f32", 4 * n)) {
        std::vector<float> pos(12 * n);
        std::vector<float> field(12 * n);
        std::vector<double> steps(n);
        detail::interleave(points, std::span<float>(pos).first(3 * n));
        for (std::size_t i = 0; i < n; ++i) {
            const double r = fixarray::norm(points[i].v);
            steps[i] = (step_re > 0.0) ? step_re : auto_step(r, DifferenceLane::Fp32Device);
            for (std::size_t j = 0; j < 3; ++j) {
                const std::size_t base = 3 * (((j + 1) * n) + i);
                pos[base + 0] = static_cast<float>(points[i].v[0]);
                pos[base + 1] = static_cast<float>(points[i].v[1]);
                pos[base + 2] = static_cast<float>(points[i].v[2]);
                pos[base + j] = static_cast<float>(points[i].v[j] + steps[i]);
            }
        }
        if (detail::igrf_on_device(model, pos, field)) {
            for (std::size_t i = 0; i < n; ++i) {
                BDerivatives& d = out[i];
                d = BDerivatives{};
                d.b = FieldVector<Frame::GEO>{fixarray::vec3d{field[(3 * i) + 0],
                                                              field[(3 * i) + 1],
                                                              field[(3 * i) + 2]}};
                d.b_mag = d.b.magnitude();
                // The step the DEVICE actually took: the fp32 position it evaluated at is not
                // exactly x + h, and dividing by the intended h instead of the realised one would
                // add a relative error of order eps32*r/h -- which at h = 1e-4*r is 6e-4, larger
                // than everything else in this routine's budget put together.
                for (std::size_t j = 0; j < 3; ++j) {
                    const std::size_t base = 3 * (((j + 1) * n) + i);
                    const double h = static_cast<double>(pos[base + j]) -
                                     static_cast<double>(pos[(3 * i) + j]);
                    const fixarray::vec3d bj{field[base + 0], field[base + 1], field[base + 2]};
                    const double hj = h != 0.0 ? h : steps[i];
                    d.grad_b_mag[j] = (fixarray::norm(bj) - d.b_mag) / hj;
                    for (std::size_t k = 0; k < 3; ++k) {
                        d.diff_b(k, j) = (bj[k] - d.b.v[k]) / hj;
                    }
                }
            }
            return {Status::Ok, true};
        }
    }
#endif

    for (std::size_t i = 0; i < n; ++i) {
        out[i] = bderivs(model, points[i], step_re).value;
    }
    return {Status::Ok, false};
}

/**
 * The guiding-centre geometry at every point of a batch — IRBEM's `COMPUTE_GRAD_CURV_CURL`.
 *
 * **A host loop, permanently, and the number that settles it:** 136 bytes cross the bus per point
 * (16 fp32 in, 18 out) for 11.6 ns of host arithmetic — measured, 2²⁰ points, `-O3 -march=native`,
 * 86 Mpts/s. At ~12 GB/s of effective PCIe bandwidth the transfer alone is ~11.3 ns/point, so
 * moving the data costs what computing it costs, *before* the ~30 µs submit floor and before a
 * single lane runs. Arithmetic intensity is ~0.4 flops/byte, which is below even the dipole
 * kernel's 0.5 — and `gpu/dispatch.hpp` already records that the dipole kernel LOSES 0.69× at
 * every batch size it was measured at. There is no size at which this one wins, so no kernel for
 * it exists in `gpu/irbem.slang` and none should be written.
 *
 * The routine is still worth having as a batch entry point: it keeps the loop in one place, it
 * vectorizes, and it lets a caller hand over an ephemeris rather than write the loop themselves.
 *
 * @param derivs the field and derivatives per point, as @ref bderivs_batch produces them.
 * @param out receives one @ref GradCurvCurl per point; same length as @p derivs.
 * @param statuses receives each point's status, so a single degenerate point reports itself instead
 *        of spoiling the batch; same length as @p derivs.
 * @return @ref Status::DomainError on a length mismatch, @ref Status::Ok when every point was
 *         computable, and the first non-`Ok` per-point status otherwise.
 * @complexity O(N) — ~50 flops and one square root per point.
 * @alloc none.
 * @test IrbemField.GradCurvCurlBatchMatchesThePointwiseRoutine
 */
inline Status grad_curv_curl_batch(std::span<const BDerivatives> derivs,
                                   std::span<GradCurvCurl> out, std::span<Status> statuses) {
    const std::size_t n = derivs.size();
    if (out.size() != n || statuses.size() != n) return Status::DomainError;
    Status worst = Status::Ok;
    for (std::size_t i = 0; i < n; ++i) {
        const Result<GradCurvCurl> r = grad_curv_curl(derivs[i]);
        out[i] = r.value;
        statuses[i] = r.status;
        if (worst == Status::Ok) worst = r.status;
    }
    return worst;
}

/**
 * The magnetic hemisphere of every point of a batch — IRBEM's `GET_HEMI_MULTI`.
 *
 * Three field evaluations per point — the base, and one step either way along `B̂` — and the same
 * one-dispatch discipline as @ref bderivs_batch, with one difference that is forced by the physics:
 * `B̂` is not known until the base evaluation has come back, so this cannot be a single `3N`
 * dispatch. It is **two**: one of `N` points to get the directions, then one of `2N` to get the
 * neighbours. Two dispatches for three evaluations is the best a data dependency of this shape
 * allows, and it is still a factor of `N` better than the per-point loop the name suggests.
 *
 * fp32 costs nothing here. Only the SIGN of the difference survives, and the sign is unambiguous
 * everywhere except within ~10⁻⁴ Re of the magnetic equator, where the two hemispheres are
 * genuinely adjacent and no precision decides the question. That is the one routine in this file
 * the device lane serves at no accuracy cost at all.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param points the points, GEO, Earth radii.
 * @param out receives one @ref Hemisphere per point; same length as @p points.
 * @param step_re the step along `B̂` in Earth radii; `0.0` selects @ref auto_step.
 * @return @ref Status::DomainError on a length mismatch or a non-finite step, @ref Status::Ok
 *         otherwise. The value is `true` when the device serviced the call.
 * @complexity O(3·N·NMAX²) field evaluations.
 * @alloc none on the host lane. The device lane stages `9N` floats each way.
 * @test IrbemField.HemisphereBatchAgreesWithTheReferenceLane
 * @test IrbemField.HemisphereBatchUsesTheDeviceWhenOneIsAvailable
 */
template <int NMAX>
[[nodiscard]] inline Result<bool> hemisphere_batch(const Igrf<NMAX>& model,
                                                   std::span<const Position<Frame::GEO>> points,
                                                   std::span<Hemisphere> out,
                                                   double step_re = 0.0) {
    const std::size_t n = points.size();
    if (out.size() != n) return {Status::DomainError, false};
    if (!std::isfinite(step_re) || step_re < 0.0) return {Status::DomainError, false};
    if (n == 0) return {Status::Ok, false};

#ifdef CHEATAH_SPACE_IRBEM_FIELD_GPU
    if (gpu::prefer_gpu("irbem_igrf_f32", 2 * n)) {
        std::vector<float> pos(3 * n);
        std::vector<float> base(3 * n);
        detail::interleave(points, pos);
        if (detail::igrf_on_device(model, pos, base)) {
            std::vector<float> probe(6 * n);
            std::vector<float> probe_b(6 * n);
            for (std::size_t i = 0; i < n; ++i) {
                const fixarray::vec3d bv{base[(3 * i) + 0], base[(3 * i) + 1], base[(3 * i) + 2]};
                const double bmag = fixarray::norm(bv);
                const double r = fixarray::norm(points[i].v);
                const double h = (step_re > 0.0) ? step_re
                                                 : auto_step(r, DifferenceLane::Fp32Device);
                const fixarray::vec3d step = bmag > 0.0 ? bv * (h / bmag) : fixarray::vec3d{};
                for (std::size_t k = 0; k < 3; ++k) {
                    probe[(3 * i) + k] = static_cast<float>(points[i].v[k] + step[k]);
                    probe[(3 * (n + i)) + k] = static_cast<float>(points[i].v[k] - step[k]);
                }
            }
            if (detail::igrf_on_device(model, probe, probe_b)) {
                for (std::size_t i = 0; i < n; ++i) {
                    const double bmag = fixarray::norm(fixarray::vec3d{
                        base[(3 * i) + 0], base[(3 * i) + 1], base[(3 * i) + 2]});
                    const double up = fixarray::norm(fixarray::vec3d{
                        probe_b[(3 * i) + 0], probe_b[(3 * i) + 1], probe_b[(3 * i) + 2]});
                    const double dn = fixarray::norm(fixarray::vec3d{probe_b[(3 * (n + i)) + 0],
                                                                     probe_b[(3 * (n + i)) + 1],
                                                                     probe_b[(3 * (n + i)) + 2]});
                    if (!(bmag > 0.0) || up == dn) {
                        out[i] = Hemisphere::Invalid;
                    } else {
                        out[i] = up > dn ? Hemisphere::North : Hemisphere::South;
                    }
                }
                return {Status::Ok, true};
            }
        }
    }
#endif

    for (std::size_t i = 0; i < n; ++i) { out[i] = hemisphere(model, points[i], step_re).value; }
    return {Status::Ok, false};
}

}  // namespace cheatah::space::irbem
