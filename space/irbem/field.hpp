#pragma once

/**
 * @file field.hpp
 * @brief space.irbem — the magnetic field at a point, its first derivatives, and everything
 *        algebraic that follows from them.
 *
 * @ref lstar.hpp is the compute-bound half of this module: one point in, ~10⁵ field evaluations,
 * six scalars out. This file is the other half — the **streaming** half, where the answer is as
 * large as the question and the work per point is a few thousand flops. Four IRBEM routines live
 * here, ordered by how much arithmetic they do per byte moved, because that ordering is the whole
 * design. Everything in this table is measured on an RTX 3070 Ti against this file's own fp64 host
 * lane, `-O3 -march=native -ffp-contract=off`, best of five, transfers included:
 *
 * | routine | evaluations / point | crossover | best measured speedup |
 * |---|---|---|---|
 * | `GET_FIELD_MULTI` (@ref field_batch) | 1 | 512 points | **21.8×** at 2¹⁴ |
 * | `GET_HEMI_MULTI` (@ref hemisphere_batch) | 3, in two dispatches | 256 points | **25.4×** at 2¹⁴ |
 * | `GET_BDERIVS` (@ref bderivs_batch) | 4, in ONE dispatch | 128 points | **26.1×** at 2¹⁴ |
 * | `COMPUTE_GRAD_CURV_CURL` (@ref grad_curv_curl_batch) | **0** | — | **host, always** |
 *
 * Those three crossovers are the same number wearing different clothes. This launcher's measured
 * floor is **~115 µs per dispatch** — read straight off the small-batch timings, where
 * @ref field_batch at 128 points costs 936 ns/point, @ref bderivs_batch at 16 points costs
 * 7.2 µs/point, and @ref hemisphere_batch at 16 points costs 14.4 µs/point, i.e. 120, 115 and
 * *230* µs respectively; the last is exactly twice the others because it is the one routine that
 * cannot avoid two dispatches. Divide that floor by the ~307 ns the host spends per field
 * evaluation and you get **~512 evaluations per dispatch** as the break-even for all three. A
 * routine that does four evaluations per point therefore reaches the device at a quarter of the
 * batch size, which is exactly what the table says.
 *
 * ## Why COMPUTE_GRAD_CURV_CURL is host-only, permanently
 *
 * It touches no field model at all — pure algebra over `GET_BDERIVS`'s outputs — and it is
 * therefore the one routine here a GPU can only make slower. **136 bytes per point would have to
 * cross the bus** (16 fp32 in, 18 out) for **19 ns of host arithmetic** (measured: 50–57 Mpts/s at
 * 2¹⁰ and 2¹⁸). The payload bandwidth of this seam is **7.15 GB/s**, measured directly by running
 * a degree-1 IGRF dispatch over 2²¹ points — 30 flops of kernel, so what is left is the round trip
 * — which puts the transfer for 136 B/point at **≥ 19.0 ns/point**. The copy alone costs exactly
 * what the whole computation costs, before the kernel runs and before the ~115 µs floor is paid.
 * Arithmetic intensity ~0.4 flops/byte, below even the dipole kernel's 0.5 — and `gpu/dispatch.hpp`
 * records that the dipole kernel LOSES 0.69× at every size it was measured at. There is no batch
 * size that recovers this; the ratio is fixed, not amortizable. No kernel for it exists in
 * `gpu/irbem.slang` and none should be written.
 *
 * ## The finite-difference step is the whole design of GET_BDERIVS
 *
 * `GET_BDERIVS` has no closed form here. IGRF's Jacobian *is* analytically available — the Legendre
 * recursion already carries `A'ⁿₘ` — but IRBEM's routine is defined as a **forward difference with
 * a caller-supplied step `dX`**, and reproducing it means differencing. That puts every derivative
 * in this file on one number, pulled in two directions:
 *
 *  - **too large** and truncation dominates. A one-sided difference has error `(h/2)·|∂²B|`, which
 *    for a field falling as `r⁻³` is a relative error of about `2h/r` on `∇|B|` — dead linear in
 *    the step, and measured as such.
 *  - **too small** and cancellation destroys it. Differencing two nearby field magnitudes loses
 *    every digit they share, so the roundoff term grows as `ε·r/h`; it *diverges* as `h` shrinks.
 *
 * The two cross at `h ∝ r·√ε`, which is why the automatic step here is **proportional to the
 * radius**: both terms scale with `r`, so one ratio serves LEO and the outer belt alike. And `ε` is
 * not the same on the two lanes — the device evaluates the whole 105-term series in fp32 and comes
 * back with a measured `7.3 × 10⁻⁷` relative error on `|B|`, not `6 × 10⁻⁸` — so the two lanes have
 * genuinely different optimal steps. Measured, absolute step, 2 000 points over `r = 1.5 … 8.5`,
 * error on `∇|B|` against a Richardson-extrapolated fp64 reference (max / median):
 *
 * | `dX` | fp32 device lane | fp64 host lane |
 * |---|---|---|
 * | 10⁻⁵ | 2.5e-01 / 4.8e-02 | 1.4e-05 / 2.9e-06 |
 * | 10⁻⁴ | 2.8e-02 / 4.8e-03 | 1.4e-04 / 2.9e-05 |
 * | 10⁻³ | **3.5e-03 / 7.1e-04** | 1.4e-03 / 2.9e-04 |
 * | 10⁻² | 1.4e-02 / 2.9e-03 | 1.4e-02 / 2.9e-03 |
 *
 * Read the bottom row: at `dX = 10⁻²` the two lanes agree to three digits, because truncation
 * swamps everything and precision is irrelevant. Read the top row: at `dX = 10⁻⁵` the device is
 * **four orders of magnitude worse** and returns 25% errors. There is no single default that
 * serves both, and this file does not pretend otherwise — @ref auto_step takes the lane, and the
 * two ratios are named after their lanes.
 *
 * Sweeping the step at fixed radius pins each lane's optimum and shows the `∝ r` scaling directly.
 * On the device, 600 points per shell:
 *
 * | shell | best `dX` | as a ratio `h/r` | max error there |
 * |---|---|---|---|
 * | r = 1.2 | 3 × 10⁻⁴ | 2.5e-04 | 1.5e-03 |
 * | r = 2.0 | 1 × 10⁻³ | 5.0e-04 | 1.5e-03 |
 * | r = 4.0 | 2 × 10⁻³ | 5.0e-04 | 1.4e-03 |
 * | r = 6.6 | 2 × 10⁻³ | 3.0e-04 | 1.2e-03 |
 * | r = 10  | 3 × 10⁻³ | 3.0e-04 | 1.4e-03 |
 *
 * A constant ratio, and a flat `~1.5 × 10⁻³` achievable accuracy across a factor of eight in
 * radius — hence @ref device_step_ratio. The same sweep in fp64 puts the host's trough at
 * `h ≈ r · 5 × 10⁻⁸` with max / median `1.1 × 10⁻⁷ / 7.3 × 10⁻⁸` over 60 points spanning
 * `r = 1.05 … 11.5` — hence @ref host_step_ratio. **The device lane is four orders of magnitude
 * less accurate on derivatives, and that is a property of the step it can afford, not of the
 * kernel.** It is stated here rather than averaged away.
 *
 * A caller who supplies `dX` explicitly gets exactly that step on both lanes, which is what a
 * differential comparison against IRBEM needs: at `dX = 10⁻³` the two implementations are compared
 * at matched resolution, the same discipline `docs/ERROR_BUDGET.md` §2(a) imposes on L\*.
 *
 * ## What the reference actually computes, established as a black box
 *
 * IRBEM is LGPL-3.0 and this repository is MIT, so its Fortran is never read. Its *behaviour* is
 * fair game, and the following were measured by running the compiled `-O2` oracle on inputs chosen
 * to separate the candidate definitions. Each is reproducible from `tools/oracle/` and each is
 * pinned by a test here:
 *
 *  - `GET_BDERIVS` is a **forward** (one-sided) difference, never central. Feeding the oracle's own
 *    `GET_FIELD_MULTI` outputs at `x` and `x + dX·ê_j` through `(B(x+dX·ê_j) − B(x))/dX` reproduces
 *    its `gradBmag` and `diffB` **bit for bit** — relative difference exactly `0.0`, at
 *    `dX = 10⁻¹, 10⁻², 10⁻³`, over 300 points. A central difference disagrees at the 10⁻³ level, so
 *    this is not an artefact of a loose tolerance.
 *  - `diffB(i,j) = ∂Bᵢ/∂xⱼ`, Fortran column-major, so component `i` is the fast index.
 *    @ref BDerivatives::diff_b indexes `(row, col) = (i, j)` — the same convention, written the way
 *    mathematics writes it.
 *  - `COMPUTE_GRAD_CURV_CURL`'s `curvature` is `Â − (Â·B̂)B̂` with `Â = (B̂·∇)B/|B|`: it projects
 *    out `Â`'s **own** parallel component and never reads `gradBmag`. For a consistent Jacobian
 *    that equals the documented `(B̂·∇)B̂` identically; for the *inconsistent* inputs a black-box
 *    probe can feed it, the candidates differ by a multiple of `B̂`, and the oracle follows this
 *    one — reproduced to `3.7 × 10⁻¹⁵` while subtracting `grad_par·B̂/|B|` instead is off by 10³.
 *  - All eight `COMPUTE_GRAD_CURV_CURL` outputs, **swept through storms**: 46 driver
 *    configurations — Kp 0…9 under T89, Dst 0…−400 nT under T96 and T01-storm, Pdyn 0.5…40 nPa,
 *    and Bz 0…−30 nT densely southward — × 300 points, **12 900 comparisons**. `grad_par`,
 *    `grad_perp`, `grad_drift`, `curl_b` and `div_b` reproduce at exactly `0.0`; `r_curv` at
 *    `5.9 × 10⁻¹⁴`, `curvature` at `7.7 × 10⁻¹²`, `curv_drift` at `9.1 × 10⁻¹²`. The sweep is what
 *    makes this meaningful for a module whose external field models are not written yet: the
 *    ALGEBRA is validated on the oracle's own storm-time `B` and Jacobian, independently of which
 *    model produced them, so it is already correct for the disturbed conditions the library exists
 *    for. What is *not* yet validated in a storm is the field itself — see the gap note below.
 *  - `GET_HEMI_MULTI` returns the sign of `d|B|/ds` along `+B̂`: the northern hemisphere is the side
 *    of the magnetic equator where the field is *rising* in the direction B points. Agreed with the
 *    oracle on **400 / 400** and again **300 / 300** random points spanning `r = 1.2 … 8.5 Re`.
 *  - `options(2) = 0` means "initialize IGRF once per year (year.5)" — the published options table
 *    says so, and it is why every oracle comparison here uses epoch **2015.5** exactly rather than
 *    a day-of-year fraction. With that epoch, @ref field_batch reproduces `GET_FIELD_MULTI` over
 *    300 points at:
 *
 * | truncation | 8 | 9 | **10** | 11 | 13 |
 * |---|---|---|---|---|---|
 * | max rel, vector | 1.5e-04 | 4.2e-05 | **1.7e-15** | 1.2e-05 | 1.5e-05 |
 *
 *    Ten orders of magnitude, in one column. That is how IRBEM's internal truncation was
 *    identified — by sweeping the degree and watching the residual collapse, not by reading a line
 *    of its source — and at degree 10 the two implementations agree to **floating-point noise**,
 *    nine orders inside the `1 × 10⁻⁶` `Bgeo` budget. Degree 13 is IGRF-14 as IAGA published it and
 *    its `1.5 × 10⁻⁵` is the *model* difference, not an error in either implementation. A
 *    differential test must say which truncation it ran; this one does.
 *  - `GET_BDERIVS` end to end: @ref bderivs_batch against the oracle's `gradBmag` and `diffB`,
 *    same 300 points, is `5.6 × 10⁻¹⁴` at `dX = 10⁻¹`, `4.4 × 10⁻¹³` at `10⁻²` and
 *    `4.5 × 10⁻¹²` at `10⁻³` — the growth is the fp64 cancellation term `ε/h` and nothing else,
 *    which is the same curve the step study measures from the other side.
 *
 * ## The gap, stated plainly
 *
 * There is no external field model in this module yet. Every *field* comparison above is at
 * `kext = 0`, internal field only. The storm sweep validates the derivative and guiding-centre
 * ALGEBRA under disturbed conditions, because that algebra is a function of `B` and its Jacobian
 * and does not care where they came from — but a claim that this module reproduces IRBEM at
 * `Kp = 7`, `Dst = −300` would today be a claim about arithmetic, not about magnetospheric physics.
 * It is not made here.
 *
 * ## Sources
 *
 *  - Alken et al., *International Geomagnetic Reference Field: the fourteenth generation*, Earth
 *    Planets Space (2025) — the field itself; evaluated by @ref Igrf.
 *  - Northrop, *The Adiabatic Motion of Charged Particles*, Interscience (1963), ch. 1 — the
 *    guiding-centre drifts these derivative products are the geometric factors of: the gradient
 *    drift `(B̂ × ∇⊥B)/B`, the curvature drift `B̂ × (B̂·∇)B̂`, and the radius of curvature. Both
 *    become velocities only after multiplication by `m v²/(qB)` — see @ref GradCurvCurl.
 *  - Roederer, *Dynamics of Geomagnetically Trapped Radiation*, Springer (1970), §1.3.
 */

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

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

// @ref GeoFieldModel moved to frames.hpp: the tracer in lstar.hpp needs it, and this header
// already sits on top of lstar.hpp, so defining it here would be a cycle. It is written purely in
// terms of Position and FieldVector, which frames.hpp owns.


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
 * Eight thousand times larger than @ref host_step_ratio, and that factor is `√(ε₃₂ᵉᶠᶠ/ε₆₄)` and
 * nothing else — where `ε₃₂ᵉᶠᶠ` is the **7.3 × 10⁻⁷** the fp32 kernel actually delivers on `|B|`,
 * an order of magnitude above a single fp32 rounding because the harmonic sum runs over 105 terms.
 * Assuming `6 × 10⁻⁸` instead lands the step three times too small and costs a factor of twenty in
 * the answer, which is how this constant was found.
 *
 * The ratio is flat across the belt: the per-shell sweep in the file brief puts the optimum between
 * `2.5 × 10⁻⁴` and `5 × 10⁻⁴` from `r = 1.2` to `r = 10`, with the achievable error a nearly
 * constant `1.2 … 1.5 × 10⁻³` throughout.
 */
inline constexpr double device_step_ratio = 4.0e-4;

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
// Lane selection — the routine's crossover, not the kernel's
// ---------------------------------------------------------------------------------------------

/**
 * The batch size at or above which @ref field_batch's device lane is measured to win.
 *
 * `gpu/dispatch.hpp`'s registry records **128** for `irbem_igrf_f32`, derived from the kernel's own
 * throughput and an assumed ~30 µs submit floor. This launcher's *measured* floor is ~115 µs — it
 * also packs two coefficient tables, acquires five buffers and does four uploads — so the routine's
 * crossover is four times the kernel's. Measured: 0.33× at 128 points, 0.62× at 256, **1.22× at
 * 512**, 2.39× at 1 024. Using the registry's number unmodified would make every batch between 128
 * and 450 points *slower* than the host, which is precisely the failure the per-kernel crossover
 * exists to prevent — so it is recorded here rather than papered over.
 */
inline constexpr std::size_t field_batch_crossover = 512;

/**
 * The batch size at or above which @ref bderivs_batch's device lane wins.
 *
 * A quarter of @ref field_batch_crossover, because each point is four field evaluations and they
 * all go in ONE dispatch: 128 points are 512 evaluations, and ~512 evaluations per dispatch is
 * where this seam breaks even. Measured: 0.64× at 64 points, **1.27× at 128**, 2.49× at 256,
 * 8.44× at 1 024, 26.1× at 2¹⁴.
 */
inline constexpr std::size_t bderivs_batch_crossover = 128;

/**
 * The batch size at or above which @ref hemisphere_batch's device lane wins.
 *
 * Between the other two, and for a reason worth naming: three evaluations per point would put it at
 * ~171, but `B̂` is not known until the first dispatch returns, so this routine pays the ~115 µs
 * floor **twice**. Measured: 0.51× at 128 points, 0.97× at 256, 3.76× at 1 024, 25.4× at 2¹⁴.
 */
inline constexpr std::size_t hemisphere_batch_crossover = 256;

namespace detail {

/**
 * Whether a routine that breaks even at @p crossover points should take the device for @p points.
 *
 * Two thresholds, and the order matters. `gpu::prefer_gpu` is consulted first because it is what
 * asks whether a device exists at all and what honours `CHEATAH_SPACE_IRBEM_GPU_CROSSOVER`; the
 * routine's own measured crossover is applied second, and **only when the operator has not set that
 * variable**. An explicit override is a decision — it is how the differential suite drives the
 * device lane at four points — and a constant compiled into this header must not silently veto it.
 *
 * @param points the batch size, in points of the calling routine.
 * @param crossover the routine's measured crossover.
 * @return true when the device lane should run.
 * @complexity O(number of registered kernels) — one `getenv` and a linear scan of
 *              @ref gpu::registered_kernels.
 * @alloc none.
 * @test IrbemField.CrossoverIsTheRoutinesNotTheKernels
 */
#ifdef CHEATAH_SPACE_IRBEM_FIELD_GPU
[[nodiscard]] inline bool prefer_device(std::size_t points, std::size_t crossover) {
    if (!gpu::prefer_gpu("irbem_igrf_f32", points)) return false;
    if (std::getenv("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER") != nullptr) return true;
    return points >= crossover;
}
#endif  // CHEATAH_SPACE_IRBEM_FIELD_GPU

}  // namespace detail

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
 * @test IrbemField.BderivsRefusesInputsItCannotAnswer
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
 * `(m v⊥²/2qB)·grad_drift` and the curvature drift is `(m v∥²/qB)·curv_drift` (Northrop 1963,
 * §1.3). The `1/B` is not decoration — `grad_drift` and `curv_drift` are both `1/Re`, so an
 * energy-over-charge multiplier alone does not make a velocity; it takes `m v²/(qB)`, which is
 * `V/T`, i.e. m²/s, over a length.
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
 * 19 ns/point on the host and rules off the device permanently.
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
 * @test IrbemField.DivergenceAndCurlAreTheDifferencingResidual
 * @test IrbemField.GradCurvCurlHandlesTheDegenerateCases
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
 * @param step_re the step along `B̂` in Earth radii; `0.0` selects @ref auto_step at
 *        @ref DifferenceLane::Fp32Device — the *device* ratio, deliberately, on the host lane too.
 *        Unlike @ref bderivs this is not accuracy-critical: only the sign of the difference
 *        survives, so the generous step costs nothing, and using one ratio on both lanes is what
 *        makes @ref hemisphere_batch's host and device answers comparable point for point instead
 *        of differing by four orders of magnitude in resolution.
 * @return the hemisphere, with @ref Status::DomainError for a non-finite or origin point and
 *         @ref Hemisphere::Invalid where the field vanishes or the derivative is exactly zero.
 * @complexity Three field evaluations.
 * @alloc none.
 * @test IrbemField.HemisphereAgreesWithTheOracleGoldens
 * @test IrbemField.HemisphereFlipsAcrossTheDipoleEquator
 * @test IrbemField.HemisphereRefusesInputsItCannotAnswer
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
 * @alloc two — `assign` sizes each vector, so each allocates once on a fresh vector. Once per
 *        DISPATCH, never per point: the host hot path measured by `valgrind --tool=memcheck` is
 *        flat at 9 allocations from 1 rep to 100.
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
 * @test IrbemField.AMissingShaderFallsBackToTheHostLane
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
 * The step a `float` lane ACTUALLY takes when asked to move `x` by `h`.
 *
 * `x + h` is not representable in `float`, so the device evaluates at `fl(fl(x) + h)` and the
 * difference it really formed is `fl(fl(x) + h) − fl(x)`, not `h`. Dividing by the intended `h`
 * would put a relative error of order `ε₃₂·|x|/h` straight onto every derivative — **measured at
 * 2.3 × 10⁻⁴ maximum over 4 096 points at the automatic device step**, against a device lane whose
 * whole achievable accuracy is ~1.5 × 10⁻³. It is a sixth of the budget, and it costs one
 * subtraction to remove, so it is removed.
 *
 * Kept as a named function rather than inlined into @ref bderivs_batch precisely so the claim is
 * checkable: from inside the device lane the correction is buried under the fp32 field noise it is
 * smaller than, and a test could not tell it had been deleted.
 *
 * @param x the base coordinate, Earth radii.
 * @param h the intended step, Earth radii.
 * @return the realised step, or @p h when the two `float` positions collapse onto each other (an
 *         `h` so small that it vanishes into `x`'s last bit), which keeps the caller's division
 *         finite instead of producing an infinity.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemField.TheRealisedStepIsTheOneTheFloatLaneCanTake
 */
[[nodiscard]] inline double realised_step(double x, double h) {
    const double taken = static_cast<double>(static_cast<float>(x + h)) -
                         static_cast<double>(static_cast<float>(x));
    return taken != 0.0 ? taken : h;
}

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
 * *least* arithmetically intense thing here worth offloading — ~1 900 flops for 24 bytes in and 24
 * out, about 20 flops/byte — which is exactly why the crossover is consulted rather than assumed.
 * Measured on an RTX 3070 Ti against this file's own fp64 host lane, `-O3 -march=native
 * -ffp-contract=off`, best of five, transfers included:
 *
 * | points | 128 | 256 | 512 | 1024 | 2048 | 2¹⁴ | 2²⁰ |
 * |---|---|---|---|---|---|---|---|
 * | speedup | 0.33× | 0.62× | **1.22×** | 2.39× | 4.56× | **21.8×** | 17.4× |
 *
 * — 14.2 ns/point on the device at 2¹⁴ against the host's 309. Below @ref field_batch_crossover the
 * ~115 µs dispatch floor alone exceeds the whole computation and the host wins; a per-point loop
 * calling @ref Igrf::evaluate can never reach the device at all, which is the reason this entry
 * point takes the whole batch. The curve turns over slightly at 2²⁰ because the routine is
 * *staging*-bound by then, not kernel-bound: the same dispatch measured in isolation costs 5.9
 * ns/point at 2²¹, so the remaining ~12 ns is this file converting fp64 positions down to fp32 and
 * fp32 fields back up, which is inherent to the typed API and not to the device.
 *
 * The device lane returns fp32. Measured maximum relative deviation against the fp64 host lane:
 * **7.3 × 10⁻⁷** over 2 000 points and **1.1 × 10⁻⁶** over 2²⁰ — so at the larger sample it
 * *just* exceeds the `1 × 10⁻⁶` `Bgeo` budget of `docs/ERROR_BUDGET.md` §4. That is stated rather
 * than rounded away: it is the cost of summing 105 harmonic terms in fp32, it grows like the tail
 * of a distribution as the sample grows, and the fix if a caller needs the budget honoured at every
 * point is the host lane, not a different kernel. Nothing here accumulates *across* points — one
 * point, one thread, no reduction — so the budget's reduction concern still bites only the
 * integrals in @ref lstar.hpp.
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
 * @test IrbemField.BatchRoutinesRefuseMismatchedSpans
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
    if (detail::prefer_device(n, field_batch_crossover)) {
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
 * points each. That pays the ~115 µs dispatch floor four times for work that has no dependency
 * between the four groups whatsoever. Building the `4N` points up front and dispatching once pays
 * it once, and at the same time quadruples the occupancy of the launch — which is what pulls the
 * crossover down to @ref bderivs_batch_crossover, a **quarter** of @ref field_batch_crossover.
 * Measured: 0.64× at 64 points, 1.27× at 128, 2.49× at 256, 8.44× at 1 024, **26.1×** at 2¹⁴
 * (49.9 ns/point against the host's 1.30 µs). Four dispatches instead of one would have moved that
 * crossover to 512 points and left everything below it on the host.
 *
 * The layout is `[all N base points][all N +x][all N +y][all N +z]` rather than four consecutive
 * points per input point. Same dispatch either way, but this way each of the four groups is a
 * contiguous run whose lanes read neighbouring coefficients in the same order, and the host-side
 * differencing walks four unit-stride streams instead of one stride-4 one.
 *
 * On the device the differenced values are fp32, and the file brief's tables are the consequence:
 * the achievable accuracy is ~1.5 × 10⁻³ relative against ~1 × 10⁻⁷ on the host, and the step that
 * achieves it is eight thousand times larger. That four-order gap is why @p step_re defaults to
 * *the lane's* step rather than to a constant, and why a caller comparing the two lanes must pass
 * an explicit step to both — at a matched `dX = 10⁻²` they agree to three digits, and at a matched
 * `dX = 10⁻⁵` the device is wrong by 25%. It is a real limitation, quantified here rather than
 * hidden behind an average.
 *
 * One detail that is worth more than it looks: the differencing divides by the step the device
 * **actually took**, recovered from the fp32 position it was handed, not by the step that was
 * intended. `x + h` is not representable in fp32, and at `h = 4 × 10⁻⁴·r` the rounding is a
 * relative error on `h` of order `ε₃₂·r/h ≈ 1.5 × 10⁻⁴` — comparable to everything else in this
 * routine's budget, and free to remove.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param points the points, GEO, Earth radii.
 * @param out receives one @ref BDerivatives per point; same length as @p points.
 * @param step_re the difference step `dX` in Earth radii; `0.0` (the default) selects
 *        @ref auto_step for whichever lane actually runs.
 * @return @ref Status::DomainError on a length mismatch or a non-finite step; @ref Status::Ok
 *         otherwise. The value is `true` when the device serviced the call. A single point the
 *         pointwise @ref bderivs would refuse — the origin, a non-finite coordinate — is zero-filled
 *         rather than spoiling the batch, **and both lanes zero-fill the same point**: the device
 *         lane applies the host's gate before differencing, because without it the device
 *         differences a field that is infinite one step inside the Earth and returns a NaN gradient
 *         where the host returns zeros.
 * @complexity O(4·N·NMAX²) field evaluations plus ~40 flops per point of differencing.
 * @alloc none on the host lane. The device lane stages `12N` floats in and `12N` out.
 * @test IrbemField.BderivsBatchAgreesWithTheReferenceLane
 * @test IrbemField.BderivsBatchMatchesTheOracleGoldens
 * @test IrbemField.BderivsBatchUsesTheDeviceWhenOneIsAvailable
 * @test IrbemField.BothDerivativeLanesRefuseTheSameDegeneratePoint
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
    if (detail::prefer_device(n, bderivs_batch_crossover)) {
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
                // The same gate the host lane applies, for the same reason and with the same
                // answer. Without it the two lanes DIVERGE at a degenerate point: the host's
                // bderivs refuses the origin and zero-fills, while the device would difference a
                // field that is infinite one step away and hand back a NaN gradient that then
                // propagates silently through grad_curv_curl into whatever the caller does next.
                // Measured before this guard existed: host b_mag 0, grad 0; device b_mag 0,
                // grad NaN, at the same point of the same batch.
                const double ri = fixarray::norm(points[i].v);
                if (!(ri > 0.0) || !std::isfinite(ri)) continue;
                d.b = FieldVector<Frame::GEO>{fixarray::vec3d{field[(3 * i) + 0],
                                                              field[(3 * i) + 1],
                                                              field[(3 * i) + 2]}};
                d.b_mag = d.b.magnitude();
                // The step the DEVICE actually took, not the one it was asked for -- see
                // detail::realised_step, whose measured worth at the automatic device step
                // (h = 4e-4*r) is 2.3e-4 relative on every derivative below.
                for (std::size_t j = 0; j < 3; ++j) {
                    const std::size_t base = 3 * (((j + 1) * n) + i);
                    const fixarray::vec3d bj{field[base + 0], field[base + 1], field[base + 2]};
                    const double hj = detail::realised_step(points[i].v[j], steps[i]);
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
 * **A host loop, permanently, and the number that settles it.** 136 bytes would have to cross the
 * bus per point (16 fp32 in, 18 out) for **19 ns of host arithmetic** — measured, 17.7 ns/point at
 * 2¹⁰ and 19.9 at 2¹⁸, i.e. 50–57 Mpts/s, `-O3 -march=native -ffp-contract=off`. This seam's
 * payload bandwidth is **7.15 GB/s**, measured directly by dispatching a *degree-1* IGRF kernel
 * over 2²¹ points — 30 flops of kernel, so what is left is the round trip — which puts 136 B/point
 * at **≥ 19.0 ns/point of transfer alone**. The copy costs exactly what the computation costs,
 * before the kernel runs and before the ~115 µs dispatch floor is paid. Arithmetic intensity is
 * ~0.4 flops/byte, below even the dipole kernel's 0.5 — and `gpu/dispatch.hpp` already records
 * that the dipole kernel LOSES 0.69× at every batch size it was measured at. There is no size at
 * which this one wins: the ratio is fixed, not amortizable. So no kernel for it exists in
 * `gpu/irbem.slang` and none should be written.
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
 * fp32 costs almost nothing here, because only the SIGN of the difference survives. Measured over
 * 2¹⁸ random points spanning `r = 1.5 … 8.5`, the device lane and the fp64 host lane disagree on
 * **11 points in 262 144** — 4 × 10⁻⁵ — and every one of them is a point sitting within a step of
 * the magnetic equator, where `d|B|/ds` passes through zero and the two hemispheres are genuinely
 * adjacent. No precision decides that question; a point *on* the equator is in neither hemisphere.
 * Measured speedup: 0.51× at 128 points, 0.97× at 256, 3.76× at 1 024, **25.4×** at 2¹⁴.
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
 * @test IrbemField.TheDeviceLaneReportsAPointWithNoField
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
    if (detail::prefer_device(n, hemisphere_batch_crossover)) {
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
