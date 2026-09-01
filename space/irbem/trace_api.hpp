#pragma once

/**
 * @file trace_api.hpp
 * @brief space.irbem — the routines that RETURN the field line: the traces and the point finders.
 *
 * `lstar.hpp` traces a field line and throws the path away. That is what makes it fast: the point
 * goes to the device once, ~10⁵ field evaluations happen on it, and six scalars come back —
 * ~9 400 flops per byte moved, which is why the device beats the host by 48.9× at 65 536 lines.
 *
 * **These five routines are the opposite shape, and the file exists to keep the two apart.** A
 * trace that returns `posit` writes four doubles per step per line — position and `|B|` — up to
 * IRBEM's `posit(3,3000)`. A single L=4 line at the module's default step is ~120 samples, ~3.8 kB;
 * the same line through @ref trace_invariant returns 48 bytes. The arithmetic per step is
 * unchanged, so the arithmetic intensity falls by two orders of magnitude — ~9 400 flops/byte to
 * ~125 — and the routines land on the *other* side of the roofline. What that costs is MEASURED in
 * @ref trace_field_line_toward_earth_batch, and it is not where the guess said it would be: the
 * device crossover came out at 256 lines rather than higher than the tracer's 512, because a
 * crossover is a ratio and this routine's HOST lane is 2.6× dearer per line too. The bandwidth
 * shows up in the CEILING instead — the invariant tracer's speedup is still climbing at 65 536
 * lines, this one flattens in the low twenties and goes noisy above 4 096.
 *
 * Two consequences are structural, not stylistic:
 *
 *  - **These never share the L\* hot path.** The `irbem_trace_i_f32` kernel stores no path
 *    precisely so that it stays compute-bound; teaching it to emit `posit` would destroy the one
 *    property the module's performance argument rests on. The device lane here is a *separate*
 *    kernel with its own registry row and its own measured crossover.
 *  - **The caller owns the buffer.** Every routine writes into a `std::span` the caller supplies
 *    and allocates nothing. A 3 000-sample path is 72 kB; a 4 096-line batch of them is 288 MB.
 *    Sizing that is the caller's decision, and a routine that made it silently would be the wrong
 *    place for it.
 *
 * ## What the five routines are
 *
 * Semantics follow `vendor/IRBEM/docs/source/api/magnetic_coordinates.rst` and were pinned against
 * the compiled oracle run as a black box (see `tests/irbem_trace_api_test.cpp`):
 *
 *  - @ref trace_field_line — the whole line, foot to foot at the reference surface `R0`, plus
 *    `Lm`, `Bmin` and `XJ` for a particle mirroring at the input point.
 *  - @ref trace_field_line_toward_earth — from the input point to the Earth at a FIXED step, in
 *    the input point's own magnetic hemisphere. For plotting.
 *  - @ref find_mirror_point — where a particle of local pitch angle α turns around.
 *  - @ref find_magequator — the minimum-`|B|` point, sub-step resolved.
 *  - @ref find_foot_point — where the line crosses a given GEODETIC altitude, in a chosen
 *    hemisphere.
 *
 * ## Three sub-step refinements, because the step is 0.08 R_E and the answers are points
 *
 * The reference implementation reports the coarse RK4 sample it happened to land on for `Bmin`
 * only after a sub-step refinement of its own, and this header does the same everywhere an answer
 * is a *location* rather than an integral:
 *
 *  - `Bmin` and the magnetic equator come from a **parabolic fit in arc length** through the three
 *    samples bracketing the minimum. The coarse minimum is wrong in the value at O(ds²) and in the
 *    position at O(ds): on an L=4 line the coarse sample sits 5.5 × 10⁻⁶ relative above the
 *    oracle's `Bmin` and 4.4 × 10⁻³ R_E away from where it is. Deleting the fit and re-running the
 *    differential sweep moves the equator position from 4.6 × 10⁻⁵ R_E of the oracle to
 *    1.03 × 10⁻² — **224× worse**, measured, which is the answer to "is this refinement earning
 *    its lines".
 *  - The mirror point, the foot point and the `R0` end caps come from **regula falsi inside the
 *    last step** — on `|B|`, on geodetic altitude and on radius respectively — evaluating the true
 *    field at each trial rather than interpolating between the bracketing samples. Four Illinois
 *    iterations reach an exact zero on the `R0` cap; the shipped setting is eight.
 *
 * @note @ref find_foot_point terminates on GEODETIC altitude, not radius, so it converts inside the
 *       step loop through `coords_geodetic.hpp`. That is not a detail: the WGS-84 ellipsoid is
 *       21.4 km flatter at the poles than at the equator, and field-line feet are at |latitude| >
 *       45°, so terminating on geocentric radius instead would miss the requested altitude by up to
 *       ~15 km — roughly 30× the ~0.5 km the oracle itself converges to.
 *
 * ## Sources
 *
 *  - PRBEM/IRBEM, *Magnetic coordinates and fields*, `docs/source/api/magnetic_coordinates.rst` —
 *    the argument lists, the units, the `hemi_flag` codes and the `posit(3,3000)` cap.
 *  - Roederer, *Dynamics of Geomagnetically Trapped Radiation*, Springer (1970), ch. 2 — the
 *    invariants; `XJ` here is `I` for a particle mirroring at the input point, which is the
 *    convention the oracle's `make_lstar` and `trace_field_line` both use (verified: they agree to
 *    2 × 10⁻³ relative at L≈4, their own discretization difference).
 *  - Hilton, *L parameter: a new approximation*, J. Geophys. Res. 76(28):6952 (1971) — `Lm`, via
 *    @ref mcilwain_l.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <span>
#include <vector>

#include "coords_geodetic.hpp"
#include "frames.hpp"
#include "igrf.hpp"
#include "lstar.hpp"
#include "status.hpp"

namespace cheatah::space::irbem {

/// IRBEM's `posit(3,3000)` cap — the largest path its Fortran entry points can return, and the
/// default ceiling here so a caller porting from it needs no new number.
inline constexpr std::size_t irbem_max_path_points = 3000;

/**
 * One sample along a traced field line: where it is, and how strong the field is there.
 *
 * Position and magnitude together, rather than IRBEM's two parallel arrays (`posit(3,N)` and
 * `Blocal(N)`), because every consumer of a path — a plot, a bounce-period quadrature, a loss-cone
 * test — reads both at the same index. Interleaving halves the number of cache lines that walk
 * touches, which for a bandwidth-bound routine is the whole game.
 */
struct PathPoint {
    /// Where, in GEO Cartesian, Earth radii.
    Position<Frame::GEO> position;
    /// `|B|` there, in nanotesla.
    double b_magnitude = 0.0;
};

/**
 * Which foot @ref find_foot_point should walk to.
 *
 * The enumerator values ARE IRBEM's `hemi_flag` codes, so a caller porting from the Fortran can
 * pass its integer through `static_cast` and a reader of either can check the other.
 *
 * @note Distinct from @ref Hemisphere, which reports the hemisphere a point IS in (IRBEM's `xHEMI`).
 *       This one SELECTS a destination and carries two values `xHEMI` has no meaning for — `Same`
 *       and `Opposite` are relative to the input point. They were both called `Hemisphere` until
 *       the umbrella header included both and the compiler rejected the pair.
 */
enum class FootTarget : int {
    Same = 0,       ///< The input point's own magnetic hemisphere — the direction of increasing `|B|`.
    North = 1,      ///< The northern foot.
    South = -1,     ///< The southern foot.
    Opposite = 2,   ///< The hemisphere the input point is NOT in.
};

/**
 * Tuning for the path-returning traces.
 *
 * Defaults match `lstar.hpp`'s @ref TraceOptions where the two overlap, so a `Bmin` from
 * @ref find_magequator and a `Bmin` from @ref trace_invariant are the same number computed on the
 * same grid rather than two resolutions of the same idea.
 */
struct PathTraceOptions {
    /// The reference surface the trace terminates on, Earth radii — IRBEM's `R0`. Values below 1
    /// are legal and meaningful: a drift-loss-cone study wants the line followed below the surface.
    double r0 = 1.0;
    /// Steps per unit L, as @ref TraceOptions::steps_per_l. `ds = L/steps_per_l`, with L the start
    /// point's dipole L — see `lstar.hpp`'s brief on why the step is proportional to L.
    double steps_per_l = 50.0;
    /// An ABSOLUTE step in Earth radii, overriding @ref steps_per_l when positive. This is what
    /// IRBEM's `TRACE_FIELD_LINE_TOWARD_EARTH` takes as its `ds`, and the only reason it exists:
    /// a caller who wants uniform sampling for a plot should not have to work backwards from L.
    double step_size = 0.0;
    /// The cap on RK4 steps in ONE direction. Reaching it is @ref Status::OpenFieldLine.
    int max_steps = 4000;
    /// Regula-falsi iterations inside the final step, for every sub-step refinement here. Eight
    /// takes the residual below the fp64 noise of the field evaluation itself; the cost is eight
    /// field evaluations against a trace's several hundred.
    int refine_iterations = 8;
};

/**
 * What @ref trace_field_line produces, alongside the path itself.
 *
 * The scalars are IRBEM's `Lm`, `Bmin` and `XJ`, computed on the SAME samples the path is made of
 * rather than by a second trace — which is what makes this one pass over the line instead of two.
 */
struct TracedLine {
    /// `|B|` at the input point, nT — IRBEM's `Blocal` at the input index.
    double b_local = 0.0;
    /// The minimum `|B|` along the line, nT, parabola-refined off the sample grid.
    double b_min = 0.0;
    /// `I` for a particle mirroring at the input point, Earth radii — IRBEM's `XJ`.
    double invariant_i = 0.0;
    /// McIlwain's `L`, Earth radii — IRBEM's `Lm`, via @ref mcilwain_l.
    double mcilwain_l = 0.0;
    /// The magnetic equator, GEO, Earth radii — the parabola-refined position of @ref b_min.
    Position<Frame::GEO> equator{};
    /// How many samples were written into the caller's span.
    std::size_t point_count = 0;
    /// The index of the input point within the path. Exact: the input point is stored verbatim.
    std::size_t start_index = 0;
    /// The index of the sample nearest @ref equator — the parabola's centre, not its vertex.
    std::size_t equator_index = 0;
    /// Whether the path ran out of caller-supplied room before the line closed; the status is then
    /// @ref Status::NotConverged rather than @ref Status::Ok.
    ///
    /// There are TWO shapes of truncation and they leave different amounts behind, so a caller
    /// must not read the scalars without checking @ref start_index:
    ///
    ///  - The BACKWARD half alone filled the buffer. The trace returns before the input point is
    ///    even stored, so @ref start_index is 0 and @ref b_min, @ref invariant_i and
    ///    @ref mcilwain_l are all **exactly zero** — nothing was integrated and nothing is
    ///    reported. Reading them as a short field line would be reading zeros.
    ///  - The backward half FITTED and the forward one did not. Then the quadrature ran over the
    ///    partial line and the scalars describe THAT, not the field line.
    bool truncated = false;
};

/// Where a particle of a given local pitch angle turns around — IRBEM's `FIND_MIRROR_POINT`.
struct MirrorPoint {
    /// `|B|` at the input point, nT.
    double b_local = 0.0;
    /// `|B|` at the mirror point, nT. Exactly `b_local / sin²α` by the first adiabatic invariant;
    /// the trace finds WHERE that value is reached, it does not recompute WHAT it is.
    double b_mirror = 0.0;
    /// The mirror point, GEO, Earth radii.
    Position<Frame::GEO> position{};
};

/// The minimum-`|B|` point on the line — IRBEM's `FIND_MAGEQUATOR`.
struct MagneticEquator {
    /// The minimum `|B|`, nT, parabola-refined off the sample grid.
    double b_min = 0.0;
    /// Where it is, GEO, Earth radii, from a partial RK4 step to the parabola's vertex.
    Position<Frame::GEO> position{};
};

/// Where the line crosses a chosen geodetic altitude — IRBEM's `FIND_FOOT_POINT`.
struct FootPoint {
    /// The foot, GEODETIC: altitude km, latitude deg, east longitude deg. IRBEM returns `XFOOT` in
    /// GDZ too, which is the one place in its API where an output frame differs from its input's.
    Position<Frame::GDZ> position{};
    /// The field vector there, GEO, nT — IRBEM's `BFOOT`.
    FieldVector<Frame::GEO> field{};
    /// Its magnitude, nT — IRBEM's `BFOOTMAG`. Stored rather than recomputed because every caller
    /// wants it and `fixarray::norm` is a square root.
    double b_magnitude = 0.0;
};

namespace detail {

/// The largest fraction of the current geocentric radius one RK4 step may cover.
///
/// `dipole_l` is `r / cos²(latitude)`, so the L-proportional step blows up toward the poles: at
/// |latitude| 77° and r = 2 R_E it is 0.79 R_E, and a single step of that length can land inside
/// the Earth or on the far side of it, at which point every subsequent quantity is fiction. This
/// binds only above |latitude| ≈ 68.5° (where `1/cos²lat > 8/steps_per_l` at the default 50) and
/// leaves every step in the belts untouched.
inline constexpr double max_step_fraction_of_radius = 1.0 / 8.0;

/// The step length a trace from @p start uses, in Earth radii.
///
/// `L/steps_per_l` — `lstar.hpp`'s rule, so the two tracers walk the same grid — with the polar
/// clamp above. **The clamp is the one place these routines and `lstar.hpp` can disagree about the
/// grid**, and only above |latitude| ≈ 68.5°, where the unclamped step would be longer than an
/// eighth of the way to the centre of the Earth.
///
/// @param start the starting position, GEO, Earth radii.
/// @param opt the options; @ref PathTraceOptions::step_size wins when positive, unclamped, because
///        a caller naming an absolute step has already decided.
/// @return the unsigned step length; never zero, so a trace cannot stall.
/// @complexity O(1). @alloc none.
/// @test IrbemTraceApi.StepSizeFollowsTheDipoleLUnlessOverridden
/// @test IrbemTraceApi.StepSizeIsClampedNearThePoles
[[nodiscard]] inline double path_step(const Position<Frame::GEO>& start,
                                      const PathTraceOptions& opt) {
    if (opt.step_size > 0.0) return opt.step_size;
    const double per_l = opt.steps_per_l > 0.0 ? opt.steps_per_l : 50.0;
    return std::min(dipole_l(start) / per_l,
                    fixarray::norm(start.v) * max_step_fraction_of_radius);
}

/// The signed step that walks toward INCREASING `|B|` from @p p — the direction of the near mirror
/// point, the near foot, and "the same magnetic hemisphere".
///
/// Probe steps rather than a gradient: the walk is along the field line, so what matters is which
/// way the LINE's field rises, and a `∇|B|` at the point answers a different question wherever the
/// line is curved.
///
/// **Both directions are probed and the LARGER wins, and that is load-bearing.** Comparing one
/// probe against the start point instead is the obvious one-probe version and it is wrong within a
/// step of the magnetic equator, where both neighbours are higher than the start: it then falls
/// back on an arbitrary tie-break and picks the wrong hemisphere half the time. `|B|` falls
/// monotonically toward the equator along a field line, so the neighbour with the LARGER `|B|` is
/// always the one further from it — which is the answer in the generic case and in the degenerate
/// one alike. This was not hypothetical: at (−7.492, −2.806, 0) in 2015 the magnetic equator sits
/// 0.006 R_E away, the one-probe rule sent `FIND_FOOT_POINT`'s "same hemisphere" 130° of latitude
/// to the wrong foot, and the differential sweep caught it. A genuine tie — a start exactly on the
/// minimum — has no right answer and returns `+ds_mag`.
///
/// @tparam NMAX the IGRF truncation degree.
/// @param model the internal field model. @param p the position, GEO, Earth radii.
/// @param b the field at @p p, nT, already known. @param ds_mag the unsigned step length.
/// @return `+ds_mag` or `-ds_mag`.
/// @complexity Eight IGRF evaluations — two RK4 steps — against a trace's several hundred.
/// @alloc none.
/// @test IrbemTraceApi.IncreasingFieldStepPointsAwayFromTheEquator
/// @test IrbemTraceApi.IncreasingFieldStepIsRightBesideTheMagneticEquator
template <int NMAX>
[[nodiscard]] inline double increasing_field_step(const Igrf<NMAX>& model,
                                                  const Position<Frame::GEO>& p,
                                                  const fixarray::vec3d& b, double ds_mag) {
    fixarray::vec3d up{};
    fixarray::vec3d down{};
    (void)rk4_step(model, p, b, ds_mag, up);
    (void)rk4_step(model, p, b, -ds_mag, down);
    return fixarray::norm(up) >= fixarray::norm(down) ? ds_mag : -ds_mag;
}

/// A point found inside one RK4 step, with the field there.
struct RefinedPoint {
    Position<Frame::GEO> position{};  ///< Where, GEO, Earth radii.
    fixarray::vec3d field{};          ///< The field there, nT.
    double arc = 0.0;                 ///< The signed arc from the step's origin, Earth radii.
};

/// Regula falsi on `f` inside a single RK4 step, evaluating the true field at every trial.
///
/// The shared refinement behind three answers that are LOCATIONS: the mirror point (`f = |B| −
/// B_m`), the foot point (`f = altitude − stop_alt`) and the `R0` end caps (`f = |r| − R0`).
/// Interpolating between the two bracketing samples instead would be one line shorter and wrong by
/// O(ds²) — 6 × 10⁻³ R_E at the module's default step, against the ~10⁻¹⁰ this leaves.
///
/// The Illinois modification is not decoration: `|B|` along a field line is strongly convex near a
/// mirror point, so plain regula falsi can retain one endpoint for every iteration and converge
/// linearly with a ratio near 1. Halving the STALE endpoint's value restores superlinearity, which
/// is what lets eight iterations be enough — measured residuals are 1e-13 R_E on the `R0` caps and
/// 3e-11 nT on a mirror field.
///
/// The halving fires only on a REPEATED retention, which is the whole of Illinois and was worth
/// getting right: halving unconditionally penalises an endpoint that was just replaced, and the
/// iteration degenerates into exact bisection — error halving every step, alternating in sign,
/// stalling at ~1e-6 after eight iterations instead of reaching round-off. That was the first
/// version of this function, and the `R0` end-cap test caught it at 2.5e-6 R_E (16 metres).
///
/// A trial is clamped into the bracket. It should never need to be — every caller here arrives
/// with a genuine sign change — but an extrapolating secant on a bracket that does not contain a
/// root walks arbitrarily far outside the step and returns a point that is not on the field line
/// at all, which is a much worse failure than returning the step's far end.
///
/// @tparam NMAX the IGRF truncation degree.
/// @tparam F the predicate type; called as `f(Position<Frame::GEO>, fixarray::vec3d) -> double`.
/// @param model the internal field model.
/// @param p the step's origin, GEO, Earth radii. @param b the field at @p p, nT.
/// @param ds the signed step length that brackets the root.
/// @param f_lo `f` evaluated at @p p; must be nonzero and of opposite sign to `f` at the far end.
/// @param f_hi `f` evaluated at the far end of the step.
/// @param f the predicate whose root is sought.
/// @param iterations how many trials to take.
/// @return the refined point and the arc at which it sits.
/// @complexity `iterations` RK4 steps — four IGRF evaluations each. @alloc none.
/// @test IrbemTraceApi.RefineInStepFindsTheRadiusToRoundoff
template <int NMAX, class F>
[[nodiscard]] inline RefinedPoint refine_in_step(const Igrf<NMAX>& model,
                                                 const Position<Frame::GEO>& p,
                                                 const fixarray::vec3d& b, double ds, double f_lo,
                                                 double f_hi, F&& f, int iterations) {
    double lo = 0.0;
    double hi = ds;
    int retained = 0;   // +1 when `lo` survived the last iteration, -1 when `hi` did
    RefinedPoint best{p, b, 0.0};
    for (int i = 0; i < iterations; ++i) {
        const double denom = f_hi - f_lo;
        if (!(std::abs(denom) > 0.0)) break;
        const double raw = lo - (f_lo * (hi - lo) / denom);
        const double t = std::clamp(raw, std::min(lo, hi), std::max(lo, hi));
        fixarray::vec3d b_t{};
        const Position<Frame::GEO> q = rk4_step(model, p, b, t, b_t);
        const double f_t = f(q, b_t);
        best = RefinedPoint{q, b_t, t};
        if (!std::isfinite(f_t)) break;
        if ((f_t < 0.0) == (f_lo < 0.0)) {
            lo = t;
            f_lo = f_t;
            if (retained > 0) f_hi *= 0.5;   // Illinois: only a STALE endpoint's weight decays
            retained = 1;
        } else {
            hi = t;
            f_hi = f_t;
            if (retained < 0) f_lo *= 0.5;
            retained = -1;
        }
    }
    return best;
}

/// The vertex of the parabola through three equally spaced samples of `|B|`.
///
/// `Bmin` is a minimum of a smooth function sampled at spacing `ds`, so the sample nearest it is
/// wrong by O(ds²) in the value and O(ds) in the position. Three samples determine the parabola
/// exactly; its vertex is the standard three-point formula
/// `δ = (ds/2)·(b₋ − b₊)/(b₋ − 2b₀ + b₊)` with value `b₀ − (b₋ − b₊)²/(8·(b₋ − 2b₀ + b₊))`.
///
/// @param b_prev the sample one step before the minimum, nT.
/// @param b_centre the smallest sample, nT.
/// @param b_next the sample one step after it, nT.
/// @param ds the sample spacing, Earth radii.
/// @param arc receives the signed arc from the centre sample to the vertex.
/// @return the fitted minimum, or @p b_centre when the three samples are collinear (a flat or
///         non-convex bracket), in which case @p arc is set to zero and the coarse sample stands.
/// @complexity O(1). @alloc none.
/// @test IrbemTraceApi.ParabolicVertexIsExactOnAParabola
[[nodiscard]] inline double parabolic_minimum(double b_prev, double b_centre, double b_next,
                                              double ds, double& arc) {
    const double curvature = b_prev - (2.0 * b_centre) + b_next;
    const double slope = b_prev - b_next;
    if (!(curvature > 0.0)) {
        arc = 0.0;
        return b_centre;
    }
    arc = 0.5 * ds * slope / curvature;
    return b_centre - ((slope * slope) / (8.0 * curvature));
}

}  // namespace detail

/**
 * The minimum-`|B|` point on the field line through @p start — IRBEM's `FIND_MAGEQUATOR`.
 *
 * Walks downhill in `|B|` until the field stops falling, then **fits a parabola in arc length**
 * through the three samples that bracket the minimum and takes a partial RK4 step to its vertex.
 * The coarse sample is not the answer and returning it would be a silent O(ds) error in the
 * position: measured on an L≈4 line at the default step, the coarse minimum sits 5.5 × 10⁻⁶
 * relative above the oracle's `Bmin` and 4.4 × 10⁻³ R_E from its position.
 *
 * A start point that is ALREADY at the minimum is not a special case in the physics and is not one
 * here: both neighbours are probed, and if both are higher the bracket is centred on the start.
 *
 * Measured against the oracle over 84 start points × 4 epochs at matched IGRF: `Bmin` to
 * 1.9 × 10⁻⁶ relative (budget 10⁻⁵) and the position to 4.6 × 10⁻⁵ R_E. The `Bmin` residual is the
 * ORACLE's, not ours — our value at that resolution is within 10⁻⁸ of our own converged value, and
 * the oracle sits the same 1.9 × 10⁻⁶ from it.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param start the starting position, GEO, Earth radii.
 * @param opt the tracing options.
 * @return the equator, with @ref Status::OpenFieldLine when the walk reached
 *         @ref PathTraceOptions::r0 or the step cap without finding a minimum, and
 *         @ref Status::DomainError for a start inside `r0` or a vanishing field. The value is
 *         populated in every case: the lowest field actually seen is still diagnostic.
 * @complexity O(steps) IGRF evaluations, ~4 per step, plus the fit; typically 50–200 steps.
 * @alloc none.
 * @test IrbemTraceApi.MagEquatorMatchesTheOracle
 * @test IrbemTraceApi.MagEquatorIsTheMinimumAlongTheLine
 * @test IrbemTraceApi.MagEquatorRefinesAStartThatIsAlreadyTheMinimum
 * @test IrbemTraceApi.MagEquatorReportsAnOpenLine
 */
template <int NMAX>
[[nodiscard]] inline Result<MagneticEquator> find_magequator(const Igrf<NMAX>& model,
                                                             const Position<Frame::GEO>& start,
                                                             const PathTraceOptions& opt = {}) {
    MagneticEquator eq{};
    const double r0 = fixarray::norm(start.v);
    if (!(r0 > opt.r0) || !std::isfinite(r0)) return {Status::DomainError, eq};

    fixarray::vec3d b = model.evaluate(start).v;
    const double b_start = fixarray::norm(b);
    if (!(b_start > 0.0)) return {Status::DomainError, eq};

    const double ds_mag = detail::path_step(start, opt);

    // Probe both ways once. Downhill is the direction to walk; if neither is downhill the start is
    // already the minimum and the two probes ARE the bracket.
    fixarray::vec3d b_plus{};
    fixarray::vec3d b_minus{};
    (void)detail::rk4_step(model, start, b, ds_mag, b_plus);
    (void)detail::rk4_step(model, start, b, -ds_mag, b_minus);
    const double f_plus = fixarray::norm(b_plus);
    const double f_minus = fixarray::norm(b_minus);

    eq.b_min = b_start;
    eq.position = start;

    if (f_plus >= b_start && f_minus >= b_start) {
        double arc = 0.0;
        eq.b_min = detail::parabolic_minimum(f_minus, b_start, f_plus, ds_mag, arc);
        fixarray::vec3d b_v{};
        eq.position = detail::rk4_step(model, start, b, arc, b_v);
        return {Status::Ok, eq};
    }

    const double ds = f_plus < f_minus ? ds_mag : -ds_mag;
    Position<Frame::GEO> prev = start;
    double b_prev = b_start;
    Position<Frame::GEO> cur = start;
    fixarray::vec3d b_cur_v = b;
    double b_cur = b_start;

    for (int i = 0; i < opt.max_steps; ++i) {
        fixarray::vec3d b_new{};
        const Position<Frame::GEO> q = detail::rk4_step(model, cur, b_cur_v, ds, b_new);
        const double b_next = fixarray::norm(b_new);
        if (fixarray::norm(q.v) < opt.r0) {
            eq.b_min = b_cur;
            eq.position = cur;
            return {Status::OpenFieldLine, eq};
        }
        if (b_next >= b_cur) {
            double arc = 0.0;
            eq.b_min = detail::parabolic_minimum(b_prev, b_cur, b_next, ds, arc);
            fixarray::vec3d b_v{};
            eq.position = detail::rk4_step(model, cur, b_cur_v, arc, b_v);
            return {Status::Ok, eq};
        }
        prev = cur;
        b_prev = b_cur;
        cur = q;
        b_cur_v = b_new;
        b_cur = b_next;
    }
    (void)prev;
    eq.b_min = b_cur;
    eq.position = cur;
    return {Status::OpenFieldLine, eq};
}

/**
 * Where a particle of local pitch angle @p alpha_deg at @p start turns around — IRBEM's
 * `FIND_MIRROR_POINT`.
 *
 * The mirror field is not searched for, it is KNOWN: the first adiabatic invariant makes it
 * `B_m = B_local / sin²α` exactly. What the trace finds is where along the line that value is
 * reached, by walking in the direction of INCREASING `|B|` — which is what puts the mirror point in
 * the particle's own magnetic hemisphere, as the reference does.
 *
 * α = 90° is returned without tracing at all: the particle mirrors where it is, so the answer is
 * the input point verbatim and `B_m = B_local`. That is not an optimisation, it is the only
 * answer that is exactly right — a trace would return a point one refinement away from where the
 * particle demonstrably is.
 *
 * A line that reaches @ref PathTraceOptions::r0 before `B_m` is a particle in the LOSS CONE: it
 * hits the atmosphere instead of mirroring. That is physics, so it is @ref Status::OpenFieldLine
 * and not an error — the reference answers `baddata` here, which cannot distinguish it from a bad
 * input.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param start the starting position, GEO, Earth radii.
 * @param alpha_deg the LOCAL pitch angle at @p start, degrees, strictly inside (0, 180). Both
 *        endpoints are refused rather than clamped: a particle with no perpendicular velocity has
 *        no mirror point at all, and `sin(180 deg)` is 1.2e-16 rather than zero in binary64, so
 *        accepting it would return a mirror field 10^32 times `B_local` and an `OpenFieldLine`
 *        that looks like a physics result instead of a bad argument.
 * @param opt the tracing options.
 * @return the mirror point; @ref Status::OpenFieldLine when the particle is in the loss cone or the
 *         step cap was reached, @ref Status::DomainError for a start inside `r0`, a vanishing
 *         field, or a pitch angle outside the open interval (0, 180).
 * @complexity O(steps) IGRF evaluations plus @ref PathTraceOptions::refine_iterations more.
 * @alloc none.
 * @test IrbemTraceApi.MirrorPointMatchesTheOracle
 * @test IrbemTraceApi.MirrorPointAtNinetyDegreesIsTheInputPoint
 * @test IrbemTraceApi.MirrorPointInTheLossConeIsReported
 */
template <int NMAX>
[[nodiscard]] inline Result<MirrorPoint> find_mirror_point(const Igrf<NMAX>& model,
                                                           const Position<Frame::GEO>& start,
                                                           double alpha_deg,
                                                           const PathTraceOptions& opt = {}) {
    MirrorPoint mp{};
    const double r_start = fixarray::norm(start.v);
    if (!(r_start > opt.r0) || !std::isfinite(r_start)) return {Status::DomainError, mp};

    fixarray::vec3d b = model.evaluate(start).v;
    mp.b_local = fixarray::norm(b);
    mp.position = start;
    if (!(mp.b_local > 0.0)) return {Status::DomainError, mp};

    if (!(alpha_deg > 0.0) || !(alpha_deg < 180.0)) return {Status::DomainError, mp};
    const double sin_a = std::sin(alpha_deg * (std::numbers::pi / 180.0));
    if (!(sin_a > 0.0)) return {Status::DomainError, mp};
    mp.b_mirror = mp.b_local / (sin_a * sin_a);

    // A locally mirroring particle mirrors HERE. Nothing to trace, and no refinement can improve
    // on the point the caller handed in.
    if (mp.b_mirror <= mp.b_local) return {Status::Ok, mp};

    const double ds_mag = detail::path_step(start, opt);
    const double ds = detail::increasing_field_step(model, start, b, ds_mag);

    Position<Frame::GEO> cur = start;
    fixarray::vec3d b_cur = b;
    double f_lo = mp.b_local - mp.b_mirror;   // negative: below the mirror field
    for (int i = 0; i < opt.max_steps; ++i) {
        fixarray::vec3d b_new{};
        const Position<Frame::GEO> q = detail::rk4_step(model, cur, b_cur, ds, b_new);
        const double b_next = fixarray::norm(b_new);
        if (b_next >= mp.b_mirror) {
            const auto f = [&](const Position<Frame::GEO>&, const fixarray::vec3d& bv) {
                return fixarray::norm(bv) - mp.b_mirror;
            };
            const detail::RefinedPoint hit = detail::refine_in_step(
                model, cur, b_cur, ds, f_lo, b_next - mp.b_mirror, f, opt.refine_iterations);
            mp.position = hit.position;
            return {Status::Ok, mp};
        }
        if (fixarray::norm(q.v) < opt.r0) {
            mp.position = q;
            return {Status::OpenFieldLine, mp};   // the loss cone: the atmosphere came first
        }
        cur = q;
        b_cur = b_new;
        f_lo = b_next - mp.b_mirror;
    }
    mp.position = cur;
    return {Status::OpenFieldLine, mp};
}

/**
 * Where the field line through @p start crosses geodetic altitude @p stop_alt_km — IRBEM's
 * `FIND_FOOT_POINT`.
 *
 * The termination is on GEODETIC altitude, which is why this converts inside the step loop rather
 * than comparing radii and converting once at the end. The WGS-84 ellipsoid's polar semi-axis is
 * 21.4 km shorter than its equatorial one and field-line feet sit at |latitude| > 45°, so a
 * geocentric stand-in would land up to ~15 km from the requested surface — some 30× the ~0.5 km the
 * oracle's own iteration converges to, and in the one place where an altitude error maps directly
 * onto an atmospheric density.
 *
 * @p hemisphere follows IRBEM's `hemi_flag`. @ref FootTarget::Same walks in the direction of
 * increasing `|B|`, which is the input point's own side of the magnetic equator;
 * @ref FootTarget::Opposite walks the other way. @ref FootTarget::North and @ref FootTarget::South
 * name a side outright: the same-hemisphere foot is traced first and, if its geodetic latitude has
 * the wrong sign, the other direction is traced instead. Classifying by the FOOT's latitude rather
 * than the start's is what makes the answer right for a start point below the magnetic equator but
 * above the geographic one — the two differ by up to 11° and the dip equator wanders further still.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param start the starting position, GEO, Earth radii.
 * @param stop_alt_km the geodetic altitude to stop at, kilometres above the WGS-84 ellipsoid.
 * @param hemisphere which foot to return.
 * @param opt the tracing options. @ref PathTraceOptions::r0 must be below the surface the altitude
 *        names, or the trace terminates before it gets there.
 * @return the foot point; @ref Status::OpenFieldLine when the line left the domain or hit the step
 *         cap before reaching the altitude, @ref Status::DomainError for a start inside `r0`, a
 *         vanishing field, a non-finite altitude, or a start already BELOW @p stop_alt_km.
 * @complexity O(steps) IGRF evaluations plus one geodetic conversion per step — the conversion is
 *             four Bowring iterations, ~40 flops, against the step's ~2 000.
 * @alloc none.
 * @test IrbemTraceApi.FootPointMatchesTheOracle
 * @test IrbemTraceApi.FootPointFootTargetFlagsSelectTheTwoFeet
 * @test IrbemTraceApi.FootPointLandsOnTheRequestedGeodeticAltitude
 */
template <int NMAX>
[[nodiscard]] inline Result<FootPoint> find_foot_point(const Igrf<NMAX>& model,
                                                       const Position<Frame::GEO>& start,
                                                       double stop_alt_km, FootTarget hemisphere,
                                                       const PathTraceOptions& opt = {}) {
    FootPoint fp{};
    const double r_start = fixarray::norm(start.v);
    if (!(r_start > opt.r0) || !std::isfinite(r_start) || !std::isfinite(stop_alt_km)) {
        return {Status::DomainError, fp};
    }
    fixarray::vec3d b = model.evaluate(start).v;
    if (!(fixarray::norm(b) > 0.0)) return {Status::DomainError, fp};
    if (geo_to_gdz(start).radius() <= stop_alt_km) return {Status::DomainError, fp};

    const double ds_mag = detail::path_step(start, opt);
    const double same = detail::increasing_field_step(model, start, b, ds_mag);
    const double opposite = -same;

    // The one walk, parameterised by direction, so North/South can run it twice without a second
    // copy of the loop.
    const auto walk = [&](double ds) -> Result<FootPoint> {
        FootPoint out{};
        Position<Frame::GEO> cur = start;
        fixarray::vec3d b_cur = b;
        double f_lo = geo_to_gdz(cur).radius() - stop_alt_km;
        for (int i = 0; i < opt.max_steps; ++i) {
            fixarray::vec3d b_new{};
            const Position<Frame::GEO> q = detail::rk4_step(model, cur, b_cur, ds, b_new);
            const double alt = geo_to_gdz(q).radius();
            if (alt <= stop_alt_km) {
                const auto f = [&](const Position<Frame::GEO>& p, const fixarray::vec3d&) {
                    return geo_to_gdz(p).radius() - stop_alt_km;
                };
                const detail::RefinedPoint hit = detail::refine_in_step(
                    model, cur, b_cur, ds, f_lo, alt - stop_alt_km, f, opt.refine_iterations);
                out.position = geo_to_gdz(hit.position);
                out.field = FieldVector<Frame::GEO>{hit.field};
                out.b_magnitude = fixarray::norm(hit.field);
                return {Status::Ok, out};
            }
            if (fixarray::norm(q.v) < opt.r0) {
                out.position = geo_to_gdz(q);
                out.field = FieldVector<Frame::GEO>{b_new};
                out.b_magnitude = fixarray::norm(b_new);
                return {Status::OpenFieldLine, out};
            }
            cur = q;
            b_cur = b_new;
            f_lo = alt - stop_alt_km;
        }
        out.position = geo_to_gdz(cur);
        out.field = FieldVector<Frame::GEO>{b_cur};
        out.b_magnitude = fixarray::norm(b_cur);
        return {Status::OpenFieldLine, out};
    };

    if (hemisphere == FootTarget::Same) return walk(same);
    if (hemisphere == FootTarget::Opposite) return walk(opposite);

    const Result<FootPoint> first = walk(same);
    const double wanted = hemisphere == FootTarget::North ? 1.0 : -1.0;
    if (first.ok() && (first.value.position.latitude() * wanted) > 0.0) return first;
    if (!first.ok()) return first;
    return walk(opposite);
}

/**
 * Trace from @p start toward the Earth at a fixed step — IRBEM's
 * `TRACE_FIELD_LINE_TOWARD_EARTH`.
 *
 * The half-line the input point sits on, sampled uniformly for a plot. The direction is the
 * input point's own magnetic hemisphere — the direction of increasing `|B|` — which is the
 * reference's behaviour and the only one that makes "toward the Earth" well defined for a point
 * off the magnetic equator. Sample 0 is the input point verbatim; the last sample is the first one
 * inside @ref PathTraceOptions::r0, so the path visibly crosses the surface rather than stopping
 * short of it, which is what the reference does and what a plot wants.
 *
 * The step is @ref PathTraceOptions::step_size when positive and `L/steps_per_l` otherwise — the
 * reference takes `ds` as a required argument here and nowhere else, and this is why.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param start the starting position, GEO, Earth radii.
 * @param path receives the samples; the caller sizes it, and the trace stops when it is full.
 * @param opt the tracing options.
 * @return how many samples were written; @ref Status::NotConverged when @p path filled before the
 *         surface was reached, @ref Status::OpenFieldLine when the step cap was, and
 *         @ref Status::DomainError for an empty span, a start inside `r0`, or a vanishing field.
 * @complexity O(samples) IGRF evaluations, ~4 per sample.
 * @alloc none — the caller's span is the only storage.
 * @test IrbemTraceApi.TowardEarthMatchesTheOracle
 * @test IrbemTraceApi.TowardEarthSamplesAreOneStepApart
 * @test IrbemTraceApi.TowardEarthReportsATruncatedPath
 */
template <int NMAX>
[[nodiscard]] inline Result<std::size_t> trace_field_line_toward_earth(
    const Igrf<NMAX>& model, const Position<Frame::GEO>& start, std::span<PathPoint> path,
    const PathTraceOptions& opt = {}) {
    if (path.empty()) return {Status::DomainError, 0};
    const double r_start = fixarray::norm(start.v);
    if (!(r_start > opt.r0) || !std::isfinite(r_start)) return {Status::DomainError, 0};

    fixarray::vec3d b = model.evaluate(start).v;
    const double b_start = fixarray::norm(b);
    if (!(b_start > 0.0)) return {Status::DomainError, 0};

    const double ds = detail::increasing_field_step(model, start, b, detail::path_step(start, opt));

    path[0] = PathPoint{start, b_start};
    std::size_t n = 1;
    Position<Frame::GEO> cur = start;
    fixarray::vec3d b_cur = b;
    for (int i = 0; i < opt.max_steps; ++i) {
        if (n == path.size()) return {Status::NotConverged, n};
        fixarray::vec3d b_new{};
        const Position<Frame::GEO> q = detail::rk4_step(model, cur, b_cur, ds, b_new);
        path[n] = PathPoint{q, fixarray::norm(b_new)};
        ++n;
        if (fixarray::norm(q.v) < opt.r0) return {Status::Ok, n};
        cur = q;
        b_cur = b_new;
    }
    return {Status::OpenFieldLine, n};
}

namespace detail {

/// Trace one half-line from @p start to the reference surface, writing samples forward from
/// `path[0]` and capping the last one exactly ON the surface.
///
/// The end cap is the reason this is not simply the loop above: `trace_field_line` reports a line
/// between two `R0` crossings, so an endpoint one whole step past the surface would make the path's
/// extent an artefact of the step size. Regula falsi on `|r| − R0` inside the last step fixes it
/// for four IGRF evaluations per iteration.
///
/// @tparam NMAX the IGRF truncation degree.
/// @param model the internal field model. @param start the origin, GEO, Earth radii, NOT written.
/// @param b_start the field at @p start, nT. @param ds the signed step.
/// @param path receives the samples, first one step away from @p start.
/// @param opt the tracing options.
/// @param status receives @ref Status::Ok when the surface was reached, @ref Status::NotConverged
///        when @p path filled first, @ref Status::OpenFieldLine when the step cap did.
/// @return how many samples were written.
/// @complexity O(samples) IGRF evaluations. @alloc none.
/// @test IrbemTraceApi.TraceFieldLineEndsOnTheReferenceSurfaceAtBothFeet
template <int NMAX>
[[nodiscard]] inline std::size_t trace_half_line(const Igrf<NMAX>& model,
                                                 const Position<Frame::GEO>& start,
                                                 const fixarray::vec3d& b_start, double ds,
                                                 std::span<PathPoint> path,
                                                 const PathTraceOptions& opt, Status& status) {
    Position<Frame::GEO> cur = start;
    fixarray::vec3d b_cur = b_start;
    double f_lo = fixarray::norm(start.v) - opt.r0;
    std::size_t n = 0;
    for (int i = 0; i < opt.max_steps; ++i) {
        if (n == path.size()) {
            status = Status::NotConverged;
            return n;
        }
        fixarray::vec3d b_new{};
        const Position<Frame::GEO> q = rk4_step(model, cur, b_cur, ds, b_new);
        const double r = fixarray::norm(q.v);
        if (r < opt.r0) {
            const auto f = [&](const Position<Frame::GEO>& p, const fixarray::vec3d&) {
                return fixarray::norm(p.v) - opt.r0;
            };
            const RefinedPoint hit =
                refine_in_step(model, cur, b_cur, ds, f_lo, r - opt.r0, f, opt.refine_iterations);
            path[n] = PathPoint{hit.position, fixarray::norm(hit.field)};
            ++n;
            status = Status::Ok;
            return n;
        }
        path[n] = PathPoint{q, fixarray::norm(b_new)};
        ++n;
        cur = q;
        b_cur = b_new;
        f_lo = r - opt.r0;
    }
    status = Status::OpenFieldLine;
    return n;
}

}  // namespace detail

/**
 * Trace the WHOLE field line through @p start, foot to foot — IRBEM's `TRACE_FIELD_LINE`.
 *
 * Samples run **along `+B̂`**: index 0 is the foot the field points away from (the southern one,
 * for the real geomagnetic field), the last index is the foot it points into, and the input point
 * sits verbatim at @ref TracedLine::start_index. The reference orders its `posit` the other way; the
 * choice is arbitrary, so it is stated rather than inherited, and a caller who cares should read
 * `start_index` instead of assuming an end.
 *
 * Both ends are capped ON the reference surface by regula falsi inside the last step, so the path's
 * extent is a property of the field line and not of the step size.
 *
 * `Lm`, `Bmin` and `XJ` come from the samples already in the caller's buffer rather than from a
 * second trace — one pass over the line, not two. `XJ` is `I` for a particle mirroring at the input
 * point, integrated by the same right-endpoint rule `lstar.hpp` uses outward from the equator in
 * both directions, so this and @ref trace_invariant are the same quadrature on the same grid rather
 * than two approximations that happen to be close. (`IrbemTraceApi.XjAgreesWithTheInvariantTracer`
 * holds them to 1 × 10⁻⁴ relative; what separates them is only that RK4 is not exactly reversible,
 * so the two grids drift apart at O(ds⁵) per step.)
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param start the starting position, GEO, Earth radii.
 * @param path receives the samples; the caller sizes it. IRBEM's own cap is
 *        @ref irbem_max_path_points, which at the default step is ~25× more room than an L=8 line
 *        needs.
 * @param opt the tracing options.
 * @return the line's scalars and the sample count. @ref Status::NotConverged when @p path filled
 *         before the line closed — @ref TracedLine::truncated is then set, and what the scalars
 *         mean depends on WHICH half filled it, which that field's brief spells out: they are all
 *         zero when the backward half alone exhausted the buffer. @ref Status::OpenFieldLine when
 *         a half hit the step cap, and @ref Status::DomainError for an empty span, a start inside
 *         `r0`, or a vanishing field.
 * @complexity O(samples) IGRF evaluations, ~4 per sample, plus O(samples) for the quadrature.
 * @alloc none — the caller's span is the only storage, and the in-place reverse that puts the two
 *        halves in order needs none.
 * @test IrbemTraceApi.TraceFieldLineMatchesTheOracle
 * @test IrbemTraceApi.TraceFieldLineEndsOnTheReferenceSurfaceAtBothFeet
 * @test IrbemTraceApi.XjAgreesWithTheInvariantTracer
 * @test IrbemTraceApi.TraceFieldLineReportsATruncatedPath
 * @test IrbemTraceApi.TraceFieldLineTruncatesInTheForwardHalfToo
 */
template <int NMAX>
[[nodiscard]] inline Result<TracedLine> trace_field_line(const Igrf<NMAX>& model,
                                                         const Position<Frame::GEO>& start,
                                                         std::span<PathPoint> path,
                                                         const PathTraceOptions& opt = {}) {
    TracedLine line{};
    if (path.empty()) return {Status::DomainError, line};
    const double r_start = fixarray::norm(start.v);
    if (!(r_start > opt.r0) || !std::isfinite(r_start)) return {Status::DomainError, line};

    const fixarray::vec3d b = model.evaluate(start).v;
    line.b_local = fixarray::norm(b);
    if (!(line.b_local > 0.0)) return {Status::DomainError, line};

    const double ds = detail::path_step(start, opt);

    // The −B̂ half first, written forward and then reversed in place: that is what puts the whole
    // path in one order with no scratch buffer and no second pass over the field.
    Status back_status = Status::Ok;
    const std::size_t back = detail::trace_half_line(model, start, b, -ds, path, opt, back_status);
    std::reverse(path.begin(), path.begin() + static_cast<std::ptrdiff_t>(back));

    if (back == path.size()) {
        line.point_count = back;
        line.truncated = true;
        return {Status::NotConverged, line};
    }
    path[back] = PathPoint{start, line.b_local};
    line.start_index = back;

    Status fwd_status = Status::Ok;
    const std::size_t fwd = detail::trace_half_line(model, start, b, ds, path.subspan(back + 1),
                                                    opt, fwd_status);
    line.point_count = back + 1 + fwd;
    line.truncated = back_status == Status::NotConverged || fwd_status == Status::NotConverged;

    // --- Bmin, parabola-refined off the sample grid ------------------------------------------
    std::size_t j = 0;
    for (std::size_t i = 1; i < line.point_count; ++i) {
        if (path[i].b_magnitude < path[j].b_magnitude) j = i;
    }
    line.equator_index = j;
    line.b_min = path[j].b_magnitude;
    line.equator = path[j].position;
    // The end caps are PARTIAL steps, so a bracket that touches one is not equally spaced and the
    // three-point formula does not apply to it. Interior brackets are — which is every real case,
    // since a minimum one step from a foot is not a magnetic equator.
    if (j > 1 && j + 2 < line.point_count) {
        double arc = 0.0;
        line.b_min = detail::parabolic_minimum(path[j - 1].b_magnitude, path[j].b_magnitude,
                                               path[j + 1].b_magnitude, ds, arc);
        if (arc != 0.0) {
            fixarray::vec3d b_v{};
            const fixarray::vec3d b_j = model.evaluate(path[j].position).v;
            line.equator = detail::rk4_step(model, path[j].position, b_j, arc, b_v);
        }
    }

    // --- XJ: I for a particle mirroring at the input point ------------------------------------
    // The same right-endpoint rule lstar.hpp integrates with, outward from the equator in both
    // directions, so the two are one quadrature and not two.
    const double b_mirror = line.b_local;
    double sum = 0.0;
    bool closed = true;
    for (const int dir : {1, -1}) {
        std::size_t i = j;
        bool reached = false;
        while (true) {
            if (dir > 0) {
                if (i + 1 >= line.point_count) break;
                ++i;
            } else {
                if (i == 0) break;
                --i;
            }
            const double b_here = path[i].b_magnitude;
            if (b_here >= b_mirror) {
                const std::size_t prev = dir > 0 ? i - 1 : i + 1;
                const double b_prev = path[prev].b_magnitude;
                if (b_here > b_prev) {
                    const double frac = (b_mirror - b_prev) / (b_here - b_prev);
                    sum += 0.5 * std::sqrt(std::max(0.0, 1.0 - (b_prev / b_mirror))) * frac;
                }
                reached = true;
                break;
            }
            sum += std::sqrt(1.0 - (b_here / b_mirror));
        }
        closed = closed && reached;
    }
    line.invariant_i = sum * ds;

    const Result<double> lm = mcilwain_l(line.invariant_i, b_mirror, dipole_moment(model));
    line.mcilwain_l = lm.value;

    if (line.truncated) return {Status::NotConverged, line};
    if (back_status != Status::Ok || fwd_status != Status::Ok || !closed) {
        return {Status::OpenFieldLine, line};
    }
    return {lm.status, line};
}

// =============================================================================================
// The device lane — a SEPARATE kernel, because these routines are on the other side of the
// roofline
// =============================================================================================

#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
namespace detail {

/// Stage one batch of fixed-step earthward traces on the device.
///
/// Everything epoch-dependent is computed on the HOST, once: the Gauss coefficients are
/// interpolated to the epoch here and the Legendre normalisation is `constexpr` here, so the
/// device never sees IGRF's 26-epoch table. That is the same split `lstar.hpp` makes and for the
/// same reason — interpolating per thread would be N redundant copies of a calculation done once.
///
/// @tparam NMAX the IGRF truncation degree.
/// @param model the internal field model. @param starts the starting positions, GEO, Earth radii.
/// @param paths receives `starts.size() × max_points` samples, line-major.
/// @param counts receives the sample count per line.
/// @param statuses receives the per-line status.
/// @param max_points the stride of @p paths, and the cap on samples per line.
/// @param opt the tracing options; @ref PathTraceOptions::step_size must be positive.
/// @return `Status::ParametersMissing` when the device could not be used after all, which the
///         caller reads as "fall through to the host" rather than as an error; otherwise
///         `Status::Ok` when every line reached the surface.
/// @complexity One dispatch; O(lines × steps) field evaluations, concurrent, and
///             `16 × max_points` bytes read back per line.
/// @alloc host staging vectors for the coefficients, positions and results, plus seven pooled
///        device buffers returned on scope exit. Per BATCH, never per line.
/// @test IrbemTraceApiGpu.PathKernelAgreesWithTheHostLane
/// @test IrbemTraceApiGpu.PathBatchHonoursTheStepCapOnEitherLane
template <int NMAX>
[[nodiscard]] inline Result<bool> trace_path_on_device(
    const Igrf<NMAX>& model, std::span<const Position<Frame::GEO>> starts,
    std::span<PathPoint> paths, std::span<std::uint32_t> counts, std::span<Status> statuses,
    std::size_t max_points, const PathTraceOptions& opt) {

    if (!gpu::available()) return {Status::ParametersMissing, false};
    if (!std::filesystem::exists(gpu::shader_path("irbem_trace_path_f32"))) {
        return {Status::ParametersMissing, false};
    }

    constexpr int kSlots = ((NMAX + 1) * (NMAX + 2)) / 2;
    const std::size_t n = starts.size();

    std::vector<float> coef(2 * kSlots, 0.0F);
    for (int deg = 1; deg <= NMAX; ++deg) {
        for (int m = 0; m <= deg; ++m) {
            const std::size_t k = ((static_cast<std::size_t>(deg) * (deg + 1)) / 2) + m;
            coef[k] = static_cast<float>(model.g(deg, m));
            coef[kSlots + k] = static_cast<float>(model.h(deg, m));
        }
    }
    constexpr auto kNorm =
        ::cheatah::space::irbem::detail::make_legendre_normalisation<NMAX, double>();
    std::vector<float> nrm((2 * kSlots) + NMAX + 1, 0.0F);
    for (int k = 0; k < kSlots; ++k) {
        nrm[static_cast<std::size_t>(k)] = static_cast<float>(kNorm.e[static_cast<std::size_t>(k)]);
        nrm[static_cast<std::size_t>(kSlots + k)] =
            static_cast<float>(kNorm.f[static_cast<std::size_t>(k)]);
    }
    for (int deg = 0; deg <= NMAX; ++deg) {
        nrm[static_cast<std::size_t>((2 * kSlots) + deg)] =
            static_cast<float>(kNorm.diagonal[static_cast<std::size_t>(deg)]);
    }

    std::vector<float> pos(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
        pos[(3 * i) + 0] = static_cast<float>(starts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(starts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(starts[i].v[2]);
    }
    // ds and r0 ride in the uint dims buffer scaled by 10^6 — the ABI has no float dims slot.
    // For `r0` that is a one-off 6 mm; for the STEP it is not, because the step is applied a few
    // hundred times and the quantum accumulates. Measured over 1 024 lines, device against host:
    // an exactly representable ds = 0.02 R_E ends 6.5e-6 R_E apart at the last sample, while
    // ds = 0.0234567891 (which 1e-6 cannot hold) ends 8.5e-5 R_E — thirteen times further, and
    // that difference is the quantisation rather than fp32. A caller who wants the two lanes to
    // agree to fp32 alone should name a step that is a whole number of microradii.
    //
    // `max_steps` rides here too, and must: the host loop is bounded by it, so a kernel that
    // ignored it would answer a capped trace differently depending on whether a device was
    // present. Measured before it was carried: at max_steps = 10 all 1 024 lines of a batch
    // disagreed — device 80 samples and `Ok`, host 11 and `OpenFieldLine`.
    const std::array<std::uint32_t, 6> dims{
        static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(NMAX),
        static_cast<std::uint32_t>(max_points),
        static_cast<std::uint32_t>(std::llround(opt.step_size * 1.0e6)),
        static_cast<std::uint32_t>(std::llround(opt.r0 * 1.0e6)),
        static_cast<std::uint32_t>(std::max(0, opt.max_steps))};

    std::vector<float> path_out(3 * max_points * n);
    std::vector<float> bmag_out(max_points * n);
    std::vector<std::uint32_t> report(2 * n);

    namespace gl = gpu::detail::gl;
    gl::detail::Context& c = gl::detail::ctx();
    gpu::detail::Leases lease;
    gl::detail::Buffer* b_pos = lease.add(c.new_data_buffer(pos.size() * sizeof(float)));
    gl::detail::Buffer* b_cf = lease.add(c.new_buffer(coef.size() * sizeof(float)));
    gl::detail::Buffer* b_nr = lease.add(c.new_buffer(nrm.size() * sizeof(float)));
    gl::detail::Buffer* b_path = lease.add(c.new_data_buffer(path_out.size() * sizeof(float)));
    gl::detail::Buffer* b_bmag = lease.add(c.new_data_buffer(bmag_out.size() * sizeof(float)));
    gl::detail::Buffer* b_rep =
        lease.add(c.new_data_buffer(report.size() * sizeof(std::uint32_t)));
    gl::detail::Buffer* b_dm = lease.add(c.new_buffer(dims.size() * sizeof(std::uint32_t)));
    c.upload(b_pos, pos.data(), pos.size() * sizeof(float));
    c.upload(b_cf, coef.data(), coef.size() * sizeof(float));
    c.upload(b_nr, nrm.data(), nrm.size() * sizeof(float));
    c.upload(b_dm, dims.data(), dims.size() * sizeof(std::uint32_t));
    {
        c.dispatch_1d(gpu::qualified("irbem_trace_path_f32").c_str(), lease.data(), 7, n);
    }
    c.download(b_path, path_out.data(), path_out.size() * sizeof(float));
    c.download(b_bmag, bmag_out.data(), bmag_out.size() * sizeof(float));
    c.download(b_rep, report.data(), report.size() * sizeof(std::uint32_t));

    bool all_ok = true;
    for (std::size_t i = 0; i < n; ++i) {
        counts[i] = report[2 * i];
        const std::uint32_t code = report[(2 * i) + 1];
        statuses[i] = code < status_count ? static_cast<Status>(code) : Status::DomainError;
        all_ok = all_ok && (statuses[i] == Status::Ok);
        for (std::uint32_t k = 0; k < counts[i]; ++k) {
            const std::size_t src = (i * max_points) + k;
            paths[src] = PathPoint{Position<Frame::GEO>{fixarray::vec3d{path_out[3 * src],
                                                                       path_out[(3 * src) + 1],
                                                                       path_out[(3 * src) + 2]}},
                                   bmag_out[src]};
        }
    }
    return {all_ok ? Status::Ok : Status::OpenFieldLine, true};
}

}  // namespace detail
#endif  // CHEATAH_SPACE_IRBEM_LSTAR_GPU

/**
 * Trace a whole batch of field lines earthward and return every path.
 *
 * The batch form of @ref trace_field_line_toward_earth, and the only shape a device can
 * accelerate: one trace is a serial RK4 chain, so no hardware makes it faster, and the parallelism
 * is entirely ACROSS lines.
 *
 * **This kernel's crossover is 256 lines — MEASURED, and lower than the invariant tracer's 512
 * rather than higher, which is not what the output shape predicts.** Both run the same RK4 chain
 * over the same field. `irbem_trace_i_f32` returns four floats per line — ~9 400 flops per byte
 * moved. This one returns four floats per STEP, so at ~250 steps a line it moves ~250× as many
 * bytes for the same arithmetic: ~125 flops/byte. The transfer term the invariant tracer never
 * pays is real, but it does not move the crossover, because a crossover is a RATIO and this
 * routine's host lane is ~2.6× dearer per line too. Where the transfer shows up instead is the
 * ceiling — see the two bullets under the table. That is why it is a separate kernel with a
 * separate registry row and a separate measured threshold, and why merging the two would be a
 * mistake in both directions.
 *
 * Measured on an RTX 3070 Ti against this header's own fp64 host lane (`-O3 -march=native
 * -ffp-contract=off`), fixed step `ds = 0.02 R_E`, 512 samples of headroom per line, starts spread
 * over L = 2…8 (`IrbemTraceApiGpu.PathKernelCrossover`):
 *
 * | lines | 64 | 128 | 256 | 512 | 1024 | 4096 | 16384 | 65536 |
 * |---|---|---|---|---|---|---|---|---|
 * | speedup | 0.33–0.39× | 0.66–0.69× | **1.05–1.15×** | 2.05–2.15× | 3.69–3.94× | 9.9–13.3× | 9.6–23.5× | 19–26.3× |
 *
 * against a host lane flat at 154–165 µs/line. Every cell is a range over FOUR full runs, not a
 * best-of: below 1 024 lines the spread is a few per cent and the crossover reproduces at 256
 * every time, while above 4 096 the same size has measured 9.6× and 23.5× on different runs. Two
 * things in that row are worth stating plainly because they contradict the obvious expectation:
 *
 *  - **The crossover came out LOWER than the invariant tracer's, not higher.** Not because this
 *    kernel is cheaper — because the HOST lane is dearer. A fixed `ds = 0.02 R_E` line is ~250
 *    steps against the invariant tracer's ~120 at `ds = L/50`, so the host costs 154 µs/line
 *    against 60, and the same submit floor is paid off in half the batch. A crossover is a ratio.
 *  - **Where the bandwidth shows is the CEILING.** The invariant tracer's speedup is still
 *    climbing at 65 536 lines (48.9×); this one flattens in the low-to-mid twenties and goes
 *    run-to-run noisy above 4 096, where a batch stages 100 MB–1 GB of results through host
 *    memory. The ranges above are quoted as ranges on purpose: at 16 384 lines the spread across
 *    four runs is a factor of 2.4, so any single number from that end of the curve is a
 *    performance claim the next run will not support.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param starts the starting positions, GEO, Earth radii.
 * @param paths receives `starts.size() × max_points` samples, line-major: line `i`'s samples are
 *        `paths[i*max_points .. i*max_points + counts[i])`. The caller sizes it, which at 3 000
 *        samples and 4 096 lines is 393 MB — a number that belongs to the caller, not here.
 * @param counts receives each line's sample count; same length as @p starts.
 * @param statuses receives each line's status; same length as @p starts.
 * @param opt the tracing options. @ref PathTraceOptions::step_size MUST be positive: the device
 *        lane takes a fixed step, and an L-proportional one would make the step a per-line value
 *        the dims buffer cannot carry.
 * @return `true` when the device lane serviced the call — asserted by a test rather than trusted,
 *         because a silent fallback is what makes a performance claim worthless. The status is
 *         @ref Status::Ok when every line reached the surface, @ref Status::DomainError on a
 *         length mismatch or a non-positive step.
 * @complexity O(lines × steps) field evaluations; on the device those run concurrently.
 * @alloc the device lane stages coefficients, positions and results per BATCH; the host lane
 *        allocates nothing.
 * @test IrbemTraceApiGpu.PathKernelAgreesWithTheHostLane
 * @test IrbemTraceApiGpu.PathBatchUsesTheDeviceWhenOneIsAvailable
 * @test IrbemTraceApiGpu.PathBatchHonoursTheStepCapOnEitherLane
 */
template <int NMAX>
[[nodiscard]] inline Result<bool> trace_field_line_toward_earth_batch(
    const Igrf<NMAX>& model, std::span<const Position<Frame::GEO>> starts,
    std::span<PathPoint> paths, std::span<std::uint32_t> counts, std::span<Status> statuses,
    const PathTraceOptions& opt) {

    const std::size_t n = starts.size();
    if (counts.size() != n || statuses.size() != n) return {Status::DomainError, false};
    if (!(opt.step_size > 0.0)) return {Status::DomainError, false};
    if (n == 0) return {Status::Ok, false};
    if (paths.size() % n != 0) return {Status::DomainError, false};
    const std::size_t max_points = paths.size() / n;
    if (max_points == 0) return {Status::DomainError, false};

#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
    if (gpu::prefer_gpu("irbem_trace_path_f32", n)) {
        if (const Result<bool> r = detail::trace_path_on_device(model, starts, paths, counts,
                                                               statuses, max_points, opt);
            r.status != Status::ParametersMissing) {
            return r;
        }
        // The device could not be used after all (no SPIR-V, no device). Fall through to the host:
        // a missing shader is a deployment problem, not a reason to refuse to compute.
    }
#endif

    bool all_ok = true;
    for (std::size_t i = 0; i < n; ++i) {
        const Result<std::size_t> r = trace_field_line_toward_earth(
            model, starts[i], paths.subspan(i * max_points, max_points), opt);
        counts[i] = static_cast<std::uint32_t>(r.value);
        statuses[i] = r.status;
        all_ok = all_ok && (r.status == Status::Ok);
    }
    return {all_ok ? Status::Ok : Status::OpenFieldLine, false};
}

}  // namespace cheatah::space::irbem
