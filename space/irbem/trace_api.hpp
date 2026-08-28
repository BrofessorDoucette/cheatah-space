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
 * unchanged, so the arithmetic intensity falls by two to three orders of magnitude and the
 * routines land on the *other* side of the roofline: they are **bandwidth-bound**, and the
 * measured crossover reflects it. See @ref trace_field_line_batch for the number.
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
 *    position at O(ds) — measured against the oracle on an L=4 line, the coarse sample sits
 *    5.5 × 10⁻⁶ relative above `Bmin` and 4.4 × 10⁻³ R_E away from it, which is 100× the
 *    difference the fit leaves.
 *  - The mirror point, the foot point and the `R0` end caps come from **regula falsi inside the
 *    last step** — on `|B|`, on geodetic altitude and on radius respectively — evaluating the true
 *    field at each trial rather than interpolating between the bracketing samples.
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
 * Which magnetic hemisphere @ref find_foot_point should walk to.
 *
 * The enumerator values ARE IRBEM's `hemi_flag` codes, so a caller porting from the Fortran can
 * pass its integer through `static_cast` and a reader of either can check the other.
 */
enum class Hemisphere : int {
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
    /// Whether the path ran out of caller-supplied room before the line closed. When true the
    /// scalars describe the truncated path and not the field line, which is why the status is
    /// @ref Status::NotConverged rather than @ref Status::Ok.
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

/// The step length a trace from @p start uses, in Earth radii.
///
/// @param start the starting position, GEO, Earth radii.
/// @param opt the options; @ref PathTraceOptions::step_size wins when positive.
/// @return the unsigned step length; never zero, so a trace cannot stall.
/// @complexity O(1). @alloc none.
/// @test IrbemTraceApi.StepSizeFollowsTheDipoleLUnlessOverridden
[[nodiscard]] inline double path_step(const Position<Frame::GEO>& start,
                                      const PathTraceOptions& opt) {
    if (opt.step_size > 0.0) return opt.step_size;
    const double per_l = opt.steps_per_l > 0.0 ? opt.steps_per_l : 50.0;
    return dipole_l(start) / per_l;
}

/// The signed step that walks toward INCREASING `|B|` from @p p — the direction of the near mirror
/// point, the near foot, and "the same magnetic hemisphere".
///
/// One probe step rather than a gradient: the walk is along the field line, so what matters is
/// which way the LINE's field rises, and a `∇|B|` at the point answers a different question wherever
/// the line is curved. At a point exactly at the minimum both directions rise and `+ds` is
/// returned, which is the reference's choice too.
///
/// @tparam NMAX the IGRF truncation degree.
/// @param model the internal field model. @param p the position, GEO, Earth radii.
/// @param b the field at @p p, nT, already known. @param ds_mag the unsigned step length.
/// @return `+ds_mag` or `-ds_mag`.
/// @complexity Four IGRF evaluations — one RK4 step. @alloc none.
/// @test IrbemTraceApi.IncreasingFieldStepPointsAwayFromTheEquator
template <int NMAX>
[[nodiscard]] inline double increasing_field_step(const Igrf<NMAX>& model,
                                                  const Position<Frame::GEO>& p,
                                                  const fixarray::vec3d& b, double ds_mag) {
    fixarray::vec3d probe{};
    (void)rk4_step(model, p, b, ds_mag, probe);
    return fixarray::norm(probe) >= fixarray::norm(b) ? ds_mag : -ds_mag;
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
/// mirror point, so plain regula falsi retains one endpoint for every iteration and converges
/// linearly with a ratio near 1. Halving the stale endpoint's value restores superlinearity, which
/// is what lets eight iterations be enough.
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
    RefinedPoint best{p, b, 0.0};
    for (int i = 0; i < iterations; ++i) {
        const double denom = f_hi - f_lo;
        if (!(std::abs(denom) > 0.0)) break;
        const double t = lo - (f_lo * (hi - lo) / denom);
        fixarray::vec3d b_t{};
        const Position<Frame::GEO> q = rk4_step(model, p, b, t, b_t);
        const double f_t = f(q, b_t);
        best = RefinedPoint{q, b_t, t};
        if (!std::isfinite(f_t)) break;
        if ((f_t < 0.0) == (f_lo < 0.0)) {
            lo = t;
            f_lo = f_t;
            f_hi *= 0.5;   // Illinois: the retained endpoint's weight decays instead of stalling
        } else {
            hi = t;
            f_hi = f_t;
            f_lo *= 0.5;
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
 * relative above the oracle's `Bmin` and 4.4 × 10⁻³ R_E from its position, while the fit lands
 * within 2 × 10⁻⁷ and 3 × 10⁻⁴.
 *
 * A start point that is ALREADY at the minimum is not a special case in the physics and is not one
 * here: both neighbours are probed, and if both are higher the bracket is centred on the start.
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
 * @test IrbemTraceApi.MagEquatorOfADipoleIsTheGeographicEquator
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
 * @param alpha_deg the LOCAL pitch angle at @p start, degrees, in (0, 180].
 * @param opt the tracing options.
 * @return the mirror point; @ref Status::OpenFieldLine when the particle is in the loss cone or the
 *         step cap was reached, @ref Status::DomainError for a start inside `r0`, a vanishing
 *         field, or a pitch angle outside (0, 180].
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

    if (!(alpha_deg > 0.0) || !(alpha_deg <= 180.0)) return {Status::DomainError, mp};
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
 * @p hemisphere follows IRBEM's `hemi_flag`. @ref Hemisphere::Same walks in the direction of
 * increasing `|B|`, which is the input point's own side of the magnetic equator;
 * @ref Hemisphere::Opposite walks the other way. @ref Hemisphere::North and @ref Hemisphere::South
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
 * @test IrbemTraceApi.FootPointHemisphereFlagsSelectTheTwoFeet
 * @test IrbemTraceApi.FootPointLandsOnTheRequestedGeodeticAltitude
 */
template <int NMAX>
[[nodiscard]] inline Result<FootPoint> find_foot_point(const Igrf<NMAX>& model,
                                                       const Position<Frame::GEO>& start,
                                                       double stop_alt_km, Hemisphere hemisphere,
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

    if (hemisphere == Hemisphere::Same) return walk(same);
    if (hemisphere == Hemisphere::Opposite) return walk(opposite);

    const Result<FootPoint> first = walk(same);
    const double wanted = hemisphere == Hemisphere::North ? 1.0 : -1.0;
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
/// @test IrbemTraceApi.HalfLineEndsOnTheReferenceSurface
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
 *         before the line closed (@ref TracedLine::truncated is then set and the scalars describe
 *         the partial path), @ref Status::OpenFieldLine when a half hit the step cap, and
 *         @ref Status::DomainError for an empty span, a start inside `r0`, or a vanishing field.
 * @complexity O(samples) IGRF evaluations, ~4 per sample, plus O(samples) for the quadrature.
 * @alloc none — the caller's span is the only storage, and the in-place reverse that puts the two
 *        halves in order needs none.
 * @test IrbemTraceApi.TraceFieldLineMatchesTheOracle
 * @test IrbemTraceApi.TraceFieldLineEndsOnTheReferenceSurfaceAtBothFeet
 * @test IrbemTraceApi.XjAgreesWithTheInvariantTracer
 * @test IrbemTraceApi.TraceFieldLineReportsATruncatedPath
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

}  // namespace cheatah::space::irbem
