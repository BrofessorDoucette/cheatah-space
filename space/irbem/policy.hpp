#pragma once

/**
 * @file policy.hpp
 * @brief space.irbem — the two compile-time axes every downstream routine is templated on:
 *        what arithmetic it computes in, and whose algorithm it computes.
 *
 * This module has to be two libraries at once. It has to be a *differential* reimplementation —
 * able to reproduce IRBEM's answer closely enough that a disagreement is a bug rather than a
 * choice — and it has to be the faster, more accurate library that was the point of rewriting it.
 * Those two goals conflict at every design decision: IRBEM's drift-shell walk is serial across
 * azimuths, its mirror-point search lands on a fixed θ grid, and its I integral is a first-order
 * rectangle sum. Reproducing that is how you check yourself; shipping it is how you stay slow.
 *
 * Both goals are met by making the difference a compile-time *policy* rather than a runtime flag.
 * A `bool improved` argument threaded through a tracer would cost a branch in the innermost loop
 * (the one that runs ~10⁵ times per L\* point), would defeat inlining across it, and would leave
 * both code paths in every binary. A policy type costs nothing: `if constexpr (C::carry_k1)`
 * compiles to one or the other, and the lane is visible in the type of the result rather than
 * buried in a caller's argument.
 *
 * ### Precision — the load-bearing invariant is `accum = double`
 *
 * `docs/ERROR_BUDGET.md` §3 measures the discretization floor at IRBEM's recommended resolution
 * as 1.2 × 10⁻³ … 4.0 × 10⁻³ relative on L\*. Against that floor:
 *
 * - one fp32 operation costs ~6 × 10⁻⁸ (the unit roundoff of binary32) — four to five orders
 *   below the floor, so the *integrand* may be fp32, which is what makes the GPU lane possible at
 *   all (SPIR-V's `GLSL.std.450` transcendentals exist only for 16- and 32-bit floats, and every
 *   Tsyganenko model calls `exp`);
 * - naively summing the ~10³ terms of a field-line quadrature in fp32 costs ~10⁻⁴ — the *same
 *   order* as the discretization floor at the tighter settings, i.e. the reduction alone would
 *   consume the whole budget it was supposed to fit inside.
 *
 * Roundoff in a sum grows with the number of terms; roundoff in an evaluation does not. So the
 * integrand type is a free choice and the accumulator type is not. Every precision policy in this
 * header sets `accum = double`, `SoundPrecision` refuses one that does not, and `Policy` is
 * constrained on it — an fp32 accumulator cannot be introduced by accident, only by editing this
 * paragraph. (ERROR_BUDGET §6: reductions are additionally *ordered* — a fixed sequence, never a
 * tree reduction whose order varies with workgroup size. That is a property of the reduction code,
 * not of a type, so it is asserted there rather than encoded here.)
 *
 * @note Nothing in this header allocates, branches at runtime, or has a representation. The policy
 *       structs are empty; they exist entirely to be looked at by the compiler.
 */

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <string_view>

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// Roundoff arithmetic — where the tolerances in docs/ERROR_BUDGET.md come from
// -------------------------------------------------------------------------------------------

/**
 * The unit roundoff of @p T: the largest relative error a single correctly-rounded operation can
 * introduce, `u = ε/2` for round-to-nearest.
 *
 * This is the `6 × 10⁻⁸` row of ERROR_BUDGET §3 for `float` (2⁻²⁴ = 5.96 × 10⁻⁸). Note the
 * distinction from `std::numeric_limits<T>::epsilon()`, which is the *spacing* 2⁻²³ — twice this.
 * Bounds below are stated in `u`, per Higham, *Accuracy and Stability of Numerical Algorithms*
 * (2nd ed., §2.2), which is the convention that makes "one operation" and "n operations"
 * comparable.
 *
 * @tparam T the floating-point format.
 */
template <std::floating_point T>
inline constexpr double unit_roundoff = std::numeric_limits<T>::epsilon() / T{2};
// The halving is done in T and then widened, which is exact either way: an epsilon is a power of
// two, so dividing it by two is exact in any binary format, and widening to double is exact for
// every format narrower than it.

/**
 * The worst-case relative error of summing @p terms values naively (left to right) in @p T.
 *
 * Higham §4.2 gives `|E| ≤ (n−1)u / (1 − (n−1)u)` for recursive summation of n terms. What is
 * returned is the numerator `(n−1)u`, so it understates that bound by the denominator — which
 * differs from one by `(n−1)u` itself, a part in ~10⁴ for the fp32 cases here and a part in ~10¹³
 * for fp64, both far under the factor-of-two the budget's own spelling already carries (see the
 * note below). It is a worst case — every rounding error the same sign and maximal — and the
 * realistic figure is the random walk of @ref random_walk_estimate.
 *
 * @param terms how many values are summed; `0` and `1` involve no addition and cost nothing.
 * @tparam T the format the *accumulator* is in — the point of the whole precision policy is that
 *           this is `double` even when the integrand is `float`.
 * @return the bound as a relative error, `0.0` for fewer than two terms.
 * @complexity O(1).
 * @alloc none.
 * @note ERROR_BUDGET §3's "naive fp32 sum of 10³ terms ~1.2 × 10⁻⁴" row states the same bound in
 *       ε rather than u, so it is 2× this function at the same n. That row is the conservative
 *       spelling; both say the reduction cannot be fp32.
 */
template <std::floating_point T>
[[nodiscard]] constexpr double naive_sum_bound(std::size_t terms) {
    if (terms < 2) return 0.0;  // n−1 additions; zero of them for a single term or none
    return static_cast<double>(terms - 1) * unit_roundoff<T>;
}

/**
 * The realistic relative error of summing @p terms values in @p T, treating the individual
 * roundings as independent and mean-zero: `√n · u`.
 *
 * This is the `√N·ε` row of ERROR_BUDGET §3 (which, as in @ref naive_sum_bound, states it in ε
 * rather than u, so its ~4 × 10⁻⁶ at n = 10³ is twice what this returns). It is an *estimate*,
 * rather than a bound — a summation whose
 * errors correlate (a monotone integrand, which a field-line quadrature very much is) can exceed
 * it, up to @ref naive_sum_bound. Both are quoted in the budget because the honest margin lies
 * between them: one to three orders above the discretization floor, not the two to three the
 * planning estimate assumed.
 *
 * @param terms how many values are summed.
 * @tparam T the accumulator format.
 * @return the estimated relative error; `0.0` for no terms.
 * @complexity O(1) — one square root.
 * @alloc none.
 * @note Not `constexpr`: `std::sqrt` is not a constant expression before C++26.
 */
template <std::floating_point T>
[[nodiscard]] double random_walk_estimate(std::size_t terms) {
    return std::sqrt(static_cast<double>(terms)) * unit_roundoff<T>;
}

/**
 * The inverse question, and the one downstream code actually asks: how many terms may be summed
 * in @p T before naive summation alone breaches @p budget?
 *
 * Answering it is what decides where a reduction may run. At the `XJ` budget of 1 × 10⁻⁴ relative
 * (ERROR_BUDGET §4) this returns ~1.7 × 10³ for `float` — the same order as the term count of a
 * single field-line quadrature, i.e. the reduction would consume its whole budget — against
 * ~9 × 10¹¹ for `double`, where the question stops mattering. That contrast is the accumulator
 * invariant, in one number.
 *
 * @param budget the largest acceptable relative error; a non-positive or NaN budget admits no
 *               terms at all rather than being treated as unlimited.
 * @tparam T the accumulator format.
 * @return the largest n for which `(n−1)u ≤ budget`, saturated at `SIZE_MAX` for a budget so loose
 *         that the count is not representable.
 * @complexity O(1).
 * @alloc none.
 */
template <std::floating_point T>
[[nodiscard]] constexpr std::size_t max_terms_within(double budget) {
    constexpr std::size_t saturated = std::numeric_limits<std::size_t>::max();
    if (!(budget > 0.0)) return 0;  // the negation also rejects NaN
    const double n = budget / unit_roundoff<T> + 1.0;
    if (n >= static_cast<double>(saturated)) return saturated;
    return static_cast<std::size_t>(n);
}

// -------------------------------------------------------------------------------------------
// Precision policy
// -------------------------------------------------------------------------------------------

/**
 * A precision policy: the arithmetic type the integrand is evaluated in, and the type the results
 * are accumulated in.
 *
 * Structural only — it says the two names exist. It deliberately does NOT check that `accum` is
 * `double`, so that `SoundPrecision` can be the thing that rejects an unsound policy and the
 * rejection can be tested rather than merely asserted.
 *
 * @tparam P the candidate policy type.
 */
template <class P>
concept Precision = requires {
    typename P::integrand;
    typename P::accum;
} && std::floating_point<typename P::integrand> && std::floating_point<typename P::accum>;

/**
 * Whether @p P accumulates in `double` — the invariant of this header, stated as a trait so a
 * downstream `static_assert` can name it directly.
 * @tparam P a `Precision` policy.
 */
template <Precision P>
inline constexpr bool accumulates_in_double = std::same_as<typename P::accum, double>;

/**
 * A precision policy that may actually be used: a `Precision` whose accumulator is `double`.
 *
 * The integrand type is free — that choice is bounded by `unit_roundoff`, which does not grow
 * with the size of the problem. The accumulator type is not free, because `naive_sum_bound`
 * does grow with it, and at the ~10³ terms of one field-line quadrature an fp32 accumulator lands
 * on the discretization floor itself (ERROR_BUDGET §3). Anything constrained on this concept
 * cannot be instantiated with such a policy.
 *
 * @tparam P the candidate policy type.
 */
template <class P>
concept SoundPrecision = Precision<P> && accumulates_in_double<P>;

/**
 * The fp64 reference lane: everything in `double`.
 *
 * What the other lanes are measured against, and the only lane whose disagreement with IRBEM is
 * attributable to algorithm rather than to arithmetic. Slower on the GPU by roughly the ratio of
 * fp64 to fp32 throughput, and on much consumer hardware not runnable there at all.
 */
struct Exact {
    /// The type the field models and the integrand are evaluated in.
    using integrand = double;
    /// The type every sum is accumulated in — `double`, always. See the file brief.
    using accum = double;
    /// Human-readable lane name, for the differential suite's per-lane report.
    static constexpr std::string_view name = "Exact";
};

/**
 * The GPU lane: fp32 integrand, fp64 accumulator.
 *
 * fp32 is not merely faster here, it is *required*: SPIR-V's `GLSL.std.450` transcendentals are
 * defined only for 16- and 32-bit floats, so a kernel that calls `exp` — which every Tsyganenko
 * external field model does — cannot be compiled in fp64 at all. ERROR_BUDGET §3 shows that
 * constraint and the error budget point the same way, with the one exception this policy encodes:
 * the reduction stays in `double`.
 */
struct Fast {
    /// The type the field models and the integrand are evaluated in.
    using integrand = float;
    /// The type every sum is accumulated in — `double`, always, even here. See the file brief.
    using accum = double;
    /// Human-readable lane name, for the differential suite's per-lane report.
    static constexpr std::string_view name = "Fast";
};

// -------------------------------------------------------------------------------------------
// Compatibility policy
// -------------------------------------------------------------------------------------------

/**
 * A compatibility policy: four independent switches saying, for each place where IRBEM's algorithm
 * and the improved one differ, which is taken.
 *
 * Each is a `static constexpr bool` so downstream code reads `if constexpr (C::carry_k1)` and
 * neither branch survives into the binary. Requiring them through `std::bool_constant` is what
 * makes "is a constant expression" part of the concept rather than a convention.
 *
 * @tparam C the candidate policy type.
 */
template <class C>
concept Compat = requires {
    typename std::bool_constant<C::carry_k1>;
    typename std::bool_constant<C::independent_azimuth_seed>;
    typename std::bool_constant<C::bracketed_root_find>;
    typename std::bool_constant<C::transformed_quadrature>;
};

/**
 * Reproduce IRBEM's algorithm, decision for decision — every switch off.
 *
 * This is the lane the differential suite runs when it compares against the oracle, because a
 * comparison is only informative when both sides discretize the same way: ERROR_BUDGET §2(a) shows
 * the 0.01 L\* target is *at or below* IRBEM's own error at its recommended resolution, so a
 * comparison across settings measures resolution rather than correctness. Not the lane to ship.
 */
struct IrbemFaithful {
    /**
     * Off — the field is evaluated a fifth time at each step's arrival point, and that value is
     * discarded rather than reused as the next step's first stage.
     *
     * The redundancy is arithmetically invisible (the same function of the same bits gives the
     * same bits), so this switch is about the *call sequence*, not the numbers: a faithful lane
     * calls the field model the same number of times in the same order, which matters while
     * cross-checking against a library whose field models carry mutable `COMMON` state.
     */
    static constexpr bool carry_k1 = false;

    /**
     * Off — the drift shell's azimuths are walked in sequence, each seeded from the previous
     * azimuth's converged result. Correct, and serial: the shell's ~25·(options(4)+1) azimuths
     * cannot be evaluated concurrently, and the answer depends on the order they were visited in.
     */
    static constexpr bool independent_azimuth_seed = false;

    /**
     * Off — the mirror point is the θ-grid point at which the field crosses `Bm`, on the fixed
     * grid `dθ = π/(720·(options(3)+1))`. The located root therefore carries O(dθ) error that
     * *jitters* as the grid changes, which is the non-monotonic convergence measured in
     * ERROR_BUDGET §2(b) (`options=4` repeatedly worse than `options=2`).
     */
    static constexpr bool bracketed_root_find = false;

    /**
     * Off — the I integral `∫ √(1 − B/Bm) ds` between mirror points (Roederer, *Dynamics of
     * Geomagnetically Trapped Radiation*, 1970, §2) is summed as first-order rectangles on the θ
     * grid. The integrand's derivative is singular at both endpoints, where `B → Bm`, so the
     * rectangle rule converges at O(h) with a large endpoint constant. ERROR_BUDGET §2(b) names
     * this the single largest accuracy lever in the original algorithm.
     */
    static constexpr bool transformed_quadrature = false;

    /// Human-readable lane name, for the differential suite's per-lane report.
    static constexpr std::string_view name = "IrbemFaithful";
};

/**
 * The improved algorithm — every switch on. The lane this module ships.
 *
 * Each switch is justified against the error budget individually in its own documentation below;
 * collectively they change L\* by up to the discretization figure of ERROR_BUDGET §2, which is
 * precisely why they are a compile-time choice and why the differential suite compares
 * @ref IrbemFaithful against the oracle instead.
 */
struct Improved {
    /**
     * On — the derivative already evaluated at a step's arrival point becomes the next step's
     * first stage (first-same-as-last), removing one field evaluation in five from the loop that
     * dominates the ~10⁵ model calls of a single L\* point.
     *
     * Bit-exact, not merely within budget: the next step starts from exactly the arrival point, so
     * the carried value is the value the discarded evaluation would have produced. The only way it
     * differs is if the field model is not a pure function of its arguments.
     */
    static constexpr bool carry_k1 = true;

    /**
     * On — each azimuth of the drift shell is seeded from its own local geometry rather than from
     * its neighbour's answer, so the azimuths are independent work items (the shape the GPU lane
     * needs) and the result no longer depends on the traversal order.
     *
     * Safe because the seed selects only which iterate the search starts from, and the bracketed
     * root-find below converges to the bracketed root regardless — to a tolerance far under the
     * `Bmirr` budget of 1e-5 relative (ERROR_BUDGET §4). The exception is a shell where the root
     * is *not* unique in the bracket (shell splitting, a field with multiple minima along the
     * line): there the two seedings can legitimately find different roots, and that case must be
     * reported as a named condition rather than silently resolved by whichever seed was used.
     */
    static constexpr bool independent_azimuth_seed = true;

    /**
     * On — the θ march is used only to *bracket* the `B = Bm` crossing, and the root is then
     * solved inside that bracket to a tolerance set by the budget rather than by the grid.
     *
     * This strictly reduces the error at a given resolution: the bracketing march is the same one,
     * and the refinement removes the O(dθ) term the march leaves behind, along with the grid
     * jitter that makes IRBEM's convergence non-monotonic (ERROR_BUDGET §2(b)). It costs a
     * handful of extra field evaluations per mirror point, against the ~10⁵ of an L\* point.
     */
    static constexpr bool bracketed_root_find = true;

    /**
     * On — the I integral is substituted so the square-root vanishing at the mirror points becomes
     * a smooth integrand (write `B/Bm = 1 − u²`, so `√(1 − B/Bm) = |u|` and the endpoint
     * singularity is absorbed into the Jacobian), and a higher-order rule is then applied to a
     * function that is actually differentiable at its endpoints.
     *
     * This is the switch that makes the 1e-4 relative budget on `XJ` (ERROR_BUDGET §4) reachable
     * at the default resolution at all; an O(h) rectangle rule is not within an order of it. It is
     * also the largest of the four changes, so a lane with this on must never be compared
     * point-for-point against the oracle and called agreement.
     */
    static constexpr bool transformed_quadrature = true;

    /// Human-readable lane name, for the differential suite's per-lane report.
    static constexpr std::string_view name = "Improved";
};

// -------------------------------------------------------------------------------------------
// The aggregate
// -------------------------------------------------------------------------------------------

/**
 * A complete lane: one precision policy and one compatibility policy, bound into the single type
 * that gets threaded through the tracer, the quadrature and the L\* driver.
 *
 * It re-exports the members both concepts require, so a `Policy` satisfies `SoundPrecision` and
 * `Compat` itself and a routine may be constrained on whichever it actually uses — a field-model
 * evaluator needs only the precision, a quadrature needs only the compat flags, and neither has to
 * name a whole lane to say so. What it does not re-export is the two halves' names, which stay
 * where they are: a report labels a row `Policy::precision::name` / `Policy::compat::name`.
 *
 * @tparam P the precision policy; constrained on `SoundPrecision`, so a policy accumulating in
 *           anything but `double` is a compile error here rather than a quiet loss of the whole
 *           error budget (see the file brief).
 * @tparam C the compatibility policy.
 */
template <SoundPrecision P, Compat C>
struct Policy {
    /// The precision half, for a routine that needs to pass the lane on.
    using precision = P;
    /// The compatibility half, likewise.
    using compat = C;

    /// The arithmetic type the integrand is evaluated in.
    using integrand = typename P::integrand;
    /// The arithmetic type every sum accumulates in — `double`, guaranteed by the constraint.
    using accum = typename P::accum;

    /// Whether the tracer carries the previous step's final derivative into the next step's first
    /// stage; see @ref IrbemFaithful::carry_k1 and @ref Improved::carry_k1.
    static constexpr bool carry_k1 = C::carry_k1;
    /// Whether each drift-shell azimuth is seeded independently of its neighbour; see
    /// @ref IrbemFaithful::independent_azimuth_seed and @ref Improved::independent_azimuth_seed.
    static constexpr bool independent_azimuth_seed = C::independent_azimuth_seed;
    /// Whether the mirror point is refined inside the bracket the theta march found; see
    /// @ref IrbemFaithful::bracketed_root_find and @ref Improved::bracketed_root_find.
    static constexpr bool bracketed_root_find = C::bracketed_root_find;
    /// Whether the I integral is substituted before quadrature; see
    /// @ref IrbemFaithful::transformed_quadrature and @ref Improved::transformed_quadrature.
    static constexpr bool transformed_quadrature = C::transformed_quadrature;
};

/// The lane the differential suite runs: IRBEM's algorithm in fp64, so a disagreement with the
/// oracle is attributable to neither arithmetic nor resolution.
using OraclePolicy = Policy<Exact, IrbemFaithful>;

/// The CPU production lane: the improved algorithm in fp64. The tightest answer this module has,
/// and what the GPU lane is measured against (ERROR_BUDGET §6 — the GPU lane cannot carry
/// committed bit-goldens, because FMA contraction is at the driver's discretion).
using ReferencePolicy = Policy<Exact, Improved>;

/// The throughput lane: the improved algorithm with an fp32 integrand and an fp64 ordered
/// reduction.
using GpuPolicy = Policy<Fast, Improved>;

}  // namespace cheatah::space::irbem
