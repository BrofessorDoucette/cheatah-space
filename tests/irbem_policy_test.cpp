// Unit tests for space.irbem's policy layer — the precision axis (which arithmetic type the
// integrand uses versus the accumulator), the compatibility axis (IRBEM's algorithm versus the
// improved one), and the Policy aggregate that binds them.
//
// Most of this surface is types, so most of the test is `static_assert`: the assertions that
// matter are the ones about what does NOT compile. The load-bearing one is that a policy
// accumulating in anything but `double` is rejected — an fp32 reduction over the ~10^3 terms of a
// field-line quadrature costs the same order as the entire discretization floor
// (docs/ERROR_BUDGET.md §3), so it must be impossible to introduce, not merely discouraged.
// A rejection is tested by asking whether the constrained construct is well-formed at all, inside
// a `requires` expression, which is the only way to assert a compile error without a compile error.
//
// The roundoff arithmetic is checked at exactly-representable points — powers of two, and term
// counts chosen so the products are exact in binary — so every assertion is `==`. A tolerance here
// would defeat the purpose: these functions ARE the tolerances everything else is judged by.
//
// The bound arguments are read out of a table at run time rather than written as literals, so the
// constexpr functions are genuinely executed (and so instrumented) instead of being folded.
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <string_view>
#include <type_traits>

#include "space/irbem/policy.hpp"

namespace ib = cheatah::space::irbem;

// ---- policies that must be REFUSED ---------------------------------------------------------

/// Structurally a precision policy, and unsound: summing in fp32 is the one thing the budget
/// forbids. This is the type SoundPrecision exists to reject.
struct SloppyFast {
    using integrand = float;
    using accum = float;
};

/// Not a precision policy at all — no member types.
struct NotAPrecision {};

/// Has both names, but `integrand` is not a floating-point type.
struct IntegerIntegrand {
    using integrand = int;
    using accum = double;
};

/// Has `integrand` and not `accum`. Distinct from @ref NotAPrecision, which is missing both: a
/// fixture missing both proves only that *some* typedef is required, and the `accum` half of the
/// requirement could then be deleted from the concept with no test noticing. (Checked by deleting
/// it — the suite still passed until this fixture existed.)
struct NoAccumType {
    using integrand = double;
};

/// Has both names, and it is the *accumulator* that is not floating-point. The mirror of
/// @ref IntegerIntegrand, and required for the same reason: with only the integrand fixture, the
/// concept's `std::floating_point<typename P::accum>` clause was deletable undetected.
struct IntegerAccum {
    using integrand = double;
    using accum = int;
};

/// Long double accumulates *more* than double, and is still refused: the invariant is a fixed,
/// reproducible accumulator type, not merely a wide one. x87 80-bit arithmetic reproduces across
/// neither platforms nor optimization levels, which ERROR_BUDGET §6 requires of the reference lane.
struct WiderThanDouble {
    using integrand = double;
    using accum = long double;
};

/// Structurally almost a compat policy — but the flags are not constant expressions, so they
/// cannot be branched on with `if constexpr`, which is the entire point of the axis.
struct RuntimeFlags {
    static bool carry_k1;
    static bool independent_azimuth_seed;
    static bool bracketed_root_find;
    static bool transformed_quadrature;
};

// A concept that requires four things is only *proven* to require all four by four separate
// refusals. With a single "missing one flag" fixture, three of the four `typename
// std::bool_constant<...>` lines could be deleted from `Compat` and the whole suite still passed —
// verified by deleting each of them in turn. So there is one fixture per flag, each omitting
// exactly that flag and nothing else.

/// Missing `transformed_quadrature`, and otherwise a complete compat policy.
struct PartialCompat {
    static constexpr bool carry_k1 = true;
    static constexpr bool independent_azimuth_seed = true;
    static constexpr bool bracketed_root_find = true;
};

/// Missing `carry_k1`.
struct NoCarryK1Flag {
    static constexpr bool independent_azimuth_seed = true;
    static constexpr bool bracketed_root_find = true;
    static constexpr bool transformed_quadrature = true;
};

/// Missing `independent_azimuth_seed`.
struct NoAzimuthSeedFlag {
    static constexpr bool carry_k1 = true;
    static constexpr bool bracketed_root_find = true;
    static constexpr bool transformed_quadrature = true;
};

/// Missing `bracketed_root_find`.
struct NoBracketedRootFindFlag {
    static constexpr bool carry_k1 = true;
    static constexpr bool independent_azimuth_seed = true;
    static constexpr bool transformed_quadrature = true;
};

// ---- the precision axis ---------------------------------------------------------------------

static_assert(ib::Precision<ib::Exact>);
static_assert(ib::Precision<ib::Fast>);
static_assert(ib::SoundPrecision<ib::Exact>);
static_assert(ib::SoundPrecision<ib::Fast>);

// The integrand type is the free choice; the accumulator is not.
static_assert(std::same_as<ib::Exact::integrand, double>);
static_assert(std::same_as<ib::Fast::integrand, float>);
static_assert(std::same_as<ib::Exact::accum, double>);
static_assert(std::same_as<ib::Fast::accum, double>);

// The rejection, three ways it can arrive.
static_assert(ib::Precision<SloppyFast>);        // structurally fine...
static_assert(!ib::SoundPrecision<SloppyFast>);  // ...and refused anyway
static_assert(!ib::accumulates_in_double<SloppyFast>);
static_assert(ib::accumulates_in_double<ib::Exact> && ib::accumulates_in_double<ib::Fast>);
static_assert(!ib::Precision<NotAPrecision> && !ib::SoundPrecision<NotAPrecision>);
static_assert(!ib::Precision<IntegerIntegrand> && !ib::SoundPrecision<IntegerIntegrand>);
static_assert(!ib::Precision<NoAccumType> && !ib::SoundPrecision<NoAccumType>);
static_assert(!ib::Precision<IntegerAccum> && !ib::SoundPrecision<IntegerAccum>);
static_assert(ib::Precision<WiderThanDouble> && !ib::SoundPrecision<WiderThanDouble>);

// A rejection is only observable through a constrained template that either is or is not
// well-formed, so the probes below stand in for "does not compile": one names the axis, one names
// the aggregate.
//
// An unsatisfied constraint on a *class* template-id IS a substitution failure in the immediate
// context, so `requires { typename ib::Policy<P, C>; }` is well-defined and yields `false` —
// verified on GCC 13.3 and Clang 18. That matters, because probing the aggregate through a
// function that *takes* a `Policy<P, C>` proves nothing about the constraint: deduction against
// a non-Policy argument fails first, so `Policy`'s `SoundPrecision` constraint could be dropped
// entirely (or its second parameter left unconstrained) with the whole suite still green. Both
// deletions were tried; only @ref LaneFormable catches them.

/// Callable only for a precision policy that accumulates in double.
template <ib::SoundPrecision P>
constexpr bool accepts_precision() {
    return ib::accumulates_in_double<P>;
}

/// Callable only for a lane whose precision half is sound.
template <ib::SoundPrecision P, ib::Compat C>
constexpr bool accepts_lane(ib::Policy<P, C> /*lane*/) {
    return true;
}

/// Whether the lane type `Policy<P, C>` can be named at all — i.e. whether both of `Policy`'s own
/// constraints are satisfied. This is the direct test of the header's load-bearing claim that an
/// unsound lane "cannot even be named".
template <class P, class C>
concept LaneFormable = requires { typename ib::Policy<P, C>; };

/// Whether @ref accepts_precision can be called with @p P at all — the probe is itself a template,
/// because a `requires` expression written outside one is evaluated eagerly and a failed
/// constraint becomes a hard error instead of `false` (both GCC 13 and Clang do this).
template <class P>
concept PrecisionIsAccepted = requires { accepts_precision<P>(); };

/// Likewise for a whole lane.
template <class L>
concept LaneIsAccepted = requires { accepts_lane(L{}); };

static_assert(PrecisionIsAccepted<ib::Fast> && PrecisionIsAccepted<ib::Exact>);
static_assert(!PrecisionIsAccepted<SloppyFast>);
static_assert(!PrecisionIsAccepted<WiderThanDouble>);
static_assert(!PrecisionIsAccepted<NotAPrecision>);

// ...and the constraint reaches the aggregate: the unsound lane cannot be formed at all, so no
// call taking it can be written.
static_assert(LaneIsAccepted<ib::GpuPolicy> && LaneIsAccepted<ib::OraclePolicy>);
static_assert(!LaneIsAccepted<SloppyFast>);

// The aggregate's own constraints, tested where they live. Every sound precision policy pairs with
// every compat policy...
static_assert(LaneFormable<ib::Exact, ib::IrbemFaithful> && LaneFormable<ib::Exact, ib::Improved>);
static_assert(LaneFormable<ib::Fast, ib::IrbemFaithful> && LaneFormable<ib::Fast, ib::Improved>);
// ...and every unsound one pairs with none: `Policy<SloppyFast, Improved>` is not a type that
// exists. This is the assertion that makes the fp32-accumulator invariant unbypassable rather
// than merely unrecommended.
static_assert(!LaneFormable<SloppyFast, ib::Improved>);
static_assert(!LaneFormable<WiderThanDouble, ib::Improved>);
static_assert(!LaneFormable<NotAPrecision, ib::Improved>);
static_assert(!LaneFormable<IntegerAccum, ib::Improved>);
// The second parameter is constrained too, so a lane cannot be built from half a compat policy or
// from a precision policy mistakenly passed in the compat slot.
static_assert(!LaneFormable<ib::Exact, PartialCompat>);
static_assert(!LaneFormable<ib::Exact, NoCarryK1Flag>);
static_assert(!LaneFormable<ib::Exact, ib::Fast>);

// ---- the compatibility axis ------------------------------------------------------------------

static_assert(ib::Compat<ib::IrbemFaithful>);
static_assert(ib::Compat<ib::Improved>);
static_assert(!ib::Compat<RuntimeFlags>);  // flags must be usable in `if constexpr`
static_assert(!ib::Compat<ib::Exact>);     // a precision policy is not a compat policy

// All four, or none — one refusal per flag, so that no single requirement in the concept is
// deletable without a test noticing.
static_assert(!ib::Compat<NoCarryK1Flag>);
static_assert(!ib::Compat<NoAzimuthSeedFlag>);
static_assert(!ib::Compat<NoBracketedRootFindFlag>);
static_assert(!ib::Compat<PartialCompat>);  // missing transformed_quadrature

// The faithful lane takes IRBEM's decision at every one of the four points; the improved lane
// takes the other. That total opposition is the invariant that makes "faithful" meaningful — a
// flag flipped in only one of the two structs would silently produce a third, unnamed algorithm.
static_assert(!ib::IrbemFaithful::carry_k1 && ib::Improved::carry_k1);
static_assert(!ib::IrbemFaithful::independent_azimuth_seed &&
              ib::Improved::independent_azimuth_seed);
static_assert(!ib::IrbemFaithful::bracketed_root_find && ib::Improved::bracketed_root_find);
static_assert(!ib::IrbemFaithful::transformed_quadrature && ib::Improved::transformed_quadrature);

// ---- the aggregate ----------------------------------------------------------------------------

// A Policy re-exports enough of both halves to be usable wherever either half is required, so a
// routine constrains on the axis it actually reads rather than on a whole lane.
static_assert(ib::Precision<ib::ReferencePolicy> && ib::SoundPrecision<ib::ReferencePolicy>);
static_assert(ib::Compat<ib::ReferencePolicy>);
static_assert(std::same_as<ib::GpuPolicy::integrand, float>);
static_assert(std::same_as<ib::GpuPolicy::accum, double>);
static_assert(std::same_as<ib::GpuPolicy::precision, ib::Fast>);
static_assert(std::same_as<ib::GpuPolicy::compat, ib::Improved>);

// The three named lanes: what each one is for.
static_assert(std::same_as<ib::OraclePolicy, ib::Policy<ib::Exact, ib::IrbemFaithful>>);
static_assert(std::same_as<ib::ReferencePolicy::integrand, double>);
static_assert(!ib::OraclePolicy::transformed_quadrature && !ib::OraclePolicy::carry_k1 &&
              !ib::OraclePolicy::bracketed_root_find &&
              !ib::OraclePolicy::independent_azimuth_seed);
static_assert(ib::ReferencePolicy::transformed_quadrature && ib::ReferencePolicy::carry_k1 &&
              ib::ReferencePolicy::bracketed_root_find &&
              ib::ReferencePolicy::independent_azimuth_seed);
static_assert(ib::GpuPolicy::transformed_quadrature);

// A policy has no representation: it is looked at by the compiler and never stored, so binding one
// into a kernel launch or a traced state costs nothing.
static_assert(std::is_empty_v<ib::Exact> && std::is_empty_v<ib::Fast>);
static_assert(std::is_empty_v<ib::IrbemFaithful> && std::is_empty_v<ib::Improved>);
static_assert(std::is_empty_v<ib::GpuPolicy>);

// The lane names, reached through the halves — what the differential suite labels a row with.
static_assert(ib::ReferencePolicy::precision::name == std::string_view{"Exact"});
static_assert(ib::ReferencePolicy::compat::name == std::string_view{"Improved"});
static_assert(ib::OraclePolicy::compat::name == std::string_view{"IrbemFaithful"});
static_assert(ib::Fast::name == std::string_view{"Fast"});

// The roundoff constants are the ones ERROR_BUDGET §3 quotes: u(binary32) = 2^-24 = 5.96e-8, the
// "6e-8, single fp32 operation" row.
static_assert(ib::unit_roundoff<float> == 0x1p-24);
static_assert(ib::unit_roundoff<double> == 0x1p-53);

// ---- run-time exercises ------------------------------------------------------------------------

// Term counts and budgets are looked up rather than written inline, so the constexpr functions
// under test are actually called (and instrumented) instead of folded at compile time.
constexpr std::array<std::size_t, 5> kTermCounts{0, 1, 2, 1001, 16385};

/// Launder @p value through a volatile so the compiler cannot constant-fold the `constexpr`
/// functions under test. Without this, calls with literal arguments are evaluated at compile time
/// and the coverage instrumentation never sees them — a green report over code that never ran.
/// @param value the value to make opaque.
/// @tparam T the value's type.
/// @return @p value, unchanged.
template <class T>
T opaque(T value) {
    volatile T sink = value;
    return sink;
}

/// The bound is (n-1)u — n-1 additions, each costing at most one unit roundoff.
/// @tparam T the accumulator format under test.
template <std::floating_point T>
void exercise_sum_bound() {
    EXPECT_EQ(0.0, ib::naive_sum_bound<T>(opaque(kTermCounts[0])));  // nothing summed, nothing lost
    EXPECT_EQ(0.0, ib::naive_sum_bound<T>(opaque(kTermCounts[1])));  // one term is not an addition
    EXPECT_EQ(ib::unit_roundoff<T>, ib::naive_sum_bound<T>(opaque(kTermCounts[2])));

    // 1001 terms = 1000 additions. 1000 = 125*8, so 1000*2^-k is exact in binary and this is ==.
    EXPECT_EQ(1000.0 * ib::unit_roundoff<T>, ib::naive_sum_bound<T>(opaque(kTermCounts[3])));

    // Monotone in the term count — the property that makes max_terms_within invertible.
    EXPECT_LT(ib::naive_sum_bound<T>(opaque(kTermCounts[3])),
              ib::naive_sum_bound<T>(opaque(kTermCounts[4])));
}

/// sqrt(n)*u, the realistic estimate; it must sit below the worst case and above zero.
/// @tparam T the accumulator format under test.
template <std::floating_point T>
void exercise_random_walk() {
    EXPECT_EQ(0.0, ib::random_walk_estimate<T>(opaque(kTermCounts[0])));

    // Neither 1001 nor 16385 is a perfect square, so pick the one the table sits next to:
    // sqrt(16384) = 128 exactly, and 128u is a power of two, so the product is exact.
    const std::size_t square = kTermCounts[4] - 1;  // 16384 = 128^2
    EXPECT_EQ(128.0 * ib::unit_roundoff<T>, ib::random_walk_estimate<T>(opaque(square)));

    // The estimate is an estimate: it must be strictly below the worst-case bound for any n > 1,
    // since sqrt(n) < n-1 there. If it ever exceeded it, one of the two formulae is wrong.
    EXPECT_LT(ib::random_walk_estimate<T>(opaque(kTermCounts[3])),
              ib::naive_sum_bound<T>(opaque(kTermCounts[3])));
}

/// max_terms_within inverts naive_sum_bound exactly: the answer fits, and one more does not.
/// @tparam T the accumulator format under test.
template <std::floating_point T>
void exercise_max_terms() {
    // 2^-10 is exactly representable and is an exact multiple of u for both formats, so the
    // round-trip lands on the boundary rather than near it.
    const std::array<double, 4> budgets{0x1p-10, 0.0, -1.0,
                                        std::numeric_limits<double>::quiet_NaN()};

    const std::size_t n = ib::max_terms_within<T>(opaque(budgets[0]));
    EXPECT_EQ(budgets[0], ib::naive_sum_bound<T>(opaque(n)))
        << "the largest admissible n must fit exactly";
    EXPECT_GT(ib::naive_sum_bound<T>(opaque(n + 1)), budgets[0]) << "and one more must not fit";

    // A budget of zero, a negative budget, and a NaN budget all admit nothing. NaN especially:
    // comparing it the naive way (`budget <= 0`) would let it through as "unlimited", which is the
    // failure mode where a corrupted tolerance silently disables the check it was guarding.
    EXPECT_EQ(0U, ib::max_terms_within<T>(opaque(budgets[1])));
    EXPECT_EQ(0U, ib::max_terms_within<T>(opaque(budgets[2])));
    EXPECT_EQ(0U, ib::max_terms_within<T>(opaque(budgets[3])));

    // A budget so loose the count is not representable saturates instead of wrapping.
    EXPECT_EQ(std::numeric_limits<std::size_t>::max(), ib::max_terms_within<T>(opaque(1e300)));
}

TEST(IrbemPolicy, RoundoffBoundsAreExactAtRepresentablePoints) {
    exercise_sum_bound<float>();
    exercise_sum_bound<double>();
    exercise_random_walk<float>();
    exercise_random_walk<double>();
    exercise_max_terms<float>();
    exercise_max_terms<double>();
}

// The reason `accum` is double, expressed as a test rather than a comment: at the ~10^3 terms of
// one field-line quadrature, an fp32 reduction lands on the discretization floor itself, while an
// fp64 reduction is ten orders below it. ERROR_BUDGET §2 measures that floor at 1.2e-3 relative
// (the tightest of the measured range) at IRBEM's recommended resolution.
TEST(IrbemPolicy, TheAccumulatorTypeIsTheOneChoiceTheBudgetForces) {
    constexpr double discretization_floor = 1.2e-3;
    const std::size_t terms = kTermCounts[3];  // 1001, i.e. 1000 additions

    const double fp32_sum = ib::naive_sum_bound<float>(opaque(terms));
    const double fp64_sum = ib::naive_sum_bound<double>(opaque(terms));

    // Within a factor of ~20 of the floor: the same order, which is the whole argument.
    EXPECT_GT(fp32_sum, discretization_floor / 100.0);
    // ...while a single fp32 operation is four orders below it, which is why the INTEGRAND may be
    // fp32 at all.
    EXPECT_LT(ib::unit_roundoff<float> * 1.0e4, discretization_floor);
    // ...and the fp64 reduction is nowhere near it.
    EXPECT_LT(fp64_sum * 1.0e9, discretization_floor);

    // The invariant, asserted at run time as well so the trait is exercised and not merely
    // compiled: every shipped policy accumulates in double, whatever its integrand.
    EXPECT_TRUE(ib::accumulates_in_double<ib::Fast>);
    EXPECT_TRUE(ib::accumulates_in_double<ib::Exact>);
    EXPECT_TRUE(ib::accumulates_in_double<ib::GpuPolicy>);
    EXPECT_FALSE(ib::accumulates_in_double<SloppyFast>);

    // Both named policies accumulate in the type that clears the budget, whatever their integrand.
    EXPECT_EQ(fp64_sum, ib::naive_sum_bound<ib::Fast::accum>(opaque(terms)));
    EXPECT_EQ(fp64_sum, ib::naive_sum_bound<ib::GpuPolicy::accum>(opaque(terms)));
    EXPECT_LT(fp64_sum, fp32_sum);
}

// The two ceilings policy.hpp quotes in prose, checked so the prose cannot rot: at the XJ budget
// of ERROR_BUDGET §4, an fp32 reduction runs out of room at the same order as the term count of a
// single field-line quadrature, while an fp64 one has eight orders to spare.
TEST(IrbemPolicy, TheDocumentedTermCeilingsHold) {
    constexpr double xj_budget = 1.0e-4;  // ERROR_BUDGET §4, "XJ (I integral) — 1e-4 rel"
    constexpr std::size_t quadrature_terms = 1000;  // the ~10^3 of ERROR_BUDGET §3

    // Exactly 1678: the arithmetic is one division and one truncation, so this is reproducible
    // rather than approximate. It is 1.7 quadratures' worth of terms — no headroom at all.
    EXPECT_EQ(1678U, ib::max_terms_within<float>(opaque(xj_budget)));
    EXPECT_LT(ib::max_terms_within<float>(opaque(xj_budget)), 2 * quadrature_terms);

    // ~9.0e11 for double. Bounded rather than pinned: the claim under test is the order of
    // magnitude, which is what makes the accumulator choice safe, not the last digit of it.
    const std::size_t fp64_ceiling = ib::max_terms_within<double>(opaque(xj_budget));
    EXPECT_GT(fp64_ceiling, 9.0e11);
    EXPECT_LT(fp64_ceiling, 9.1e11);
    EXPECT_GT(fp64_ceiling, quadrature_terms * 100000000U);
}

// How downstream code reads the compat axis. Both branches must exist and be selected by the
// policy alone — no runtime flag reaches the inner loop.
/// @tparam C the compat policy driving the branch.
/// @return the number of field evaluations per RK step the policy implies.
template <ib::Compat C>
constexpr int field_evals_per_step() {
    if constexpr (C::carry_k1) {
        return 4;  // first-same-as-last: the arrival derivative becomes the next step's k1
    } else {
        return 5;  // IRBEM re-evaluates at the arrival point and discards it
    }
}

TEST(IrbemPolicy, CompatFlagsSelectBranchesAtCompileTime) {
    EXPECT_EQ(5, field_evals_per_step<ib::IrbemFaithful>());
    EXPECT_EQ(4, field_evals_per_step<ib::Improved>());
    // ...and a whole lane can be handed to a routine that only wants the compat half.
    EXPECT_EQ(4, field_evals_per_step<ib::GpuPolicy>());
    EXPECT_EQ(5, field_evals_per_step<ib::OraclePolicy>());

    // The saving is one evaluation in five of the ~10^5 model calls a single L* point costs.
    EXPECT_EQ(1, field_evals_per_step<ib::IrbemFaithful>() - field_evals_per_step<ib::Improved>());
}
