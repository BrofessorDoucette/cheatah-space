#pragma once

/**
 * @file igrf.hpp
 * @brief space.irbem — the IGRF internal geomagnetic field, evaluated without a trig call.
 *
 * Everything else in this module stands on this kernel. An L\* evaluation costs on the order of
 * 10⁵ field-model calls, and the internal field is in every one of them, so its cost sets the cost
 * of the library. Two things dominate a naive implementation, and this file removes both.
 *
 * **The normalisation is computed at compile time.** The Schmidt semi-normalised Legendre
 * recursion needs the factors `(2n-1)/√(n²-m²)` and `√((n-1)²-m²)/√(n²-m²)` for every degree/order
 * pair, plus the constant value of each diagonal term `Aⁿₙ`. IRBEM rebuilds those square roots on
 * every call — 224 `sqrt`s per field evaluation that depend on nothing but the loop indices. Here
 * they are a `static constexpr` table built by @ref detail::make_legendre_normalisation using a
 * `constexpr` Newton square root, so the emitted kernel contains **no `sqrt` at all** beyond the
 * single one that turns a position into a radius.
 *
 * **The kernel is written in Cartesian, so there is no trigonometry and nothing to divide by
 * zero.** The textbook formulation evaluates `Pⁿₘ(cos θ)`, `cos mφ` and `sin mφ`, and then divides
 * the eastward component by `sin θ` — which is a `sin`/`cos` pair, a chain of angle recursions and
 * a pole singularity that every implementation special-cases with a branch. Factor the latitude
 * dependence out instead: write `Pⁿₘ(cos θ) = Aⁿₘ(u)·sinᵐθ` with `u = cos θ`, and carry
 * `cₘ = sinᵐθ·cos mφ` and `sₘ = sinᵐθ·sin mφ` rather than the angles themselves. `Aⁿₘ` is a
 * polynomial in `u = z/r`; `cₘ`, `sₘ` obey `cₘ = x̂·cₘ₋₁ - ŷ·sₘ₋₁`, `sₘ = x̂·sₘ₋₁ + ŷ·cₘ₋₁` with
 * `x̂ = x/r`, `ŷ = y/r`. Every `sinᵐθ` and every `1/sin θ` then cancels algebraically before any
 * arithmetic happens (see @ref Igrf::evaluate for the three closed forms), leaving a kernel that is
 * **one square root, zero transcendentals, zero divisions by a quantity that can vanish, and zero
 * branches** — the pole is an ordinary point, not a case.
 *
 * The physics is the published IGRF definition and nothing else. The field is the negative gradient
 * of the scalar potential
 *
 *     V(r, θ, φ) = a Σₙ (a/r)ⁿ⁺¹ Σₘ [gⁿₘ cos mφ + hⁿₘ sin mφ] Pⁿₘ(cos θ)
 *
 * with `a = 6371.2 km` the IGRF reference radius, `θ` geocentric colatitude, `φ` east longitude,
 * and `Pⁿₘ` Schmidt semi-normalised — the convention stated in the header of the IAGA coefficient
 * file itself and in Alken et al., *International Geomagnetic Reference Field: the fourteenth
 * generation*, Earth Planets Space (2025). The Legendre recursions are the standard ones (Winch,
 * Ivers, Turner & Stening, *Geomagnetism and Schmidt quasi-normalization*, Geophys. J. Int. 160
 * (2005) 487-504; also the recursion set given in the NOAA/BGS World Magnetic Model technical
 * reports). The coefficients live in @ref igrf14.hpp, generated from IAGA's own release file.
 *
 * @note Degree. IGRF-14 publishes degree 13, and @ref Igrf defaults to it. IRBEM's internal field
 *       truncates at degree 10, which is a real difference of order 10 nT near the surface, not a
 *       rounding difference — `Igrf<10>` reproduces IRBEM's truncation for differential testing,
 *       and `Igrf<13>` is the model as IAGA published it. Neither is "compatibility mode": they are
 *       different truncations of the same series, and a comparison must say which it used.
 *
 * @note Units and the meaning of "Earth radii". @ref Position carries Cartesian components in Earth
 *       radii, and this file fixes that radius to the IGRF reference radius
 *       @ref Igrf::reference_radius_km = 6371.2 km — which is what makes `a/r` simply `1/r`. That
 *       is IRBEM's convention too. A position given in kilometres must be divided by that constant
 *       and not by an ellipsoid semi-axis.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>

#include "space/irbem/frames.hpp"
#include "space/irbem/policy.hpp"
#include "space/irbem/tables/igrf14.hpp"

namespace cheatah::space::irbem {

namespace detail {

/**
 * The number of (degree, order) slots in a triangular coefficient array up to degree @p n.
 * @param n the highest degree, `>= 0`.
 * @return `(n+1)(n+2)/2` — slots for every `(degree, order)` with `0 <= order <= degree <= n`.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::size_t triangular_slots(int n) {
    return static_cast<std::size_t>((n + 1) * (n + 2) / 2);
}

/**
 * The flat index of the coefficient of degree @p n and order @p m.
 *
 * The triangular packing `n(n+1)/2 + m` is what lets the degree-`n-2` term of the Legendre
 * recursion be read without a bounds test: at `m = n-1` it aliases the (already computed,
 * always finite) slot `(n-1, 0)`, and the recursion multiplies that read by an exact zero.
 *
 * @param n the degree, `>= 0`.
 * @param m the order, `0 <= m <= n`.
 * @return the slot index.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::size_t slot_index(int n, int m) {
    return static_cast<std::size_t>((n * (n + 1) / 2) + m);
}

/**
 * A square root usable in a constant expression.
 *
 * `std::sqrt` is not `constexpr` before C++26, and the whole point of the normalisation table is
 * that its square roots never run. Three steps, and the third is what makes the answer exactly the
 * one `std::sqrt` would give rather than merely close to it:
 *
 *  1. Scale the radicand into `[0.25, 4)` by exact powers of four, so the seed `2·scale` is within
 *     a factor of two of the root for any exponent. (Seeding with `x` itself would need hundreds of
 *     iterations for a tiny or huge argument, and no fixed budget would be honest.)
 *  2. Newton-Raphson from above, which descends monotonically, until the step stops decreasing —
 *     that lands within one ulp.
 *  3. One more Newton step taken with the *exact* residual `g² - x`, recovered by Dekker's
 *     two-product splitting (the `2²⁷+1` constant), since `g*g` alone loses precisely the bits that
 *     decide the last ulp. The correction is smaller than an ulp, so the subtraction rounds to the
 *     correctly rounded root. Verified against `std::sqrt` over a wide sweep in the test.
 *
 * @param x the radicand. Negative values and NaN return `0.0` — this is an internal helper whose
 *          callers pass `n² - m² >= 0` by construction, and a silent zero is preferable to a
 *          non-constant-evaluable throw in a table initializer.
 * @return the correctly rounded square root of @p x.
 * @complexity O(1) — bounded scaling plus at most 200 Newton steps, all at compile time.
 * @alloc none.
 */
constexpr double const_sqrt(double x) {
    if (!(x > 0.0)) return 0.0;
    double scale = 1.0;
    double reduced = x;
    while (reduced >= 4.0) {
        reduced *= 0.25;
        scale *= 2.0;
    }
    while (reduced < 0.25) {
        reduced *= 4.0;
        scale *= 0.5;
    }
    double guess = 2.0 * scale;  // >= sqrt(x), because sqrt(reduced) <= 2
    for (int i = 0; i < 200; ++i) {
        const double next = 0.5 * (guess + (x / guess));
        if (next >= guess) break;
        guess = next;
    }
    // Dekker two-product: the error term of guess*guess, exactly.
    constexpr double kSplit = 134217729.0;  // 2^27 + 1
    const double squared = guess * guess;
    const double t = kSplit * guess;
    const double hi = t - (t - guess);
    const double lo = guess - hi;
    const double product_error = (((hi * hi) - squared) + (2.0 * hi * lo)) + (lo * lo);
    const double residual = (squared - x) + product_error;
    return guess - (residual / (2.0 * guess));
}

/**
 * The Schmidt-normalisation constants of the associated Legendre recursion, for degrees up to
 * @p NMAX, in the working type @p T.
 *
 * Recurrences (in `Aⁿₘ(u) = Pⁿₘ(u) / (1-u²)^{m/2}`, which is a polynomial in `u`):
 *  - diagonal: `Aⁿₙ = kₙ · Aⁿ⁻¹ₙ₋₁` with `k₁ = 1` and `kₙ = √((2n-1)/(2n))` for `n >= 2` — a
 *    *constant*, since the `sin θ` that the classic diagonal recursion carries is exactly the
 *    factor that was divided out. So @ref diagonal stores the finished values.
 *  - off-diagonal: `Aⁿₘ = eₙₘ·u·Aⁿ⁻¹ₘ - fₙₘ·Aⁿ⁻²ₘ`, and by differentiating,
 *    `A'ⁿₘ = eₙₘ·(Aⁿ⁻¹ₘ + u·A'ⁿ⁻¹ₘ) - fₙₘ·A'ⁿ⁻²ₘ`.
 *
 * The `k₁ = 1` exception is not a fudge: Schmidt semi-normalisation carries an extra `√2` for
 * every order `m >= 1`, so the single step that crosses from `m = 0` to `m = 1` has no `√((2n-1)/2n)`
 * factor. `P¹₁ = sin θ`, so `A¹₁ = 1`.
 *
 * @tparam NMAX the highest degree the table covers.
 * @tparam T the floating-point type the field kernel works in.
 */
template <int NMAX, class T>
struct LegendreNormalisation {
    /// `Aⁿₙ`, the (constant) diagonal term, indexed by degree `n`; `A⁰₀ = 1`.
    std::array<T, static_cast<std::size_t>(NMAX) + 1> diagonal{};
    /// `eₙₘ = (2n-1)/√(n²-m²)`, at @ref slot_index(n, m); zero where the recursion is unused.
    std::array<T, triangular_slots(NMAX)> e{};
    /// `fₙₘ = √((n-1)²-m²)/√(n²-m²)`, at @ref slot_index(n, m); exactly zero at `m = n-1`, which
    /// is what makes the aliased degree-`n-2` read harmless.
    std::array<T, triangular_slots(NMAX)> f{};
};

/**
 * Build @ref LegendreNormalisation — the whole point being that this runs at compile time.
 * @tparam NMAX the highest degree.
 * @tparam T the working floating-point type.
 * @return the filled table.
 * @complexity O(NMAX²) — 105 slots at degree 13, every one a compile-time Newton square root.
 * @alloc none (the table is returned by value into a `static constexpr` member).
 */
template <int NMAX, class T>
constexpr LegendreNormalisation<NMAX, T> make_legendre_normalisation() {
    LegendreNormalisation<NMAX, T> table;
    table.diagonal[0] = static_cast<T>(1.0);
    double diag = 1.0;
    for (int n = 1; n <= NMAX; ++n) {
        // n == 1 is the m=0 -> m=1 step, where Schmidt's extra sqrt(2) cancels the recursion's.
        const double k =
            (n == 1) ? 1.0
                     : const_sqrt(static_cast<double>((2 * n) - 1) / static_cast<double>(2 * n));
        diag *= k;
        table.diagonal[static_cast<std::size_t>(n)] = static_cast<T>(diag);
        for (int m = 0; m < n; ++m) {
            const double nn = static_cast<double>(n);
            const double mm = static_cast<double>(m);
            const double denom = const_sqrt((nn * nn) - (mm * mm));
            const double numer = const_sqrt(((nn - 1.0) * (nn - 1.0)) - (mm * mm));
            table.e[slot_index(n, m)] = static_cast<T>(((2.0 * nn) - 1.0) / denom);
            table.f[slot_index(n, m)] = static_cast<T>(numer / denom);
        }
    }
    return table;
}

/// Degrees per radian — the one place the spherical entry point converts.
inline constexpr double radians_per_degree = 3.14159265358979323846 / 180.0;

}  // namespace detail

/**
 * The IGRF internal geomagnetic field at one instant, ready to evaluate.
 *
 * The date is bound once, at construction, by @ref at — which is where the two things that *can*
 * fail live: an out-of-range year, and the choice between interpolating between tabulated epochs
 * and extrapolating with the published secular variation. @ref evaluate is therefore total,
 * branch-free and allocation-free, and an ephemeris of a million points pays the coefficient
 * interpolation once instead of a million times. `at` is `constexpr`, so a model for a fixed epoch
 * costs nothing at run time at all.
 *
 * @tparam NMAX the truncation degree, `1..13`. 13 is IGRF-14 in full; 10 is what IRBEM uses
 *              internally. See the file-level note on degree.
 * @tparam P the precision policy (@ref SoundPrecision). `P::integrand` is what the Legendre and
 *           longitude recursions are carried in; `P::accum` — which the concept pins to `double` —
 *           is what the ~200 harmonic terms are summed in. That split is exactly the one
 *           `docs/ERROR_BUDGET.md` §3 argues for: one fp32 evaluation costs ~6e-8 and is free
 *           against the discretization floor, but a naive fp32 sum of 10³ terms costs ~1.2e-4 and
 *           lands *on* it. Roundoff in an evaluation does not grow with the size of the problem;
 *           roundoff in a reduction does.
 */
template <int NMAX = tables::igrf14_max_degree, SoundPrecision P = Exact>
class Igrf;

/// Whether @p M is an @ref Igrf instantiation — the trait the drift-shell machinery uses to decide
/// whether its DEVICE fast-paths apply. The device kernels stage one internal model's coefficients;
/// a composed model (internal plus external) must take the host lane for the stages that have no
/// combined kernel yet, and `if constexpr (is_igrf_v<M>)` is what routes that at compile time
/// instead of a runtime branch that would still have to compile the staging for types it cannot
/// stage.
template <class>
inline constexpr bool is_igrf_v = false;
/// The specialization that answers `true`: exactly the @ref Igrf instantiations, at any degree
/// and precision policy, and nothing else — a composed model must not match, since matching is
/// what routes a model into device paths that stage only an internal field.
template <int NMAX, SoundPrecision P>
inline constexpr bool is_igrf_v<Igrf<NMAX, P>> = true;

template <int NMAX, SoundPrecision P>
class Igrf {
    static_assert(NMAX >= 1 && NMAX <= tables::igrf14_max_degree,
                  "IGRF-14 publishes degrees 1 through 13");

    /// Flat coefficient slots at this truncation; the packing matches @ref tables::igrf14_g, so
    /// truncating the model is a prefix copy rather than a re-index.
    static constexpr std::size_t kSlots = detail::triangular_slots(NMAX);

  public:
    /// The type the recursions and the stored coefficients are carried in.
    using integrand = typename P::integrand;

    /// The type the harmonic sums are accumulated in — `double`, by @ref SoundPrecision.
    using accum = typename P::accum;

  private:
    /// The compile-time normalisation table — 224 square roots that never run.
    static constexpr auto kNorm = detail::make_legendre_normalisation<NMAX, integrand>();

  public:
    /// The truncation degree this instantiation evaluates.
    static constexpr int degree = NMAX;

    /// The IGRF reference radius `a`, in kilometres — and hence the Earth radius that
    /// @ref Position's Cartesian components are counted in, throughout this module.
    static constexpr double reference_radius_km = 6371.2;

    /// The first tabulated epoch. Before this the model is not defined and @ref at reports so.
    static constexpr double earliest_year = 1900.0;

    /// The last tabulated main-field epoch; beyond it the published secular variation is used.
    static constexpr double latest_epoch_year = 2025.0;

    /// The end of the secular-variation prediction. Past this the model is not defined; @ref at
    /// reports so rather than extrapolating a five-year linear prediction indefinitely.
    static constexpr double latest_year = 2030.0;

    /**
     * The model at a decimal year.
     *
     * Between tabulated epochs the coefficients are linear in time — the IGRF's own defined
     * interpolation, not an approximation to something smoother. From @ref latest_epoch_year to
     * @ref latest_year they are the last main field plus the published secular variation times the
     * elapsed years. At a tabulated epoch the result is bit-identical to the table.
     *
     * @param decimal_year the epoch, e.g. `2007.5` for mid-2007.
     * @return the model, or `std::nullopt` when @p decimal_year lies outside
     *         `[earliest_year, latest_year]` or is NaN. The empty case is never quietly clamped:
     *         a 2035 ephemeris is a question this model cannot answer, and saying so is the answer.
     * @complexity O(NMAX²) — one pass over the 105 coefficient slots.
     * @alloc none.
     */
    static constexpr std::optional<Igrf> at(double decimal_year) {
        if (!(decimal_year >= earliest_year) || !(decimal_year <= latest_year)) {
            return std::nullopt;
        }
        Igrf model;
        model.year_ = decimal_year;

        constexpr std::size_t last = tables::igrf14_epoch_count - 1;
        if (decimal_year >= latest_epoch_year) {
            const double dt = decimal_year - latest_epoch_year;
            for (std::size_t i = 0; i < kSlots; ++i) {
                model.g_[i] = static_cast<integrand>(tables::igrf14_g[last][i] +
                                                     (dt * tables::igrf14_g_sv[i]));
                model.h_[i] = static_cast<integrand>(tables::igrf14_h[last][i] +
                                                     (dt * tables::igrf14_h_sv[i]));
            }
            return model;
        }
        // The epochs are exactly 5 years apart from 1900.0, so the bracket is arithmetic rather
        // than a search, and (year - epoch)/5 is exact at every tabulated epoch.
        constexpr double kSpacing = 5.0;
        const auto lower =
            static_cast<std::size_t>((decimal_year - earliest_year) / kSpacing);
        const double w = (decimal_year - tables::igrf14_epochs[lower]) / kSpacing;
        for (std::size_t i = 0; i < kSlots; ++i) {
            model.g_[i] = static_cast<integrand>(
                tables::igrf14_g[lower][i] +
                (w * (tables::igrf14_g[lower + 1][i] - tables::igrf14_g[lower][i])));
            model.h_[i] = static_cast<integrand>(
                tables::igrf14_h[lower][i] +
                (w * (tables::igrf14_h[lower + 1][i] - tables::igrf14_h[lower][i])));
        }
        return model;
    }

    /**
     * The epoch this model was built for.
     * @return the decimal year passed to @ref at.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] constexpr double year() const { return year_; }

    /**
     * A time-interpolated `g` Gauss coefficient.
     * @param n the degree. @param m the order.
     * @return `gⁿₘ` in nT, or exactly zero when `(n, m)` is outside this truncation — the
     *         coefficients above @ref degree, and every `n < 1`, are structurally absent from the
     *         model rather than merely unknown.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] constexpr integrand g(int n, int m) const {
        if (n < 1 || n > NMAX || m < 0 || m > n) return static_cast<integrand>(0);
        return g_[detail::slot_index(n, m)];
    }

    /**
     * A time-interpolated `h` Gauss coefficient.
     * @param n the degree. @param m the order.
     * @return `hⁿₘ` in nT, or exactly zero outside the truncation. `hⁿ₀` is zero by definition —
     *         `sin 0φ` is identically zero, so the coefficient has no meaning.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] constexpr integrand h(int n, int m) const {
        if (n < 1 || n > NMAX || m < 0 || m > n) return static_cast<integrand>(0);
        return h_[detail::slot_index(n, m)];
    }

    /**
     * The internal field at a geocentric Cartesian position — the hot kernel.
     *
     * Writing `u = z/r`, `x̂ = x/r`, `ŷ = y/r`, `Rₙ = r^-(n+2)`, `Gⁿₘ = gⁿₘcₘ + hⁿₘsₘ`, and
     * summing over the model,
     *
     *     T0 = Σ (n+1)Rₙ Gⁿₘ Aⁿₘ            T1 = Σ Rₙ Gⁿₘ A'ⁿₘ         T2 = Σ m Rₙ Gⁿₘ Aⁿₘ
     *     U  = Σ m Rₙ Aⁿₘ (gⁿₘcₘ₋₁ + hⁿₘsₘ₋₁)                          V  = Σ m Rₙ Aⁿₘ (hⁿₘcₘ₋₁ - gⁿₘsₘ₋₁)
     *
     * the Cartesian field is exactly
     *
     *     Bx = x̂(T0 + u·T1 + T2) - U      By = ŷ(T0 + u·T1 + T2) - V      Bz = u·T0 - (1-u²)T1 + u·T2
     *
     * The `1/sin θ` of the eastward component and the `sinᵐθ` of the Legendre functions have both
     * cancelled analytically — `cₘcos φ + sₘsin φ = sin θ·cₘ₋₁` is the identity that does it — so
     * the pole needs no test and no epsilon. The `m = 0` term reads a padded zero slot for `cₘ₋₁`
     * and is multiplied by `m`, so it contributes exactly nothing without a branch.
     *
     * @param p the position, geocentric Cartesian, in units of @ref reference_radius_km. Inside the
     *          Earth is arithmetically fine and physically meaningless — the series diverges below
     *          the core; callers that care enforce `r >= 1` themselves. `r == 0` is the one input
     *          with no answer and yields infinities rather than a silent number.
     * @return the field in nT, in the same geographic frame.
     * @complexity O(NMAX²) — 105 coefficient slots at degree 13, ~1900 flops, one `sqrt`, no
     *             transcendentals, no branches in the summation.
     * @alloc none — three fixed `std::array`s on the stack, 2 KiB at degree 13.
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::GEO>& p) const {
        const double x = p.v[0];
        const double y = p.v[1];
        const double z = p.v[2];
        const double r = std::sqrt((x * x) + (y * y) + (z * z));  // the one square root
        const double inv_r = 1.0 / r;
        const integrand ux = static_cast<integrand>(x * inv_r);
        const integrand uy = static_cast<integrand>(y * inv_r);
        const integrand uz = static_cast<integrand>(z * inv_r);

        // Aⁿₘ(u) and its u-derivative, triangular. Seeded at n = 0 and n = 1 so the general
        // recursion below never has to name a degree below zero.
        std::array<integrand, kSlots> a{};
        std::array<integrand, kSlots> da{};
        a[0] = static_cast<integrand>(1);
        a[1] = uz;
        da[1] = static_cast<integrand>(1);
        a[2] = kNorm.diagonal[1];
        for (int n = 2; n <= NMAX; ++n) {
            const std::size_t base = detail::slot_index(n, 0);
            const std::size_t prev = base - static_cast<std::size_t>(n);
            const std::size_t prev2 = prev - static_cast<std::size_t>(n - 1);
            for (int m = 0; m < n; ++m) {
                const std::size_t k = base + static_cast<std::size_t>(m);
                const std::size_t k1 = prev + static_cast<std::size_t>(m);
                const std::size_t k2 = prev2 + static_cast<std::size_t>(m);
                a[k] = (kNorm.e[k] * uz * a[k1]) - (kNorm.f[k] * a[k2]);
                da[k] = (kNorm.e[k] * (a[k1] + (uz * da[k1]))) - (kNorm.f[k] * da[k2]);
            }
            a[base + static_cast<std::size_t>(n)] =
                kNorm.diagonal[static_cast<std::size_t>(n)];  // derivative is zero: no u in it
        }

        // cs[m+1] = sinᵐθ·cos mφ, ss[m+1] = sinᵐθ·sin mφ. Index 0 is the pad that lets the m = 0
        // term of the U/V sums read "cₘ₋₁" without a branch.
        std::array<integrand, static_cast<std::size_t>(NMAX) + 2> cs{};
        std::array<integrand, static_cast<std::size_t>(NMAX) + 2> ss{};
        cs[1] = static_cast<integrand>(1);
        for (int m = 1; m <= NMAX; ++m) {
            const auto i = static_cast<std::size_t>(m);
            cs[i + 1] = (ux * cs[i]) - (uy * ss[i]);
            ss[i + 1] = (ux * ss[i]) + (uy * cs[i]);
        }

        accum t0 = 0.0;
        accum t1 = 0.0;
        accum t2 = 0.0;
        accum su = 0.0;
        accum sv = 0.0;
        accum rn = inv_r * inv_r * inv_r;  // (a/r)^(n+2) at n = 1, with a = 1 by choice of unit
        for (int n = 1; n <= NMAX; ++n) {
            const std::size_t base = detail::slot_index(n, 0);
            integrand p0 = static_cast<integrand>(0);
            integrand p1 = static_cast<integrand>(0);
            integrand p2 = static_cast<integrand>(0);
            integrand pu = static_cast<integrand>(0);
            integrand pv = static_cast<integrand>(0);
            for (int m = 0; m <= n; ++m) {
                const std::size_t k = base + static_cast<std::size_t>(m);
                const std::size_t i = static_cast<std::size_t>(m);
                const integrand gc = g_[k];
                const integrand hc = h_[k];
                const integrand bigg = (gc * cs[i + 1]) + (hc * ss[i + 1]);
                const integrand ga = bigg * a[k];
                const integrand fm = static_cast<integrand>(m);
                p0 += ga;
                p2 += fm * ga;
                p1 += bigg * da[k];
                const integrand ma = fm * a[k];
                pu += ma * ((gc * cs[i]) + (hc * ss[i]));
                pv += ma * ((hc * cs[i]) - (gc * ss[i]));
            }
            // Per-degree partials in P::integrand, totals in P::accum: the ordered fp64 reduction
            // the error budget requires, without paying a conversion inside the inner loop.
            t0 += rn * static_cast<accum>(n + 1) * static_cast<accum>(p0);
            t1 += rn * static_cast<accum>(p1);
            t2 += rn * static_cast<accum>(p2);
            su += rn * static_cast<accum>(pu);
            sv += rn * static_cast<accum>(pv);
            rn *= inv_r;
        }

        const accum w = static_cast<accum>(uz);
        const accum radial = t0 + (w * t1) + t2;
        return FieldVector<Frame::GEO>{
            fixarray::vec3d{(static_cast<double>(ux) * radial) - su,
                            (static_cast<double>(uy) * radial) - sv,
                            (w * t0) - ((1.0 - (w * w)) * t1) + (w * t2)}};
    }

    /**
     * The internal field at a geocentric *spherical* position.
     *
     * The convenience entry point, and the only one that costs trigonometry: it turns
     * (radius, latitude, east longitude) into Cartesian and calls the kernel. A loop over an
     * ephemeris should convert once and stay Cartesian rather than call this per point.
     *
     * @param p the position — radius in units of @ref reference_radius_km, geocentric latitude and
     *          east longitude in degrees.
     * @return the field in nT, in geographic Cartesian components (**not** local north/east/down).
     * @complexity O(NMAX²), plus two sine/cosine pairs.
     * @alloc none.
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::SPH>& p) const {
        const double lat = p.v[1] * detail::radians_per_degree;
        const double lon = p.v[2] * detail::radians_per_degree;
        const double rho = p.v[0] * std::cos(lat);
        return evaluate(Position<Frame::GEO>{fixarray::vec3d{rho * std::cos(lon),
                                                             rho * std::sin(lon),
                                                             p.v[0] * std::sin(lat)}});
    }

  private:
    /// Only @ref at builds one, because only @ref at can decide whether the date is answerable.
    constexpr Igrf() = default;

    /// `g` coefficients in nT at @ref year_, triangular packing.
    std::array<integrand, kSlots> g_{};
    /// `h` coefficients in nT at @ref year_, triangular packing.
    std::array<integrand, kSlots> h_{};
    /// The epoch the coefficients were interpolated to.
    double year_{};
};

}  // namespace cheatah::space::irbem
