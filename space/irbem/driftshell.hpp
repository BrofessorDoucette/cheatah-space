#pragma once

/**
 * @file driftshell.hpp
 * @brief space.irbem — Roederer's L\*: the drift shell, its flux, and the parallel that IRBEM's
 *        algorithm cannot express.
 *
 * `lstar.hpp` produces `I`, `B_min` and `B_mirr` for ONE field line, on the device, at 1.24 µs a
 * line. This header is the step that turns those into the number the module exists for:
 *
 * ```
 * for each of Nder azimuths around the drift shell:
 *     root-find the field line whose (I, B_m) matches the starting line's
 *     -> its ionospheric footpoint
 * Phi = ∬ B·dA over the polar cap bounded by those footpoints
 * L* = 2πk₀ / (Phi R_E)
 * ```
 *
 * ## Why this is the whole argument for a device lane
 *
 * IRBEM computes one L\* in a measured 15.5 ms (`docs/ERROR_BUDGET.md` §5, `-O2` build) — about 26
 * minutes for a 100 000-point ephemeris, single-threaded, with no threading available. Essentially
 * all of it is field-line tracing: ~10⁵ field evaluations per point. Those traces are INDEPENDENT
 * of each other, and every one of them is a serial RK4 chain that no hardware makes faster. The
 * only parallelism in the problem is *across* traces — and IRBEM cannot take it, for a reason worth
 * stating precisely.
 *
 * **IRBEM seeds azimuth *i*'s root-find by extrapolating from azimuths *i−1* and *i−2*.** That is a
 * good serial heuristic: consecutive azimuths of a drift shell are close, so the extrapolation
 * lands near the answer and the fixed-step θ march that follows is short. It is also the ONLY
 * reason its azimuth loop must run in order. Break that dependency — seed every azimuth
 * independently from the starting line's dipole `L_m` — and all `Nder` azimuths of all `ntime`
 * points become one batch. That is what this header does, and it is the entire structural
 * difference between the two implementations.
 *
 * The seed is worse in isolation, and it does not matter: an independent seed lands within ~7% of
 * the root (measured over L = 2…8), a bracketed Brent root-find (Brent, *Algorithms for
 * Minimization without Derivatives*, Prentice-Hall 1973, ch. 4) closes that in a bounded ~8 traces,
 * and the batch is `ntime × Nder` wide instead of 1. IRBEM's fixed-step θ march, by contrast, has
 * no bound at all — it steps until it arrives.
 *
 * ## The three stages, and where each runs
 *
 * | stage | shape | lane |
 * |---|---|---|
 * | reference line: `(B_m, I₀)` for the starting point | `ntime` traces | `trace_invariant_batch` |
 * | root-find: which line at each azimuth carries `(B_m, I₀)` | `ntime × Nder × trials` traces | `trace_invariant_batch` |
 * | footpoints: where those lines meet `r = 1` | `ntime × Nder` walks | `irbem_shell_foot_f32` |
 * | flux: `∬ B·dA` over the ragged polar cap | `ntime × Nder × ~96` cells | `irbem_igrf_f32` |
 *
 * Every stage is ONE batch. Nothing loops over points calling a scalar routine — that shape cannot
 * be accelerated no matter what hardware is present, which is exactly why `lstar.hpp` exposes
 * `trace_invariant_batch` rather than only `trace_invariant`.
 *
 * ## The residual is continuous where the physics is degenerate
 *
 * At each azimuth the unknown is the radius `r` at which the shell's field line crosses the
 * magnetic-equatorial (MAG `z = 0`) plane. The mirror field is conserved BY CONSTRUCTION — the
 * trace is seeded at that crossing with the local pitch angle `α = asin √(B(r)/B_m)`, so
 * `lstar.hpp` computes `B_m/sin²α = B_m` exactly — which leaves one equation, `I(r) = I₀`, in one
 * unknown. `I(r)` increases monotonically with `r`: a larger shell is a weaker field, a lower
 * `B_min`, and a longer bounce path between the same two mirror fields.
 *
 * That construction needs `B(r) ≤ B_m`, and below the shell it is not. Rather than declare those
 * radii invalid — which would leave the root sitting on the edge of a region with no residual at
 * all and make Brent return anywhere in it — the residual continues onto them analytically:
 *
 * ```
 * f(r) = I(r) − I₀                                  where B(r) ≤ B_m   (a trace)
 * f(r) = −I₀ − w·(B(r)/B_m − 1)                     where B(r) >  B_m   (no trace)
 * ```
 *
 * Both branches increase with `r`, and the analytic one is negative everywhere, so `f` has exactly
 * one sign change and Brent's bracket is safe. `w` is a positive scale in Earth radii; it changes
 * how fast the bracket closes and cannot move the root.
 *
 * **The two branches do not meet, and the size of the step is a known error term.** At the junction
 * the analytic branch reaches `−I₀`, but the traced branch does not start at `I = 0`: the seed is
 * the field line's MAGNETIC-EQUATORIAL-PLANE crossing, and that is not quite its `B_min`, so a
 * particle mirroring there still bounces through a short but non-zero path. The step is therefore
 * `I(B_min(r_junction))`, and it matters only for a particle whose own `I₀` is smaller than it —
 * i.e. one mirroring essentially at the field minimum. MEASURED over 40 (L, azimuth) pairs with
 * IGRF-14, expressed as the relative shell-radius error it can cause:
 *
 * | L | 2 | 3 | 4 | 6 | 8 |
 * |---|---|---|---|---|---|
 * | worst `δr/r` over 8 azimuths | 4.7e-3 | 1.4e-3 | 5.7e-4 | 9.4e-5 | 1.9e-5 |
 *
 * So it is bounded by ~0.009 in L\* at L = 2 and negligible above L = 4, it applies only to the
 * equatorially-mirroring limit, and it is the largest remaining term in this implementation's error
 * budget after the quadrature and the trace step. Closing it means seeding on the minimum-`B`
 * surface rather than the MAG `z = 0` plane — a per-candidate latitude minimisation, not a change
 * of algorithm — and it has not been done because it is not the term that is binding: of the 41
 * oracle comparisons measured (L = 1.5…9, both hemispheres, all longitudes, pitch angles 15…90°,
 * epochs 1900…2029) the one that leaves the 0.01 budget is an off-equatorial start, where the
 * dominant term is `lstar.hpp`'s fixed trace step and not this junction. See **Measured** below.
 *
 * ## Two things that are quietly wrong if you do the obvious thing
 *
 * **The cap boundary lives at the FOOTPOINT's longitude, not the seed's.** A field line does not
 * stay at one magnetic longitude on its way down: measured against IGRF-14, the footpoint of a
 * shell seeded at MAG azimuth φ arrives up to 11° away from φ — comparable to the 14.4° spacing of
 * the 25 default azimuths. Building the polar cap as `θ_fp(φ_seed)` rather than `θ_fp(φ_foot)`
 * distorts the boundary in a way that CORRELATES with `dθ/dφ`, so it does not average out. It cost
 * a systematic 1.7 % at L = 2 rising to 4.3 % at L = 8 — every value low, every value plausible —
 * until the contour was re-indexed onto the footpoints' own longitudes and the deviation from the
 * oracle fell to 0.002. The azimuths are therefore SORTED by footpoint longitude and integrated on
 * that non-uniform grid.
 *
 * **`k₀` is the epoch's dipole moment, not a constant.** @ref dipole_moment reads it from the same
 * IGRF coefficients the traces use. A stale moment shows up as a CONSTANT relative offset at every
 * shell, which reads like a units error and is not one; the 1960s-era `0.311653e5` is 4.3 % high
 * against a present-day epoch and buys a uniform 1.4 % error in every L\* the library returns.
 *
 * ## The quadrature is deliberately not IRBEM's
 *
 * `ERROR_BUDGET.md` §2(b) measures IRBEM's L\* converging NON-monotonically — `options = 4` is
 * repeatedly worse than `options = 2` — and identifies the cause as a first-order rectangle rule
 * over a fixed-step θ march whose last cell jitters as the root moves across a grid point. Here the
 * θ integral is a midpoint rule whose FINAL cell is truncated exactly at `θ_fp`, so the boundary
 * is not quantised to the grid, and the azimuthal integral is a trapezoid over the sorted
 * footpoint longitudes. Same cell count, second order instead of first, and no jitter.
 *
 * ## Measured
 *
 * Measured by [`tools/oracle/lstar_diff.cpp`](../../tools/oracle/lstar_diff.cpp) against
 * `libirbem-O2.so` (IRBEM `e7cecb0`, gfortran 13.3.0, `-O2 -ffp-contract=off -fno-fast-math`), run
 * as a black box, and by [`tests/irbem_driftshell_test.cpp`](../../tests/irbem_driftshell_test.cpp)
 * against the goldens that harness produced. IGRF internal field only (`kext = 0`),
 * 2015-180 12:00 UT, twelve geometries spanning L = 2…8.1 including two off the equator, plus the
 * eleven awkward ones the accuracy tables below add; `get_igrf_version_` reports 14, the generation
 * `tables/igrf14.hpp` carries, so both libraries evaluate the same field.
 *
 * **Accuracy.** Deviations in L\*, absolute:
 *
 * | matched `options(3,4)` | vs IRBEM at the same setting | vs IRBEM's converged (9,9) answer | IRBEM's own distance from (9,9) |
 * |---|---|---|---|
 * | 0 (`Nder` 25, `dθ` 0.25°) | 0.0354 max | **0.0099 max, 0.0044 mean** | 0.0268 max, 0.0109 mean |
 * | 9 (`Nder` 250, `dθ` 0.025°) | **0.0066 max, 0.0018 mean** | 0.0066 max | — |
 *
 * Read the first row carefully, because the largest number in it is the least interesting.
 * `ERROR_BUDGET.md` §2 measures IRBEM's default-resolution discretization at 0.010–0.017, and the
 * third column reproduces it. Two implementations each ~0.01 from converged, erring in opposite
 * directions, differ by ~0.02–0.035 — which says nothing about either. The column that means
 * something is the second: at IRBEM's own recommended resolution this implementation sits **2.5×
 * closer to the converged answer than IRBEM does** (0.0044 mean against 0.0109), which is the
 * second-order midpoint quadrature with an exactly-truncated final cell earning its keep against a
 * first-order rectangle rule on a grid the boundary jitters across. At matched `options = 9` both
 * are near converged and the agreement is 0.0066 on these twelve — inside the budget's 0.01, but
 * see the off-equatorial row below for where that stops being true and why.
 *
 * Coverage of the parameter space, all at matched `options = 9`:
 *
 * | sweep | range | worst deviation |
 * |---|---|---|
 * | pitch angle (shell splitting) | 15…90°, `I` from 0.002 to 11.6 Re | 0.0079 |
 * | epoch | 1900…2029, `k₀` moving 7.8 % | 0.0050 |
 * | **off-equator and off-eastern** | 4 southern, 3 at negative x/y, L = 1.5 and 9, one tilted | **0.0105** |
 * | fp32 device lane vs fp64 host lane | 32 shells | 7.6 × 10⁻⁶ (RTX 3070 Ti) |
 * | centred dipole, `L*` against `L_m` (closed form, no oracle) | L = 2…6.6 | 0.0038 |
 *
 * **The third row is outside the 0.01 budget, and the reason is the trace step rather than the
 * shell.** Ten of the twelve geometries in the table above are equatorial and the other two are
 * northern; that sample is where this implementation looks best, and it is not the domain. At
 * (0, 4, −1.5) Re the deviation is 0.0105. For a 90° particle at an OFF-EQUATORIAL start the
 * mirror point IS the start, so the integrand of `I` leaves its endpoint with a square-root
 * singularity that `lstar.hpp`'s fixed-step rule resolves slowly, and
 * @ref DriftShellOptions::trace defaults to its `steps_per_l = 50`. MEASURED at that geometry:
 *
 * | `trace.steps_per_l` | 50 (default) | 100 | 200 | 400 |
 * |---|---|---|---|---|
 * | \|ours − IRBEM\| | 0.0105 | 0.0075 | 0.0059 | 0.0053 |
 *
 * Monotone, converging to a NON-ZERO ~0.005 that is a real algorithmic difference from IRBEM.
 * One doubling of the trace step buys the budget off the equator and costs a factor of two in
 * trace time, which is why it is a field on the options struct and not the default. The footpoint
 * resolution is irrelevant here: 100 against 400 moves the answer by 5 × 10⁻⁶.
 *
 * **The diagnostic fields are looser than L\*, and by more than the budget.** `DriftShell` also
 * returns `lm`, `invariant_i` and `b_local` as IRBEM's `Lm`, `XJ` and `Blocal`. Over the same
 * eleven geometries, relative: `Blocal` 1.7 × 10⁻⁵ against a 10⁻⁶ budget (upstream, in `igrf.hpp`),
 * `Lm` 5.4 × 10⁻³ against 10⁻³, `XJ` 9.2 × 10⁻² against 10⁻⁴. All of the last two is the same
 * fixed trace step. L\* survives it because the root-find matches `I(r)` against an `I₀` computed
 * THE SAME WAY, so the two errors are correlated and largely divide out — which is why the shell
 * is good to 10⁻³ while the invariant it is built on is good to 10⁻². Stated here rather than
 * left to be discovered, and asserted at the measured level by
 * `IrbemDriftShell.HeavyDifferentialReportsTheDiagnosticFieldsIrbemDoesAndSaysHowCloselyTheyAgree`.
 *
 * **Cost**, milliseconds per L\* point at IRBEM's default resolution, same machine, same run,
 * ~220 traces per point either way (RTX 3070 Ti; host lane `-O3 -march=native -ffp-contract=off`):
 *
 * | | IRBEM `-O2` | this, host lane | this, device, batch 64 | batch 512 | batch 2048 |
 * |---|---|---|---|---|---|
 * | ms per L\* point | 14.3 | 12.8 | 1.45 | 0.435 | **0.356** |
 * | speedup vs IRBEM | 1.0× | 1.1× | 9.9× | 33× | **40×** |
 *
 * Re-measured independently on a later day, same machine: IRBEM 14.0–14.6, host lane 11.6–13.3,
 * batch 64 1.40–1.41, batch 512 0.434–0.441, batch 2048 0.346–0.350. The row reproduces.
 * One caveat that is easy to trip over: with `irbem_trace_i_f32.spv` missing from the shader
 * directory the same batch runs at 7.9 ms/point — the flux and footpoint kernels still fire, so
 * the seam looks alive while the stage that matters has fallen back. That is what
 * @ref make_lstar_batch's `Result<bool>` is narrowed to report.
 *
 * **A single point is 10…33 ms with a device present, worse than the oracle and NOISY**, and that
 * is the honest shape of the problem rather than a shortfall: `Nder = 25` shells is an order of
 * magnitude below the ~512-line crossover `gpu/dispatch.hpp` measures, so the traces run on the
 * host anyway while the flux stage still pays a dispatch. Repeated runs of the same call have
 * measured 10.6, 11.6, 32.9 and once 128 ms — the spread is submit-and-fence jitter on a dispatch
 * with nothing in it, and it is the reason no small-batch number here is quoted to three figures.
 * Forcing the host lane (`CHEATAH_SPACE_IRBEM_NO_GPU=1`) gives a flat 11.6…13.3 and is the better
 * choice for a single point. The parallelism lives across POINTS, and it is
 * @ref make_lstar_batch that takes it — which is why a 100 000-point ephemeris falls from IRBEM's
 * ~24 minutes to ~36 seconds.
 *
 * **The device answer is reproducible and portable.** Run to run on the same device the 32-point
 * cross-lane deviation repeats to every printed digit (7.64 × 10⁻⁶ on consecutive runs), and the
 * same batch on llvmpipe (9.9 × 10⁻⁶) and an Intel iGPU (2.04 × 10⁻⁵) lands within a factor of
 * three, so the agreement is the fp32 arithmetic rather than one driver's — and all three are two
 * to three orders inside the 0.01 L\* budget. The footpoint colatitude agrees to 0.0002° on all
 * three, which is the number the flux quadrature actually consumes.
 *
 * **Allocation.** Measured with `valgrind --tool=memcheck`: 117 allocations for one L\* point,
 * 147 for a batch of eight and 147 for a batch of SIXTY-FOUR — 0 leaks, 0 memcheck errors. The
 * count is O(rounds), never O(points), never per trace, and never per field evaluation.
 *
 * ## Sources
 *
 *  - Roederer, *Dynamics of Geomagnetically Trapped Radiation*, Springer (1970), ch. 2 — the third
 *    invariant `Φ`, the drift shell, and `L* = 2πk₀/(Φ R_E)`.
 *  - Schulz & Lanzerotti, *Particle Diffusion in the Radiation Belts*, Springer (1974), §1.3.
 *  - Brent, *Algorithms for Minimization without Derivatives*, Prentice-Hall (1973), ch. 4 — the
 *    bracketed root-find that replaces IRBEM's serial θ march.
 *  - IRBEM's published option table (`docs/source/api/general_information.rst`): `options(3)` sets
 *    `dθ = π/(720·(options(3)+1))` and `options(4)` sets `Nder = 25·(options(4)+1)`. Those two
 *    numbers are what @ref DriftShellOptions::from_irbem reproduces, so a differential comparison
 *    is comparing algorithms rather than resolutions.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <span>
#include <vector>

#include "coords_rotations.hpp"
#include "frames.hpp"
#include "igrf.hpp"
#include "lstar.hpp"
#include "status.hpp"

namespace cheatah::space::irbem {

/**
 * How finely a drift shell is resolved, and how hard the root-find works.
 *
 * The two resolutions carry IRBEM's meaning exactly, so @ref from_irbem can translate its
 * `options(3)` and `options(4)` and a differential comparison can be made at MATCHED settings —
 * which `docs/ERROR_BUDGET.md` §2(a) shows is the only setting at which a 0.01 L\* target means
 * anything, IRBEM's own default-resolution error being 0.010–0.017 at L ≈ 6.
 */
struct DriftShellOptions {
    /// How many azimuths the shell is sampled at — IRBEM's `Nder`. The default 25 is its
    /// `options(4) = 0`.
    int azimuths = 25;
    /// The colatitude step of the polar-cap quadrature, degrees — IRBEM's `dθ`. The default 0.25°
    /// is its `options(3) = 0`, i.e. `π/720`.
    double colatitude_step_deg = 0.25;
    /// How many radii the bracketing sweep tries, spread over
    /// `[bracket_low, bracket_high] × L_m`. These all go into ONE batch, which is what makes the
    /// bracketing phase `ntime × Nder × trials` wide rather than a serial scan.
    int bracket_trials = 6;
    /// The low end of the bracketing sweep, as a fraction of the dipole seed `L_m`.
    double bracket_low = 0.82;
    /// The high end of the bracketing sweep, as a fraction of the dipole seed `L_m`. Measured over
    /// L = 2…8 with IGRF-14 the true root sits within ±7 % of `L_m`; the sweep is widened to
    /// −18 %/+32 % so a disturbed shell still brackets.
    double bracket_high = 1.32;
    /// How many Brent iterations follow the sweep. Bounded on purpose: IRBEM's θ march has no
    /// bound, and an unbounded refinement inside a lock-step batch would make every shell wait for
    /// the worst one.
    int refine_iterations = 8;
    /// The radius the root-find is converged to, Earth radii. 5 × 10⁻⁴ Re is ~1/20 of a
    /// quadrature cell in the footpoint colatitude it feeds, so the root-find is not the
    /// error term.
    double radius_tolerance = 5.0e-4;
    /// The scale of the analytic continuation below the shell, Earth radii. Positive; it sets how
    /// fast an out-of-range bracket closes and cannot move the root. See the file brief.
    double inaccessible_scale = 1.0;
    /// Tracing options for the invariant `I` — handed straight to @ref trace_invariant_batch.
    TraceOptions trace{};
    /// Tracing options for the footpoint walk. Longer than @ref trace — it runs from the magnetic
    /// equator all the way to `r = 1` rather than stopping at the mirror points — and finer,
    /// because its colatitude is what the flux quadrature's boundary is. MEASURED at IRBEM's
    /// default resolution, worst deviation from the converged oracle over 12 points:
    /// `steps_per_l` = 25 → 0.0121, 50 → 0.0101, 100 → 0.0100, 200 → 0.0099. The knee is at 50 and
    /// the default sits one doubling past it, which costs ~10 % of a host-lane L\* and nothing
    /// measurable on the device.
    TraceOptions footpoint{100.0, 8000, 1.0};

    /**
     * The settings IRBEM's `options(3)` and `options(4)` select.
     *
     * `dθ = π/(720·(t_resolution+1))` and `Nder = 25·(r_resolution+1)`, both from its published
     * option table. Used by the differential suite so the two libraries are compared at the same
     * resolution rather than across resolutions — the distinction `ERROR_BUDGET.md` §2(a) insists
     * on.
     *
     * @param t_resolution IRBEM's `options(3)`, 0–9; the colatitude resolution.
     * @param r_resolution IRBEM's `options(4)`, 0–9; the azimuth count.
     * @return the options, with everything else left at its default.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemDriftShell.IrbemOptionsTranslateToTheDocumentedResolutions
     */
    [[nodiscard]] static DriftShellOptions from_irbem(int t_resolution, int r_resolution) {
        DriftShellOptions o;
        o.azimuths = 25 * (r_resolution + 1);
        o.colatitude_step_deg = 180.0 / (720.0 * (t_resolution + 1));
        return o;
    }
};

/// Everything one L\* evaluation produces. Fixed size, trivially copyable — the drift shell itself
/// is not returned, because a caller computing an ephemeris wants six scalars per point and not
/// `Nder × steps` positions.
struct DriftShell {
    double lstar = 0.0;        ///< Roederer's L\*, Earth radii.
    double phi = 0.0;          ///< The third invariant `Φ`, nT·Re²; always positive.
    double lm = 0.0;           ///< McIlwain's L of the starting line, Earth radii.
    double invariant_i = 0.0;  ///< `I` of the starting line, Earth radii — IRBEM's `XJ`.
    double b_mirror = 0.0;     ///< The mirror field the shell conserves, nT.
    double b_min = 0.0;        ///< `B_min` on the starting line, nT.
    double b_local = 0.0;      ///< `|B|` at the starting point, nT.
    int azimuths = 0;          ///< How many azimuths converged; equals `DriftShellOptions::azimuths`
                               ///< when the status is @ref Status::Ok.
    int traces = 0;            ///< Field lines traced for this point — the cost, reported rather
                               ///< than assumed.
};

namespace detail {

/// Two π, the full turn the azimuths sweep.
inline constexpr double kTwoPi = 2.0 * std::numbers::pi;

/**
 * One bracketed root-find, held as data so `ntime × Nder` of them can advance in LOCK STEP.
 *
 * This is the whole reason the shell parallelises. A root-find written as a loop that calls a
 * function owns its own control flow, and thousands of them cannot share a device dispatch. Split
 * into "state plus one step", every shell proposes its next radius, all the proposals go into one
 * batch, and the batch comes back and advances every state. The algorithm is unchanged; only who
 * owns the loop is.
 *
 * The method is Brent's (1973, ch. 4): inverse quadratic interpolation through the three most
 * recent points, falling back to a secant when only two are distinct, and to bisection whenever
 * the interpolated step would leave the bracket or fail to shrink it fast enough. The bracket is
 * never given up, so a pathological residual costs iterations and never divergence.
 */
struct RootState {
    double a = 0.0;   ///< The PREVIOUS best estimate — the third point of the interpolation.
                      ///< (At @ref root_begin, before the first step has rotated anything, it is
                      ///< instead the low end of the bracket; @ref root_step's first action moves
                      ///< the contra-point out of it and into @ref c.)
    double b = 0.0;   ///< The current best estimate, and the radius to probe next.
    double c = 0.0;   ///< The CONTRA-POINT: `f(c)` has the opposite sign to `f(b)`, so the root is
                      ///< between them and `c - b` is the bracket that has to close.
    double fa = 0.0;  ///< The residual at @ref a.
    double fb = 0.0;  ///< The residual at @ref b.
    double fc = 0.0;  ///< The residual at @ref c.
    double step = 0.0;       ///< The step just taken.
    double prev_step = 0.0;  ///< The step before it — Brent's test for "is interpolation working".
    bool done = false;       ///< Whether the bracket has closed to the requested tolerance.
    bool valid = false;      ///< Whether a bracket was found at all. A shell with no bracket is
                             ///< @ref Status::NotConverged, never a silently bisected guess.
};

/**
 * Begin a root-find over a bracket that is already known to straddle the root.
 *
 * @param lo the low end. @param f_lo its residual, which must be `≤ 0`.
 * @param hi the high end. @param f_hi its residual, which must be `> 0`.
 * @return the state, with @ref RootState::b holding the first radius to probe.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemDriftShell.BrentFindsAKnownRootOfAnAnalyticFunction
 */
[[nodiscard]] inline RootState root_begin(double lo, double f_lo, double hi, double f_hi) {
    RootState s;
    s.a = lo;
    s.fa = f_lo;
    s.b = hi;
    s.fb = f_hi;
    s.c = hi;
    s.fc = f_hi;
    s.valid = true;
    return s;
}

/**
 * Absorb the residual at @ref RootState::b and propose the next probe.
 *
 * Brent's iteration, written as one transition so a batch of root-finds can share a dispatch. The
 * order of the three cases is the method: re-bracket if the new point landed on the same side as
 * the contra-point, rotate so `b` is the best estimate, then either interpolate or bisect.
 *
 * @param s the state, advanced in place. @ref RootState::done is set when the bracket has closed.
 * @param f_at_b the residual at the radius @ref RootState::b that was just probed.
 * @param tol the bracket width to converge to, in the same units as the abscissa.
 * @complexity O(1) — a fixed handful of arithmetic operations and no evaluation of the residual.
 * @alloc none.
 * @test IrbemDriftShell.BrentFindsAKnownRootOfAnAnalyticFunction
 * @test IrbemDriftShell.BrentBisectsWhenInterpolationWouldLeaveTheBracket
 */
inline void root_step(RootState& s, double f_at_b, double tol) {
    s.fb = f_at_b;
    if ((s.fb > 0.0) == (s.fc > 0.0)) {
        // b and c fell on the same side of the root, so c is no longer a contra-point; a — which
        // holds the previous iterate, or the low end of the bracket on the first step — is.
        s.c = s.a;
        s.fc = s.fa;
        s.step = s.b - s.a;
        s.prev_step = s.step;
    }
    if (std::abs(s.fc) < std::abs(s.fb)) {
        // c is the better estimate: rotate it into b and push the old b out to both a and c, which
        // re-establishes "b is best, c is the contra-point" with the bracket unchanged.
        const double t = s.b;
        s.a = t;
        s.b = s.c;
        s.c = t;
        const double ft = s.fb;
        s.fa = ft;
        s.fb = s.fc;
        s.fc = ft;
    }
    // The floor under the step: a bracket cannot be closed below the abscissa's own resolution.
    const double band = (2.0 * 2.22e-16 * std::abs(s.b)) + (0.5 * tol);
    const double half = 0.5 * (s.c - s.b);
    if (std::abs(half) <= band || s.fb == 0.0) {
        s.done = true;
        return;
    }
    if (std::abs(s.prev_step) >= band && std::abs(s.fa) > std::abs(s.fb)) {
        const double r1 = s.fb / s.fa;
        double p = 0.0;
        double q = 0.0;
        if (s.a == s.c) {
            // Only two distinct points: the secant.
            p = 2.0 * half * r1;
            q = 1.0 - r1;
        } else {
            // Three distinct points: inverse quadratic interpolation.
            const double r2 = s.fa / s.fc;
            const double r3 = s.fb / s.fc;
            p = r1 * ((2.0 * half * r2 * (r2 - r3)) - ((s.b - s.a) * (r3 - 1.0)));
            q = (r2 - 1.0) * (r3 - 1.0) * (r1 - 1.0);
        }
        if (p > 0.0) q = -q;
        p = std::abs(p);
        // Accept the interpolation only if it stays inside the bracket AND shrinks at least as
        // fast as the step before it; otherwise bisect. This is what bounds the iteration count.
        const double limit =
            std::min((3.0 * half * q) - std::abs(band * q), std::abs(s.prev_step * q));
        if (2.0 * p < limit) {
            s.prev_step = s.step;
            s.step = p / q;
        } else {
            s.step = half;
            s.prev_step = half;
        }
    } else {
        s.step = half;
        s.prev_step = half;
    }
    s.a = s.b;
    s.fa = s.fb;
    s.b += (std::abs(s.step) > band) ? s.step : std::copysign(band, half);
}

/**
 * Stage a model's Gauss coefficients and Legendre normalisation as the `float` buffers the kernels
 * bind.
 *
 * The interpolation to the epoch happens HERE, once per batch, exactly as `lstar.hpp`'s trace lane
 * does it: the device never sees IGRF's 26-epoch table, because interpolating it per thread would
 * be one redundant copy per point of a calculation the host does once.
 *
 * @tparam NMAX the truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param coef receives `g[slots]` then `h[slots]`.
 * @param norm receives `e[slots]`, `f[slots]`, `diagonal[NMAX+1]`.
 * @complexity O(NMAX²).
 * @alloc none — both spans are the caller's.
 * @test IrbemDriftShell.FluxCellsAgreeBetweenLanes
 */
template <GeoFieldModel M>
inline void stage_model(const M& model, std::span<float> coef, std::span<float> norm) {
    constexpr int NMAX = M::degree;
    constexpr std::size_t kSlots = ((NMAX + 1) * (NMAX + 2)) / 2;
    std::fill(coef.begin(), coef.end(), 0.0F);
    for (int deg = 1; deg <= NMAX; ++deg) {
        for (int m = 0; m <= deg; ++m) {
            const std::size_t k = (static_cast<std::size_t>(deg) * (deg + 1)) / 2 + m;
            coef[k] = static_cast<float>(model.g(deg, m));
            coef[kSlots + k] = static_cast<float>(model.h(deg, m));
        }
    }
    constexpr auto kNorm = make_legendre_normalisation<NMAX, double>();
    for (std::size_t k = 0; k < kSlots; ++k) {
        norm[k] = static_cast<float>(kNorm.e[k]);
        norm[kSlots + k] = static_cast<float>(kNorm.f[k]);
    }
    for (std::size_t deg = 0; deg <= static_cast<std::size_t>(NMAX); ++deg) {
        norm[(2 * kSlots) + deg] = static_cast<float>(kNorm.diagonal[deg]);
    }
}

/// How many `float` slots @ref stage_model's coefficient buffer needs at degree @p nmax.
/// @param nmax the truncation degree.
/// @return `2·(nmax+1)(nmax+2)/2`.
/// @complexity O(1). @alloc none.
/// @test IrbemDriftShell.FluxCellsAgreeBetweenLanes
[[nodiscard]] inline constexpr std::size_t coefficient_slots(int nmax) {
    return static_cast<std::size_t>((nmax + 1) * (nmax + 2));
}

/// How many `float` slots @ref stage_model's normalisation buffer needs at degree @p nmax.
/// @param nmax the truncation degree.
/// @return `2·(nmax+1)(nmax+2)/2 + nmax + 1`.
/// @complexity O(1). @alloc none.
/// @test IrbemDriftShell.FluxCellsAgreeBetweenLanes
[[nodiscard]] inline constexpr std::size_t normalisation_slots(int nmax) {
    return coefficient_slots(nmax) + static_cast<std::size_t>(nmax) + 1;
}

/**
 * The internal field at a whole batch of GEO points — the device lane where it pays, the host lane
 * where it does not.
 *
 * The one field-evaluation primitive this header needs: the root-find wants `|B|` at each candidate
 * seed to turn a mirror field into a pitch angle, the footpoint stage wants the field VECTOR to
 * pick a hemisphere, and the flux quadrature wants the radial component at every polar-cap cell.
 * All three are the same batched evaluation, so there is one of it.
 *
 * @tparam NMAX the truncation degree.
 * @param model the internal field model. @param coef the staged coefficients.
 * @param norm the staged normalisation.
 * @param pos the points, xyz-interleaved, `3N` floats, GEO, Earth radii.
 * @param out receives the field, xyz-interleaved, `3N` floats, nT.
 * @return `true` when the device serviced the call.
 * @complexity O(N) field evaluations; concurrent on the device.
 * @alloc none here; the device lane's buffers come from the context's pool.
 * @test IrbemDriftShell.FluxCellsAgreeBetweenLanes
 */
template <GeoFieldModel M>
[[nodiscard]] inline bool field_batch(const M& model, std::span<const float> coef,
                                      std::span<const float> norm, std::span<const float> pos,
                                      std::span<float> out) {
    const std::size_t n = pos.size() / 3;
#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
    // The device fast-path stages ONE internal model's coefficients, so it applies only when the
    // model IS that internal model. A composed model (internal plus external) takes the host loop
    // below, whose `model.evaluate` is the total field — which is also the correct physics for the
    // flux quadrature: IRBEM's own cap integral evaluates the total field, and evaluating only the
    // internal part there would bias Phi by the external field's ~0.1-0.3% at r = 1, comparable to
    // the whole L* budget during a storm.
    if constexpr (is_igrf_v<M>) {
        if (gpu::prefer_gpu("irbem_igrf_f32", n)) {
            const std::array<std::uint32_t, 2> dims{static_cast<std::uint32_t>(n),
                                                    static_cast<std::uint32_t>(M::degree)};
            if (gpu::launch_igrf(pos, coef, norm, dims, out)) return true;
        }
    }
#endif
    (void)coef;
    (void)norm;
    for (std::size_t i = 0; i < n; ++i) {
        const Position<Frame::GEO> p{fixarray::vec3d{pos[(3 * i) + 0], pos[(3 * i) + 1],
                                                     pos[(3 * i) + 2]}};
        const fixarray::vec3d b = model.evaluate(p).v;
        out[(3 * i) + 0] = static_cast<float>(b[0]);
        out[(3 * i) + 1] = static_cast<float>(b[1]);
        out[(3 * i) + 2] = static_cast<float>(b[2]);
    }
    return false;
}

/// Where the shell's field line crosses the magnetic-equatorial plane, in GEO.
/// @param mag_to_geo the MAG→GEO rotation for the epoch.
/// @param azimuth_rad the MAG azimuth. @param radius the crossing radius, Earth radii.
/// @return the point in GEO.
/// @complexity O(1) — two trigonometric evaluations and one matrix-vector product.
/// @alloc none.
/// @test IrbemDriftShell.SeedsLieInTheMagneticEquatorialPlane
[[nodiscard]] inline Position<Frame::GEO> equatorial_seed(const fixarray::mat3d& mag_to_geo,
                                                          double azimuth_rad, double radius) {
    const fixarray::vec3d m{radius * std::cos(azimuth_rad), radius * std::sin(azimuth_rad), 0.0};
    return Position<Frame::GEO>{mag_to_geo * m};
}

/**
 * Walk a field line from @p start down to `r = 1` — the fp64 reference the `irbem_shell_foot_f32`
 * kernel is verified against, and the lane a machine without a device runs.
 *
 * The final step is HALVED onto the sphere rather than interpolated across it. Interpolating the
 * crossing linearly along the last chord costs O(ds²) in the footpoint's colatitude, which at
 * `ds = L/50` is a tenth of a degree — an order worse than the 0.25° quadrature cell it feeds, and
 * a bias rather than noise. Halving costs ~2 steps per level and leaves the RK4 truncation as the
 * only error.
 *
 * @tparam NMAX the truncation degree.
 * @param model the internal field model.
 * @param start the magnetic-equator crossing, GEO, Earth radii.
 * @param northward whether to walk along `+B` (true) or `−B`.
 * @param opt the tracing options; @ref TraceOptions::steps_per_l and
 *        @ref TraceOptions::max_steps are read.
 * @return the footpoint ON the unit sphere, with @ref Status::OpenFieldLine when the walk never
 *         reached the surface within @ref TraceOptions::max_steps and @ref Status::DomainError for
 *         a start inside the atmosphere or in a null field.
 * @complexity O(steps) field evaluations, ~4 per step; ~250 steps typically.
 * @alloc none — the walk stores no path.
 * @test IrbemDriftShell.FootpointsAgreeBetweenLanes
 * @test IrbemDriftShell.FootpointOfADipoleLineIsTheAnalyticColatitude
 */
template <GeoFieldModel M>
[[nodiscard]] inline Result<Position<Frame::GEO>> walk_to_surface(const M& model,
                                                                  const Position<Frame::GEO>& start,
                                                                  bool northward,
                                                                  const TraceOptions& opt) {
    /// The step below which the walk counts as having landed, Earth radii — ~64 m, two orders
    /// below the finest quadrature cell and far above fp32's resolution near 1.0. The kernel
    /// carries the same constant, so the two lanes stop in the same place.
    constexpr double kFootTol = 1.0e-5;
    const double r0 = fixarray::norm(start.v);
    if (!(r0 > opt.min_radius)) return {Status::DomainError, start};
    fixarray::vec3d b = model.evaluate(start).v;
    if (!(fixarray::norm(b) > 0.0)) return {Status::DomainError, start};

    double ds = (northward ? 1.0 : -1.0) * dipole_l(start) / opt.steps_per_l;
    Position<Frame::GEO> p = start;
    bool landed = false;
    for (int k = 0; k < opt.max_steps; ++k) {
        fixarray::vec3d b_new{};
        const Position<Frame::GEO> q = rk4_step(model, p, b, ds, b_new);
        if (fixarray::norm(q.v) < 1.0) {
            if (std::abs(ds) < kFootTol) {
                p = q;
                landed = true;
                break;
            }
            ds *= 0.5;
            continue;
        }
        p = q;
        b = b_new;
    }
    const double rp = fixarray::norm(p.v);
    const Position<Frame::GEO> foot{rp > 0.0 ? p.v * (1.0 / rp) : fixarray::vec3d{0.0, 0.0, 1.0}};
    return {landed ? Status::Ok : Status::OpenFieldLine, foot};
}

}  // namespace detail

namespace detail {

/// The largest number of polar-cap cells evaluated in one dispatch. Bounds the flux stage's
/// staging memory at ~24 MB no matter how many points the batch holds, and cannot change the
/// answer: cells are produced and consumed in ascending index order, so the fp64 accumulation
/// sequence is the same whatever the chunk size.
inline constexpr std::size_t kFluxChunk = 1U << 20U;

/**
 * Evaluate the drift-shell residual `f(r)` for a whole batch of candidate radii — the step that
 * makes the root-find parallel.
 *
 * Every candidate that needs a trace goes into ONE call to @ref trace_invariant_batch. Candidates
 * that do not need one — the radii where `|B|` on the magnetic-equatorial plane already exceeds
 * the mirror field, so no particle mirroring at `B_m` reaches them — are answered analytically and
 * never enter the batch. That is not an optimisation for its own sake: it is what keeps the
 * residual defined and continuous below the shell, which is what lets Brent keep a bracket when
 * `I₀ ≈ 0`. See the file brief.
 *
 * An open field line is reported as `−I₀ − w`: a line that reaches the atmosphere before it
 * reaches its mirror point means the candidate shell is too SMALL, so a negative residual pushes
 * the bracket outward. That is IRBEM's "drift shell intersects the atmosphere" case, and here it
 * steers the search rather than aborting it.
 *
 * @tparam NMAX the truncation degree.
 * @param model the internal field model. @param mag_to_geo the epoch's MAG→GEO rotation.
 * @param coef the staged coefficients. @param norm the staged normalisation.
 * @param azimuth each candidate's MAG azimuth, radians.
 * @param b_mirror each candidate's conserved mirror field, nT.
 * @param invariant0 each candidate's target invariant `I₀`, Earth radii.
 * @param radius the radii to probe, Earth radii — same length as the three above.
 * @param residual receives `f(r)`, same length.
 * @param opt the drift-shell options; @ref DriftShellOptions::trace and
 *        @ref DriftShellOptions::inaccessible_scale are read.
 * @param traces incremented by however many field lines were actually traced.
 * @return `true` when the device serviced the trace batch.
 * @complexity O(candidates × steps) field evaluations, concurrent on the device.
 * @alloc a fixed number of vectors per ROUND — never per candidate and never per trace.
 * @test IrbemDriftShell.ResidualIsMonotoneAcrossTheAccessibilityBoundary
 * @test IrbemDriftShell.ResidualIsMonotoneInTheShellRadius
 */
template <GeoFieldModel M>
[[nodiscard]] inline bool residual_round(const M& model, const fixarray::mat3d& mag_to_geo,
                                         std::span<const float> coef, std::span<const float> norm,
                                         std::span<const double> azimuth,
                                         std::span<const double> b_mirror,
                                         std::span<const double> invariant0,
                                         std::span<const double> radius, std::span<double> residual,
                                         const DriftShellOptions& opt, int& traces) {
    const std::size_t m = radius.size();
    std::vector<float> pos(3 * m);
    std::vector<float> fld(3 * m);
    std::vector<Position<Frame::GEO>> seed(m);
    for (std::size_t j = 0; j < m; ++j) {
        seed[j] = equatorial_seed(mag_to_geo, azimuth[j], radius[j]);
        pos[(3 * j) + 0] = static_cast<float>(seed[j].v[0]);
        pos[(3 * j) + 1] = static_cast<float>(seed[j].v[1]);
        pos[(3 * j) + 2] = static_cast<float>(seed[j].v[2]);
    }
    (void)field_batch(model, coef, norm, pos, fld);

    std::vector<std::size_t> work;
    std::vector<Position<Frame::GEO>> starts;
    std::vector<double> pitch;
    work.reserve(m);
    starts.reserve(m);
    pitch.reserve(m);
    const double w = opt.inaccessible_scale;
    for (std::size_t j = 0; j < m; ++j) {
        const double bp = std::hypot(std::hypot(static_cast<double>(fld[(3 * j) + 0]),
                                                static_cast<double>(fld[(3 * j) + 1])),
                                     static_cast<double>(fld[(3 * j) + 2]));
        if (!(bp < b_mirror[j])) {
            // Below the shell: no particle mirroring at B_m reaches this radius. The analytic
            // continuation keeps f strictly increasing through the degenerate I₀ = 0 root.
            residual[j] = -invariant0[j] - (w * ((bp / b_mirror[j]) - 1.0));
            continue;
        }
        work.push_back(j);
        starts.push_back(seed[j]);
        // sin α = √(B/B_m) makes the trace's own B_local/sin²α exactly B_m, so the mirror field is
        // conserved by construction and the root-find is one equation in one unknown.
        pitch.push_back(std::asin(std::sqrt(bp / b_mirror[j])) * (180.0 / std::numbers::pi));
    }

    std::vector<FieldLine> lines(work.size());
    std::vector<Status> st(work.size());
    const Result<bool> r =
        trace_invariant_batch(model, starts, pitch, lines, st, opt.trace);
    traces += static_cast<int>(work.size());
    for (std::size_t q = 0; q < work.size(); ++q) {
        const std::size_t j = work[q];
        residual[j] = st[q] == Status::Ok ? lines[q].invariant_i - invariant0[j]
                                          : -invariant0[j] - w;
    }
    return r.value;
}

}  // namespace detail

/**
 * Roederer's L\* for a whole batch of points — **this is the routine to call**.
 *
 * The batch form is not a convenience wrapper. It is the only shape in which the problem is
 * parallel: one L\* point is `Nder` independent root-finds of a few traces each, and a device
 * dispatch does not pay for `Nder = 25` lines (`gpu/dispatch.hpp` measures the crossover at ~512).
 * Hand over `ntime` points at once and every stage becomes `ntime × Nder` wide, which is where the
 * measured 48.9× on the trace kernel actually lands. A loop calling @ref make_lstar per point
 * cannot be accelerated — the same fact `lstar.hpp` states about `trace_invariant`.
 *
 * Each point is treated as a LOCALLY MIRRORING particle at @p pitch_angles_deg, matching IRBEM's
 * `make_lstar` (90°) and `make_lstar_shell_splitting` (arbitrary). The shell that particle drifts
 * on is the set of field lines carrying its `(B_m, I)`, and `Φ` is the flux through the polar cap
 * those lines' footpoints enclose.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param rotations the epoch's rotations — only `geo_to_mag` is read, so the drift shell is
 *        organised about the DIPOLE axis rather than the geographic one.
 * @param starts the points, GEO, Earth radii.
 * @param pitch_angles_deg the local pitch angle at each point; same length as @p starts.
 * @param out receives one @ref DriftShell per point; same length as @p starts.
 * @param statuses receives each point's status; same length as @p starts.
 * @param opt the resolution and root-find settings.
 * @return @ref Status::Ok when every point closed; @ref Status::NotConverged when any shell failed
 *         to bracket or any footpoint failed to reach the surface; @ref Status::OpenFieldLine when
 *         a starting line did not close; @ref Status::DomainError on a length mismatch or a
 *         non-physical model. The value is `true` when the device serviced THE TRACES, and it is
 *         deliberately not an "any stage ran on the device" flag: the footpoint walk and the flux
 *         quadrature are folded out of it, because both take the device lane at batch sizes far
 *         below the ~512-line trace crossover and an OR across all three would read `true` on a
 *         call whose traces — 95 % of the cost, and the whole performance claim — ran on the host.
 *         MEASURED: with `irbem_trace_i_f32.spv` removed from the shader directory and every other
 *         kernel present, the OR'd flag still read `true` while the batch ran at 7.9 ms/point
 *         instead of 0.35, i.e. it reported success for exactly the silent fallback it exists to
 *         catch. `IrbemDriftShell.UsesTheDeviceWhenOneIsAvailable` asserts this narrower flag,
 *         so that fallback now fails the suite.
 * @complexity O(points × Nder × (trials + iterations)) traces, plus O(points × Nder × 180/dθ)
 *             field evaluations for the flux. All of it concurrent on the device.
 * @alloc O(rounds) vectors for the root-find and O(1) per flux chunk; nothing per trace and
 *        nothing per field evaluation.
 * @test IrbemDriftShell.BatchMatchesThePointAtATimeCall
 * @test IrbemDriftShell.MatchesTheOracleAtIrbemDefaultResolution
 * @test IrbemDriftShell.UsesTheDeviceWhenOneIsAvailable
 */
template <GeoFieldModel M>
[[nodiscard]] inline Result<bool> make_lstar_batch(const M& model,
                                                   const Rotations& rotations,
                                                   std::span<const Position<Frame::GEO>> starts,
                                                   std::span<const double> pitch_angles_deg,
                                                   std::span<DriftShell> out,
                                                   std::span<Status> statuses,
                                                   const DriftShellOptions& opt = {}) {
    const std::size_t n = starts.size();
    if (pitch_angles_deg.size() != n || out.size() != n || statuses.size() != n) {
        return {Status::DomainError, false};
    }
    if (n == 0) return {Status::Ok, false};
    const double k0 = dipole_moment(model);
    const int nd = opt.azimuths;
    if (nd < 3 || opt.bracket_trials < 2 || !(opt.colatitude_step_deg > 0.0) || !(k0 > 0.0)) {
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = DriftShell{};
            statuses[i] = Status::DomainError;
        }
        return {Status::DomainError, false};
    }
    const std::size_t m = n * static_cast<std::size_t>(nd);
    const fixarray::mat3d geo_to_mag = rotation_matrix<Frame::MAG, Frame::GEO>(rotations);
    const fixarray::mat3d mag_to_geo = rotation_matrix<Frame::GEO, Frame::MAG>(rotations);

    std::vector<float> coef(detail::coefficient_slots(M::degree));
    std::vector<float> norm(detail::normalisation_slots(M::degree));
    detail::stage_model(model, coef, norm);

    // ---- 1. the reference lines: one trace per point, in one batch ---------------------------
    std::vector<FieldLine> ref(n);
    std::vector<Status> ref_status(n);
    bool device = trace_invariant_batch(model, starts, pitch_angles_deg, ref, ref_status,
                                        opt.trace)
                      .value;
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = DriftShell{};
        out[i].b_local = ref[i].b_local;
        out[i].b_min = ref[i].b_min;
        out[i].b_mirror = ref[i].b_mirror;
        out[i].invariant_i = ref[i].invariant_i;
        out[i].lm = mcilwain_l(ref[i].invariant_i, ref[i].b_mirror, k0).value;
        statuses[i] = ref_status[i];
    }

    // ---- 2. the per-azimuth unknowns, flattened so every shell of every point is one row ------
    int traced = static_cast<int>(n);   // the reference lines, one per point
    std::vector<double> azimuth(m);
    std::vector<double> bm(m);
    std::vector<double> i0(m);
    std::vector<double> radius(m);
    std::vector<double> resid(m);
    std::vector<detail::RootState> root(m);
    std::vector<char> live(m, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const Position<Frame::MAG> pm{geo_to_mag * starts[i].v};
        const double phi0 = std::atan2(pm.v[1], pm.v[0]);
        for (int k = 0; k < nd; ++k) {
            const std::size_t j = (i * static_cast<std::size_t>(nd)) + static_cast<std::size_t>(k);
            azimuth[j] = phi0 + (detail::kTwoPi * k / nd);
            bm[j] = ref[i].b_mirror;
            i0[j] = ref[i].invariant_i;
            live[j] = ref_status[i] == Status::Ok && out[i].lm > 1.0 ? 1 : 0;
        }
    }

    // ---- 3. bracketing: every azimuth of every point, every trial radius, ONE batch -----------
    // The independent seed is what breaks IRBEM's serial azimuth chain. `L_m` alone lands within
    // ~7% of the root over L = 2..8; the sweep is wider than that so a disturbed shell still
    // brackets, and the whole sweep is a single dispatch rather than a scan that stops early.
    const int trials = opt.bracket_trials;
    const std::size_t sweep = m * static_cast<std::size_t>(trials);
    std::vector<double> sweep_az(sweep);
    std::vector<double> sweep_bm(sweep);
    std::vector<double> sweep_i0(sweep);
    std::vector<double> sweep_r(sweep);
    std::vector<double> sweep_f(sweep);
    for (std::size_t j = 0; j < m; ++j) {
        const double seed_l = out[j / static_cast<std::size_t>(nd)].lm;
        for (int t = 0; t < trials; ++t) {
            const std::size_t q = (j * static_cast<std::size_t>(trials)) + static_cast<std::size_t>(t);
            const double frac = opt.bracket_low + ((opt.bracket_high - opt.bracket_low) * t /
                                                   (trials - 1));
            sweep_az[q] = azimuth[j];
            sweep_bm[q] = bm[j];
            sweep_i0[q] = i0[j];
            sweep_r[q] = std::max(1.05, seed_l * frac);
        }
    }
    device = detail::residual_round(model, mag_to_geo, coef, norm, sweep_az, sweep_bm, sweep_i0,
                                    sweep_r, sweep_f, opt, traced) ||
             device;
    for (std::size_t j = 0; j < m; ++j) {
        if (live[j] == 0) continue;
        const std::size_t base = j * static_cast<std::size_t>(trials);
        bool bracketed = false;
        for (int t = 0; t + 1 < trials && !bracketed; ++t) {
            if (sweep_f[base + t] <= 0.0 && sweep_f[base + t + 1] > 0.0) {
                root[j] = detail::root_begin(sweep_r[base + t], sweep_f[base + t],
                                             sweep_r[base + t + 1], sweep_f[base + t + 1]);
                bracketed = true;
            }
        }
        if (!bracketed) live[j] = 0;
    }

    // ---- 4. refinement: all live shells advance in lock step, one batch per iteration ---------
    std::vector<std::size_t> active;
    active.reserve(m);
    for (int it = 0; it < opt.refine_iterations; ++it) {
        active.clear();
        for (std::size_t j = 0; j < m; ++j) {
            if (live[j] != 0 && !root[j].done) active.push_back(j);
        }
        if (active.empty()) break;
        std::vector<double> az(active.size());
        std::vector<double> mb(active.size());
        std::vector<double> ii(active.size());
        std::vector<double> rr(active.size());
        std::vector<double> ff(active.size());
        for (std::size_t q = 0; q < active.size(); ++q) {
            az[q] = azimuth[active[q]];
            mb[q] = bm[active[q]];
            ii[q] = i0[active[q]];
            rr[q] = root[active[q]].b;
        }
        device = detail::residual_round(model, mag_to_geo, coef, norm, az, mb, ii, rr, ff, opt,
                                        traced) ||
                 device;
        for (std::size_t q = 0; q < active.size(); ++q) {
            detail::root_step(root[active[q]], ff[q], opt.radius_tolerance);
        }
    }
    for (std::size_t j = 0; j < m; ++j) radius[j] = root[j].b;

    // ---- 5. footpoints: where each converged line meets r = 1 ---------------------------------
    std::vector<float> fpos(3 * m);
    std::vector<float> ffld(3 * m);
    std::vector<float> fdir(m);
    std::vector<float> foot(3 * m);
    std::vector<std::uint32_t> fstat(m, static_cast<std::uint32_t>(Status::DomainError));
    for (std::size_t j = 0; j < m; ++j) {
        const Position<Frame::GEO> s = detail::equatorial_seed(mag_to_geo, azimuth[j], radius[j]);
        fpos[(3 * j) + 0] = static_cast<float>(s.v[0]);
        fpos[(3 * j) + 1] = static_cast<float>(s.v[1]);
        fpos[(3 * j) + 2] = static_cast<float>(s.v[2]);
    }
    (void)detail::field_batch(model, coef, norm, fpos, ffld);
    for (std::size_t j = 0; j < m; ++j) {
        // Which way along B leads to the NORTHERN cap: the kernel works in GEO and has no dipole
        // axis, so the hemisphere choice is made once, here, where the axis is known.
        const double along = (static_cast<double>(ffld[(3 * j) + 0]) * rotations.dipole_geo[0]) +
                             (static_cast<double>(ffld[(3 * j) + 1]) * rotations.dipole_geo[1]) +
                             (static_cast<double>(ffld[(3 * j) + 2]) * rotations.dipole_geo[2]);
        fdir[j] = along >= 0.0 ? 1.0F : -1.0F;
    }
    bool foot_on_device = false;
#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
    // Igrf-only: irbem_shell_foot_f32 walks the INTERNAL field it was staged with. A total-field
    // footpoint must walk the total field — the external part bends the last few steps near the
    // cusp — so a composed model takes the generic host walk below, which evaluates the model it
    // was given. When a combined footpoint kernel exists this guard widens; until then routing a
    // total field through the internal kernel would be fast and wrong.
    if constexpr (is_igrf_v<M>) {
        if (gpu::prefer_gpu("irbem_shell_foot_f32", m)) {
            const std::array<std::uint32_t, 4> dims{
                static_cast<std::uint32_t>(m), static_cast<std::uint32_t>(M::degree),
                static_cast<std::uint32_t>(opt.footpoint.max_steps),
                static_cast<std::uint32_t>(opt.footpoint.steps_per_l * 1000.0)};
            foot_on_device = gpu::launch_shell_foot(fpos, fdir, coef, norm, dims, foot, fstat);
        }
    }
#endif
    if (!foot_on_device) {
        for (std::size_t j = 0; j < m; ++j) {
            const Position<Frame::GEO> s = detail::equatorial_seed(mag_to_geo, azimuth[j],
                                                                   radius[j]);
            const Result<Position<Frame::GEO>> f =
                detail::walk_to_surface(model, s, fdir[j] > 0.0F, opt.footpoint);
            foot[(3 * j) + 0] = static_cast<float>(f.value.v[0]);
            foot[(3 * j) + 1] = static_cast<float>(f.value.v[1]);
            foot[(3 * j) + 2] = static_cast<float>(f.value.v[2]);
            fstat[j] = static_cast<std::uint32_t>(f.status);
        }
    }

    // ---- 6. the ragged polar cap, and the ordered fp64 reduction ------------------------------
    const double dth = opt.colatitude_step_deg * (std::numbers::pi / 180.0);
    std::vector<double> theta(m, 0.0);
    std::vector<double> phi_f(m, 0.0);
    for (std::size_t j = 0; j < m; ++j) {
        if (live[j] != 0 && fstat[j] != static_cast<std::uint32_t>(Status::Ok)) live[j] = 0;
        const fixarray::vec3d fm =
            geo_to_mag * fixarray::vec3d{static_cast<double>(foot[(3 * j) + 0]),
                                         static_cast<double>(foot[(3 * j) + 1]),
                                         static_cast<double>(foot[(3 * j) + 2])};
        const double rn = fixarray::norm(fm);
        theta[j] = rn > 0.0 ? std::acos(std::clamp(fm[2] / rn, -1.0, 1.0)) : 0.0;
        phi_f[j] = std::atan2(fm[1], fm[0]);
    }

    // The work list is RAGGED — every azimuth has its own cell count — so it is flattened to one
    // dimension and dispatched flat. The reduction that follows runs on the HOST, in fp64, in
    // ascending cell order: a tree reduction whose order varies with workgroup size would not be
    // reproducible, and ERROR_BUDGET.md §6 makes ordered reductions a standing requirement.
    std::vector<std::size_t> cell_begin(m + 1, 0);
    for (std::size_t j = 0; j < m; ++j) {
        const std::size_t cells =
            live[j] != 0 ? static_cast<std::size_t>(std::ceil(theta[j] / dth)) : 0;
        cell_begin[j + 1] = cell_begin[j] + cells;
    }
    const std::size_t total_cells = cell_begin[m];
    std::vector<double> column(m, 0.0);
    std::vector<float> cpos(3 * std::min(total_cells, detail::kFluxChunk));
    std::vector<float> cfld(3 * std::min(total_cells, detail::kFluxChunk));
    for (std::size_t base = 0; base < total_cells; base += detail::kFluxChunk) {
        const std::size_t count = std::min(detail::kFluxChunk, total_cells - base);
        std::size_t j = static_cast<std::size_t>(
            std::upper_bound(cell_begin.begin(), cell_begin.end(), base) - cell_begin.begin() - 1);
        for (std::size_t c = 0; c < count; ++c) {
            const std::size_t g = base + c;
            while (g >= cell_begin[j + 1]) ++j;
            const std::size_t i_cell = g - cell_begin[j];
            const double t0 = static_cast<double>(i_cell) * dth;
            const double t1 = std::min(t0 + dth, theta[j]);
            const double tm = 0.5 * (t0 + t1);
            const fixarray::vec3d u{std::sin(tm) * std::cos(phi_f[j]),
                                    std::sin(tm) * std::sin(phi_f[j]), std::cos(tm)};
            const fixarray::vec3d g_geo = mag_to_geo * u;
            cpos[(3 * c) + 0] = static_cast<float>(g_geo[0]);
            cpos[(3 * c) + 1] = static_cast<float>(g_geo[1]);
            cpos[(3 * c) + 2] = static_cast<float>(g_geo[2]);
        }
        // Deliberately NOT folded into `device`: see the return's documentation. The flux cells
        // are a cheap IGRF batch that takes the device lane at 128 points, so folding it in would
        // make the flag read `true` on a call whose traces — 95 % of the cost, and the whole
        // performance claim — ran on the host.
        (void)detail::field_batch(model, coef, norm,
                                  std::span<const float>(cpos).first(3 * count),
                                  std::span<float>(cfld).first(3 * count));
        j = static_cast<std::size_t>(
            std::upper_bound(cell_begin.begin(), cell_begin.end(), base) - cell_begin.begin() - 1);
        for (std::size_t c = 0; c < count; ++c) {
            const std::size_t g = base + c;
            while (g >= cell_begin[j + 1]) ++j;
            const std::size_t i_cell = g - cell_begin[j];
            const double t0 = static_cast<double>(i_cell) * dth;
            const double t1 = std::min(t0 + dth, theta[j]);
            const double tm = 0.5 * (t0 + t1);
            const fixarray::vec3d u{std::sin(tm) * std::cos(phi_f[j]),
                                    std::sin(tm) * std::sin(phi_f[j]), std::cos(tm)};
            const fixarray::vec3d g_geo = mag_to_geo * u;
            // B_r on the unit sphere: the radial component is frame-independent, so the dot
            // product is taken in GEO where the field was evaluated.
            const double br = (static_cast<double>(cfld[(3 * c) + 0]) * g_geo[0]) +
                              (static_cast<double>(cfld[(3 * c) + 1]) * g_geo[1]) +
                              (static_cast<double>(cfld[(3 * c) + 2]) * g_geo[2]);
            column[j] += br * std::sin(tm) * (t1 - t0);
        }
    }

    // ---- 7. Phi on the FOOTPOINTS' own longitudes, then L* ------------------------------------
    std::vector<std::size_t> order(static_cast<std::size_t>(nd));
    for (std::size_t i = 0; i < n; ++i) {
        if (statuses[i] != Status::Ok) continue;   // the starting line already failed
        const std::size_t base = i * static_cast<std::size_t>(nd);
        int good = 0;
        for (int k = 0; k < nd; ++k) good += live[base + static_cast<std::size_t>(k)] != 0 ? 1 : 0;
        out[i].azimuths = good;
        if (good != nd) {
            statuses[i] = Status::NotConverged;   // a gap in the contour is not a cap
            continue;
        }
        for (int k = 0; k < nd; ++k) order[static_cast<std::size_t>(k)] = base + static_cast<std::size_t>(k);
        std::sort(order.begin(), order.end(),
                  [&](std::size_t x, std::size_t y) { return phi_f[x] < phi_f[y]; });
        double flux = 0.0;
        for (int k = 0; k < nd; ++k) {
            const std::size_t j = order[static_cast<std::size_t>(k)];
            const double before = phi_f[order[static_cast<std::size_t>((k + nd - 1) % nd)]] +
                                  (k == 0 ? -detail::kTwoPi : 0.0);
            const double after = phi_f[order[static_cast<std::size_t>((k + 1) % nd)]] +
                                 (k == nd - 1 ? detail::kTwoPi : 0.0);
            flux += column[j] * 0.5 * (after - before);
        }
        out[i].phi = std::abs(flux);
        // A cap enclosing no flux is not a drift shell. Unreachable with a closed contour on a
        // real field — which is exactly why it is a guard and not an assumption.
        out[i].lstar = out[i].phi > 0.0 ? detail::kTwoPi * k0 / out[i].phi : 0.0;
        if (!(out[i].lstar > 0.0)) statuses[i] = Status::NotConverged;
    }
    // The batch's own status is the FIRST thing that went wrong, reported once rather than
    // threaded through every branch above.
    Status worst = Status::Ok;
    for (std::size_t i = 0; i < n; ++i) {
        if (worst == Status::Ok) worst = statuses[i];
    }
    // The trace count is accumulated over the whole batch; report it per point, which is the
    // number that compares against IRBEM's per-point cost.
    const int per_point = traced / static_cast<int>(n);
    for (std::size_t i = 0; i < n; ++i) out[i].traces = per_point;
    return {worst, device};
}

/**
 * Roederer's L\* for one point.
 *
 * The reference lane, and the one to reach for when there is a single point to compute — but NOT
 * the fast one. `Nder = 25` root-finds is a batch of 25, an order of magnitude below the ~512-line
 * crossover `gpu/dispatch.hpp` measures, so a single L\* runs almost entirely on the host however
 * much hardware is present. That is a property of the problem, not of this implementation: the
 * parallelism in L\* lives across POINTS, and @ref make_lstar_batch is where it is taken.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param rotations the epoch's rotations; only `geo_to_mag` and the dipole axis are read.
 * @param start the point, GEO, Earth radii.
 * @param pitch_angle_deg the local pitch angle; 90° is IRBEM's `make_lstar` convention.
 * @param opt the resolution and root-find settings.
 * @return the shell, with the same statuses @ref make_lstar_batch reports. The value is populated
 *         in every case — a failed shell still carries `L_m`, `I` and the fields, which is what
 *         makes the failure diagnosable.
 * @complexity O(Nder × (trials + iterations)) traces plus O(Nder × 180/dθ) field evaluations.
 * @alloc as @ref make_lstar_batch, at `n = 1`.
 * @test IrbemDriftShell.MatchesTheOracleAtIrbemDefaultResolution
 * @test IrbemDriftShell.BatchMatchesThePointAtATimeCall
 */
template <GeoFieldModel M>
[[nodiscard]] inline Result<DriftShell> make_lstar(const M& model,
                                                   const Rotations& rotations,
                                                   const Position<Frame::GEO>& start,
                                                   double pitch_angle_deg = 90.0,
                                                   const DriftShellOptions& opt = {}) {
    const std::array<Position<Frame::GEO>, 1> starts{start};
    const std::array<double, 1> pitch{pitch_angle_deg};
    std::array<DriftShell, 1> shell{};
    std::array<Status, 1> st{};
    const Result<bool> r =
        make_lstar_batch(model, rotations, starts, pitch, shell, st, opt);
    return {r.status == Status::DomainError ? Status::DomainError : st[0], shell[0]};
}

}  // namespace cheatah::space::irbem
