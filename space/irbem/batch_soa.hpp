#pragma once

/**
 * @file batch_soa.hpp
 * @brief space.irbem — the IGRF harmonic sum evaluated for a STRIP of points at once, with the
 *        point index as the SIMD lane.
 *
 * @ref igrf.hpp's kernel cannot vectorize within one evaluation, and this is structural, not a
 * missed optimisation: the Legendre recursion `Aⁿₘ = eₙₘ·u·Aⁿ⁻¹ₘ - fₙₘ·Aⁿ⁻²ₘ` is loop-carried in
 * `n`, its inner `m` loop is at most `NMAX`-wide with a different trip count every iteration, and
 * the compiler proves it: at `-O3 -march=native` the scalar kernel's hot path retires **1 819
 * scalar mul/add against 348 packed** — 84 % scalar. The dependency-free axis is ACROSS points:
 * no term of one point's sum reads anything of another's. So this file transposes the batch loop.
 * The `n`/`m` recursion loops stay outer, and the innermost, contiguous axis is a strip of
 * @ref igrf_strip_points points evaluated in lockstep — each vector lane carries one point.
 *
 * **This is a scheduling change, not an arithmetic change.** Every lane performs bit-for-bit the
 * scalar kernel's operation sequence: same expressions, same associativity, same
 * `integrand`/`accum` casts (@ref policy.hpp), no reassociation, no fused multiply-add introduced
 * (the build is `-ffp-contract=off` and vector-extension arithmetic obeys it). The result is
 * asserted `memcmp`-identical to @ref Igrf::evaluate per point over the whole test corpus —
 * poles, axis points, all local times and both truncations — by
 * IrbemBatchSimd.BatchLaneIsBitIdenticalToTheScalarLane. Because the arithmetic is unchanged,
 * everything proven of the scalar kernel (the div-B identity, the oracle differentials, the
 * ERROR_BUDGET numbers) transfers to this lane verbatim.
 *
 * @note Bit identity is a property of the BUILD, not of this source alone. It holds under the
 *       `-ffp-contract=off` this repository sets globally (CMakeLists.txt, bench/CMakeLists.txt)
 *       — the same flag the scalar kernel's pinned goldens already require. Under GCC's default
 *       contraction (`-std=gnu++20` with no `-ffp-contract` flag) the compiler fuses different
 *       multiply-adds in the scalar loop and the strip, and the suite catches it: 2 of 7 tests
 *       fail (BatchLaneIsBitIdenticalToTheScalarLane and FastPolicyStripIsBitIdenticalToo —
 *       sparse per-point mismatches over 4099 points, measured 2026-08-28). A consumer compiling
 *       these headers outside the repository's flags gets a lane that is still within the error
 *       budget but NOT bit-identical; the memcmp tests are the sentinel, not the prose.
 *
 * The lanes are spelled with the GCC/Clang vector extension (`__attribute__((vector_size)))`,
 * @ref detail::LaneVec) rather than left to autovectorisation, because autovectorisation was
 * tried first and lost: with the strip written as ordinary arrays indexed `[slot][lane]`, GCC 13
 * fully unrolls the 8-iteration lane loops and then fails to SLP-vectorise the unrolled recursion
 * (`-fopt-info-vec-missed`: "complicated access pattern"), measuring **2.0×** — the win capped by
 * the still-scalar recursion. Typed vectors make the packed code the compiler's only choice; the
 * emitted `mulpd`/`fmadd…pd` presence is itself a regression-gated assertion
 * (IrbemBatchSimd.StripCodegenCarriesPackedVectorOps), so this cannot silently decay back to
 * scalar. One 512-bit vector per row was also tried and measured **0.65×** on AVX2 (GCC's
 * double-pumped emulation spills); two native 4-lane vectors per strip — @ref detail::kStripUnroll
 * — both fit the register file and double the independent dependency chains, which is what hides
 * the recursion's multiply-add latency.
 *
 * Measured by [`bench/irbem_bench.cpp`](../../bench/irbem_bench.cpp) (`BM_cpu_igrf_batch` vs
 * `BM_cpu_igrf_batch_soa`, 2¹⁶ points, degree 13, `-O3 -march=native -ffp-contract=off`, seven
 * interleaved repetitions pinned to one i7-12700H P-core, cv ≈ 3 %): scalar **23.7 ms = 362
 * ns/point**, strip **6.57 ms = 100 ns/point** — **3.6×**, the residual scalar work being the
 * per-strip gather/scatter, the square-root head and `|B|`. Pinning is part of the method, not a
 * flourish: on that hybrid part an unpinned run lands on an E-core often enough to read anywhere
 * from 107 to 200 ns/point, and an E-core's narrower vector units measure only 1.9× — the win is a
 * property of the core's SIMD width, and the number quoted names the core. The GPU lanes of
 * @ref field.hpp are untouched: this file raises the floor the device crossover is measured
 * against, nothing else.
 */

#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>

#include "frames.hpp"
#include "igrf.hpp"
#include "policy.hpp"
#include "status.hpp"

namespace cheatah::space::irbem {

namespace detail {

/**
 * The native SIMD row type: @p V lanes of @p T, one point per lane.
 *
 * A fully specialised trait rather than `T __attribute__((vector_size(sizeof(T) * V)))` written
 * inline, because GCC silently ignores `vector_size` on a type that is still dependent — inside a
 * template the attribute must land on a concrete type, and a trait specialisation is where one
 * exists. Only the widths the strip evaluator actually uses are specialised: 32 bytes, the widest
 * vector AVX2 executes natively (a 64-byte row was measured at 0.65× — see the file brief).
 *
 * @tparam T the lane type — a policy's `integrand` or `accum` (@ref policy.hpp).
 * @tparam V lanes per row.
 */
template <std::floating_point T, std::size_t V>
struct LaneVec;

/// The 4-lane fp64 row — one AVX2 `ymm` of doubles. @copydetails LaneVec
template <>
struct LaneVec<double, 4> {
    /// Four doubles, elementwise arithmetic, per-lane bit-identical to scalar.
    using type = double __attribute__((vector_size(32)));
};

/// The 8-lane fp32 row — one AVX2 `ymm` of floats, for the `Fast` integrand.
template <>
struct LaneVec<float, 8> {
    /// Eight floats, elementwise arithmetic, per-lane bit-identical to scalar.
    using type = float __attribute__((vector_size(32)));
};

/// Lanes in one native 32-byte row of @p T: 4 doubles or 8 floats.
/// @tparam T the lane type.
template <std::floating_point T>
inline constexpr std::size_t kNativeLanes = 32UL / sizeof(T);

/// Independent rows evaluated per strip. Two, because the Legendre recursion is a multiply-add
/// dependency chain: a second row costs registers but doubles the in-flight chains, measured worth
/// it (3.6× vs ~2.9× single-row); four rows measured WORSE (3.3× and falling — register spills).
inline constexpr std::size_t kStripUnroll = 2;

}  // namespace detail

/**
 * Points evaluated per strip by @ref igrf_batch_host — the natural blocking for a caller that
 * wants its batch length to divide evenly (any length works; a remainder runs through the scalar
 * kernel, which is bit-identical anyway).
 *
 * 8 for the fp64 `Exact` policy (2 rows × 4 lanes), 16 for the fp32 `Fast` integrand.
 *
 * @tparam P the precision policy (@ref policy.hpp).
 */
template <SoundPrecision P>
inline constexpr std::size_t igrf_strip_points =
    detail::kStripUnroll * detail::kNativeLanes<typename P::integrand>;

namespace detail {

/**
 * Evaluate the IGRF kernel for one strip of @ref igrf_strip_points points in SoA form — the
 * vector twin of @ref Igrf::evaluate, lane `l` of every row carrying point `l`.
 *
 * The four stages are the scalar kernel's, transposed (equations and provenance in
 * @ref igrf.hpp's brief — Alken et al., *IGRF: the fourteenth generation*, Earth Planets Space
 * (2025); recursions per Winch et al., Geophys. J. Int. 160 (2005) 487-504):
 *
 *  1. head: `r`, `1/r`, `x̂ ŷ û` per lane — scalar `std::sqrt`/divide per point, exactly the
 *     scalar kernel's;
 *  2. the `Aⁿₘ`/`A′ⁿₘ` Legendre recursion, `n`/`m` loops outer, rows inner;
 *  3. the `cₘ`/`sₘ` azimuthal recursion, likewise;
 *  4. the five harmonic sums, per-degree partials in `integrand`, totals in `accum` — the same
 *     ordered fp64 reduction, per lane.
 *
 * Deliberately `noinline`: IrbemBatchSimd.StripCodegenCarriesPackedVectorOps disassembles this
 * exact symbol to assert the packed instructions are still there, and an inlined copy would leave
 * it nothing to find. The call costs one `call`/`ret` per @ref igrf_strip_points points.
 *
 * @tparam NMAX the truncation degree.
 * @tparam P the precision policy; lanes are `P::integrand`, sums `P::accum`.
 * @param model the model, already built for the epoch.
 * @param px @ref igrf_strip_points x-coordinates, GEO, Earth radii.
 * @param py the y-coordinates, likewise.
 * @param pz the z-coordinates, likewise.
 * @param bx receives @ref igrf_strip_points x-components, nT.
 * @param by the y-components, likewise.
 * @param bz the z-components, likewise.
 * @complexity O(NMAX²) vector operations per strip — the scalar kernel's flop count once, for
 *             @ref igrf_strip_points points.
 * @alloc none — every stage is a fixed `std::array` of vector rows on the stack, ~17 KiB at
 *        degree 13.
 * @test IrbemBatchSimd.BatchLaneIsBitIdenticalToTheScalarLane
 * @test IrbemBatchSimd.StripCodegenCarriesPackedVectorOps
 */
template <int NMAX, SoundPrecision P>
[[gnu::noinline]] inline void igrf_evaluate_strip(const Igrf<NMAX, P>& model,
                                                  const double* px, const double* py,
                                                  const double* pz, double* bx, double* by,
                                                  double* bz) {
    using integrand = typename P::integrand;
    using accum = typename P::accum;
    constexpr std::size_t V = kNativeLanes<integrand>;
    constexpr std::size_t U = kStripUnroll;
    using ivec = typename LaneVec<integrand, V>::type;
    using avec = typename LaneVec<accum, kNativeLanes<accum>>::type;
    // One accum row per integrand row: with a float integrand each 8-lane ivec is carried into
    // TWO 4-lane fp64 avecs, so the accum side is indexed [U][V/kNativeLanes<accum>] flattened.
    constexpr std::size_t AV = kNativeLanes<accum>;
    constexpr std::size_t AU = (U * V) / AV;
    using irow = std::array<ivec, U>;
    using arow = std::array<avec, AU>;
    constexpr std::size_t kSlots = triangular_slots(NMAX);
    static constexpr auto kNorm = make_legendre_normalisation<NMAX, integrand>();

    // 1. Head: radius, inverse radius, direction cosines — per lane, exactly the scalar kernel's
    // one square root and one divide. `inv_r` stays fp64 (the scalar kernel's `1.0 / r`); the
    // direction cosines are cast to the integrand as the scalar kernel casts them.
    std::array<double, U * V> inv_r;
    irow ux;
    irow uy;
    irow uz;
    for (std::size_t j = 0; j < U; ++j) {
        for (std::size_t l = 0; l < V; ++l) {
            const std::size_t q = (j * V) + l;
            const double x = px[q];
            const double y = py[q];
            const double z = pz[q];
            const double r = std::sqrt((x * x) + (y * y) + (z * z));
            inv_r[q] = 1.0 / r;
            ux[j][l] = static_cast<integrand>(x * inv_r[q]);
            uy[j][l] = static_cast<integrand>(y * inv_r[q]);
            uz[j][l] = static_cast<integrand>(z * inv_r[q]);
        }
    }

    // 2. Aⁿₘ(u) and its u-derivative — the scalar kernel's recursion with the rows innermost.
    // Zero-initialised for the same reason the scalar arrays are: the k2 read at m = n-1 aliases
    // an always-finite slot and is annihilated by an exact-zero fₙₘ, and the diagonal slots of
    // `da` stay exactly zero.
    std::array<irow, kSlots> a{};
    std::array<irow, kSlots> da{};
    for (std::size_t j = 0; j < U; ++j) {
        a[0][j] = ivec{} + static_cast<integrand>(1);
        a[1][j] = uz[j];
        da[1][j] = ivec{} + static_cast<integrand>(1);
        a[2][j] = ivec{} + kNorm.diagonal[1];
    }
    for (int n = 2; n <= NMAX; ++n) {
        const std::size_t base = slot_index(n, 0);
        const std::size_t prev = base - static_cast<std::size_t>(n);
        const std::size_t prev2 = prev - static_cast<std::size_t>(n - 1);
        for (int m = 0; m < n; ++m) {
            const std::size_t k = base + static_cast<std::size_t>(m);
            const std::size_t k1 = prev + static_cast<std::size_t>(m);
            const std::size_t k2 = prev2 + static_cast<std::size_t>(m);
            const integrand e = kNorm.e[k];
            const integrand f = kNorm.f[k];
            for (std::size_t j = 0; j < U; ++j) {
                a[k][j] = (e * uz[j] * a[k1][j]) - (f * a[k2][j]);
                da[k][j] = (e * (a[k1][j] + (uz[j] * da[k1][j]))) - (f * da[k2][j]);
            }
        }
        for (std::size_t j = 0; j < U; ++j) {
            a[base + static_cast<std::size_t>(n)][j] =
                ivec{} + kNorm.diagonal[static_cast<std::size_t>(n)];
        }
    }

    // 3. cₘ = sinᵐθ·cos mφ, sₘ = sinᵐθ·sin mφ — index 0 is the scalar kernel's pad slot that lets
    // the m = 0 term of the U/V sums read "cₘ₋₁" branch-free.
    std::array<irow, static_cast<std::size_t>(NMAX) + 2> cs{};
    std::array<irow, static_cast<std::size_t>(NMAX) + 2> ss{};
    for (std::size_t j = 0; j < U; ++j) cs[1][j] = ivec{} + static_cast<integrand>(1);
    for (int m = 1; m <= NMAX; ++m) {
        const auto i = static_cast<std::size_t>(m);
        for (std::size_t j = 0; j < U; ++j) {
            cs[i + 1][j] = (ux[j] * cs[i][j]) - (uy[j] * ss[i][j]);
            ss[i + 1][j] = (ux[j] * ss[i][j]) + (uy[j] * cs[i][j]);
        }
    }

    // 4. The five harmonic sums — per-degree partials in the integrand, totals in the accumulator,
    // the scalar kernel's ordered reduction per lane. The coefficient reads go through the public
    // accessors, hoisted to one scalar broadcast per (n, m) — in range by construction, so they
    // return exactly the slot the scalar kernel indexes.
    arow t0{};
    arow t1{};
    arow t2{};
    arow su{};
    arow sv{};
    arow rpow;
    arow irv;  // inv_r re-rowed to the accumulator width, so `rpow *= 1/r` is one vector op
    for (std::size_t j = 0; j < AU; ++j) {
        for (std::size_t l = 0; l < AV; ++l) {
            const double ir = inv_r[(j * AV) + l];
            irv[j][l] = ir;
            rpow[j][l] = ir * ir * ir;
        }
    }
    for (int n = 1; n <= NMAX; ++n) {
        irow p0{};
        irow p1{};
        irow p2{};
        irow pu{};
        irow pv{};
        const std::size_t base = slot_index(n, 0);
        for (int m = 0; m <= n; ++m) {
            const std::size_t k = base + static_cast<std::size_t>(m);
            const auto i = static_cast<std::size_t>(m);
            const integrand gc = model.g(n, m);
            const integrand hc = model.h(n, m);
            const auto fm = static_cast<integrand>(m);
            for (std::size_t j = 0; j < U; ++j) {
                const ivec bigg = (gc * cs[i + 1][j]) + (hc * ss[i + 1][j]);
                const ivec ga = bigg * a[k][j];
                p0[j] += ga;
                p2[j] += fm * ga;
                p1[j] += bigg * da[k][j];
                const ivec ma = fm * a[k][j];
                pu[j] += ma * ((gc * cs[i][j]) + (hc * ss[i][j]));
                pv[j] += ma * ((hc * cs[i][j]) - (gc * ss[i][j]));
            }
        }
        for (std::size_t j = 0; j < AU; ++j) {
            avec ap0;
            avec ap1;
            avec ap2;
            avec apu;
            avec apv;
            for (std::size_t l = 0; l < AV; ++l) {
                const std::size_t q = (j * AV) + l;
                ap0[l] = static_cast<accum>(p0[q / V][q % V]);
                ap1[l] = static_cast<accum>(p1[q / V][q % V]);
                ap2[l] = static_cast<accum>(p2[q / V][q % V]);
                apu[l] = static_cast<accum>(pu[q / V][q % V]);
                apv[l] = static_cast<accum>(pv[q / V][q % V]);
            }
            t0[j] += rpow[j] * static_cast<accum>(n + 1) * ap0;
            t1[j] += rpow[j] * ap1;
            t2[j] += rpow[j] * ap2;
            su[j] += rpow[j] * apu;
            sv[j] += rpow[j] * apv;
            rpow[j] *= irv[j];
        }
    }

    // Recombine — the scalar kernel's closing three lines, per lane.
    for (std::size_t j = 0; j < AU; ++j) {
        for (std::size_t l = 0; l < AV; ++l) {
            const std::size_t q = (j * AV) + l;
            const auto w = static_cast<accum>(uz[q / V][q % V]);
            const accum radial = t0[j][l] + (w * t1[j][l]) + t2[j][l];
            bx[q] = (static_cast<double>(ux[q / V][q % V]) * radial) - su[j][l];
            by[q] = (static_cast<double>(uy[q / V][q % V]) * radial) - sv[j][l];
            bz[q] = (w * t0[j][l]) - ((1.0 - (w * w)) * t1[j][l]) + (w * t2[j][l]);
        }
    }
}

}  // namespace detail

/**
 * The internal field at every point of a batch, on the host, through the SIMD strip evaluator —
 * the vectorised twin of the `Igrf::evaluate`-per-point loop, bit-identical to it.
 *
 * Whole strips of @ref igrf_strip_points go through @ref detail::igrf_evaluate_strip; the
 * remainder (at most a strip minus one) runs the scalar kernel, which computes the identical
 * bits, so a caller can neither observe nor need to care where the boundary fell. This is what
 * @ref field_batch's host lane calls; it is public so a caller who KNOWS the batch is host-bound
 * can skip the device-crossover reasoning entirely.
 *
 * Measured by [`bench/irbem_bench.cpp`](../../bench/irbem_bench.cpp): 100 ns/point against the
 * scalar loop's 362 at 2¹⁶ points, degree 13, pinned to a P-core — 3.6×; see the file brief for
 * the methodology and the variants that measured worse.
 *
 * @tparam NMAX the truncation degree.
 * @tparam P the precision policy (@ref policy.hpp).
 * @param model the model, already built for the epoch.
 * @param points the points, GEO, Earth radii.
 * @param b receives one field vector per point, GEO, nT; same length as @p points.
 * @param b_mag receives `|B|` per point, nT, computed from the returned vector exactly as the
 *        scalar loop computes it; same length as @p points.
 * @return @ref Status::DomainError on a length mismatch (nothing written), @ref Status::Ok
 *         otherwise.
 * @complexity O(N·NMAX²) — the scalar loop's flop count, retired ~3.6× faster.
 * @alloc none — the strip workspace is a fixed stack array; asserted by
 *        IrbemBatchSimd.BatchLaneNeverTouchesTheHeap.
 * @test IrbemBatchSimd.BatchLaneIsBitIdenticalToTheScalarLane
 * @test IrbemBatchSimd.EveryTailLengthIsBitIdentical
 * @test IrbemBatchSimd.BatchRefusesMismatchedSpans
 * @test IrbemBatchSimd.BatchLaneNeverTouchesTheHeap
 */
template <int NMAX, SoundPrecision P>
[[nodiscard]] inline Status igrf_batch_host(const Igrf<NMAX, P>& model,
                                            std::span<const Position<Frame::GEO>> points,
                                            std::span<FieldVector<Frame::GEO>> b,
                                            std::span<double> b_mag) {
    const std::size_t n = points.size();
    if (b.size() != n || b_mag.size() != n) return Status::DomainError;

    constexpr std::size_t W = igrf_strip_points<P>;
    std::array<double, W> px;
    std::array<double, W> py;
    std::array<double, W> pz;
    std::array<double, W> bx;
    std::array<double, W> by;
    std::array<double, W> bz;
    std::size_t i = 0;
    for (; i + W <= n; i += W) {
        for (std::size_t l = 0; l < W; ++l) {
            px[l] = points[i + l].v[0];
            py[l] = points[i + l].v[1];
            pz[l] = points[i + l].v[2];
        }
        detail::igrf_evaluate_strip(model, px.data(), py.data(), pz.data(), bx.data(), by.data(),
                                    bz.data());
        for (std::size_t l = 0; l < W; ++l) {
            b[i + l] = FieldVector<Frame::GEO>{fixarray::vec3d{bx[l], by[l], bz[l]}};
            b_mag[i + l] = b[i + l].magnitude();
        }
    }
    for (; i < n; ++i) {
        b[i] = model.evaluate(points[i]);
        b_mag[i] = b[i].magnitude();
    }
    return Status::Ok;
}

}  // namespace cheatah::space::irbem
