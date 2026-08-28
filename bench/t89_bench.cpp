// t89_bench.cpp — the harness behind every performance number quoted for Tsyganenko (1989).
//
// A SEPARATE translation unit from irbem_bench.cpp on purpose: that suite carries a routine
// manifest and a baseline CSV that the whole module's regression tracking hangs off, and a kernel
// landing in the middle of it would have to land in three places at once. This one measures one
// thing — the T89 host and device lanes against each other across batch size — which is exactly
// what the crossover recorded in `space/irbem/gpu/dispatch.hpp`'s registry row is derived from.
//
// It is OFF by default like the rest of bench/, and it needs no CMake. Compile it directly, the way
// scripts/bench_run.sh compiles the main suite (from the repository root, with $B a build tree that
// already fetched google/benchmark and $CHEATAH_DIR a cheatah checkout):
//
//   g++ -std=c++20 -O3 -march=native -ffp-contract=off bench/t89_bench.cpp -I. \
//     -I$CHEATAH_DIR/stdlib/ndarray -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray \
//     -I$CHEATAH_GPU_LINALG_DIR/include -I$CHEATAH_GPU_LINALG_DIR/../cheatah-gpu \
//     -I$CHEATAH_GPU_LINALG_DIR/../cheatah-gpu/build/vk/_deps/volk-src \
//     -I$VULKAN_SDK/include -I$B/_deps/benchmark-src/include \
//     -DCHEATAH_SPACE_IRBEM_SPV_DIR='"'"'"$PWD/build/gpu/shaders"'"'"' \
//     -DCHEATAH_GPU_LINALG_SPV_DIR='"'"'"$PWD/build/gpu/shaders"'"'"' \
//     $B/libbenchmark.a $CHEATAH_GPU_LINALG_DIR/build/libcheatah_gpu_linalg_volk.a -ldl -pthread \
//     -o /tmp/t89_bench && CHEATAH_GPU_LINALG_VK_DEVICE=NVIDIA /tmp/t89_bench
//
// Without cheatah-gpu-linalg on the include path the device benchmarks simply are not compiled and
// the host ones still run: the same __has_include seam the header itself uses.
//
// HOW TO READ IT. Google Benchmark reports a MEAN, and a synchronous GPU dispatch has a long tail
// (driver scheduling, another process on the device) that a mean tracks and a user never sees. The
// crossover recorded in the registry is a BEST-OF, which is the convention the other kernel rows
// use, so run this with `--benchmark_repetitions=5` and read the `_median` aggregate rows (the
// pinned google/benchmark reports mean, median, stddev and cv, not min). The host lane is steady to
// within a percent either way; the device's coefficient of variation is ~1.4% at 2^22 and ~6% at
// 2^18, which is why the mean and the best-of differ at all.
//
// WHAT THE NUMBERS MEAN. T89 is closed-form straight-line arithmetic — no loop, no data-dependent
// branch — at roughly 400 flops for 24 bytes in and 12 out, i.e. ~11 flops/byte. That is the whole
// prediction being tested: the streaming dipole kernel at 0.5 flops/byte LOSES to the host on this
// seam, IGRF at ~20 wins by ~9x, and T89 should land on IGRF's side. The crossover is where the
// device dispatch floor stops dominating, and it is a measurement, not a policy.

#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

#include "space/irbem/ext_t89.hpp"

namespace {

using cheatah::space::irbem::t89_field_host;
using cheatah::space::irbem::t89_param_block;
using cheatah::space::irbem::t89_param_count;
using cheatah::space::irbem::t89_parameters;
namespace ir = cheatah::space::irbem;

/// The dipole tilt these runs use — a realistic mid-range value, not zero, so every `sin(psi)` and
/// `tan(psi)` term in the model is actually evaluated. Timing a model at the one tilt where half its
/// terms vanish would flatter it.
constexpr double kTilt = 0.35;

/// The Kp bin. Every bin costs exactly the same — the model has no data-dependent branch — so the
/// choice is arbitrary and stated rather than swept.
constexpr int kBin = 4;

/// A deterministic scatter of GSM points over 2..20 R_E: the inner magnetosphere and near tail,
/// which is where a real drift-shell batch lives. A 64-bit LCG, so the array is reproducible and is
/// never a lattice — a lattice would put whole workgroups on identical geometry.
const std::vector<float>& points(std::size_t n) {
    static std::vector<float> buf;
    if (buf.size() < 3 * n) {
        buf.resize(3 * n);
        std::uint64_t s = 0x9E3779B97F4A7C15ULL;
        const auto next = [&s] {
            s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
            return static_cast<double>(s >> 11) / 9007199254740992.0;
        };
        for (std::size_t i = 0; i < n; ++i) {
            const double r = 2.0 + (18.0 * next());
            const double th = std::acos(1.0 - (2.0 * next()));
            const double ph = 6.283185307179586 * next();
            buf[(3 * i) + 0] = static_cast<float>(r * std::sin(th) * std::cos(ph));
            buf[(3 * i) + 1] = static_cast<float>(r * std::sin(th) * std::sin(ph));
            buf[(3 * i) + 2] = static_cast<float>(r * std::cos(th));
        }
    }
    return buf;
}

/// The batch sizes both lanes are measured at, straddling the crossover so the curve — and not just
/// its endpoint — is on the record.
void crossover_sizes(::benchmark::Benchmark* b) {
    for (std::int64_t n : {64, 256, 1024, 2048, 4096, 16384, 65536, 262144, 1048576, 4194304}) {
        b->Arg(n);
    }
}

/// The fp32 host lane — the same expressions the kernel evaluates, in the same precision.
/// @param state the benchmark state; `range(0)` is the batch size.
/// @complexity O(N), ~400 flops per point.
/// @alloc none inside the timed loop.
void BM_cpu_t89_batch(benchmark::State& state) {
    const auto n = static_cast<std::size_t>(state.range(0));
    const std::vector<float>& pos = points(n);
    std::vector<float> out(3 * n);
    const auto sp = static_cast<float>(std::sin(kTilt));
    const auto cp = static_cast<float>(std::cos(kTilt));
    for (auto _ : state) {
        benchmark::DoNotOptimize(
            t89_field_host(std::span<const float>(pos.data(), 3 * n), out, sp, cp, kBin));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(n));
}
BENCHMARK(BM_cpu_t89_batch)->Apply(crossover_sizes)->Unit(benchmark::kMicrosecond)->UseRealTime();

/// The fp64 scalar entry point, one point per call — the lane a caller who has not batched gets.
/// @param state the benchmark state.
/// @complexity O(1) per call.
/// @alloc none.
void BM_cpu_t89_field(benchmark::State& state) {
    const std::vector<float>& pos = points(1024);
    std::size_t i = 0;
    for (auto _ : state) {
        const ir::Position<ir::Frame::GSM> p{cheatah::fixarray::vec3d(
            pos[(3 * i) + 0], pos[(3 * i) + 1], pos[(3 * i) + 2])};
        benchmark::DoNotOptimize(ir::t89_field(p, kTilt, 35.0).value.v[2]);
        i = (i + 1) % 1024;
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_cpu_t89_field);

#if CHEATAH_SPACE_IRBEM_T89_GPU
/// The device lane, transfers included — the honest comparison, since a caller cannot use a kernel
/// without paying to get its data there and back.
/// @param state the benchmark state; `range(0)` is the batch size.
/// @complexity One dispatch over `ceil(N/256)` workgroups, plus `5N` floats over the bus.
/// @alloc `dispatch_batch`'s four pooled device buffers, released before each iteration returns.
void BM_gpu_t89_batch(benchmark::State& state) {
    if (!ir::gpu::available()) {
        state.SkipWithError(ir::gpu::unavailable_reason().c_str());
        return;
    }
    const auto n = static_cast<std::size_t>(state.range(0));
    const std::vector<float>& pos = points(n);
    std::vector<float> out(3 * n);
    const std::array<float, t89_param_count> block = t89_param_block(
        static_cast<float>(std::sin(kTilt)), static_cast<float>(std::cos(kTilt)), kBin);
    // One warm-up outside the timed loop: the first launch of a kernel reads and compiles its
    // pipeline, and timing that once as though it were per-batch work would be a lie about the
    // steady state every real caller sees.
    ir::gpu::dispatch_batch("irbem_t89_f32", std::span<const float>(pos.data(), 3 * n), out,
                            std::span<const float>(block));
    for (auto _ : state) {
        ir::gpu::dispatch_batch("irbem_t89_f32", std::span<const float>(pos.data(), 3 * n), out,
                                std::span<const float>(block));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(n));
}
BENCHMARK(BM_gpu_t89_batch)->Apply(crossover_sizes)->Unit(benchmark::kMicrosecond)->UseRealTime();
#endif

}  // namespace

BENCHMARK_MAIN();
