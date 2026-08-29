// om97_bench.cpp — the harness behind the performance numbers quoted for Ostapenko & Maltsev (1997).
//
// A SEPARATE translation unit from irbem_bench.cpp on purpose, for the reason t89_bench.cpp gives:
// that suite carries a routine manifest and a baseline CSV that the whole module's regression
// tracking hangs off. This one measures one thing — the OM97 host and device lanes against each
// other across batch size — which is what the crossover recorded in `space/irbem/gpu/dispatch.hpp`'s
// `irbem_om97_f32` row is derived from.
//
// It is OFF by default like the rest of bench/, and it needs no CMake. Compile it directly, the way
// t89_bench.cpp is compiled (from the repository root, with $B a build tree that already fetched
// google/benchmark and $CHEATAH_DIR a cheatah checkout):
//
//   g++ -std=c++20 -O3 -march=native -ffp-contract=off bench/om97_bench.cpp -I. \
//     -I$CHEATAH_DIR/stdlib/ndarray -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray \
//     -I$CHEATAH_GPU_LINALG_DIR/include -I$CHEATAH_GPU_LINALG_DIR/../cheatah-gpu \
//     -I$CHEATAH_GPU_LINALG_DIR/../cheatah-gpu/build/vk/_deps/volk-src \
//     -I$VULKAN_SDK/include -I$B/_deps/benchmark-src/include \
//     -DCHEATAH_SPACE_IRBEM_SPV_DIR='"'"'"$PWD/build/gpu/shaders"'"'"' \
//     -DCHEATAH_GPU_LINALG_SPV_DIR='"'"'"$PWD/build/gpu/shaders"'"'"' \
//     $B/libbenchmark.a $CHEATAH_GPU_LINALG_DIR/build/libcheatah_gpu_linalg_volk.a -ldl -pthread \
//     -o /tmp/om97_bench && CHEATAH_GPU_LINALG_VK_DEVICE=NVIDIA /tmp/om97_bench
//
// Without cheatah-gpu-linalg on the include path the device benchmarks simply are not compiled and
// the host ones still run: the same __has_include seam the header itself uses.
//
// HOW TO READ IT. Run with `--benchmark_repetitions=5` and read the `_median` aggregate rows; the
// registry's crossover is a best-of, per the convention the other kernel rows use.
//
// WHAT THE NUMBERS MEAN. OM97 is 17 polynomial harmonics under one rotation — ~250 flops for 24
// bytes in and 12 out, ~10 flops/byte, a little under T89's ~11. The prediction under test is that
// it lands beside T89 on the winning side of the ridge, with a crossover of the same order.

#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

#include "space/irbem/ext_ostapenko.hpp"

namespace {

namespace ir = cheatah::space::irbem;
using ir::om97_amplitudes;
using ir::om97_field_host;
using ir::om97_harmonic_count;
using ir::om97_param_block;
using ir::om97_param_count;
using ir::Om97Drivers;

/// A realistic mid-range tilt, so the five tilt harmonics are actually evaluated.
constexpr double kTilt = 0.31;

/// The corpus's moderate drivers.
constexpr Om97Drivers kDrivers{-42.0, 3.2, 35.0, -5.0};

/// A deterministic scatter of `n` points at 3..10 Re, xyz-interleaved in float.
std::vector<float> scatter(std::size_t n) {
    std::vector<float> out(3 * n);
    std::uint64_t s = 0x9E3779B97F4A7C15ULL;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) / 9007199254740992.0;
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double r = 3.0 + (7.0 * next());
        const double th = std::acos(1.0 - (2.0 * next()));
        const double ph = 6.283185307179586 * next();
        out[(3 * i) + 0] = static_cast<float>(r * std::sin(th) * std::cos(ph));
        out[(3 * i) + 1] = static_cast<float>(r * std::sin(th) * std::sin(ph));
        out[(3 * i) + 2] = static_cast<float>(r * std::cos(th));
    }
    return out;
}

/// The fp32 host lane at batch size `state.range(0)`.
void BM_om97_host_f32(benchmark::State& state) {
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const std::vector<float> pos = scatter(n);
    std::vector<float> out(3 * n);
    const std::array<float, om97_harmonic_count> amp = om97_amplitudes<float>(kDrivers);
    const float sp = static_cast<float>(std::sin(kTilt));
    const float cp = static_cast<float>(std::cos(kTilt));
    for (auto _ : state) {
        benchmark::DoNotOptimize(om97_field_host(pos, out, sp, cp, amp));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(n) * static_cast<std::int64_t>(state.iterations()));
}
BENCHMARK(BM_om97_host_f32)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12)->Arg(1 << 14)->Arg(1 << 16)->Arg(1 << 20)->Arg(1 << 22)->UseRealTime();

#if CHEATAH_SPACE_IRBEM_OM97_GPU
/// The device lane at batch size `state.range(0)`, transfers included, one warm-up dispatch first.
void BM_om97_device_f32(benchmark::State& state) {
    if (!ir::gpu::available()) {
        state.SkipWithError(ir::gpu::unavailable_reason().c_str());
        return;
    }
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const std::vector<float> pos = scatter(n);
    std::vector<float> out(3 * n);
    const std::array<float, om97_param_count> block = om97_param_block(
        static_cast<float>(std::sin(kTilt)), static_cast<float>(std::cos(kTilt)),
        om97_amplitudes<float>(kDrivers));
    ir::gpu::dispatch_batch("irbem_om97_f32", std::span<const float>(pos.data(), 3 * n), out,
                            std::span<const float>(block));
    for (auto _ : state) {
        ir::gpu::dispatch_batch("irbem_om97_f32", std::span<const float>(pos.data(), 3 * n), out,
                                std::span<const float>(block));
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(n) * static_cast<std::int64_t>(state.iterations()));
}
BENCHMARK(BM_om97_device_f32)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12)->Arg(1 << 14)->Arg(1 << 16)->Arg(1 << 20)->Arg(1 << 22)->UseRealTime();
#endif

}  // namespace

BENCHMARK_MAIN();
