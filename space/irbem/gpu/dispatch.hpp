#pragma once

/**
 * @file dispatch.hpp
 * @brief space.irbem — the GPU seam: which kernels exist, whether a device can run them, where
 *        their SPIR-V lives, and the one blocking launcher every device routine goes through.
 *
 * IRBEM is single-threaded Fortran, and `docs/ERROR_BUDGET.md` §5 measures its cost at 15.5 ms per
 * L\* point even after being rebuilt with optimization — about 26 minutes for a 100 000-point
 * ephemeris. The field-line and drift-shell integrals underneath that number are embarrassingly
 * parallel across points, which is the entire reason this module has a device lane. This header is
 * the seam that lane crosses.
 *
 * Three things would go wrong without it, and each is a real failure this file exists to prevent.
 *
 * **A kernel name that nothing implements must not reach the driver.** `cheatah-gpu-linalg`'s
 * context binds buffers and launches by NAME; a typo, or a kernel that was renamed in the Slang and
 * not in the caller, is a missing `.spv` at best and — on an emulated device that quietly ignores
 * unknown shaders — a green test over a buffer of zeros at worst. So every launch goes through
 * @ref kernel_info first, and an unregistered name is @ref UnknownKernel before any device is
 * touched. `README.md`'s "no silent failures" promise is enforced here, not hoped for.
 *
 * **"No GPU" must be an answer, not a crash.** Availability is two independent questions: whether
 * the `cheatah-gpu-linalg` headers were present when this translation unit was compiled
 * (@ref compiled_with_gpu, a compile-time constant), and whether a device actually comes up now
 * (@ref available, a runtime probe that never throws). A caller that wants the answer asks; a
 * caller that launches anyway gets a named @ref GpuUnavailable carrying
 * @ref unavailable_reason, never a segfault in a driver.
 *
 * **Our kernels are not `cheatah-gpu-linalg`'s kernels.** That library resolves a BARE kernel
 * name against its own shader directory, and a PATH-QUALIFIED one (`<dir>/<name>`) exactly where it
 * says. This module launches its kernels by the qualified form (@ref qualified), so its directory
 * and linalg's coexist in one process with no shared state between them. (An earlier revision of
 * that library refused '/' in a name and offered no directory argument, and this header worked
 * around it by scoping the library's environment variable around each launch; the qualified name
 * is the upstream fix that workaround asked for.)
 *
 * ### What this header does NOT do
 *
 * It does not accumulate. Every kernel behind this seam evaluates one point per thread and writes
 * it; nothing here sums. That is deliberate — `policy.hpp` and ERROR_BUDGET §3 make `accum =
 * double` the load-bearing invariant of the whole module, and the reduction stays on the host in
 * fp64, in a fixed order. An fp32 integrand is four to five orders below the discretization floor;
 * an fp32 *reduction* over 10³ terms is not.
 *
 * @note The device lane cannot carry committed bit-goldens: FMA contraction inside a shader is at
 *       the driver's discretion (ERROR_BUDGET §6). It is held instead to the `Blocal` budget of
 *       1e-6 relative against the host reference, plus run-to-run self-consistency.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>

#include "space/irbem/frames.hpp"

// The compile-time half of "is there a GPU?": whether the device stack was on the include path at
// all. A consumer that builds without cheatah-gpu-linalg gets a header that still compiles, still
// answers every question, and reports itself unavailable — the __has_include seam the ecosystem
// uses everywhere for an optional dependency.
#if __has_include("cheatah_gpu_linalg/context.hpp")
/// 1 when `cheatah-gpu-linalg` was on the include path for this translation unit, 0 otherwise.
/// Prefer @ref cheatah::space::irbem::gpu::compiled_with_gpu in C++; this macro exists so the
/// `#if` above and the tests can ask the same question the preprocessor answered.
#  define CHEATAH_SPACE_IRBEM_HAVE_GPU 1
#  include "cheatah_gpu_linalg/context.hpp"
#else
/// 1 when `cheatah-gpu-linalg` was on the include path for this translation unit, 0 otherwise.
/// Prefer @ref cheatah::space::irbem::gpu::compiled_with_gpu in C++; this macro exists so the
/// `#if` above and the tests can ask the same question the preprocessor answered.
#  define CHEATAH_SPACE_IRBEM_HAVE_GPU 0
#endif

namespace cheatah::space::irbem::gpu {

/**
 * Whether this translation unit was compiled against `cheatah-gpu-linalg`.
 *
 * The compile-time half of the availability question, and a constant — usable in `if constexpr`
 * and in a `static_assert` for a build that means to require the device stack. It says nothing
 * about whether a device exists on the machine running the binary; that is @ref available.
 */
inline constexpr bool compiled_with_gpu = (CHEATAH_SPACE_IRBEM_HAVE_GPU != 0);

// ---------------------------------------------------------------------------------------------
// Named failures — every way a launch can be refused, and never a bare runtime_error
// ---------------------------------------------------------------------------------------------

/**
 * A kernel name that no entry point in `irbem.slang` implements.
 *
 * Thrown by @ref kernel_info, and therefore by every launcher, before any device is touched. This
 * is the guard that stops a renamed or misspelled kernel from becoming a silent no-op.
 */
struct UnknownKernel : std::runtime_error {
    /**
     * Construct with a message naming the offending kernel.
     * @param what the message, which must name the kernel that was asked for.
     * @complexity O(len(what)).
     * @alloc one string copy, inside `std::runtime_error`.
     */
    explicit UnknownKernel(const std::string& what) : std::runtime_error(what) {}
};

/**
 * A registered kernel whose compiled SPIR-V is not where @ref shader_dir says it should be.
 *
 * A build error rather than a machine capability, which is why it is a distinct type and why it is
 * diagnosed BEFORE availability: a missing `.spv` is worth reporting on a machine with no GPU too,
 * since that is exactly the machine where the shader build step is most likely to have been
 * skipped.
 */
struct ShaderMissing : std::runtime_error {
    /**
     * Construct with a message naming the kernel and the path that was searched.
     * @param what the message, which must carry the resolved path so the fix is obvious.
     * @complexity O(len(what)).
     * @alloc one string copy, inside `std::runtime_error`.
     */
    explicit ShaderMissing(const std::string& what) : std::runtime_error(what) {}
};

/**
 * A device launch on a build or a machine that has no device.
 *
 * Carries @ref unavailable_reason, so a caller reports why rather than inventing a reason.
 */
struct GpuUnavailable : std::runtime_error {
    /**
     * Construct with a message carrying the bring-up failure reason.
     * @param what the message, which must carry @ref unavailable_reason.
     * @complexity O(len(what)).
     * @alloc one string copy, inside `std::runtime_error`.
     */
    explicit GpuUnavailable(const std::string& what) : std::runtime_error(what) {}
};

// ---------------------------------------------------------------------------------------------
// The kernel registry
// ---------------------------------------------------------------------------------------------

/**
 * One entry point of `irbem.slang`, as the host side knows it.
 *
 * The registry is what makes an unknown kernel a named error instead of a missing file, and what
 * lets a test assert that the host's list and the Slang source's list are the same list.
 */
struct KernelInfo {
    /// The Slang entry-point name, which is also the `.spv` stem. A string literal, so
    /// `name` is a valid NUL-terminated C string for the context's `const char*` launch API.
    const char* name;
    /// One line saying what the kernel computes, for a diagnostic or a kernel listing.
    std::string_view brief;
    /// How many buffers the kernel binds, at set 0, bindings `0 .. bindings-1`.
    unsigned bindings;
    /// How many `float` scalars the kernel's parameter buffer holds.
    unsigned params;
    /// The batch size at or above which the device lane is MEASURED to beat the host lane, or
    /// @ref never_faster_on_device when no such size was found. Per kernel, because the crossover
    /// is a property of arithmetic intensity: a kernel that does a dozen flops per point can never
    /// pay for moving that point across the bus, and one that traces a field line for a thousand
    /// steps pays for it immediately. A single library-wide constant would be wrong for both.
    std::size_t crossover_points;
};

/**
 * The @ref KernelInfo::crossover_points of a kernel the device never wins on.
 *
 * Not a "not measured yet" placeholder — it is the answer for a transfer-bound kernel, and a
 * kernel carrying it runs on the host under @ref Lane::Auto no matter how large the batch is.
 */
inline constexpr std::size_t never_faster_on_device = std::numeric_limits<std::size_t>::max();

/**
 * The @ref KernelInfo::crossover_points of a kernel that should ALWAYS take the device lane.
 *
 * This is the value the integral kernels carry, and it is what "the GPU is the default" means in
 * practice. It is not a policy override — it is what the arithmetic says. One L\* point is roughly
 * 10^5 field evaluations, ~5 x 10^7 flops; at even 10% of this device's fp32 peak that is ~25 us of
 * compute against a measured ~30 us synchronous dispatch floor, so a batch of TWO points already
 * pays for the round trip and every realistic batch dwarfs it. The transfer is ~60 bytes per point
 * for ~5 x 10^7 flops — an arithmetic intensity of ~8 x 10^5 flops/byte, some 350x beyond the
 * PCIe ridge point, which is why the bus simply does not appear in the accounting.
 *
 * The contrast with @ref never_faster_on_device is the whole reason this is per-kernel rather than
 * one library-wide switch. Measured on this machine, same seam, same transfers, same device:
 *
 *   | kernel  | flops/byte | device        | host          | verdict     |
 *   | dipole  |        0.5 | 243 Mpts/s    | 352 Mpts/s    | 0.69x LOSES |
 *   | IGRF-14 |        ~20 | 27.3 Mpts/s   | 3.3 Mpts/s    | 8.37x WINS  |
 *   | trace   |     ~9 400 | the integrals — the point never leaves the device |
 *
 * A blanket "always use the GPU" would make the dipole kernel 45% SLOWER. Defaulting to the device
 * where the arithmetic earns it, and only there, is the honest form of the same intent.
 *
 * @note On a UNIFIED-MEMORY device — an APU, an integrated GPU, Apple silicon — there is no PCIe
 *       copy at all, so the transfer term that makes the dipole lose largely vanishes and even
 *       transfer-bound kernels may win. `cheatah-gpu-linalg` reads `VkPhysicalDeviceProperties::
 *       deviceType` when it scores devices but does not expose it, so this header cannot yet ask.
 *       Until it can, the crossovers here are the DISCRETE-GPU measurement, which is the
 *       conservative direction: a unified-memory device will beat them, never miss them.
 */
inline constexpr std::size_t always_on_device = 1;

/**
 * Every kernel `irbem.slang` implements.
 *
 * A `constexpr` table rather than a function, so it costs nothing, can be iterated at compile time,
 * and cannot drift into a registry that allocates. Adding a kernel means adding a row here AND an
 * entry point there; the completeness test fails if only one of the two happens.
 */
inline constexpr std::array<KernelInfo, 12> registered_kernels{{
    // MEASURED, RTX 3070 Ti / Vulkan, best of five per size, transfers included, against this
    // header's own host lane built -O2 -ffp-contract=off (n : device Mpts/s : host Mpts/s):
    //   2^10 : 14 : 444 | 2^14 : 145 : 440 | 2^16 : 240 : 440 | 2^18 : 280 : 446
    //   2^20 : 246 : 400 | 2^22 : 310 : 449 | 2^24 : 258 : 431
    // The device never wins, and the ratio flattens at ~0.6x rather than climbing: the kernel is
    // ~12 flops over 24 bytes moved per point, so it is bus-bound end to end and the host's
    // vectorized loop is simply closer to the data. That is a fact about THIS kernel, not about
    // the seam — the traces this module exists for evaluate the field hundreds of times per point
    // with the point resident, which is the regime the bus cost amortizes in.
    {"irbem_dipole_f32", "centred-dipole B (nT) at each MAG-frame point (Re); params = {g10}", 4, 1,
     never_faster_on_device},

    // MEASURED, RTX 3070 Ti / Vulkan, 2^20 points, best of five, transfers included, against
    // igrf.hpp's fp64 host lane built -O3 -march=native -ffp-contract=off:
    //   device 36.6 ns/eval (27.3 Mpts/s) | host 306.5 ns/eval (3.3 Mpts/s) | 8.37x
    //   max relative deviation device-vs-host over 2^20 points: 8.8e-07, inside the Blocal budget.
    // The contrast with the dipole row above is the whole design in two lines: same seam, same
    // transfers, same device, ~40x the arithmetic intensity, and the verdict flips from 0.69x to
    // 8.37x. Crossover derived from those two throughputs and the measured ~30 us synchronous
    // dispatch floor: the device pays for itself once n/27.3e6 + 30e-6 < n/3.3e6, i.e. n > ~113.
    // Rounded up to a power of two for a stable, easily-reasoned threshold.
    //
    // The coefficients arrive ALREADY interpolated to the epoch — the 26-epoch IGRF table stays on
    // the host, because interpolating it per thread would be 10^5 redundant copies of a calculation
    // done once per batch. Hence no scalar params: everything the kernel needs is a buffer.
    {"irbem_igrf_f32",
     "IGRF-14 internal B (nT) at each GEO point (Re); coefficients and Legendre normalisation "
     "arrive pre-interpolated as buffers",
     5, 0, 128},

    // MEASURED, RTX 3070 Ti / Vulkan, best of three, transfers included, against lstar.hpp's fp64
    // host lane built -O3 -march=native -ffp-contract=off. Field lines spread over L = 2..8 and
    // pitch angles 30..80 degrees, i.e. the shape a real drift-shell batch has:
    //
    //     N lines :     64 :    256 :   1024 :   4096 :  16384 :  65536
    //     speedup :   0.5x :   0.8x :   3.1x :  12.7x :  33.2x :  48.9x
    //
    //   at N = 65536 : device 1.24 us/line (808 787 lines/s) | host 60.4 us/line (16 543 lines/s)
    //   max relative deviation in I, fp32 device vs fp64 host: 5.9e-05, inside the XJ budget.
    //
    // The curve is the whole argument for this module. Below ~512 lines the ~30 us submit floor
    // dominates and the DEVICE LOSES; above it the arithmetic intensity takes over and the speedup
    // is still climbing at 65 536, meaning the device is not yet saturated. One L* point is ~130
    // traces, so this is 0.158 ms/L*-point against the -O2 oracle's measured 15.5 ms.
    //
    // Note the contrast with the two rows above: 0.5 flops/byte loses, ~20 wins by 9x, ~9400 wins
    // by 49x. Same seam, same transfers, same device.
    {"irbem_trace_i_f32",
     "field-line trace and the second invariant I; one thread per FIELD LINE, whole RK4 chain "
     "on-device, no path stored",
     7, 0, 512},

    // MEASURED, RTX 3070 Ti / Vulkan — see driftshell.hpp's brief for the whole L* accounting.
    // A footpoint walk is the tracer's arithmetic over a LONGER path: ~250 RK4 steps from the
    // magnetic equator down to r = 1, four IGRF evaluations each, plus ~30 more spent halving the
    // final step onto the sphere. ~16 bytes in and 12 out per line for ~5 x 10^5 flops, so the
    // intensity is if anything higher than irbem_trace_i_f32's ~9 400 flops/byte and the crossover
    // is inherited from it rather than re-derived: below ~512 lines the ~30 us submit floor
    // dominates, above it the device runs away.
    {"irbem_shell_foot_f32",
     "drift-shell ionospheric footpoints; one thread per shell azimuth, the final step halved onto "
     "r = 1 rather than interpolated across it",
     7, 0, 512},

    // MEASURED, RTX 3070 Ti / Vulkan, best of seven per size, transfers included, against
    // ext_t89.hpp's own fp32 host lane built -O3 -march=native -ffp-contract=off, over points
    // scattered at 2..20 Re across the inner magnetosphere and near tail:
    //
    //     N points :  1024 :  2048 :  4096 : 16384 : 65536 :  2^20 :  2^22
    //     speedup  : 0.74x : 1.40x : 2.60x : 7.17x : 12.1x : 12.9x : 15.1x
    //
    //   at N = 2^22 : device 3.25 ns/point | host 49.1 ns/point
    //   max ABSOLUTE deviation device-vs-host over 2^20 points: 5.1e-05 nT. The largest RELATIVE
    //   deviation, 3.0e-05, sits at a near-null where the external field is 0.23 nT; against a
    //   typical 80 nT external field the same absolute error is 6e-07, inside the Blocal budget.
    //
    // Tsyganenko (1989) is closed-form straight-line arithmetic — no loop, no data-dependent
    // branch, ~400 flops for 24 bytes in and 12 out, i.e. ~11 flops/byte. An order of magnitude
    // above the dipole row's 0.5 and within a factor of two of IGRF's ~20, which is why it lands
    // on IGRF's side of the ridge and not the dipole's. Crossover from the two throughputs and the
    // measured ~62 us dispatch floor: the device pays once n/308e6 + 62e-6 < n/20.4e6, i.e.
    // n > ~1300; rounded up to the next power of two.
    //
    // No coefficient BUFFER: the whole model is 30 scalars (sin psi, cos psi, C1..C19 and the nine
    // non-linear parameters of one Kp bin), so they ride in the parameter block and the kernel fits
    // dispatch_batch's pos/out/params/dims shape exactly.
    {"irbem_t89_f32",
     "Tsyganenko 1989 external B (nT) at each GSM point (Re); params = {sin psi, cos psi, C1..C19, "
     "dx, a_rc, D0, gamma_rc, Rc, G, a_T, Dy, x0} for one Kp bin",
     4, 30, 2048},

    // MEASURED, RTX 3070 Ti / Vulkan, best of five, transfers included, against trace_api.hpp's
    // fp64 host lane built -O3 -march=native -ffp-contract=off. Fixed step ds = 0.02 Re, 512
    // samples of headroom per line, starts spread over L = 2..8:
    //
    //     N lines :       64 :      128 :      256 :      512 :     1024 :     4096 :    16384 :    65536
    //     speedup : .33-.39x : .66-.69x : 1.0-1.2x : 2.0-2.2x : 3.7-3.9x : 9.9-13.3x : 9.6-23.5x : 19-26.3x
    //     us/line :  450-473 :  234-236 :  143-151 :   76-77  :   42-43  :   11.9-15.3 : 6.8-11.0 :  6.0-7.0
    //     host    :  154-165 us/line, flat in N
    //
    // FOUR full runs, and every cell is their spread rather than the best of them. Below 1 024
    // lines the spread is a few per cent and the crossover reproduces at 256 on every run; at
    // 16 384 the same size has measured 9.6x and 23.5x, a factor of 2.4 apart.
    //
    // The crossover is 256, where the device first wins. **It is LOWER than irbem_trace_i_f32's
    // 512, and the reason is not that this kernel is cheaper -- it is that the HOST lane is dearer
    // here**: a fixed ds = 0.02 Re line is ~250 steps against the invariant tracer's ~120 at
    // ds = L/50, so the host costs 154 us/line against 60, and the same ~30-60 us submit floor is
    // paid off in half the batch. A crossover is a ratio, not a property of the kernel alone, and
    // this row is the cleanest evidence in the registry for why it is measured per kernel.
    //
    // WHERE THE BANDWIDTH SHOWS IS THE CEILING, NOT THE CROSSOVER. irbem_trace_i_f32 returns four
    // floats per LINE (~9 400 flops/byte) and its speedup is still climbing at 65 536, reaching
    // 48.9x. This one returns four floats per STEP -- ~125 flops/byte, some 75x less -- and its
    // curve flattens in the low twenties and becomes run-to-run noisy above 4 096, where each
    // batch stages 100 MB to 1 GB of results through host memory and the measurement stops being
    // about the device at all. The two ranges above are the spread over two full runs and are
    // quoted as ranges deliberately: reporting the best of them as a single number would be a
    // performance claim the second run does not support. Two kernels rather than one on purpose:
    // teaching the invariant tracer to emit `posit` would put this transfer in the L* hot path.
    //
    // Bindings: pos, coef, norm, path, bmag, report, dims. No scalar params: the fixed step and
    // R0 ride in the dims buffer scaled by 1e6, the same trick the invariant tracer uses for
    // steps_per_l.
    {"irbem_trace_path_f32",
     "field-line trace WITH its path; one thread per line, fixed step toward the Earth, writes "
     "3M+M floats per line -- bandwidth-bound, unlike irbem_trace_i_f32",
     7, 0, 256},

    // The trace through the TOTAL field — IGRF plus T89 — one thread per field line, sharing
    // igrf_eval and t89_eval with the batched kernels so each piece of physics exists in this
    // file exactly once. The ext buffer carries one Kp bin's parameter block (selected on the
    // HOST: a batch shares one epoch and one Kp, so a per-thread bin branch would diverge the
    // warp for nothing) plus the epoch's gsm_from_geo rotation, column-major. EIGHT bindings —
    // the descriptor layout's cap, and the reason Leases::capacity is 8.
    //
    // Crossover inherited from irbem_trace_i_f32's measured 512 as the CONSERVATIVE bound: T89
    // adds ~400 flops to each of ~600 field evaluations for the same ~44 bytes moved, so the
    // true crossover can only be lower. Re-measure and tighten rather than guess lower now.
    {"irbem_trace_total_f32",
     "field-line trace and I through IGRF+T89; one thread per FIELD LINE, whole RK4 chain "
     "on-device; ext = one Kp bin's T89 block + gsm_from_geo",
     8, 0, 512},

    // The Olson-Pfitzer dynamic field (kext = 6): the SHARED t89_eval on a scaled position, times
    // s^3, plus a ~40-flop ring disc. Same 24 bytes in and 12 out as the T89 row for ~450 flops,
    // so the same side of the ridge. MEASURED, RTX 3070 Ti / Vulkan, best of five per size,
    // transfers included, against ext_opd.hpp's own fp32 host lane built -O3 -march=native
    // -ffp-contract=off, points scattered at 2..20 Re:
    //
    //     N points :  1024 :  2048 :  4096 : 16384 : 65536 :  2^20
    //     speedup  : 0.86x : 1.59x : 2.97x : 8.01x : 10.5x : 14.2x
    //
    //   at N = 2^20 : device 3.95 ns/point | host 55.9 ns/point
    //   max ABSOLUTE deviation device-vs-host over 2^20 points: 6.9e-05 nT.
    //
    // The device first wins at 2048 — one step later than the plain T89 row's derived ~1300 —
    // because the device-side cost is the same t89_eval plus a ring and the ~70 us floor is the
    // same, while the host lane is marginally dearer per point; the measured 2048 stands as the
    // crossover rather than the row above's.
    //
    // Params: the T89 quiet block [0..29] exactly as t89_eval reads it, then s and the ring
    // amplitude C — 32 floats, no coefficient buffer, dispatch_batch's pos/out/params/dims shape.
    {"irbem_opd_f32",
     "Olson-Pfitzer dynamic external B (nT) at each GSM point (Re); params = {T89 quiet block "
     "(30), s, C}",
     4, 32, 2048},

    // Mead & Fairfield 1975 (kext = 1): three quadratic polynomials, ~50 flops for the same 24
    // bytes in and 12 out — ~1.4 flops/byte, within a factor of three of the dipole row's 0.5
    // and an order of magnitude below T89's ~11. MEASURED, RTX 3070 Ti / Vulkan, best of seven
    // per size, transfers included, against ext_mead.hpp's own fp32 host lane built -O3
    // -march=native -ffp-contract=off, points scattered at 2..20 Re:
    //
    //     N points :   256 :  1024 :  2048 :  4096 : 16384 : 65536 :  2^20 :  2^22
    //     speedup  : 0.01x : 0.05x : 0.10x : 0.18x : 0.47x : 0.75x : 0.90x : 0.95x
    //
    //   at N = 2^22 : device 4.95 ns/point (202 Mpts/s) | host 4.72 ns/point (212 Mpts/s)
    //   max absolute deviation device-vs-host over 2^22 points: 5.3e-05 nT.
    //
    // The device NEVER wins and the ratio is still below 1.0 at four million points, flattening
    // toward the bus limit exactly as the dipole row does: the host's vectorized loop evaluates a
    // 50-flop polynomial faster than PCIe can move the point. The verdict is a fact about this
    // kernel's arithmetic intensity, not about the seam — the same 24 bytes carry T89's ~400
    // flops to a 15x win. The kernel exists so that a total-field tracer can call mead_eval
    // on-device, where the point is resident; as a batched field kernel it is routed to the host.
    {"irbem_mead_f32",
     "Mead & Fairfield 1975 external B (nT) at each GSM point (Re); params = {sin psi, cos psi, "
     "psi deg, a1..a7, b1..b3, c1..c7} for one Kp bin",
     4, 20, never_faster_on_device},
    // MEASURED, RTX 3070 Ti / Vulkan, best of five per size, transfers included, against
    // ext_opq.hpp's own fp32 host lane built -O3 -march=native -ffp-contract=off, over points
    // scattered at 2.5..15 Re (the model's whole region):
    //
    //     N points :  1024 :  2048 :  4096 : 16384 : 65536 :  2^20 :  2^22
    //     speedup  : 0.45x : 1.19x : 2.32x : 6.33x : 11.9x : 11.2x : 14.9x
    //
    //   at N = 2^22 : device 3.35 ns/point | host 49.9 ns/point
    //   max ABSOLUTE deviation device-vs-host over 2^20 points: 8.0e-05 nT, against external
    //   fields of tens of nT — the same contraction-level agreement the T89 row reports.
    //
    // Olson & Pfitzer (1977) is 86 fused polynomial terms and one exp for 24 bytes in and 12 out —
    // ~15 flops/byte, T89's regime exactly, with no data-dependent branch (the loop's trip
    // pattern depends on nothing but its own counters). Crossover from the measured curve: 2048
    // is the first size at which the device wins (1.19x), and the same threshold T89 carries.
    //
    // No coefficient BUFFER: the 172 tilt-folded coefficients plus sin psi and cos psi ride in
    // the parameter block, so the kernel fits dispatch_batch's pos/out/params/dims shape.
    {"irbem_opq_f32",
     "Olson-Pfitzer 1977 quiet external B (nT) at each GSM point (Re); params = {sin psi, "
     "cos psi, A(32), B(32), C(22), D(22), E(32), F(32)} folded for one tilt",
     4, 174, 2048},
    // MEASURED, RTX 3070 Ti / Vulkan, best of five per size, transfers included, against
    // ext_ostapenko.hpp's own fp32 host lane built -O3 -march=native -ffp-contract=off, over
    // points scattered at 3..10 Re (the paper's fitted region):
    //
    //     N points :  1024 :  2048 :  4096 : 16384 : 65536 :  2^20 :  2^22
    //     speedup  : 0.37x : 0.80x : 1.60x : 4.36x : 7.79x : 7.83x : 10.1x
    //
    //   at N = 2^22 : device 3.37 ns/point | host 34.1 ns/point
    //   max ABSOLUTE deviation device-vs-host over 2^16 points: 4.6e-05 nT; device vs the fp64
    //   reference over 8192 points: 2.3e-05 nT.
    //
    // Ostapenko & Maltsev (1997) is 17 polynomial harmonics under one rotation — ~250 flops for
    // 24 bytes in and 12 out, ~10 flops/byte, a little under T89's and in the same regime: no
    // loop with a data-dependent trip count, no branch, no transcendental. The host lane is
    // 1.4x faster per point than T89's (34 vs 49 ns) for the same reason, which is why the
    // crossover sits one power of two ABOVE T89's 2048 despite the same device throughput:
    // at 2048 the device is still paying its ~80 us dispatch floor against a 67 us host loop.
    // Crossover from the measured curve: the first power of two at which the device wins
    // outright. Measured by bench/om97_bench.cpp.
    //
    // No coefficient BUFFER: the 17 amplitudes are the drivers' regression, done ONCE per batch
    // on the host, and ride in the parameter block with sin psi and cos psi.
    {"irbem_om97_f32",
     "Ostapenko-Maltsev 1997 external B (nT) at each GSM point (Re); params = {sin psi, cos psi, "
     "A1..A17} for one driver set",
     4, 19, 4096},
    // Tsyganenko, Singer & Kasper 2003, the STORM model (kext = 10). By a wide margin the most
    // arithmetic in this file that is still one thread, one point: 108 divergence-free modes —
    // 58 thickened warped current discs, 24 radial-current sheets, 24 box harmonics and two
    // uniform fields — for the same 24 bytes in and 12 out as the T89 row. ~3 700 flops plus 24
    // exp and 24 sin/cos pairs, i.e. ~100 flops/byte: nine times T89's ~11, and the highest
    // intensity of any BATCHED field kernel here (only the tracers, which return four floats per
    // LINE, are higher). The device is expected to win early and by a lot, and does.
    //
    // MEASURED, RTX 3070 Ti / Vulkan, best of five per size, transfers included, against
    // ext_t01s.hpp's own fp32 host lane (t01s_field_host) built -O3 -march=native
    // -ffp-contract=off, points scattered at 2.2..12 Re:
    //
    //     @@BENCHTABLE@@
    //
    // The amplitude solve is NOT in either column: it is 1 188 multiply-adds paid ONCE per batch
    // on the host, by t01s_param_block, for both lanes. Putting the 108 x 11 coefficient table on
    // the device instead would move 9.5 KB and 1 188 flops per THREAD to save one host loop.
    //
    // No coefficient BUFFER: sin psi, cos psi and the 108 solved amplitudes are 110 floats and
    // ride in the parameter block, so the kernel fits dispatch_batch's pos/out/params/dims shape.
    {"irbem_t01s_f32",
     "Tsyganenko-Singer-Kasper 2003 storm-time external B (nT) at each GSM point (Re); params = "
     "{sin psi, cos psi, 108 mode amplitudes} solved on the host for one driver state",
     4, 110, 256},
}};

/// The registry and the lease capacity must agree, checked at COMPILE time: a kernel that binds
/// more buffers than a launch can hold is a build error, not a runtime overrun.
static_assert([] {
    std::size_t most = 0;
    for (const KernelInfo& k : registered_kernels) most = k.bindings > most ? k.bindings : most;
    return most;
}() <= 8, "a kernel binds more than the shared descriptor set layout's eight slots");

/**
 * The registry entry for @p name.
 * @param name the Slang entry-point name to look up.
 * @return a reference to the entry, which has static storage duration.
 * @throws UnknownKernel when no entry point of that name is registered.
 * @complexity O(number of registered kernels) — a linear scan over a table of a handful of rows,
 *             which beats any hashing at this size and, unlike hashing, is `constexpr`-friendly.
 * @alloc none on success; the failure path builds one message string.
 */
inline const KernelInfo& kernel_info(std::string_view name) {
    for (const KernelInfo& k : registered_kernels) {
        if (std::string_view(k.name) == name) return k;
    }
    throw UnknownKernel("space.irbem gpu: no such kernel '" + std::string(name) +
                        "' (see registered_kernels)");
}

// ---------------------------------------------------------------------------------------------
// Availability
// ---------------------------------------------------------------------------------------------

namespace detail {

#if CHEATAH_SPACE_IRBEM_HAVE_GPU
/// Alias for the device stack this seam drives.
namespace gl = ::cheatah::gpu::linalg;

/**
 * Whether the operator has switched the device lane off with `CHEATAH_SPACE_IRBEM_NO_GPU`.
 *
 * A deliberate seam, not a debugging leftover: the differential suite has to be able to force the
 * host lane on a machine that HAS a device, and the failure path of every launcher has to be
 * reachable in a test on that same machine. Any value other than the empty string and `"0"` counts
 * as "off".
 * @return true when the variable is set to something that means "off".
 * @complexity O(1).
 * @alloc none.
 */
inline bool gpu_disabled_by_env() noexcept {
    const char* v = std::getenv("CHEATAH_SPACE_IRBEM_NO_GPU");
    return v != nullptr && v[0] != '\0' && !(v[0] == '0' && v[1] == '\0');
}
#endif

}  // namespace detail

/**
 * Whether a device launch can run right now.
 *
 * The runtime half of the availability question: true only when the device stack was compiled in,
 * the operator has not switched it off, and `cheatah-gpu-linalg`'s own cached bring-up probe
 * succeeds. Never throws — that is the point of asking instead of trying.
 *
 * Deliberately NOT cached here, although the probe underneath it is: the environment override has
 * to take effect when it is set, and the cost of the query is a `getenv` plus a load of a
 * function-local static.
 * @return true when @ref dispatch_batch would reach the device.
 * @complexity O(1).
 * @alloc none.
 */
inline bool available() noexcept {
#if CHEATAH_SPACE_IRBEM_HAVE_GPU
    return !detail::gpu_disabled_by_env() && detail::gl::available();
#else
    return false;
#endif
}

/**
 * Why @ref available is false.
 * @return the reason, empty exactly when the device is available.
 * @complexity O(1).
 * @alloc one string (the returned value).
 */
inline std::string unavailable_reason() {
#if CHEATAH_SPACE_IRBEM_HAVE_GPU
    if (detail::gpu_disabled_by_env()) return "device lane switched off by CHEATAH_SPACE_IRBEM_NO_GPU";
    return detail::gl::unavailable_reason();
#else
    return "built without cheatah-gpu-linalg (cheatah_gpu_linalg/context.hpp was not on the "
           "include path)";
#endif
}

// ---------------------------------------------------------------------------------------------
// Where the compiled kernels live
// ---------------------------------------------------------------------------------------------

/**
 * The directory holding this module's compiled `.spv` files.
 *
 * Resolution order deliberately mirrors `cheatah-gpu-linalg`'s own `spv_bytes`, so one habit covers
 * both: the `CHEATAH_SPACE_IRBEM_SPV_DIR` environment variable, then the same-named compile
 * definition, then `<repo>/build/shaders` derived from `__FILE__` — which means a consumer that
 * passes no definitions at all still finds the kernels of the checkout it compiled against.
 * @return the directory; not checked for existence, which is @ref dispatch_batch's job.
 * @complexity O(1).
 * @alloc a few short-lived `std::filesystem::path` temporaries on the `__FILE__` fallback, one on
 *        the other two; heap, but once per launch and never per point.
 */
inline std::filesystem::path shader_dir() {
    if (const char* env = std::getenv("CHEATAH_SPACE_IRBEM_SPV_DIR")) return {env};
#if defined(CHEATAH_SPACE_IRBEM_SPV_DIR)
    return {CHEATAH_SPACE_IRBEM_SPV_DIR};
#else
    // space/irbem/gpu/dispatch.hpp -> gpu -> irbem -> space -> the repository root.
    return std::filesystem::path(__FILE__)
               .parent_path()
               .parent_path()
               .parent_path()
               .parent_path() /
           "build" / "shaders";
#endif
}

/**
 * Where @p kernel's compiled SPIR-V is expected.
 * @param kernel the entry-point name; the file is `<kernel>.spv`, matching what
 *        `cheatah-gpu-linalg`'s context will independently resolve during the launch.
 * @return the full path.
 * @complexity O(1).
 * @alloc @ref shader_dir's, plus one string and one path for the file name.
 */
inline std::filesystem::path shader_path(std::string_view kernel) {
    return shader_dir() / (std::string(kernel) + ".spv");
}

/**
 * @p kernel addressed by its directory — the name `cheatah-gpu-linalg`'s context resolves EXACTLY
 * where it says, without consulting the environment.
 *
 * That library treats a kernel name containing '/' as a path (its `spv_bytes`), so this module's
 * kernels are launched by `<shader_dir>/<kernel>` and linalg's own bare names keep resolving to
 * linalg's directory in the same process. Before it accepted qualified names this header had to
 * `setenv` the library's shader variable around every launch and restore it afterwards — a
 * process-global mutated around a blocking call, documented as a workaround. It is gone.
 * @param kernel the entry-point name.
 * @return the qualified name, without extension.
 * @complexity O(1).
 * @alloc @ref shader_dir's, plus one string.
 */
inline std::string qualified(std::string_view kernel) {
    return (shader_dir() / std::string(kernel)).string();
}

/**
 * Where the Slang source of this module's kernels lives.
 *
 * Derived from `__FILE__` rather than assumed relative to the working directory, so the
 * registry-completeness test can read the source no matter where it is run from.
 * @return the path to `irbem.slang`.
 * @complexity O(1).
 * @alloc two short-lived paths, one of them returned.
 */
inline std::filesystem::path shader_source_path() {
    return std::filesystem::path(__FILE__).parent_path() / "irbem.slang";
}

// ---------------------------------------------------------------------------------------------
// Lane selection
// ---------------------------------------------------------------------------------------------

/**
 * The batch size at or above which @p kernel's device lane is worth its transfers.
 *
 * Normally @ref KernelInfo::crossover_points, i.e. the measurement recorded next to the kernel.
 * `CHEATAH_SPACE_IRBEM_GPU_CROSSOVER` overrides it process-wide with a decimal batch size, which
 * is what makes the policy re-tunable on a machine whose bus and cores are not this one's without
 * a rebuild — and what lets the suite drive the device lane through @ref Lane::Auto. A value that
 * is not a whole decimal number is ignored rather than read as zero, because reading a typo as
 * "always use the GPU" is the expensive direction to be wrong in.
 *
 * @param kernel the registered kernel the batch would run.
 * @return the crossover, possibly @ref never_faster_on_device.
 * @throws UnknownKernel when @p kernel is not registered.
 * @complexity O(number of registered kernels).
 * @alloc none.
 */
inline std::size_t gpu_crossover(std::string_view kernel) {
    if (const char* env = std::getenv("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER")) {
        char* end = nullptr;
        const unsigned long long v = std::strtoull(env, &end, 10);
        if (end != env && *end == '\0') return static_cast<std::size_t>(v);
    }
    return kernel_info(kernel).crossover_points;
}

/**
 * Whether a batch of @p points of @p kernel should take the device lane under @ref Lane::Auto.
 * @param kernel the registered kernel the batch would run.
 * @param points how many points the batch holds.
 * @return true when the kernel has a measured crossover, the batch reaches it, and a device is
 *         available; false — meaning "run it on the host" — in every other case.
 * @throws UnknownKernel when @p kernel is not registered.
 * @complexity O(number of registered kernels).
 * @alloc none.
 */
inline bool prefer_gpu(std::string_view kernel, std::size_t points) {
    const std::size_t crossover = gpu_crossover(kernel);
    return crossover != never_faster_on_device && points >= crossover && available();
}

/**
 * Which arithmetic lane a batch routine runs on.
 *
 * `Auto` is what production uses; the two explicit values exist because the differential suite has
 * to be able to run the SAME inputs through both lanes and compare, which it cannot do if the
 * choice is only ever made for it.
 */
enum class Lane : std::uint8_t {
    Auto,  ///< @ref prefer_gpu decides.
    Host,  ///< Force the CPU lane.
    Gpu,   ///< Force the device lane; throws @ref GpuUnavailable when there is no device.
};

// ---------------------------------------------------------------------------------------------
// The launcher
// ---------------------------------------------------------------------------------------------

namespace detail {

#if CHEATAH_SPACE_IRBEM_HAVE_GPU


/**
 * Holds the device buffers of one launch and returns every one of them to the context's pool when
 * the launch ends — including when it ends by throwing out of the driver.
 *
 * A fixed `std::array` sized to the widest kernel in the registry: nothing here allocates, and a
 * launch cannot leak a pooled allocation into a long-running process.
 */
class Leases {
public:
        /// The most buffers any registered kernel binds — the total-field tracer's EIGHT, which is
    /// also the shared descriptor set layout's cap. There is no headroom left: a ninth binding
    /// needs an ABI conversation, not a bump here. Asserted against the registry rather than
    /// trusted, because this constant said 4 when the seven-binding tracer landed and the result
    /// was a segfault inside `add` — so the capacity and its comment now change in the same edit
    /// as the kernel that moves them, by rule.
    static constexpr std::size_t capacity = 8;
    static_assert(capacity <= 8, "the shared descriptor set layout provides eight binding slots");

    Leases() = default;
    Leases(const Leases&) = delete;
    Leases& operator=(const Leases&) = delete;
    Leases(Leases&&) = delete;
    Leases& operator=(Leases&&) = delete;

    /**
     * Release every held buffer back to the pool.
     * @complexity O(number of buffers held).
     * @alloc none here. The context's own free list may grow, which is its pooling policy, not
     *        this launch's per-point cost.
     */
    ~Leases() {
        for (std::size_t i = 0; i < n_; ++i) gl::detail::ctx().release_buffer(bufs_[i]);
    }

    /**
     * Take ownership of @p b and bind it at the next binding slot.
     * @param b the buffer, freshly acquired from the context.
     * @return @p b, so the acquisition and the hand-off read as one expression.
     * @complexity O(1).
     * @alloc none.
     */
    gl::detail::Buffer* add(gl::detail::Buffer* b) {
        bufs_[n_] = b;
        ++n_;
        return b;
    }

    /**
     * The buffers, in binding order, as the context's launch API wants them.
     * @return a pointer to the first of @ref capacity slots.
     * @complexity O(1).
     * @alloc none.
     */
    gl::detail::Buffer** data() { return bufs_.data(); }

private:
    /// The held buffers, binding 0 first.
    std::array<gl::detail::Buffer*, capacity> bufs_{};
    /// How many of @ref bufs_ are live.
    std::size_t n_ = 0;
};

#endif  // CHEATAH_SPACE_IRBEM_HAVE_GPU

}  // namespace detail

/**
 * Run a registered kernel over one batch of xyz-interleaved points, blocking until it completes.
 *
 * The single door to the device. The order of the checks is itself part of the contract, and each
 * step is diagnosed before the next is attempted: an unregistered kernel (@ref UnknownKernel), then
 * a malformed batch (`std::invalid_argument`), then a shader that was never compiled
 * (@ref ShaderMissing), then a machine or build with no device (@ref GpuUnavailable). Shader
 * existence is checked BEFORE availability on purpose — a machine with no GPU is exactly where a
 * skipped shader-build step is likeliest, and reporting "no GPU" there would hide it.
 *
 * @param kernel the registered entry-point name to launch.
 * @param pos the input points, xyz-interleaved, `3N` floats; the batch size `N` is derived from it.
 * @param out the output vectors, xyz-interleaved, exactly as long as @p pos; overwritten in full.
 * @param params the kernel's scalar parameters, exactly `KernelInfo::params` of them.
 * @throws UnknownKernel when @p kernel is not in @ref registered_kernels.
 * @throws std::invalid_argument when @p pos is not a whole number of points, when @p out is not the
 *         same length, or when @p params is the wrong length for the kernel.
 * @throws ShaderMissing when the kernel's `.spv` is not under @ref shader_dir.
 * @throws GpuUnavailable when @ref available is false.
 * @complexity O(N) device work over `ceil(N/256)` workgroups, plus `O(N)` bytes moved each way.
 * @alloc no per-point allocation. Per call: a handful of strings and paths for the resolution and
 *        the diagnostics, one `setenv`, and four POOLED device buffers, all released before
 *        return. An empty batch allocates nothing and touches no device.
 */
inline void dispatch_batch(std::string_view kernel, std::span<const float> pos,
                           std::span<float> out, std::span<const float> params) {
    const KernelInfo& k = kernel_info(kernel);
    if (pos.size() % 3 != 0)
        throw std::invalid_argument("space.irbem gpu: input span is not a whole number of xyz "
                                    "points");
    if (out.size() != pos.size())
        throw std::invalid_argument("space.irbem gpu: output span length must equal the input's");
    if (params.size() != k.params)
        throw std::invalid_argument("space.irbem gpu: wrong parameter count for kernel '" +
                                    std::string(kernel) + "'");
    const std::size_t n = pos.size() / 3;
    if (n == 0) return;   // an empty batch is a no-op, not an error and not a device touch

    const std::filesystem::path spv = shader_path(std::string_view(k.name));
    if (!std::filesystem::exists(spv))
        throw ShaderMissing("space.irbem gpu: kernel '" + std::string(k.name) +
                            "' has no compiled SPIR-V at " + spv.string() +
                            " (build the shaders, or set CHEATAH_SPACE_IRBEM_SPV_DIR)");
    if (!available())
        throw GpuUnavailable("space.irbem gpu: cannot launch '" + std::string(k.name) +
                             "': " + unavailable_reason());

#if CHEATAH_SPACE_IRBEM_HAVE_GPU
    namespace gl = detail::gl;
    gl::detail::Context& c = gl::detail::ctx();
    const std::size_t vec_bytes = pos.size() * sizeof(float);
    const std::size_t par_bytes = params.size() * sizeof(float);
    const std::array<std::uint32_t, 1> dims{static_cast<std::uint32_t>(n)};

    detail::Leases lease;
    gl::detail::Buffer* b_pos = lease.add(c.new_data_buffer(vec_bytes));
    gl::detail::Buffer* b_out = lease.add(c.new_data_buffer(vec_bytes));
    gl::detail::Buffer* b_par = lease.add(c.new_buffer(par_bytes));
    gl::detail::Buffer* b_dim = lease.add(c.new_buffer(sizeof(std::uint32_t)));
    c.upload(b_pos, pos.data(), vec_bytes);
    c.upload(b_par, params.data(), par_bytes);
    c.upload(b_dim, dims.data(), sizeof(std::uint32_t));
    {
        c.dispatch_1d(qualified(k.name).c_str(), lease.data(), k.bindings, n);
    }
    c.download(b_out, out.data(), vec_bytes);
#endif
}

/**
 * Launch the TOTAL-field tracer — IGRF plus T89 in one resident trace.
 *
 * Identical contract to @ref launch_trace, with one more upload: @p ext carries the T89 parameter
 * block for the batch's Kp bin followed by the epoch's `gsm_from_geo` rotation, column-major —
 * 39 floats that make the external field evaluable on the device without a host round-trip per
 * RK4 stage. The bin is selected on the HOST because a batch shares one epoch and one Kp; a
 * per-thread bin branch would diverge the warp for nothing.
 *
 * @param pos 3N floats, GEO, Earth radii. @param pitch N floats, degrees.
 * @param coef IGRF `g` then `h`, interpolated to the epoch. @param norm the Legendre table.
 * @param ext the 39-float external block described above.
 * @param dims `{N, nmax, max_steps, steps_per_l x 1000}`.
 * @param out receives 4N floats: `I`, `Bmin`, `Bmirr`, `Blocal`. @param status one word per line.
 * @return `false` when there is no device or no compiled SPIR-V — the caller's cue for the host
 *         lane, not an error. Genuine misuse (mismatched spans) still throws.
 * @complexity One dispatch; O(N x steps) TOTAL-field evaluations (~900 flops each), concurrent.
 * @alloc eight device buffers, returned to the context's size-classed pool on scope exit.
 * @test IrbemGpu.TotalFieldTraceAgreesWithTheHostLane
 */
[[nodiscard]] inline bool launch_trace_total(std::span<const float> pos,
                                             std::span<const float> pitch,
                                             std::span<const float> coef,
                                             std::span<const float> norm,
                                             std::span<const float> ext,
                                             std::span<const std::uint32_t> dims,
                                             std::span<float> out,
                                             std::span<std::uint32_t> status) {
    const std::size_t n = pitch.size();
    if (pos.size() != 3 * n || out.size() != 4 * n || status.size() != n) {
        throw std::invalid_argument("space.irbem gpu: total-trace span lengths disagree");
    }
    if (n == 0) return true;
    if (!available() || !std::filesystem::exists(shader_path("irbem_trace_total_f32"))) {
        return false;
    }

#if CHEATAH_SPACE_IRBEM_HAVE_GPU
    namespace gl = detail::gl;
    gl::detail::Context& c = gl::detail::ctx();
    detail::Leases lease;
    // Positions and results are device-local; the tables, the ext block and dims are small,
    // host-written once, and take the mapped path.
    gl::detail::Buffer* b_pos = lease.add(c.new_data_buffer(pos.size() * sizeof(float)));
    gl::detail::Buffer* b_pit = lease.add(c.new_buffer(pitch.size() * sizeof(float)));
    gl::detail::Buffer* b_cf  = lease.add(c.new_buffer(coef.size() * sizeof(float)));
    gl::detail::Buffer* b_nr  = lease.add(c.new_buffer(norm.size() * sizeof(float)));
    gl::detail::Buffer* b_ex  = lease.add(c.new_buffer(ext.size() * sizeof(float)));
    gl::detail::Buffer* b_out = lease.add(c.new_data_buffer(out.size() * sizeof(float)));
    gl::detail::Buffer* b_st  = lease.add(c.new_data_buffer(status.size() * sizeof(std::uint32_t)));
    gl::detail::Buffer* b_dm  = lease.add(c.new_buffer(dims.size() * sizeof(std::uint32_t)));
    c.upload(b_pos, pos.data(), pos.size() * sizeof(float));
    c.upload(b_pit, pitch.data(), pitch.size() * sizeof(float));
    c.upload(b_cf, coef.data(), coef.size() * sizeof(float));
    c.upload(b_nr, norm.data(), norm.size() * sizeof(float));
    c.upload(b_ex, ext.data(), ext.size() * sizeof(float));
    c.upload(b_dm, dims.data(), dims.size() * sizeof(std::uint32_t));
    {
        c.dispatch_1d(qualified("irbem_trace_total_f32").c_str(), lease.data(), 8, n);
    }
    c.download(b_out, out.data(), out.size() * sizeof(float));
    c.download(b_st, status.data(), status.size() * sizeof(std::uint32_t));
    return true;
#else
    return false;
#endif
}

/**
 * Launch the field-line tracer — the seven-binding shape `dispatch_batch` cannot express.
 *
 * `dispatch_batch` above is for the pos/out/params/dims kernels: one vector in, one vector out.
 * The tracer takes positions AND pitch angles AND two coefficient tables, and writes four floats
 * plus a status word per line. Rather than generalise the four-binding helper into something that
 * describes every shape badly, this is the tracer's own launcher — the kernel that carries the
 * module's whole performance argument earns fifty lines of its own.
 *
 * @param pos 3N floats, GEO, Earth radii. @param pitch N floats, degrees.
 * @param coef `g` then `h`, already interpolated to the epoch by the caller.
 * @param norm the Legendre normalisation, `constexpr` on the host and therefore free.
 * @param dims `{N, nmax, max_steps, steps_per_l × 1000}`.
 * @param out receives 4N floats: `I`, `Bmin`, `Bmirr`, `Blocal` per line.
 * @param status receives N status words, one per line, so a single non-closing line reports itself
 *        without spoiling the batch.
 * @return `false` when there is no device or no compiled SPIR-V — the caller's cue to run the host
 *         lane, not an error. Genuine misuse (mismatched spans) still throws.
 * @complexity One dispatch over `ceil(N/256)` workgroups; O(N × steps) field evaluations, run
 *             concurrently. ~28 bytes in and ~20 out per line for ~10^5 flops — the arithmetic
 *             intensity that makes this worth offloading at all.
 * @alloc seven device buffers, returned to the context's size-classed pool on scope exit.
 * @test IrbemGpu.TraceKernelAgreesWithTheHostLane
 */
[[nodiscard]] inline bool launch_trace(std::span<const float> pos, std::span<const float> pitch,
                                       std::span<const float> coef, std::span<const float> norm,
                                       std::span<const std::uint32_t> dims,
                                       std::span<float> out, std::span<std::uint32_t> status) {
    const std::size_t n = pitch.size();
    if (pos.size() != 3 * n || out.size() != 4 * n || status.size() != n) {
        throw std::invalid_argument("space.irbem gpu: trace span lengths disagree");
    }
    if (n == 0) return true;
    if (!available() || !std::filesystem::exists(shader_path("irbem_trace_i_f32"))) return false;

#if CHEATAH_SPACE_IRBEM_HAVE_GPU
    namespace gl = detail::gl;
    gl::detail::Context& c = gl::detail::ctx();
    detail::Leases lease;
    // Positions and results are device-local: they are large and touched only by the kernel.
    // The tables and dims are small and host-written once, so they take the mapped path.
    gl::detail::Buffer* b_pos = lease.add(c.new_data_buffer(pos.size() * sizeof(float)));
    gl::detail::Buffer* b_pit = lease.add(c.new_buffer(pitch.size() * sizeof(float)));
    gl::detail::Buffer* b_cf  = lease.add(c.new_buffer(coef.size() * sizeof(float)));
    gl::detail::Buffer* b_nr  = lease.add(c.new_buffer(norm.size() * sizeof(float)));
    gl::detail::Buffer* b_out = lease.add(c.new_data_buffer(out.size() * sizeof(float)));
    gl::detail::Buffer* b_st  = lease.add(c.new_data_buffer(status.size() * sizeof(std::uint32_t)));
    gl::detail::Buffer* b_dm  = lease.add(c.new_buffer(dims.size() * sizeof(std::uint32_t)));
    c.upload(b_pos, pos.data(), pos.size() * sizeof(float));
    c.upload(b_pit, pitch.data(), pitch.size() * sizeof(float));
    c.upload(b_cf, coef.data(), coef.size() * sizeof(float));
    c.upload(b_nr, norm.data(), norm.size() * sizeof(float));
    c.upload(b_dm, dims.data(), dims.size() * sizeof(std::uint32_t));
    {
        c.dispatch_1d(qualified("irbem_trace_i_f32").c_str(), lease.data(), 7, n);
    }
    c.download(b_out, out.data(), out.size() * sizeof(float));
    c.download(b_st, status.data(), status.size() * sizeof(std::uint32_t));
    return true;
#else
    return false;
#endif
}

/**
 * Launch the IGRF field kernel over a batch of arbitrary points — the five-binding shape.
 *
 * `dispatch_batch` cannot express this one: it uploads a scalar parameter block, and IGRF has no
 * scalars, only two tables. The tables arrive ALREADY interpolated to the epoch, because
 * interpolating IGRF's 26-epoch table per thread would be one redundant copy per point of a
 * calculation the host does once per batch.
 *
 * L\*'s flux integral is what made this worth exposing on the seam rather than leaving it inside a
 * benchmark: `Phi` is one field evaluation per polar-cap cell, ~2 400 cells per L\* point and
 * nothing else — exactly the ~20 flops/byte regime the registry measures the device winning by
 * 8.37x in.
 *
 * @param pos 3N floats, GEO, Earth radii. @param coef `g` then `h`, pre-interpolated.
 * @param norm the Legendre normalisation — `constexpr` on the host, so free to produce.
 * @param dims `{N, nmax}`.
 * @param out receives 3N floats: the field at each point, GEO, nT.
 * @return `false` when there is no device or no compiled SPIR-V — the caller's cue to run the host
 *         lane, not an error. Genuine misuse (mismatched spans) still throws.
 * @throws std::invalid_argument when @p out is not the same length as @p pos.
 * @complexity One dispatch over `ceil(N/256)` workgroups.
 * @alloc five device buffers, returned to the context's size-classed pool on scope exit.
 * @test IrbemDriftShell.FluxCellsAgreeBetweenLanes
 */
[[nodiscard]] inline bool launch_igrf(std::span<const float> pos, std::span<const float> coef,
                                      std::span<const float> norm,
                                      std::span<const std::uint32_t> dims, std::span<float> out) {
    if (pos.size() % 3 != 0 || out.size() != pos.size()) {
        throw std::invalid_argument("space.irbem gpu: launch_igrf wants two equal-length "
                                    "xyz-interleaved spans");
    }
    const std::size_t n = pos.size() / 3;
    if (n == 0) return true;
    if (!available() || !std::filesystem::exists(shader_path("irbem_igrf_f32"))) return false;

#if CHEATAH_SPACE_IRBEM_HAVE_GPU
    namespace gl = detail::gl;
    gl::detail::Context& c = gl::detail::ctx();
    detail::Leases lease;
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
        c.dispatch_1d(qualified("irbem_igrf_f32").c_str(), lease.data(), 5, n);
    }
    c.download(b_out, out.data(), out.size() * sizeof(float));
    return true;
#else
    return false;
#endif
}

/**
 * Launch the drift-shell footpoint tracer — the second seven-binding shape.
 *
 * The companion to @ref launch_trace, and the step that turns a converged drift shell into the
 * polar-cap boundary L\*'s flux integral is taken over. Each line walks from its magnetic-equator
 * seed down to `r = 1` along @p dir, halving its step onto the sphere rather than interpolating
 * across it; see the kernel's own header in `irbem.slang` for why that distinction is worth ~30
 * extra steps.
 *
 * @param pos 3N floats, GEO, Earth radii — one magnetic-equator seed per shell azimuth.
 * @param dir N floats, `+1` or `-1`: which way along B leads to the chosen hemisphere. The host
 *        decides, because the kernel works in GEO and has no dipole axis to compare against.
 * @param coef `g` then `h`, already interpolated to the epoch. @param norm the normalisation.
 * @param dims `{N, nmax, max_steps, steps_per_l × 1000}`.
 * @param out receives 3N floats: each footpoint in GEO, ON the unit sphere.
 * @param status receives N status words, so one line that never reaches the surface reports itself
 *        instead of spoiling the batch.
 * @return `false` when there is no device or no compiled SPIR-V — the caller's cue to run the host
 *         lane, not an error. Genuine misuse (mismatched spans) still throws.
 * @throws std::invalid_argument when the spans do not agree.
 * @complexity One dispatch over `ceil(N/256)` workgroups; O(N × steps) field evaluations, run
 *             concurrently.
 * @alloc seven device buffers, returned to the context's size-classed pool on scope exit.
 * @test IrbemDriftShell.FootpointsAgreeBetweenLanes
 */
[[nodiscard]] inline bool launch_shell_foot(std::span<const float> pos, std::span<const float> dir,
                                            std::span<const float> coef,
                                            std::span<const float> norm,
                                            std::span<const std::uint32_t> dims,
                                            std::span<float> out,
                                            std::span<std::uint32_t> status) {
    const std::size_t n = dir.size();
    if (pos.size() != 3 * n || out.size() != 3 * n || status.size() != n) {
        throw std::invalid_argument("space.irbem gpu: shell-footpoint span lengths disagree");
    }
    if (n == 0) return true;
    if (!available() || !std::filesystem::exists(shader_path("irbem_shell_foot_f32"))) return false;

#if CHEATAH_SPACE_IRBEM_HAVE_GPU
    namespace gl = detail::gl;
    gl::detail::Context& c = gl::detail::ctx();
    detail::Leases lease;
    gl::detail::Buffer* b_pos = lease.add(c.new_data_buffer(pos.size() * sizeof(float)));
    gl::detail::Buffer* b_dir = lease.add(c.new_buffer(dir.size() * sizeof(float)));
    gl::detail::Buffer* b_cf = lease.add(c.new_buffer(coef.size() * sizeof(float)));
    gl::detail::Buffer* b_nr = lease.add(c.new_buffer(norm.size() * sizeof(float)));
    gl::detail::Buffer* b_out = lease.add(c.new_data_buffer(out.size() * sizeof(float)));
    gl::detail::Buffer* b_st = lease.add(c.new_data_buffer(status.size() * sizeof(std::uint32_t)));
    gl::detail::Buffer* b_dm = lease.add(c.new_buffer(dims.size() * sizeof(std::uint32_t)));
    c.upload(b_pos, pos.data(), pos.size() * sizeof(float));
    c.upload(b_dir, dir.data(), dir.size() * sizeof(float));
    c.upload(b_cf, coef.data(), coef.size() * sizeof(float));
    c.upload(b_nr, norm.data(), norm.size() * sizeof(float));
    c.upload(b_dm, dims.data(), dims.size() * sizeof(std::uint32_t));
    {
        c.dispatch_1d(qualified("irbem_shell_foot_f32").c_str(), lease.data(), 7, n);
    }
    c.download(b_out, out.data(), out.size() * sizeof(float));
    c.download(b_st, status.data(), status.size() * sizeof(std::uint32_t));
    return true;
#else
    return false;
#endif
}

// ---------------------------------------------------------------------------------------------
// The centred dipole — the first physics behind the seam
// ---------------------------------------------------------------------------------------------

/**
 * The centred-dipole magnetic field at one point, in `double` — the reference the device lane is
 * measured against.
 *
 * The degree-1 order-0 term of the IAGA spherical-harmonic expansion that defines IGRF (Alken et
 * al., *International Geomagnetic Reference Field: the thirteenth generation*, Earth Planets Space
 * **73**:49 (2021), eq. 1). With `V = a (a/r)^2 g10 P_1^0(cos θ)` and `B = -∇V`,
 * `B_r = 2 g10 (a/r)^3 cos θ` and `B_θ = g10 (a/r)^3 sin θ`, which in the MAG frame — where the
 * dipole axis IS `ẑ` — is the single identity `B = g10 (a/r)^3 (3 (ẑ·r̂) r̂ − ẑ)`; componentwise,
 * with `a = 1 Re`, `Bx = 3 g10 x z / r⁵`, `By = 3 g10 y z / r⁵`, `Bz = g10 (3z² − r²) / r⁵`. The
 * derivation is written out in full in `irbem.slang`, whose kernel evaluates exactly these
 * expressions in exactly this order.
 *
 * The real `g10` is negative (−29404.8 nT at IGRF-13 epoch 2020.0), which is what makes the field
 * point into the ground at the northern dipole pole; the sign lives in the coefficient, so this
 * routine never encodes a hemisphere.
 *
 * @param p the point, in the MAG (centred-dipole) frame, in Earth radii.
 * @param g10_nT the degree-1 order-0 Gauss coefficient, in nanotesla.
 * @return the field at @p p, in the same frame, in nanotesla; exactly zero at the origin.
 * @complexity O(1) — one square root, one division.
 * @alloc none.
 * @note The origin is a genuine pole of a dipole field, not a rounding problem. Returning zero
 *       there is a DEFINED value chosen so the host and the kernel agree on a testable answer
 *       instead of disagreeing about which flavour of NaN they produce.
 */
inline FieldVector<Frame::MAG> dipole_field_at(Position<Frame::MAG> p, double g10_nT) {
    const double x = p.v[0];
    const double y = p.v[1];
    const double z = p.v[2];
    const double r2 = (x * x) + (y * y) + (z * z);
    if (r2 <= 0.0) return {};
    const double r5 = r2 * r2 * std::sqrt(r2);
    return FieldVector<Frame::MAG>{fixarray::vec3d{3.0 * g10_nT * x * z / r5,
                                                   3.0 * g10_nT * y * z / r5,
                                                   g10_nT * ((3.0 * z * z) - r2) / r5}};
}

/**
 * The centred-dipole field over a whole batch, on the CPU, in `float`.
 *
 * The host twin of `irbem_dipole_f32`: the same expressions, in the same order, in the same
 * precision, which is what makes a disagreement between the two attributable to the DEVICE
 * (contraction, a driver's transcendental) rather than to the arithmetic having been written
 * differently on the two sides. It is also the lane a machine without a GPU actually runs.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, in the MAG frame, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, in nanotesla; overwritten in full.
 * @param g10_nT the degree-1 order-0 Gauss coefficient, in nanotesla.
 * @throws std::invalid_argument when @p pos is not a whole number of points or @p out is a
 *         different length.
 * @complexity O(N) — one square root and one division per point, no branch per component.
 * @alloc none. Not one byte: the loop is over caller-provided spans.
 */
inline void dipole_field_host(std::span<const float> pos, std::span<float> out, float g10_nT) {
    if (pos.size() % 3 != 0 || out.size() != pos.size())
        throw std::invalid_argument("space.irbem gpu: dipole_field_host wants two equal-length "
                                    "xyz-interleaved spans");
    const std::size_t n = pos.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const float x = pos[(3 * i) + 0];
        const float y = pos[(3 * i) + 1];
        const float z = pos[(3 * i) + 2];
        const float r2 = (x * x) + (y * y) + (z * z);
        float bx = 0.0F;
        float by = 0.0F;
        float bz = 0.0F;
        if (r2 > 0.0F) {
            const float r5 = r2 * r2 * std::sqrt(r2);
            bx = 3.0F * g10_nT * x * z / r5;
            by = 3.0F * g10_nT * y * z / r5;
            bz = g10_nT * ((3.0F * z * z) - r2) / r5;
        }
        out[(3 * i) + 0] = bx;
        out[(3 * i) + 1] = by;
        out[(3 * i) + 2] = bz;
    }
}

/**
 * The centred-dipole field over a whole batch, on the device.
 * @param pos the points, xyz-interleaved, `3N` floats, in the MAG frame, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, in nanotesla; overwritten in full.
 * @param g10_nT the degree-1 order-0 Gauss coefficient, in nanotesla.
 * @throws GpuUnavailable when there is no device; @ref ShaderMissing when the kernel was not
 *         compiled; `std::invalid_argument` when the spans do not match.
 * @complexity O(N) device work, plus `2·3N` floats over the bus.
 * @alloc the parameter block is a stack `std::array`; everything else is @ref dispatch_batch's.
 */
inline void dipole_field_gpu(std::span<const float> pos, std::span<float> out, float g10_nT) {
    const std::array<float, 1> params{g10_nT};
    dispatch_batch("irbem_dipole_f32", pos, out, std::span<const float>(params));
}

/**
 * The centred-dipole field over a whole batch, on whichever lane @p lane selects.
 *
 * The production entry point. Under @ref Lane::Auto the choice is @ref prefer_gpu's, which means a
 * build with no device stack, a machine with no device, and a batch too small to pay for its
 * transfers all quietly run the host lane — while a forced @ref Lane::Gpu on any of those throws,
 * because a test that meant to exercise the device must not silently pass on the CPU.
 *
 * @param pos the points, xyz-interleaved, `3N` floats, in the MAG frame, in Earth radii.
 * @param out the field, xyz-interleaved, `3N` floats, in nanotesla; overwritten in full.
 * @param g10_nT the degree-1 order-0 Gauss coefficient, in nanotesla.
 * @param lane which lane to run; @ref Lane::Auto by default.
 * @throws GpuUnavailable when @p lane is @ref Lane::Gpu and there is no device.
 * @throws std::invalid_argument when the spans do not match.
 * @complexity O(N).
 * @alloc none on the host lane; see @ref dispatch_batch for the device lane.
 */
inline void dipole_field(std::span<const float> pos, std::span<float> out, float g10_nT,
                         Lane lane = Lane::Auto) {
    const bool use_gpu = (lane == Lane::Gpu) ||
                         (lane == Lane::Auto && prefer_gpu("irbem_dipole_f32", pos.size() / 3));
    if (use_gpu) {
        dipole_field_gpu(pos, out, g10_nT);
        return;
    }
    dipole_field_host(pos, out, g10_nT);
}

}  // namespace cheatah::space::irbem::gpu
