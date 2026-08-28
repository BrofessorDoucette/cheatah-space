# space.irbem benchmarks

The generated benchmark tables for the [`space.irbem`](README.md) module: every routine that
exists, on every lane it has, against the vendored IRBEM library run as a black box. The module
page describes what the routines compute; this page holds the numbers.

Three things are being measured, and they answer three different questions:

1. **[Per routine](#every-routine-on-every-lane-it-has)** — how fast is each routine on the host,
   on the device, and in IRBEM? The last column of that table is the one to read first: it says
   whether a routine can reach the GPU **at all**, and it is read out of `gpu/dispatch.hpp` — the
   registry and its launchers — rather than asserted, so the porting gap is visible instead of
   implied.
2. **[The crossover](#the-crossover-where-the-device-starts-winning)** — at what batch size does the
   device start to pay? This is the module's central argument and the reason the batch entry points
   exist at all.
3. **[Arithmetic intensity](#arithmetic-intensity-the-number-that-predicts-the-verdict)** — *why*
   the answer to (2) differs so wildly between kernels. Same seam, same transfers, same device;
   only the flops per byte changed, and the verdict flips.

## The ground rules

**The IRBEM lane quotes the `-O2` rebuild, never the shipped binary.** IRBEM's own makefile passes
no `-O` at all, and `docs/ERROR_BUDGET.md` measures the shipped library at 2.7× the `-O2` one.
Benchmarking against the shipped artefact would inflate every ratio in these tables by that factor,
and the resulting claims would be about GNU Fortran's default optimization level rather than about
this code. The oracle is `dlopen`ed and called through the documented C bindings — run as a black
box, never read, never linked.

**Every rate is real time, never CPU time.** A device lane spends nearly all of its wall clock
blocked rather than running, so Google Benchmark's default CPU-time accounting would credit the
GPU with several times its actual throughput. Every batch benchmark is marked `UseRealTime`, so the
host and device columns of a ratio are the same quantity.

**Transfers are inside the timed region.** A kernel time that excludes the bus is a number about
the ALU, not about an operation a caller can perform. The device columns include staging the
inputs, the dispatch, and reading the results back.

**Nothing here is hand-typed.** [`scripts/bench_run.sh`](../../scripts/bench_run.sh) runs the
binary and regenerates the region between each suite's `BENCH` markers from that run's JSON. The
"on GPU?" column and the flops/byte column are printed by the benchmark binary itself, from the
kernel registry and from a measured step count — so they cannot drift away from the code they
describe.

**A ratio is only quoted once its device lane has been checked against its host lane.** For a long
time this page did not do that, and could not: a benchmark that calls `launch_igrf` and gets `true`
back has learned that a dispatch was submitted, not that the kernel unpacked the coefficient buffer
the way the host packed it — and a kernel that reads the normalisation table at the wrong stride is
exactly as fast as one that reads it correctly. `irbem_bench --verify` now runs every device lane
against its own host lane and reports the largest relative deviation against
[`docs/ERROR_BUDGET.md`](docs/ERROR_BUDGET.md) §5, and `bench_run.sh check` runs it *before* it
grades a single row. Measured here, over the same input sets the tables are measured on:

| lane | quantity | n | RTX 3070 Ti | Intel Xe | llvmpipe | budget |
|---|---|--:|--:|--:|--:|--:|
| `irbem_dipole_f32` | \|B\| vs the fp32 host twin | 262 144 | 5.32e-07 | 4.59e-07 | 0 | 1e-06 |
| `irbem_igrf_f32` | \|B\| vs the fp64 host lane | 65 536 | 9.22e-07 | 8.21e-07 | 9.22e-07 | 1e-06 |
| `irbem_trace_i_f32` | `I` vs the fp64 host lane | 4 096 | 5.87e-05 | 2.46e-05 | 7.12e-05 | 1e-04 |
| `irbem_trace_i_f32` | `Bmin` vs the fp64 host lane | 4 096 | 8.91e-07 | 8.91e-07 | 9.62e-07 | 1e-05 |

Every lane is inside budget on all three backends, and the NVIDIA figures repeat digit for digit
between runs. Two of them are not comfortably inside it: the IGRF magnitude sits at 92% of the
`Blocal` budget and llvmpipe's `I` at 71% of the `XJ` budget. That is what an fp32 kernel measured
against an fp64 reference looks like, and it is why those budgets are where they are rather than an
order tighter — but it also means a kernel change that costs half an ulp will be caught by this
gate rather than by a user.

**Never `-ffast-math`.** The suite builds `-O3 -march=native -ffp-contract=off`, which is the
arithmetic the correctness tests pinned. `-ffast-math` would let the compiler reassociate the
spherical-harmonic sums and the RK4 accumulation, and the binary would then be measuring different
mathematics from the one that was verified. `bench_run.sh check` fails the build if the flag ever
appears.

## Every routine, on every lane it has

Read the **on GPU?** column, not just the timings. Three of its four answers are failure states,
and the interesting one is `registered, NO LAUNCHER`: a kernel can be written, compiled, correct
and even measurable in isolation while remaining unreachable through the seam every caller goes
through. That is a porting gap, not a success.

The column is *read out of* `gpu/dispatch.hpp` — the registry for the binding count, and the
header's own `launch_*` functions for whether anything actually reaches the kernel. It used to be
read out of the registry plus a rule written into the benchmark ("four bindings means
`dispatch_batch`, `irbem_trace_i_f32` means `launch_trace`, everything else is unreachable"), which
is not a derivation at all: it is a second list, in a second file, and it went stale the way second
lists do. `launch_igrf` landed in the seam and this column went on publishing `registered, NO
LAUNCHER` about a kernel the seam had learned to launch. It now goes and looks.

Routines with an epoch-dependent geometry appear twice. `(hot)` reuses one `Rotations` across the
ephemeris — the shape this module's API is built for — while `(cold)` rebuilds it per point, which
is the like-for-like against IRBEM, whose entry points take a date with every call and cannot
amortize. Reporting only the hot number would flatter this module; reporting only the cold one
would hide its actual design advantage.

Measured by [`bench/irbem_bench.cpp`](../../bench/irbem_bench.cpp); reproduce with `scripts/bench_run.sh publish irbem-per-routine`.

<!-- BENCH:irbem-per-routine begin -->
<!-- cheatah-bench-stamp v1
     suite:        irbem-per-routine
     generated:    2026-08-28
     commit:       950e2e1
     host:         pop-os, 20 CPUs @ 400 MHz
     cpu-scaling:  powersave
     build:        g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0, Google Benchmark @ 192ef10025eb2c4cdd392bc502f0c852196baa48
     device:       NVIDIA
     competitors:  IRBEM /tmp/irbem-builds/libirbem-O2.so (the -O2 REBUILD, never the shipped no-`-O` library)
     harness:      reps=5, min_time=0.2s, random-interleaving=on
     statistic:    median over repetitions; every rate is REAL time, never CPU time — a device
                   lane is blocked, not running, for nearly all of its wall clock
     publishable:  true
     layout:       routines
     watch:        space/irbem/, bench/irbem_bench.cpp

     PRODUCED BY:
       SPACE_IRBEM_ORACLE='/tmp/irbem-builds/libirbem-O2.so' CHEATAH_GPU_LINALG_VK_DEVICE='NVIDIA' \
           scripts/bench_run.sh publish irbem-per-routine
       which runs: build/bench/irbem_bench --benchmark_filter='^BM_' \
           --benchmark_repetitions=5 --benchmark_min_time=0.2s \
           --benchmark_enable_random_interleaving=true \
           --benchmark_out_format=json --benchmark_out=docs/bench/irbem-per-routine.json
-->

| routine | CPU per call | GPU per call | GPU speedup | IRBEM `-O2` per call | IRBEM / best lane | on GPU? |
|---|--:|--:|--:|--:|--:|---|
| datetime: julian_day_number | 4.75 ns | — | — | 2.69 ns | 0.57× | no - host only |
| datetime: calendar_date | 6.27 ns | — | — | 10.10 ns | 1.61× | no - host only |
| datetime: day_of_year | 6.30 ns | — | — | 6.17 ns | 0.98× | no - host only |
| datetime: decimal_year | 10.32 ns | — | — | 10.14 ns | 0.98× | no - host only |
| datetime: date_and_time_from_decimal_year | 28.21 ns | — | — | 13.84 ns | 0.49× | no - host only |
| datetime: date_and_time_from_doy_and_ut | 20.18 ns | — | — | 17.58 ns | 0.87× | no - host only |
| rotations: Rotations::at (per epoch) | 202.66 ns | — | — | — | — | no - host only |
| transform: GEO->GSM (hot) | 1.10 ns | — | — | 445.64 ns | 406.30× | no - host only |
| transform: GEO->GSM (cold) | 202.89 ns | — | — | 445.64 ns | 2.20× | no - host only |
| transform: GEO->GSE (hot) | 1.10 ns | — | — | 445.17 ns | 403.65× | no - host only |
| transform: GEO->GSE (cold) | 202.24 ns | — | — | 445.17 ns | 2.20× | no - host only |
| transform: GEO->SM (hot) | 1.10 ns | — | — | 447.53 ns | 408.49× | no - host only |
| transform: GEO->SM (cold) | 201.93 ns | — | — | 447.53 ns | 2.22× | no - host only |
| transform: GEO->MAG (hot) | 1.15 ns | — | — | 312.08 ns | 271.98× | no - host only |
| transform: GEO->MAG (cold) | 203.51 ns | — | — | 312.08 ns | 1.53× | no - host only |
| transform: GEO->GEI (hot) | 1.26 ns ±9% | — | — | 445.53 ns | 353.58× | no - host only |
| transform: GEO->GEI (cold) | 201.24 ns | — | — | 445.53 ns | 2.21× | no - host only |
| transform: GSM->SM (hot) | 1.11 ns | — | — | 444.89 ns | 399.36× | no - host only |
| transform: GSM->SM (cold) | 203.86 ns | — | — | 444.89 ns | 2.18× | no - host only |
| transform: GSM->GEO, transposed (hot) | 1.25 ns | — | — | 445.54 ns | 356.84× | no - host only |
| transform: GSM->GEO, transposed (cold) | 203.11 ns | — | — | 445.54 ns | 2.19× | no - host only |
| geodetic: gdz_to_geo | 13.98 ns | — | — | 16.26 ns | 1.16× | no - host only |
| geodetic: geo_to_gdz (Bowring) | 84.98 ns | — | — | 120.32 ns | 1.42× | no - host only |
| geodetic: sph_to_car | 11.85 ns | — | — | 16.13 ns | 1.36× | no - host only |
| geodetic: car_to_sph | 24.85 ns | — | — | 21.32 ns | 0.86× | no - host only |
| geodetic: rll_to_gdz | 9.19 ns | — | — | 12.08 ns | 1.31× | no - host only |
| helio: helio_geometry (per epoch) | 120.73 ns | — | — | — | — | no - host only |
| helio: HAE->HEE (hot) | 1.51 ns | — | — | 163.04 ns | 107.80× | no - host only |
| helio: HAE->HEE (cold) | 95.20 ns | — | — | 163.04 ns | 1.71× | no - host only |
| helio: HEE->HAE (hot) | 1.53 ns ±6% | — | — | 163.11 ns | 106.50× | no - host only |
| helio: HEE->HAE (cold) | 95.11 ns | — | — | 163.11 ns | 1.71× | no - host only |
| helio: HAE->HEEQ (hot) | 1.52 ns ±6% | — | — | 160.92 ns | 105.83× | no - host only |
| helio: HAE->HEEQ (cold) | 120.81 ns | — | — | 160.92 ns | 1.33× | no - host only |
| helio: HEEQ->HAE (hot) | 1.52 ns | — | — | 160.94 ns | 105.74× | no - host only |
| helio: HEEQ->HAE (cold) | 122.95 ns | — | — | 160.94 ns | 1.31× | no - host only |
| helio: GSE->HEE, position (hot) | 0.49 ns | — | — | 162.43 ns | 334.55× | no - host only |
| helio: HEE->GSE, position (hot) | 0.50 ns | — | — | 162.69 ns | 326.99× | no - host only |
| igrf: Igrf<13>::at (per epoch) | 48.40 ns | — | — | — | — | no - host only |
| igrf: evaluate, degree 13, GEO | 307.63 ns | — | — | 778.71 ns | 2.53× | yes - irbem_igrf_f32 via launch_igrf |
| igrf: evaluate, batch of 65536 | 308.83 ns | 48.67 ns ±116% | 6.35× | 778.71 ns | 16.00× | yes - irbem_igrf_f32 via launch_igrf |
| dipole: batch, fp32 | 2.09 ns | 3.54 ns ±44% | 0.59× | 641.41 ns | 306.51× | yes - irbem_dipole_f32 via dispatch_batch |
| lstar: trace_invariant (one line) | 58.51 µs | — | — | 551.36 µs | 9.42× | yes - irbem_trace_i_f32 via launch_trace |
| lstar: trace batch of 65536 | 58.92 µs | 1.26 µs | 46.89× | 551.36 µs | 438.79× | yes - irbem_trace_i_f32 via launch_trace |
| T89: t89_field, Kp = 0 (quiet) | 63.69 ns | — | — | 893.50 ns | 14.03× | no - host only |
| T89: t89_field, Kp = 9- (extreme storm) | 64.53 ns | — | — | 876.27 ns | 13.58× | no - host only |
| lstar: mcilwain_l (Hilton) | 38.48 ns | — | — | — | — | no - host only |
| lstar: dipole_moment | 1.32 ns | — | — | — | — | no - host only |

<details><summary><b>Also measured</b> — variants, and routines this module has not ported</summary>

| routine | CPU per call | GPU per call | GPU speedup | IRBEM `-O2` per call | IRBEM / best lane | on GPU? |
|---|--:|--:|--:|--:|--:|---|
| igrf: evaluate, degree 10, GEO | 177.00 ns | — | — | — | — | yes - irbem_igrf_f32 via launch_igrf |
| igrf: evaluate, degree 13, spherical | 329.78 ns | — | — | — | — | yes - irbem_igrf_f32 via launch_igrf |
| dipole: scalar, fp64 | 3.08 ns | — | — | 641.41 ns | 208.37× | yes - irbem_dipole_f32 via dispatch_batch |
| helio: GSE->HEE, position (cold) | 83.10 ns | — | — | 162.43 ns | 1.95× | no - host only |
| T89: t89_field, Kp = 1.0 | 65.14 ns | — | — | 888.84 ns | 13.64× | no - host only |
| T89: t89_field, Kp = 2.0 | 65.48 ns | — | — | 877.24 ns | 13.40× | no - host only |
| T89: t89_field, Kp = 3.0 | 65.60 ns | — | — | 882.92 ns | 13.46× | no - host only |
| T89: t89_field, Kp = 4.0 | 65.74 ns | — | — | 879.82 ns | 13.38× | no - host only |
| T89: t89_field, Kp = 5.0 | 65.16 ns | — | — | 878.52 ns | 13.48× | no - host only |
| T89: t89_field, Kp = 6.0 | 64.80 ns | — | — | 877.42 ns | 13.54× | no - host only |
| T89: t89_field, Kp = 7.0 | 64.29 ns | — | — | 882.90 ns | 13.73× | no - host only |
| T89: t89_field, Kp = 8.0 | 64.31 ns | — | — | 877.86 ns | 13.65× | no - host only |
| lstar: trace at alpha = 90 deg only (locally mirroring) | — | — | — | 15.60 µs | — | no - host only |
| lstar: L* drift shell (ported; NOT YET BENCHMARKED here) | — | — | — | 13874 µs | — | no - host only |

</details>
<!-- BENCH:irbem-per-routine end -->

### Four rows that need reading carefully

**The trace comparison is against `MAKE_LSTAR_SHELL_SPLITTING`, and that is not a detail.** IRBEM's
`MAKE_LSTAR` is the obvious counterpart and is the one this page quoted until it was checked. It is
the wrong one. IRBEM's own documentation says `MAKE_LSTAR` "computes the L\* parameter for locally
mirroring particles (local pitch angle of 90 degrees)" and points the reader at
`MAKE_LSTAR_SHELL_SPLITTING` "to compute L\* for arbitrary pitch angles"
(`docs/source/api/magnetic_coordinates.rst`). A locally mirroring particle mirrors *at* the
spacecraft, so its `I` spans the short arc between the point and its conjugate. This suite's lines
mirror at `B_local / sin²α` for α = 30…80° — up to four times the local field, and a far longer arc
to integrate. Called at the same eight start points, the two entry points return `XJ` = 1.85 Re and
5.5 × 10⁻³ Re at L = 2: two different integrals, not two implementations of one. They differ by 35×
in time, and quoting the 90° one turned this module's trace win into a 0.27× loss. The 90° entry
point is still measured — in **Also measured**, labelled for what it is — and at α = 90° the
shell-splitting entry point reproduces its `Lm` and `XJ` to every printed digit while costing about
twice as much, so roughly a 2× constant of its own sits inside the row this page now quotes.

**What survives that correction is the smaller caveat: the two sides still trace at their own
default step sizes.** This module's `TraceOptions::steps_per_l = 50`, i.e. `ds = L/50`, against
whatever `make_lstar_shell_splitting1` chooses for itself. So the trace rows are a statement about
each library's defaults as much as about its inner loop, and matching the two step counts is still
owed before either is quoted as an inner-loop comparison.

**`L* drift shell` has no CPU cell at all.** `driftshell.hpp::make_lstar` landed while this suite
was being written, so the row no longer says NOT PORTED — but nothing here times it yet, and an
empty CPU cell beside IRBEM's cost is the honest rendering of "ported, not yet measured". A CPU and
a device lane for it are owed.

**Two device cells carry a wide spread.** The dipole and IGRF device lanes finish in microseconds,
so they sit close to the dispatch floor and contend with whatever else holds the queue. Their
medians are directionally right — the dipole loses, the IGRF wins — but neither is a number to build
a budget on. Across four independent runs of this suite on this machine the dipole device lane
measured 3.50, 4.18, 6.42 and 14.01 ns/point and the IGRF device lane 7.09, 7.65, 12.06 and 48.67
(the last at ±116%), and the published speedup for the IGRF batch has read anywhere between 6.35×
and 43×. Neither verdict has ever changed sign. Read the verdict and the spread; do not read the
median as a budget. This is also why `bench_run.sh check` refuses to grade any row whose own
repetitions disagree by more than 15% — the gate threshold would be inside that spread, and a gate
that fires on the driver's mood is a coin flip.

### Storm conditions, not only quiet ones

A radiation-belt library exists for the disturbed magnetosphere, so a benchmark suite that only ever
measures `maginput = 0` is measuring the one condition it is not for. Until `ext_t89.hpp` landed
there was nothing here parameterized by activity to sweep; there is now, and the T89 rows are swept
across the whole published Kp envelope — Kp × 10 from 0 to 90, IRBEM's `maginput` slot-1 scaling
(`docs/source/api/general_information.rst`) — on both sides, at a real ~20° dipole tilt rather than
the ψ = 0 case where Tsyganenko (1989)'s `sin ψ` terms drop out and the model evaluates a strictly
cheaper expression than it ever does in flight.

The quiet and extreme-storm ends are in the table above; the eight bins between them are in **Also
measured**. The sweep is flat in time on both sides, which is the expected answer and is now a
measured one: T89 selects one of seven coefficient sets by Kp and then evaluates the same closed
form, so activity changes the numbers and not the work. A model whose cost *does* move with activity
— a T96 or a T01 root-find — would show it in exactly this sweep, and that is what these rows are
for once one lands. **Dst, Pdyn and southward Bz are still not swept anywhere**, because no model in
this module reads them yet; that sweep is owed with the first model that does.

## The crossover — where the device starts winning

This is the whole argument for the module's shape. A single field-line trace is a serial RK4 chain,
so no hardware makes it faster; the parallelism is entirely **across** lines, and it only becomes
available if the caller hands over the whole batch at once. A loop calling `trace_invariant` per
point cannot be accelerated no matter what hardware is present, which is why `trace_invariant_batch`
exists and why it is the entry point to call.

Below the crossover the synchronous dispatch floor dominates and **the host wins** — which is why
`prefer_gpu` consults the per-kernel crossover instead of assuming the device is always right. A
blanket "always use the GPU" would make small batches slower and the dipole kernel slower at every
size.

Every point of this curve traces the same *distribution* of field lines. The input set is indexed by
the base-2 radical inverse rather than swept linearly, so the first 64 lines already cover L = 2…8
evenly and so do the first 65 536; the only thing that changes along the curve is the batch size.

Measured by [`bench/irbem_bench.cpp`](../../bench/irbem_bench.cpp); reproduce with `scripts/bench_run.sh publish irbem-crossover`.

<!-- BENCH:irbem-crossover begin -->
<!-- cheatah-bench-stamp v1
     suite:        irbem-crossover
     generated:    2026-08-28
     commit:       950e2e1
     host:         pop-os, 20 CPUs @ 400 MHz
     cpu-scaling:  powersave
     build:        g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0, Google Benchmark @ 192ef10025eb2c4cdd392bc502f0c852196baa48
     device:       NVIDIA
     competitors:  IRBEM /tmp/irbem-builds/libirbem-O2.so (the -O2 REBUILD, never the shipped no-`-O` library)
     harness:      reps=5, min_time=0.2s, random-interleaving=on
     statistic:    median over repetitions; every rate is REAL time, never CPU time — a device
                   lane is blocked, not running, for nearly all of its wall clock
     publishable:  true
     layout:       crossover
     watch:        space/irbem/, bench/irbem_bench.cpp

     PRODUCED BY:
       SPACE_IRBEM_ORACLE='/tmp/irbem-builds/libirbem-O2.so' CHEATAH_GPU_LINALG_VK_DEVICE='NVIDIA' \
           scripts/bench_run.sh publish irbem-crossover
       which runs: build/bench/irbem_bench --benchmark_filter='trace_batch' \
           --benchmark_repetitions=5 --benchmark_min_time=0.2s \
           --benchmark_enable_random_interleaving=true \
           --benchmark_out_format=json --benchmark_out=docs/bench/irbem-crossover.json
-->

| field lines | host per line | device per line | speedup | verdict |
|--:|--:|--:|--:|---|
| 64 | 58.58 µs | 219.13 µs ±16% | 0.27× | **host wins** |
| 256 | 58.26 µs | 66.05 µs ±8% | 0.88× | **host wins** |
| 1024 | 57.91 µs | 25.68 µs ±20% | 2.26× | device wins |
| 4096 | 57.87 µs | 4.93 µs ±25% | 11.74× | device wins |
| 16384 | 57.71 µs | 1.52 µs ±10% | 37.87× | device wins |
| 65536 | 57.73 µs | 1.31 µs ±5% | 44.17× | device wins |
<!-- BENCH:irbem-crossover end -->

## Arithmetic intensity — the number that predicts the verdict

Same seam. Same transfers. Same device. The three kernels differ in one thing — how much arithmetic
they do per byte moved — and that one thing decides the verdict, including the sign of it.

The flop counts are a **term model**, stated rather than hidden: one spherical-harmonic coefficient
slot is charged five flops (`G = g·cos mφ + h·sin mφ`, the multiply by the Legendre value, and the
accumulation), which is the conventional way a harmonic model's cost is quoted. A full instruction
count of the Slang inner loop — which also carries the Legendre recursion and the two horizontal
components — comes out roughly six times higher. The argument these rows make is about the *ratio*
between them, and a ratio is unchanged by a constant factor applied to every row.

The trace row's flop count is the only one that is measured rather than counted, because it depends
on how many RK4 steps a line actually takes. It is also a **lower bound**: `FieldLine::steps`
reports the two mirror-half integrations, and the initial walk to the magnetic equator is on top of
that.

That is why the trace row here reads lower than the `~9 400 flops/byte` quoted in
[`gpu/dispatch.hpp`](gpu/dispatch.hpp)'s registry comment. Both are the same term model; they differ
in the step count they assume, and this one's is measured over *this* input set — field lines
started on the dipole equator of L = 2…8 at 30…80° pitch, where the walk to the true minimum-B point
is short. A set of lines started off-equator would trace further and this number would rise. The
figure is a property of the workload as much as of the kernel, which is exactly why it is measured
here rather than written down once.

Byte counts are exact, and are the compulsory traffic per point. The IGRF coefficient and Legendre
tables are ~1.7 KB for the *whole batch* — 0.03 bytes/point at 65 536 — so they do not appear.

All three rows are measured at 65 536 points (65 536 field lines for the tracer), which is past
every one of their crossovers, so each row is the verdict at a batch size where the device has
whatever advantage it is going to get. The two small kernels' device lanes are the noisy ones: they
are over in microseconds and contend with the driver, so read their spreads and not only their
medians. The IGRF row in particular has published a speedup between 6× and 43× across runs while
never once losing; the row's claim is the *sign* and the order of magnitude, not the digits.

Measured by [`bench/irbem_bench.cpp`](../../bench/irbem_bench.cpp); reproduce with `scripts/bench_run.sh publish irbem-intensity`.

<!-- BENCH:irbem-intensity begin -->
<!-- cheatah-bench-stamp v1
     suite:        irbem-intensity
     generated:    2026-08-28
     commit:       950e2e1
     host:         pop-os, 20 CPUs @ 836 MHz
     cpu-scaling:  powersave
     build:        g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0, Google Benchmark @ 192ef10025eb2c4cdd392bc502f0c852196baa48
     device:       NVIDIA
     competitors:  IRBEM /tmp/irbem-builds/libirbem-O2.so (the -O2 REBUILD, never the shipped no-`-O` library)
     harness:      reps=5, min_time=0.2s, random-interleaving=on
     statistic:    median over repetitions; every rate is REAL time, never CPU time — a device
                   lane is blocked, not running, for nearly all of its wall clock
     publishable:  true
     layout:       intensity
     watch:        space/irbem/, bench/irbem_bench.cpp

     PRODUCED BY:
       SPACE_IRBEM_ORACLE='/tmp/irbem-builds/libirbem-O2.so' CHEATAH_GPU_LINALG_VK_DEVICE='NVIDIA' \
           scripts/bench_run.sh publish irbem-intensity
       which runs: build/bench/irbem_bench --benchmark_filter='(dipole_field_(host_)?batch|igrf_batch|trace_batch/65536)' \
           --benchmark_repetitions=5 --benchmark_min_time=0.2s \
           --benchmark_enable_random_interleaving=true \
           --benchmark_out_format=json --benchmark_out=docs/bench/irbem-intensity.json
-->

| kernel | computes | bytes/point | flops/point | flops/byte | host per point | device per point | speedup | verdict |
|---|---|--:|--:|--:|--:|--:|--:|---|
| `irbem_dipole_f32` | centred dipole | 24 | 12 | **0.5** | 2.07 ns | 4.18 ns ±40% | 0.50× | GPU **loses** |
| `irbem_igrf_f32` | IGRF-14, degree 13 | 24 | 525 | **21.9** | 301.87 ns | 12.06 ns ±43% | 25.02× | GPU **wins** |
| `irbem_trace_i_f32` | trace + second invariant | 36 | 132530 | **3681.4** | 58.76 µs | 1.27 µs | 46.42× | GPU **wins** |

> mean RK4 steps per traced line, measured over 256 lines of the input set: 63.1
<!-- BENCH:irbem-intensity end -->

## Reproducing, and the regression gate

```sh
scripts/bench_run.sh build                    # Google Benchmark at the pin, then the suite
scripts/bench_run.sh verify                   # every device lane against its own host lane
scripts/bench_run.sh publish all              # run everything and regenerate all three tables
scripts/bench_run.sh check                    # verify, then gate against bench/baseline.csv
```

`bench/baseline.csv` is a committed row per measured lane, in nanoseconds per point.
`bench_run.sh check` re-measures and fails on a regression that clears **both** a 1.40× ratio and a
2 ns absolute floor — a routine measured in single nanoseconds moves 30% on machine load alone, and
gating on ratio alone would make the gate a coin flip. A row that trips the gate is re-measured once
before it is called a regression. The gate also fails if `-ffast-math` ever appears in the benchmark
build, if `bench/CMakeLists.txt` and `scripts/bench_run.sh` stop agreeing about the flags, or if any
device lane leaves its error budget — that check runs first, because a throughput gate cannot tell a
fast kernel from a fast wrong one.

All four failure modes have been exercised deliberately and then reverted: a doctored `baseline.csv`
row (fails, `21.01×`), a deleted `-march=native` (fails, "the two build paths have drifted"), an
added `-ffast-math` (fails, "it is banned"), and a device IGRF launch truncated from degree 13 to
12 (fails, `max_rel_dev 1.02` against a `1e-06` budget). A gate nobody has watched fail is a gate
nobody knows works.

Every table is stamped with the CPU governor it was measured under. These numbers were taken with
scaling enabled and the governor in `powersave`, which is the honest default state of the machine
rather than a tuned one; a pinned-frequency run moves the host columns and leaves the ratios
roughly where they are, but the stamp is there so a future run that disagrees can be told apart
from a regression.

The device is selected with `CHEATAH_GPU_LINALG_VK_DEVICE` (`NVIDIA`, `Intel`, `llvmpipe`). On a
machine with no device, or a build with no `cheatah-gpu-linalg` on the include path, the GPU rows
report themselves skipped and the table shows dashes: a truthful table, not a broken build. The
same is true of the IRBEM columns when `SPACE_IRBEM_ORACLE` points at nothing.
