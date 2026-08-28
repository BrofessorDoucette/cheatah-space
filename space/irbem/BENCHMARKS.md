# space.irbem benchmarks

The generated benchmark tables for the [`space.irbem`](README.md) module: every routine that
exists, on every lane it has, against the vendored IRBEM library run as a black box. The module
page describes what the routines compute; this page holds the numbers.

Three things are being measured, and they answer three different questions:

1. **[Per routine](#every-routine-on-every-lane-it-has)** — how fast is each routine on the host,
   on the device, and in IRBEM? The last column of that table is the one to read first: it says
   whether a routine can reach the GPU **at all**, and it is derived from the kernel registry
   rather than asserted, so the porting gap is visible instead of implied.
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

**Never `-ffast-math`.** The suite builds `-O3 -march=native -ffp-contract=off`, which is the
arithmetic the correctness tests pinned. `-ffast-math` would let the compiler reassociate the
spherical-harmonic sums and the RK4 accumulation, and the binary would then be measuring different
mathematics from the one that was verified. `bench_run.sh check` fails the build if the flag ever
appears.

## Every routine, on every lane it has

Read the **on GPU?** column, not just the timings. Three of its four answers are failure states,
and the interesting one is `registered, NO LAUNCHER`: a kernel can be written, compiled, correct
and even measurable in isolation while remaining unreachable through the seam every caller goes
through. That is a porting gap, not a success, and this column is generated from
`gpu/dispatch.hpp`'s `registered_kernels` so it says so out loud.

Routines with an epoch-dependent geometry appear twice. `(hot)` reuses one `Rotations` across the
ephemeris — the shape this module's API is built for — while `(cold)` rebuilds it per point, which
is the like-for-like against IRBEM, whose entry points take a date with every call and cannot
amortize. Reporting only the hot number would flatter this module; reporting only the cold one
would hide its actual design advantage.

Measured by [`bench/irbem_bench.cpp`](../../bench/irbem_bench.cpp); reproduce with `scripts/bench_run.sh publish irbem-per-routine`.

<!-- BENCH:irbem-per-routine begin -->
<!-- BENCH:irbem-per-routine end -->

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
     commit:       5692c3b
     host:         pop-os, 20 CPUs @ 1982 MHz
     cpu-scaling:  powersave
     build:        g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0, Google Benchmark @ 192ef10025eb2c4cdd392bc502f0c852196baa48
     device:       NVIDIA
     competitors:  IRBEM /tmp/irbem-builds/libirbem-O2.so (the -O2 REBUILD, never the shipped no-`-O` library)
     harness:      reps=3, min_time=0.1s, random-interleaving=on
     statistic:    median over repetitions; every rate is REAL time, never CPU time — a device
                   lane is blocked, not running, for nearly all of its wall clock
     publishable:  true
     layout:       crossover
     watch:        space/irbem/, bench/irbem_bench.cpp

     PRODUCED BY:
       SPACE_IRBEM_ORACLE='/tmp/irbem-builds/libirbem-O2.so' CHEATAH_GPU_LINALG_VK_DEVICE='NVIDIA' \
           scripts/bench_run.sh publish irbem-crossover
       which runs: build/bench/irbem_bench --benchmark_filter='trace_batch' \
           --benchmark_repetitions=3 --benchmark_min_time=0.1s \
           --benchmark_enable_random_interleaving=true \
           --benchmark_out_format=json --benchmark_out=docs/bench/irbem-crossover.json
-->

| field lines | host per line | device per line | speedup | verdict |
|--:|--:|--:|--:|---|
| 64 | 60.21 µs | 216.68 µs | 0.28× | **host wins** |
| 256 | 59.19 µs | 66.10 µs | 0.90× | **host wins** |
| 1024 | 59.23 µs | 17.18 µs | 3.45× | device wins |
| 4096 | 59.51 µs | 4.26 µs | 13.97× | device wins |
| 16384 | 59.35 µs | 1.53 µs | 38.87× | device wins |
| 65536 | 59.39 µs | 1.21 µs | 48.95× | device wins |
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

Byte counts are exact, and are the compulsory traffic per point. The IGRF coefficient and Legendre
tables are ~1.7 KB for the *whole batch* — 0.03 bytes/point at 65 536 — so they do not appear.

Measured by [`bench/irbem_bench.cpp`](../../bench/irbem_bench.cpp); reproduce with `scripts/bench_run.sh publish irbem-intensity`.

<!-- BENCH:irbem-intensity begin -->
<!-- BENCH:irbem-intensity end -->

## Reproducing, and the regression gate

```sh
scripts/bench_run.sh build                    # Google Benchmark at the pin, then the suite
scripts/bench_run.sh publish all              # run everything and regenerate all three tables
scripts/bench_run.sh check                    # gate against bench/baseline.csv
```

`bench/baseline.csv` is a committed row per measured lane, in nanoseconds per point.
`bench_run.sh check` re-measures and fails on a regression that clears **both** a 1.40× ratio and a
2 ns absolute floor — a routine measured in single nanoseconds moves 30% on machine load alone, and
gating on ratio alone would make the gate a coin flip. A row that trips the gate is re-measured once
before it is called a regression. The gate also fails if `-ffast-math` ever appears in the benchmark
build, or if `bench/CMakeLists.txt` and `scripts/bench_run.sh` stop agreeing about the flags.

The device is selected with `CHEATAH_GPU_LINALG_VK_DEVICE` (`NVIDIA`, `Intel`, `llvmpipe`). On a
machine with no device, or a build with no `cheatah-gpu-linalg` on the include path, the GPU rows
report themselves skipped and the table shows dashes: a truthful table, not a broken build. The
same is true of the IRBEM columns when `SPACE_IRBEM_ORACLE` points at nothing.
