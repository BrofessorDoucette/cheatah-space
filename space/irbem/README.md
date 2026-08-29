# space.irbem

Radiation-belt & magnetic-field models for cheatah — a **from-scratch reimplementation** of
**[PRBEM/IRBEM](https://github.com/PRBEM/IRBEM)** (the canonical library, ~96k lines of Fortran 77).
Written in C++20 with zero external dependencies: no Fortran runtime, no IRBEM link, nothing
required to build, install, or use `space.irbem`.

The reference library is single-threaded, unvectorized, ships compiled with no `-O` flag at all,
and carries ~20 mutable `COMMON` blocks that make it actively hostile to threading. One L\*
evaluation costs **~10⁵ magnetic-field model calls** — the [LANL\*
model](https://gmd.copernicus.org/articles/2/113/2009/) resorted to a neural-network surrogate
rather than pay it. That cost is the opportunity.

## Status

> **Implemented is not the same claim as verified.** Which models have been checked against the
> Fortran reference, to what tolerance, **over which part of the domain** (which activity regimes,
> which real storm events, which epochs, tilts, radii and local times), and which carry a measured
> and permanent gap rather than parity, is recorded in
> **[docs/VERIFICATION.md](docs/VERIFICATION.md)** — including what is *not* established. Read it
> before trusting or extending any routine here.

| piece | state |
|---|---|
| [`frames.hpp`](frames.hpp) — frames, `Position<F>` / `FieldVector<F>`, the `sysaxes` boundary | **in, gated** |
| coordinate transforms — geodetic (WGS84/Bowring), the Hapgood rotations, heliospheric | **in, gated** |
| [`igrf.hpp`](igrf.hpp) — IGRF-14 to degree 13, coefficients sourced from IAGA | **in, gated** |
| external field models — [T89](ext_t89.hpp) (4), [Mead](ext_mead.hpp) (1), [OP-quiet](ext_opq.hpp) (5), [OP-dynamic](ext_opd.hpp) (6), [Ostapenko](ext_ostapenko.hpp) (8) | **in, gated** |
| tracing, the bounce integral, L\*, drift shells, shell splitting | **in, gated** |
| the GPU lane — IGRF, T89, the total field, the line integral, the flux cap | **in, gated** |
| [`batch_soa.hpp`](batch_soa.hpp) — the CPU batch lane, point index as the SIMD lane (3.61×) | **in, gated** |
| T96 (7), T01 (9), T01-storm (10), T87 (2/3) | in flight |
| TS07D (13/14), the belt/atmosphere/effects models, TA15/TA16/GEO | not started |

Two models carry a **measured, documented gap** rather than parity, because the clean room cannot
reach what was never published: [T89](ext_t89.hpp) (IRBEM's `kext=4` is the unpublished 1992 T89c
revision) and [OP-dynamic](ext_opd.hpp) (neither its form nor its coefficients were ever published;
~50% RMS-relative, structure floor 67.8%). Each header states its own number. Everything else
targets oracle parity — Mead reaches 2.1e-9 relative.

## Licensing — why this is written to the papers

IRBEM is **LGPL-3.0**. Tsyganenko's own reference code, on
[his site](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/), is **GPL-3.0** — stricter
still. cheatah-space is MIT, so neither is a safe thing to derive from.

So the models are implemented from the **published definitions**: IGRF/IAGA; Tsyganenko
Planet. Space Sci. 1987 & 1989, JGR 1995/2002a,b/2005/2007 (the 2002 pair is titled *Mathematical
structure* and *Parameterization and fitting* — the formulation is published in full); Olson &
Pfitzer 1977/1988; Mead & Fairfield 1975; Roederer 1970 for L\*; Vallado for SGP4; NRL for MSIS;
NIST for SHIELDOSE-2.

IRBEM is used **only as an optional, dev-only oracle**: fetched on demand into `vendor/`
(git-ignored), built separately, and run as a black box to cross-check numbers and benchmark. A
normal `biome add cheatah-space` never touches it, and no copyleft-licensed code enters this repo.
Where a numeric result diverges and the papers cannot explain why, the specific fact that resolved
it — a constant, a threshold, a summation order — is cited in a code comment naming its source, so
the derivation trail is auditable.

Large fitted parameter files (TS07D's coefficients, the RBF models') are **never vendored and
never committed**; they are provisioned by a documented fetch step, and a model whose parameters
are absent reports that plainly rather than degrading silently.

## What it computes

Magnetic coordinates (L-shell, L\*, B), drift shells, coordinate transforms (GEO/GSM/SM/GEI…),
field-line tracing, and the radiation-belt model evaluations — over single points and, batched,
over whole ephemerides.

## Accuracy is a budget, not an accident

`docs/ERROR_BUDGET.md` is the load-bearing document: for each output quantity, what accuracy the
physics needs, what the discretization already costs, and therefore what arithmetic precision is
justified.

The headline is that there is far more headroom than intuition suggests. `L* = 2πk₀/(Φ R_E)`, so a
0.01 absolute budget on L\* is a *relative* budget of ~1.7 × 10⁻³ at L\*=6 — while IRBEM's own
recommended resolution carries **~2 × 10⁻² at L=6**, by its own documentation. Discretization
dominates roundoff by two to three orders of magnitude. Spending double precision on the integrand
buys nothing measurable; spending it on the *accumulation* buys a clean order.

Every number in that budget is measured by a convergence study, not asserted.

## Performance is measured, not asserted

[`BENCHMARKS.md`](BENCHMARKS.md) holds every performance number this module claims, generated by
[`scripts/bench_run.sh`](../../scripts/bench_run.sh) from [`bench/irbem_bench.cpp`](../../bench/irbem_bench.cpp)
— per routine against the IRBEM `-O2` rebuild, the trace-throughput crossover curve, and the
arithmetic-intensity study that explains why the device wins some kernels and loses others. Nothing
in those tables is hand-typed, and the "on GPU?" column is derived from the kernel registry, so a
routine that cannot reach the device says so rather than being quietly omitted.

## Lanes

| lane | precision | role |
|---|---|---|
| `reference` | fp64, one point at a time | the tightest check; what the others are measured against |
| `batch` | fp64, structure-of-arrays over N points | the CPU production lane |
| `gpu` | fp32 integrand, **fp64 ordered** reduction | the throughput lane |

The GPU lane is fp32 because it has to be, not merely because it is faster: SPIR-V's
`GLSL.std.450` transcendentals exist only for 16- and 32-bit floats, so an fp64 kernel calling
`exp()` — which every Tsyganenko model does — cannot be compiled at all. That constraint and the
error budget point the same way, and the result runs on every device including Apple and Intel,
which an fp64 kernel never could.

Reductions are **ordered**: Φ and I accumulate in a fixed sequence in fp64, never a tree reduction
whose order varies with workgroup size.

## Design notes

- **The geometry is [`cheatah::fixarray`](../../../cheatah/stdlib/fixarray/)** — `vec3d`/`mat3d`
  with the shape in the type, inline storage, nothing allocated, and column-major matrices so a
  frame rotation `m * v` is a sum of contiguous scaled columns. We add no container types of our
  own.
- **The frame lives in the type.** `Position<Frame::GEO>` and `Position<Frame::GSM>` are different
  types, so the classic defect of this domain — handing one frame's numbers to a routine expecting
  another, where every component is plausible and only the answer is wrong — is a compile error.
  IRBEM's runtime `sysaxes` integer enters through exactly one boundary conversion.
- **No silent failures.** A bad `sysaxes`, an unregistered GPU kernel, a missing coefficient file
  and an open field line all produce a named error. IRBEM's `baddata = -1e31` sentinel is kept at
  the API boundary for compatibility, but internally every failure carries a reason.
- **Nothing on the heap in a hot path**, asserted by a test that counts global `operator new`.
- **SIMD is verified, not assumed** — a codegen test disassembles the hot kernels and fails the
  gate if one came out scalar.

## House rules

Every function documents `@complexity` and `@alloc`, and the purr-facing surface adds `@systest`
and a compiling `@par Example`. Models cite their source paper and coefficients in a comment.
