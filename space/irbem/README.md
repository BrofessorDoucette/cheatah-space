# space.irbem (roadmap)

Radiation-belt & magnetic-field models for cheatah — a clean-room reimplementation of
**[PRBEM/IRBEM](https://github.com/PRBEM/IRBEM)** (the canonical library, ~93% Fortran).
Written **from scratch in C++ with zero dependencies**: no Fortran runtime, no IRBEM link,
nothing required to build, install, or use `space.irbem`. Correctness parity first, then beat
it on speed.

> IRBEM is **LGPL-3.0**. We do **not** copy or link its code — this is a clean-room reimpl from
> the published model definitions (IGRF, the Tsyganenko T89/T96/T01/T04 external fields, etc.).
> IRBEM is used **only as an optional, dev-only oracle**: fetched on demand into
> `space/irbem/vendor/` (git-ignored), built separately, and used to cross-check outputs and
> benchmark. A normal `biome add cheatah-space` never touches it, and no IRBEM-licensed code
> enters this repo.

## What it computes

Magnetic coordinates (L-shell, L*, B), drift shells, coordinate transforms (GEO/GSM/SM/GEI…),
field-line tracing, and the radiation-belt model evaluations — over scalars and, vectorized,
over `ndarray` time/position batches (sharing the `space.time` epoch scales).

## Plan

1. **Coordinate transforms + time** — the GEO/GSM/SM/GEI/MAG rotations, driven by `space.time`
   (sidereal time, J2000). Pure, vectorized, exhaustively tested first.
2. **Internal field: IGRF** — spherical-harmonic geomagnetic field; the foundation everything
   else stands on.
3. **External fields: Tsyganenko** T89 → T96 → T01 → T04.
4. **Field-line tracing + L/L\*/drift shells** — the headline `make_lstar`/`get_field`-class
   routines, as concept-templated cheatah-callable functions.

## Verification & benchmark (optional)

A dev harness runs identical inputs through IRBEM (the oracle) and `space.irbem`, asserts
agreement within model tolerance over a fixed grid of epochs/positions, and records throughput.
The measured speedup vs IRBEM becomes each function's `@perf`; we ship only once we match its
numbers and beat its speed. Skipped automatically when IRBEM isn't present.

## House rules

Every function documents `@complexity`, `@alloc`, `@test`/`@systest`, and `@perf` (speedup vs
IRBEM). Models cite their source paper/coefficients in a brief comment.
