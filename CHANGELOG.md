# Changelog

All notable changes to cheatah-space. This project is **alpha** — expect breaking changes
between releases. It is a cheatah standard-library extension and joins the Biome Standard
alongside the other cheatah extensions.

## Unreleased

### space.irbem — the frame layer

Work has begun on `space.irbem`, a from-scratch reimplementation of
[PRBEM/IRBEM](https://github.com/PRBEM/IRBEM) written to the published papers. This first piece is
the vocabulary the rest of the module is built from; no field model ships yet.

- `space/irbem/frames.hpp`: the twelve reference frames, `FrameKind` (a frame's three components
  are Cartesian, spherical or geodetic — that distinction is load-bearing and now typed), and the
  frame-tagged `Position<F>` / `FieldVector<F>`.
- **The frame lives in the type.** `Position<Frame::GEO>` and `Position<Frame::GSM>` are distinct
  types, so the classic defect of this domain — passing one frame's numbers to a routine expecting
  another, where every component is plausible and only the answer is wrong — is a compile error.
  Both wrappers are the size of a `vec3d` and trivially copyable: the tag costs nothing at runtime.
- IRBEM's runtime `sysaxes` integer enters the typed world through exactly one boundary conversion
  (`frame_from_sysaxes`), which reports an out-of-range code as absent rather than defaulting to a
  frame. The heliospheric frames correctly report that `sysaxes` cannot name them.
- The geometry is `cheatah::fixarray` (`vec3d`, `norm`, column-major `mat3d`) — the module adds no
  container types of its own.

### The IRBEM oracle, and a measured error budget

- `tools/oracle/` — a dev-only harness that builds the vendored IRBEM checkout and runs it as a
  **black box** through its documented C entry points (`dlopen`, never linked into anything we
  ship). It is outside every QA-gate scope and never runs in the gate.
- `space/irbem/docs/ERROR_BUDGET.md` now carries **measured** numbers rather than estimates. The
  two that matter:
  - **Discretization at IRBEM's own recommended resolution is 1.2e-3 to 4.0e-3 relative** in L\*
    (0.010–0.017 absolute at L≈6), across L≈2–8 and two external fields. That is *at or above* a
    0.01 absolute target — so 0.01 is only meaningful at matched `options`, and the differential
    suite must record which resolution it ran at.
  - Convergence is **not monotonic**: `options=4` is repeatedly worse than `options=2`. The
    first-order rectangle quadrature plus a fixed-step theta march is the largest accuracy lever
    in the original algorithm.
- The oracle is built **twice**, as-shipped and `-O2`, because **2.7x of IRBEM's cost is its
  missing `-O` flag alone** (42.4 ms vs 15.5 ms per L\* evaluation). Benchmarks quote the `-O2`
  build; anything else inflates our numbers before we write a line of physics.

### Repo

- `-ffp-contract=off` is now set **globally**. FMA contraction silently changes the last bit of
  `a*b+c`, and `space.irbem`'s self-goldens are byte comparisons; a flag that must be remembered at
  each target is one that will be forgotten at a target.
- `scripts/cppcheck.sh` no longer analyses `space/*/vendor/`, the git-ignored dev-only reference
  implementations. They are third-party, never built and never shipped, so their findings are not
  ours to fix — the Doxyfile already excluded them for the same reason.
- Coverage now reports **branches** (100%), which `space.time` alone could not exercise: it is
  straight-line arithmetic with no branches at all.

## v0.1.0-alpha (2026-08-14) — the time module, house-gated

The first release-shaped state of the repo: one working module, `space.time` (Julian Date,
Modified Julian Date, the J2000 epoch offsets, and the NASA CDF_EPOCH bridge), hand-authored
as header-only, concept-templated C++20 — scalar and ndarray-vectorized through the same
functions — plus the full extension-grade QA harness around it. `space.cdf` and `space.irbem`
remain roadmaps (design notes in their subdirectory READMEs); nothing of them ships here.

### The module

- `space.time`: 3 reference-epoch functions (`jd_unix_epoch`, `jd_j2000`,
  `cdf_epoch_unix_offset_ms`) + 10 conversions (`unix_to_jd`/`jd_to_unix`,
  `jd_to_mjd`/`mjd_to_jd`, `unix_to_mjd`/`mjd_to_unix`, `jd_to_j2000_seconds`/
  `jd_to_j2000_centuries`, `unix_to_cdf_epoch`/`cdf_epoch_to_unix`), each accepting a
  numeric scalar or a numeric `ndarray` through one `TimeInput`-constrained template.
- Every public function's docs carry `@complexity`, `@alloc`, truthful `@systest` tags naming
  the systests that cover it, and an `@par Example` with a compiling `@code{.purr}` block
  (the doc-examples convention; the gate compiles all 13 blocks).

### Biome packaging

- `cheatah.toml` manifest (pure cheatah: no extension or system dependencies; targets the
  1.10.0-alpha toolchain).
- `scripts/sign-modules.sh` writes the `.sha512` sidecars that make `import space` /
  `import space.time` resolve as VERIFIED modules on the extension path.
- `scripts/test-biome-install.sh` sandboxes the real consumer flow: the `space/` package
  copied to a throwaway dir, a fresh project compiled with cheatah env vars cleared and only
  `CHEATAH_MODULE_PATH` wired — a Julian-date round trip must print `RESULT: PASS`.

### QA gate grown to the extension bar

`scripts/qa_gate.sh` (pre-push enforced) now runs: 100% unit coverage (clang source-based,
lines + functions over the space headers, README table drift-checked) → 100% Javadoc
(`scripts/doc_coverage.sh`) + compiling doc examples → module-sidecar verification →
debug build → the `.purr` systests (moved from `tests/` to `systests/`) → the biome-install
sandbox → GoogleTest unit suite (`tests/space_time_test.cpp`, new) → ASan + UBSan →
Valgrind memcheck (`security/run-valgrind.sh`) → cppcheck → the private-reference scan
(also guarding commit messages via the commit-msg and pre-push hooks).

### Removed

- `cmake/CPM.cmake` — dead: nothing referenced it (the toolchain is the sibling checkout /
  `CHEATAH_DIR`, and biome consumers get headers only).

### Roadmap (unchanged)

- `space.cdf` — native NASA Common Data Format reader/writer, from scratch, zero
  dependencies; brings CDF_TT2000 + EPOCH16 with its leap-second table.
- `space.irbem` — clean-room reimplementation of the radiation-belt / magnetic-field model
  library (coordinate transforms → IGRF → Tsyganenko fields → L*/drift shells).
