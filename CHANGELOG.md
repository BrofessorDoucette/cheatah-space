# Changelog

All notable changes to cheatah-space. This project is **alpha** — expect breaking changes
between releases. It is a cheatah standard-library extension and joins the Biome Standard
alongside the other cheatah extensions.

## Unreleased

### space.cdf — reads real NASA CDF files into ndarrays

`space.cdf` opens a CDF 3.x file and hands its variables back as cheatah `ndarray`s, GZIP
included. That is the whole of this tranche's goal: the flow from a file on disk to a plot.

```purr
let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
let b = cdf.values(f, "F")                  # ndarray[float64], IMF magnitude in nT
let fig = figure.line(figure.new_figure(), hours, b)
plot.save(fig, "out/omni_imf.png")
```

- **The reader** — `open`, `var_names`, `record_count`, `data_type_name`, `shape`, `values`,
  `values_i64`. A `File` is immutable after open, so it copies cheaply and any number of threads
  may read one. Records are decoded in a single pass straight into the ndarray's buffer.
- **Nine headers**: the format vocabulary, one bounds-checked byte reader, mmap behind a testable
  syscall seam, the internal records, the VXR index walk, the decode kernels, DEFLATE, and the
  leap-second table. All at 100% lines and functions.
- **DEFLATE from scratch** (RFC 1951), no zlib. GZIP is the only compression that occurs in
  NASA's public archive, so this is what turns "opens OMNI" into "opens RBSP, MMS and THEMIS".
- **`examples/purr_space/`** — three runnable programs ending in a rendered PNG, plus a runner
  that skips cleanly without the toolchain, cheatah-plot, or the corpus.

**Verified against bytes NASA produced, not against our own opinion of them.** The expected
values were decoded by hand from the OMNI byte stream before the reader existed: `F` opens
6.92 / 5.84 / 5.71 nT and `Epoch` at 63587289600000.0 ms. And `test_alltypes.cdf` carries
`Longitude` (GZIP level 9) beside `longitude_copy` (uncompressed) holding the same values — they
now agree value for value, which says our inflate agrees with the encoder NASA's writer used.

Four things learned the hard way, recorded so they are not rediscovered:

- **The VXR index describes allocation, not truth.** OMNI's `Epoch` index ends at record 45055
  while `maxRec` is 44639, because VVRs are allocated at blocking-factor granularity. A reader
  that trusts the index hands back 416 records of uninitialised slack.
- **The index is a multi-level tree with sibling chains** in the very first real file opened —
  a naive VXR→VVR walk fails immediately. The walk is iterative and capped on depth and node
  count, because one byte makes `VXRnext` point at itself.
- **CDF's fields are unaligned by construction**, so every multi-byte load is `memcpy` + `bswap`
  rather than a cast to a wider pointer, which is undefined behaviour with UBSan in the gate.
- **A guard before a switch can make the arms it guards unreachable** — and an arm no
  instantiation can enter is a line no test can cover. The lossless-conversion refusals live
  inside the decode arms for exactly that reason.

Deferred, and refused by name rather than guessed at: CDF 2.x, whole-file compression,
multi-file CDFs, VAX floats, sparse records with gaps, RLE/Huffman/adaptive Huffman, column-major
N-D variables, the attributes API, the writer, the checksum, signing, and the byte-for-byte
differential harness. The three legacy codecs are last on purpose: the specification names them
but never documents their bitstreams, and no file in the archive uses any of them.

### space.irbem — four external field models, and a vectorised CPU lane

`kext=1` Mead & Fairfield (1975), `kext=5` Olson-Pfitzer quiet (1977), `kext=6` Olson-Pfitzer
dynamic (1988) and `kext=8` Ostapenko & Maltsev (1997) join T89, each written to its published
paper and verified against the vendored Fortran oracle across the four corpus regimes and the four
real storm events. Each ships the same shape: a host evaluator templated on the float type, per
model validity through `status.hpp`, a **shared** Slang function outside every kernel guard so a
tracer calls the same physics as a direct evaluation, a guarded kernel, and a registry row that a
completeness test checks against the Slang entry points and the CMake shader table.

Mead reaches **oracle parity at 2.1e-9 relative RMS** (870 grid points x 4 Kp bins x 3 tilts, worst
single point 1.2e-7). Establishing that took a black-box provenance experiment rather than a
transcription: a free 20-term quadratic fit to the isolated oracle field converges only in SM
coordinates at a 4-degree-aberrated position, and refitting on the paper's own y-symmetric basis
recovers all 68 published coefficients as 3-4 significant-figure decimals to 1e-9 relative.

**Olson-Pfitzer dynamic carries a documented gap, not parity.** Its 1988 citation is a conference
abstract with no equations and the OP77 parent is a McDonnell Douglas report: neither the
functional form nor the coefficients were ever published. What ships is the documented structure --
`B = s^3 B_q(s r) + dDst* R(r)` with `s = (P/P0)^(1/6)` -- with every constant from a citable
source (CODATA, Chapman-Ferraro/Mead 1964, O'Brien & McPherron 2000, Dessler-Parker-Sckopke). It
sits **~50% RMS-relative** from IRBEM's `kext=6` in the belts, with a measured structure floor of
67.8%. The oracle's tail does follow the published pressure law (fitted `s` = 1.110/1.260/1.420
/1.475 against a predicted 1.122/1.260/1.414/1.468), but inside 6.5 Re it is an amplitude scaling
that no published relation predicts. Parity is unreachable without reading the LGPL source, which
this clean room does not do. `div B = 0` holds regardless: the stencil residual falls as h^2.

Every model is checked for divergence-free-ness by a second-order stencil whose residual must fall
as h^2 -- the one correctness check that needs no oracle and cannot be satisfied by agreeing with a
wrong reference.

### space.irbem — the CPU batch lane vectorises across points

The Legendre recursion cannot vectorise *within* one field evaluation (loop-carried, ~7 wide), so
the batch lane was 84% scalar. `batch_soa.hpp` makes the **point index** the SIMD lane and keeps
the n/m recursion loops outer: **362 -> 100 ns/point, 3.61x**, bit-identical to the scalar lane by
`memcmp` across truncations, epochs, policies and tail lengths. An objdump test counts packed
double ops in the strip symbol and fails if the lane ever decays to scalar (a scalar-row control
takes it from 390 packed to 0). Bit identity is conditional on `-ffp-contract=off`, which the repo
sets globally and the header now documents for consumers who build outside it.

Three variants measured flat or worse and were reverted rather than shipped: autovectorised plain
arrays (2.0x), 512-bit vector rows (0.65x -- double-pumped on AVX2), and 16-point strips (3.3x).

### space.cdf — the format layer

The three headers everything above the bytes stands on. No records are parsed yet; this is what
makes parsing them safe and testable.

- **`types.hpp`** — the format's closed vocabulary: 17 data types, 21 encodings, 14 record types,
  the CDR flag bits and the error set. The 21 encodings collapse to **four** decode classes, which
  is the simplification the decoder rests on; NETWORK is plain big-endian IEEE, not XDR despite the
  name. `CDF_REAL4`/`CDF_FLOAT` and `CDF_REAL8`/`CDF_DOUBLE` are distinct format codes with
  identical layouts and both occur in real files, so both are carried rather than normalized away.
  Every `switch` omits `default:` so an added enumerator is a compiler warning, not a silent gap.
- **`bytes.hpp`** — the one bounds-checked reader. Every multi-byte load is `memcpy` + `bswap`,
  never a cast to a wider pointer: CDF's fields are **unaligned by construction** (`GDRoffset` is
  an 8-byte read at CDR offset 12; a VXR's `Last[]` starts at `28 + 4*Nentries`, odd-aligned for
  odd counts), and a wide-pointer read of an unaligned address is undefined behaviour with UBSan
  in the gate. All control integers are big-endian *regardless of the file's data encoding* — the
  encoding governs values only.
- **`require()`** — the chokepoint that makes 100% line coverage reachable on a parser. Each
  validation is one line executed on every successful parse, with the `throw` in a single
  function. Scattered `if (bad) throw` would make every throw its own uncovered line needing a
  bespoke malformed file; that is why parsers famously stall short of full coverage.
- **`mapping.hpp`** — the file is mapped once rather than read through a buffer pool (NASA's
  library seeks and copies per record, which is what `CDFsetCacheSize` tunes). The syscalls sit
  behind a `SysOps` policy because `mmap` failing is unreachable from any craftable input; in
  production the policy is a stateless empty struct that inlines away.

**Verified against a real NASA file, not against our own parser's opinion of one.** The expected
values were decoded by hand from the OMNI `hro2_1min` byte stream before this code existed: CDR
size 312, GDR at 320 = 8 + 312, encoding NETWORK, flags row-major + single-file, `zVDRhead` 21601,
the `Epoch` zVDR at `CDF_EPOCH`/maxRec 44639, `NumElems` at field offset 64 — and `eof` equal to
the file size exactly, the invariant that holds only when no checksum is present.

Two traps found and recorded rather than patched over. `mapping.hpp` first measured 38% because
the real-file tests were **skipping**: ctest runs the suite from the repo root and
`scripts/coverage.sh` from `build/cov`, so a bare relative path resolved in one and not the other.
A skipped test reads green, so the coverage run had quietly stopped exercising the real syscalls
at all — the fake was covered and the production path was not. The corpus path is now searched
upward from the working directory. And reaching the real `mmap`-failure branch needed something
that opens and stats plausibly but cannot be mapped: a directory does exactly that.

### space.irbem — the TOTAL field runs on the GPU, and storms are verified against the reference

The stated scope limit — "a TotalField batch runs on the CPU today" — is fixed, and the loop the
storm corpus was built for is closed against the Fortran oracle.

- **`irbem_trace_total_f32`**: the field-line trace through IGRF ⊕ T89, one thread per field line,
  both models evaluated per RK4 stage on the device. T89's physics was factored into a shared
  `t89_eval` (as `igrf_eval` already was), so each piece of physics exists on the device exactly
  once. Eight bindings — the shared descriptor layout's cap; `Leases::capacity` moves 7 → 8 in the
  same edit, by the rule the capacity-said-4 segfault taught. Verified device-run in all four
  activity regimes, ~2 µs/line steady-state, fp32-vs-fp64 agreement 1.3–2.7e-05 (budget 1e-4).
- **`make_lstar` is generic over `GeoFieldModel`.** The drift-shell machinery traces through the
  combined kernel; the flux quadrature and footpoint walks evaluate the TOTAL field on the host
  (IRBEM's own cap integral does the same — evaluating only the internal part there would bias Φ
  by the external field's ~0.1–0.3% at r = 1, comparable to the whole L\* budget in a storm). The
  internal-only device fast-paths are guarded by `if constexpr (is_igrf_v<M>)`: routing a total
  field through a kernel staged with half its physics is a compile-time impossibility, not a
  runtime surprise.
- **The storm differential vs the oracle, `kext=4`, all seven Kp bins, matched `options=9`**:
  worst |dL\*| 0.028–0.088 (bins 1–6), 0.123 (bin 7) — against 0.0066 for the internal field,
  where both sides evaluate the same model. The order-of-magnitude contrast is the finding: the
  numerics agree, the MODEL FAMILIES differ. IRBEM's `kext=4` is T89c, the 1992 revision that
  added two tilt-modulation tail terms and refit with ISEE-1/2 data (Oulu reference
  documentation), distributed only as GPL code and never published as equations — and the
  free-parameter experiment already in `ext_t89.hpp` proves no refit of the published form can
  close the gap. A clean room implements what was published; the harness reports the difference
  instead of hiding it.
- **`Status::NotConverged` is now covered by real physics**: a midnight L=9 shell under Kp ≥ 6
  runs into T89's stretched tail and cannot close. The branch was unreachable with a pure
  internal field (measured: every shell closes out to L=40) and untestable by mock (the machinery
  was Igrf-typed); the first caller able to produce it was the storm itself.
- The real-syscall lane of `space/cdf/mapping.hpp` got its first tests — the SysOps seam's fake
  was covered while the production `PosixSysOps` sat at 0%, a seam inverted. Six float-induction
  loops in the T89 divergence test became integer induction: accumulated `+=` rounding decided
  whether endpoint samples ran, which in a divergence VERIFIER means silently weakening the check.

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
