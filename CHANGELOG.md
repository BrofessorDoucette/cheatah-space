# Changelog

All notable changes to cheatah-space. This project is **alpha** — expect breaking changes
between releases. It is a cheatah standard-library extension and joins the Biome Standard
alongside the other cheatah extensions.

## v0.1.0-alpha (unreleased) — the time module, house-gated

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
