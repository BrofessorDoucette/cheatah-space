# cheatah-space

> ⚠️ **Work in progress.** Early alpha. All three modules are in and QA-gated at 100% line and
> function coverage: `space.time`, `space.cdf` (readers over real NASA files), and `space.irbem`
> (IGRF-14, five external field models, L\*, drift shells, and a GPU lane for the integrals).
> APIs, layout, and namespaces may change until the first tagged release.

The first cheatah standard-library **extension**: the `space` package — astronomy and
space-physics. Hand-authored, header-only C++20 (templated with concepts), surfaced in the
`cheatah::space::*` namespaces and consumed from cheatah as `import space.<module>`.

**Zero dependencies**, like cheatah itself. A researcher with a working biome installs it with
one line — nothing else to fetch, build, or link:

```sh
biome add cheatah-space
```

```purr
import space.time as st
import ndarray

let jd  = st.unix_to_jd(1700000000.0)                   # one timestamp -> Julian Date
let jds = st.unix_to_jd(ndarray.array([0.0, 86400.0]))  # a whole array, vectorized
```

Every module header ships with a `.sha512` sidecar (scripts/sign-modules.sh), so purrc
resolves `import space.time` as a VERIFIED module on the extension path
(`CHEATAH_MODULE_PATH`, which `biome add` sets) and the runtime checks the header on load —
no git checkout, no `--import-root`. The QA gate sandboxes that exact consumer flow on every
push (scripts/test-biome-install.sh).

## Modules

| module | what it is | status |
|---|---|---|
| **`space.time`** | Julian Date, Modified Julian Date, J2000, and the NASA **CDF_EPOCH** bridge. Templated with concepts (scalar `Value` \| numeric `ndarray` \| future datetime-struct `ndarray`). The pun (spacetime) is intended. | **working** → [space/time](space/time/) |
| **`space.cdf`** | NASA **Common Data Format** I/O, written from scratch in C++ — the record layer, the variable index tree, and a reader verified against real mission files. | **working** → [space/cdf](space/cdf/) |
| **`space.irbem`** | Radiation-belt & magnetic-field models — a from-scratch reimplementation of **[PRBEM/IRBEM](https://github.com/PRBEM/IRBEM)**, written to the published papers, vectorized on the CPU and with the field-line and drift-shell integrals evaluated in parallel on the GPU. | **working** → [space/irbem](space/irbem/) |

The through-line: each module has a canonical reference (NASA CDF C lib, IRBEM Fortran) we
verify against and then outperform. Those references are **optional, dev-only** oracles — never
required and never linked. Performance wins are documented per function with `@perf`.

## Design principles

- **Zero dependencies.** Build, install, and use with only the cheatah toolchain. No NASA CDF,
  no IRBEM, no zlib — nothing external is ever required.
- **Cross-platform.** Portable C++20: no platform headers, no global state, no I/O in the math.
- **Performance.** Concept-templated and vectorized over `ndarray` (SIMD); benchmarked against
  the references, shipped only when we win.
- **Security.** Explicit concepts reject bad inputs at compile time; the modules are pure and
  allocation-free for scalars; every module header is signed and verified on load. Threat
  model + standing review: [SECURITY.md](SECURITY.md).
- **Honest docs.** Every function carries `@complexity`, `@alloc`, truthful `@systest` tags,
  and an `@par Example` whose `@code{.purr}` block the gate COMPILES; `@perf` where
  benchmarked. The doc gate holds the public API at 100% Javadoc.

## Layout

```
cheatah-space/
├── space/                  # the `space` package (import root; headers + .sha512 sidecars)
│   ├── space.hpp           # package umbrella — includes the submodules
│   ├── time/time.hpp       # space.time  (C++20, concepts, vectorized)   [working]
│   ├── cdf/                # space.cdf    (records, index tree, reader)      [working]
│   └── irbem/              # space.irbem  (IGRF, external models, L*, GPU)   [working]
├── systests/               # cheatah (.purr) system tests, importing space.* end to end
├── tests/                  # C++ unit tests (GoogleTest) — the coverage + memcheck harness
├── scripts/                # qa_gate.sh + the individual checks it runs (see Developing)
├── security/               # run-valgrind.sh + suppressions
├── .githooks/              # commit-msg (private-refs scan) + pre-push (the qa_gate)
├── cheatah.toml            # biome manifest (pure cheatah — no dependencies)
├── CMakeLists.txt          # header-only INTERFACE target + standalone test/hook wiring
├── CMakePresets.json       # debug / release / asan presets
└── Doxyfile                # the doc parser config (XML out; doc_coverage.sh strict mode)
```

## Developing

The extension is header-only, so there is nothing to compile to *use* it. To develop it you
need the cheatah toolchain: the **local sibling checkout** at `../cheatah` (override with
`CHEATAH_DIR`) provides `purrc`, the `cheatah` runtime, and the stdlib sources the C++ unit
tests compile against. The toolchain is **not** vendored here.

```sh
./scripts/setup-hooks.sh    # once: point git at .githooks (also done by cmake configure)
./scripts/qa_gate.sh        # the full gate the pre-push hook enforces
```

The QA gate runs, in order — every stage a hard gate:

1. **Unit coverage** — clang source-based, **100% lines + functions** over the `space/`
   headers via `cheatah_space_tests` (scripts/coverage.sh; the table below is drift-checked).
2. **Doc coverage** — 100% Javadoc (scripts/doc_coverage.sh) and every `@par Example`
   block compiles (scripts/check_doc_examples.sh).
3. **Module sidecars** — every header verifies against its `.sha512` (sign-modules.sh).
4. **Build** (debug preset).
5. **System tests** — each `systests/test_*.purr` compiled by purrc against `space/` and run
   under the cheatah runtime (`RESULT: PASS`), which also gates that the templates
   instantiate cleanly on the scalar AND ndarray paths. Then the **biome-install sandbox**.
6. **Unit tests** (ctest), then **ASan + UBSan**, then **Valgrind memcheck**, then
   **cppcheck**, then the **private-reference scan**.

### Unit-test coverage

Measured by scripts/coverage.sh (clang source-based) over the space package headers;
committed by the gate:

<!-- coverage:start -->
| Metric | space package |
|--------|---------------|
| **Lines** | 100.00% (3447/3447) |
| **Functions** | 100.00% (365/365) |
| Regions | 98.74% |
| Branches | 93.73% |
<!-- coverage:end -->

(`space.time` is straight-line arithmetic and has no branches of its own; the Regions and Branches
metrics come from `space.cdf` and `space.irbem`. Lines and functions are the gated pair — the gate
fails below 100% on either.)

## License

MIT — © 2026 BigBrain LLC (Joshua Doucette, on its behalf). See [LICENSE](LICENSE) and [NOTICE](NOTICE).
