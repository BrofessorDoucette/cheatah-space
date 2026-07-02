# cheatah-space

> ⚠️ **Work in progress.** Early scaffold. `space.time` works and is tested; `space.cdf` and
> `space.irbem` are roadmaps only. End-to-end `biome add cheatah-space` also depends on a
> cheatah toolchain increment that is still WIP (wiring an extension's headers into purrc's
> link line) — today the modules build and the tests pass via `purrc --import-root`. APIs,
> layout, and namespaces may change.

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

## Modules

| module | what it is | status |
|---|---|---|
| **`space.time`** | Julian Date, Modified Julian Date, J2000, and the NASA **CDF_EPOCH** bridge. Templated with concepts (scalar `Value` \| numeric `ndarray` \| future datetime-struct `ndarray`). The pun (spacetime) is intended. | **working** → [space/time](space/time/) |
| **`space.cdf`** | NASA **Common Data Format** I/O, written from scratch in C++. | roadmap → [space/cdf](space/cdf/) |
| **`space.irbem`** | Radiation-belt & magnetic-field models — a clean-room reimplementation of **[PRBEM/IRBEM](https://github.com/PRBEM/IRBEM)**. | roadmap → [space/irbem](space/irbem/) |

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
  allocation-free for scalars; the toolchain signs and verifies module headers.
- **Honest docs.** Every function carries `@complexity`, `@alloc`, `@test`/`@systest`, and
  `@perf` where benchmarked. Brief everywhere else.

## Layout

```
cheatah-space/
├── space/                  # the `space` package (import root)
│   ├── space.hpp           # package header — includes the submodules
│   ├── time/time.hpp       # space.time  (C++20, concepts, vectorized)   [working]
│   ├── cdf/                # space.cdf    (roadmap)
│   └── irbem/              # space.irbem  (roadmap)
├── tests/                  # the test suite, written in cheatah (.purr), importing space.*
│   ├── check.purr          # tiny assertion helper
│   └── test_*.purr
├── scripts/qa_gate.sh      # compiles + runs every test through purrc/cheatah
├── .githooks/pre-push      # blocks pushes unless the qa_gate passes
├── cmake/CPM.cmake         # CPM bootstrap (only for fetching the toolchain when standalone)
├── CMakeLists.txt          # header-only INTERFACE target + standalone test/hook wiring
└── .vscode/                # IntelliSense → Ctrl-click into ../cheatah; .purr highlighting
```

## Developing

The extension is header-only, so there is nothing to compile to *use* it. To run the tests you
need the cheatah toolchain; while it is mid-release we build against the **local sibling
checkout** at `../cheatah`.

```sh
./scripts/setup-hooks.sh    # once: point git at .githooks (also done by cmake configure)
./scripts/qa_gate.sh        # compile every tests/*.purr against space.* and run it
```

The qa_gate locates a built `purrc`/`cheatah` under `../cheatah/build/*` (override with
`CHEATAH_DIR`), builds the toolchain if needed, then compiles and runs each cheatah test —
which also gates that the C++ templates instantiate cleanly on both the scalar and ndarray
paths. The **pre-push hook** runs this gate and blocks the push on failure (`--no-verify`
overrides, in a real emergency).

The full cheatah toolchain is **not** vendored here; CPM (or the local checkout) provides it.

## License

MIT — © 2026 BigBrain LLC (Joshua Doucette, on its behalf). See [LICENSE](LICENSE) and [NOTICE](NOTICE).
