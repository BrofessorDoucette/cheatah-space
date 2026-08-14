# Security

## Threat model

Single-trust, the cheatah household standard: the process that compiles against this
extension is trusted; its *inputs* may still be hostile or wrong. cheatah-space's attack
surface is deliberately tiny — `space.time` is pure, header-only C++20 arithmetic with no
I/O, no platform headers, no global state, and no allocations on the scalar path — so the
security posture is mostly about keeping it that way as `space.cdf` (file parsing!) and
`space.irbem` land.

## Standing security review

| Area | Finding | Status |
|------|---------|--------|
| input types | a wrong argument type reaching a template body | **By design** — the explicit `TimeInput` concepts reject non-time inputs at compile time |
| module integrity | a tampered module header on the extension path | **By design** — every header ships with a `.sha512` sidecar; purrc/the runtime verify it on resolve/load (scripts/sign-modules.sh, gated) |
| memory safety | leaks / UB in the conversions or the vectorized ndarray path | **Checked every push** — the QA gate runs the unit suite under ASan + UBSan and Valgrind memcheck |
| numerical range | a Julian Date double keeps ~50 µs resolution near the present epoch | **Documented** — round-trips are specified and tested with an absolute tolerance, not sold as exact |
| timing | conversions are not constant-time | **By design** — numerics, no secrecy claims |
| future file I/O | `space.cdf` will parse untrusted files | **Roadmap rule** — the reader lands with fuzz-style hostile-input tests before it ships; nothing is parsed today |

## Reporting

Open a GitHub issue for anything without exploitation impact. For something sensitive,
email the maintainer (see NOTICE) before filing publicly.
