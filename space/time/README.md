# space.time

Julian dates & epoch conversions — an extension *for* the cheatah standard time library.
The pun (spacetime) is intended. Hand-authored, header-only C++20, **zero dependencies**.

```purr
import space.time as st
import ndarray

let jd  = st.unix_to_jd(st.jd_unix_epoch())      # scalar
let jds = st.unix_to_jd(ndarray.array([0.0, 86400.0, 946728000.0]))   # whole array at once
let ms  = st.unix_to_cdf_epoch(0.0)              # 62167219200000.0 — feeds space.cdf
```

## Templated with explicit concepts

Every conversion is one concept-constrained template (`time/time.hpp`) that accepts:

- a **sensible scalar `Value`** (`Numeric` — a number of seconds/days/ms), or
- a **numeric `ndarray`** (`NumericArray` — vectorized over the ndarray SIMD path), or
- an **ndarray of datetime structs** (`DatetimeArray`) — lit up via the `is_datetime_v`
  trait hook once `cheatah::datetime` defines its struct.

The union is the `TimeInput` concept, so passing the wrong thing yields a crisp compiler
error instead of a template spew. ndarray support is opt-in by include path (`__has_include`),
so a scalar-only program never needs `ndarray`.

## Functions

Reference epochs (functions, so callers never hard-code magic numbers): `jd_unix_epoch`,
`jd_j2000`, `cdf_epoch_unix_offset_ms`.

Conversions (each O(n), allocation-free for scalars):

- `unix_to_jd` / `jd_to_unix` — Unix seconds ⇄ Julian Date.
- `jd_to_mjd` / `mjd_to_jd`, `unix_to_mjd` / `mjd_to_unix` — Modified Julian Date.
- `jd_to_j2000_seconds` / `jd_to_j2000_centuries` — time since J2000.0.
- `unix_to_cdf_epoch` / `cdf_epoch_to_unix` — the NASA **CDF_EPOCH** (ms since year 0) bridge.

## Numerical note

A Julian Date stored as a `double` carries a ~2.44e6-day offset, leaving ~50 µs of resolution
near the present epoch; round-tripping a small value through a JD therefore lands within ~µs,
not bit-exactly (the tests assert this with an absolute tolerance). Known instants that land on
exact half-days (Unix epoch, J2000) convert exactly.

## Roadmap

- **CDF_TT2000** (ns since J2000 with leap seconds) — needs a leap-second table; lands with
  `space.cdf`. **CDF_EPOCH16** (picosecond, two doubles).
- GMST / sidereal time once a consumer needs it.

## Tests

Exercised from cheatah in [../../tests/](../../tests/) (`import space.time`), run by the
[qa_gate](../../scripts/qa_gate.sh). Every public function carries `@complexity`, `@alloc`,
and `@test`; `@perf` is reserved for the benchmarked modules (cdf/irbem).
