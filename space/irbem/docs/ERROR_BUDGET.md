# space.irbem — the error budget

For each output quantity: what accuracy the physics needs, what the discretization already costs,
and therefore what arithmetic precision is justified. The precision policy, the GPU decomposition
and every test tolerance are *derived* from this document rather than chosen independently.

**Every number below with a measurement note is observed, not asserted.** The measurements come
from [`tools/oracle/convergence.cpp`](../../../tools/oracle/convergence.cpp), which runs the IRBEM
oracle as a black box across its two resolution knobs. Rows still marked *pending* are estimates
and are labelled as such.

---

## 1. Why there is headroom: the L\* identity

L\* is defined through the third adiabatic invariant, the magnetic flux Φ enclosed by the drift
shell:

```
L* = 2πk₀ / (Φ R_E)        so        δL*/L* = δΦ/Φ        exactly.
```

An absolute budget on L\* is therefore a *relative* budget on Φ. A 0.01 absolute target at L\*=6
is a relative target of 1.7 × 10⁻³; at L\*=2 it is 5 × 10⁻³.

## 2. Measured: what IRBEM's own discretization costs

IRBEM exposes two resolution knobs, and its documentation describes `options(3)=0` as the
recommended setting, "an error of ~2% at L=6".

- `options(3)` = θ resolution: `dθ = π / (720·(options(3)+1))` — the step of the drift-shell
  root-find *and* of the flux quadrature.
- `options(4)` = azimuth count: `Nder = 25·(options(4)+1)`.

Sweeping both together and measuring how far L\* moves from its most-converged value:

> **Setup** — IRBEM `e7cecb0`, gfortran 13.3.0, built `-O2 -ffp-contract=off -fno-fast-math`.
> 2015-180 12:00 UT, GEO Cartesian, IGRF internal (`options(5)=0`), `options(1)=1`.
> 7 points from L≈2 to L≈8 including geosynchronous and one off-equator; two external fields
> (`kext=0` internal-only, `kext=5` Olson–Pfitzer quiet). Reference is `options(3)=options(4)=9`.

| `options(3,4)` | L\* move vs (9,9), relative | at L≈6, absolute | mean wall time / point |
|---|---|---|---|
| **0, 0** (recommended) | **1.2 × 10⁻³ … 4.0 × 10⁻³** | **0.010 … 0.017** | **15.5 ms** |
| 1, 1 | 1.6 × 10⁻⁵ … 3.4 × 10⁻³ | | 28 ms |
| 2, 2 | 1.1 × 10⁻⁴ … 1.1 × 10⁻³ | | 44 ms |
| 4, 4 | 2.0 × 10⁻⁵ … 9.6 × 10⁻⁴ | | 79 ms |
| 9, 9 | (reference) | | 183 ms |

### Two findings that change the design

**(a) The 0.01 L\* target is at or below IRBEM's own default-resolution error.** At L≈6 the
default setting sits 0.010–0.017 from its own converged value. So a 0.01 agreement target is only
meaningful **at matched settings**, where both implementations are deterministic and we are
comparing algorithms rather than resolutions. Claiming 0.01 *accuracy* — closeness to the
converged answer — requires `options(3,4) ≥ 2`. The differential suite must therefore record the
options it ran at, and never compare across settings.

**(b) Convergence is not monotonic.** `options=4` is repeatedly *worse* than `options=2`
(kext=0 L≈3: 1.9 × 10⁻⁴ → 5.6 × 10⁻⁴; L≈8: 3.4 × 10⁻⁴ → 9.6 × 10⁻⁴; kext=5 L≈6.6: 1.8 × 10⁻³ at
res 0 → 3.4 × 10⁻³ at res 1). That is the signature of a first-order rectangle quadrature combined
with a fixed-step θ march whose root lands on a grid point that jitters as `dθ` changes. Two
consequences: the (9,9) reference is not "truth" — its own residual is plausibly 10⁻⁴–10⁻³ — and
the quadrature is confirmed as the single largest accuracy lever in the original algorithm, which
is what the `Compat::Improved` policy exists to address.

The documentation's "~2% at L=6" is roughly 7× larger than the ~0.3% self-convergence measured
here; the doc figure most likely describes accuracy against truth rather than against its own
converged value. We quote our measurement, and say which one it is.

## 3. The precision decision

| error source | magnitude | basis |
|---|---|---|
| **Discretization at IRBEM's default** | **1.2 × 10⁻³ … 4.0 × 10⁻³** | measured, §2 |
| fp32 roundoff over a 10³-step trace, worst case | ~1.2 × 10⁻⁴ | `N·ε`, ε = 6 × 10⁻⁸ — *pending* |
| fp32 roundoff over a 10³-step trace, random walk | ~4 × 10⁻⁶ | `√N·ε` — *pending* |
| naive fp32 sum of 10³ terms | ~1.2 × 10⁻⁴ | **avoided: accumulate in fp64** — *pending* |
| fp32 roundoff, single operation | 6 × 10⁻⁸ | `ε(binary32)/2` |

**Discretization exceeds fp32 trace roundoff by roughly 24× (against the pessimistic worst-case
roundoff bound) to 700× (against the random-walk estimate).** That is comfortable margin, and it
justifies an fp32 integrand — but it is *not* the "two to three orders of magnitude" the planning
estimate assumed from the documentation's 2% figure. The margin is real; it is one to three orders,
not two to three. The roundoff rows are analytic bounds and remain **pending** direct measurement
against the fp64 reference lane once that lane exists.

The one place fp32 would breach the budget is the **accumulation**: 10³ terms summed naively in
fp32 costs ~1.2 × 10⁻⁴, comparable to the discretization floor at the tighter settings. So every
policy keeps `accum = double`, and reductions run in a fixed order. That is precisely the
"integrals on the GPU, summed on the CPU" split.

## 4. Per-quantity tolerances

Provisional except where noted; each is replaced by a measured distribution as its wave lands. The
differential suite reports **max / p99 / RMS** per quantity, not just pass/fail — a regression that
moves p99 without breaching the cap is still a signal.

| quantity | budget | basis |
|---|---|---|
| coordinate transforms | 1e-10 rel | pure geometry, cheap, stays fp64 — *pending* |
| MLT | 1e-8 h | pure geometry — *pending* |
| `Bgeo`, `Blocal` | 1e-6 rel | fp32-limited; the models are empirical to ~10% — *pending* |
| `Bmin`, `Bmirr` | 1e-5 rel | a root-find on top of B — *pending* |
| `XJ` (I integral) | 1e-4 rel | quadrature-limited — *pending* |
| `Lm` | 1e-3 abs | Hilton's approximation carries ~1e-4 of its own — *pending* |
| **`L*`** | **0.01 abs, at matched `options`** | measured, §2(a) |
| footpoints | 1e-3 deg (~100 m) | trace-termination limited — *pending* |
| flux / dose models | 1e-4 rel | table interpolation — *pending* |

## 5. Measured: the performance baseline

Same setup, mean over 14 (point, kext) pairs at default resolution, single L\* evaluation:

| oracle build | mean per L\* point |
|---|---|
| **as IRBEM ships it** (`FFLAGS` with no `-O` at all) | **42.4 ms** |
| rebuilt `-O2 -ffp-contract=off -fno-fast-math` | **15.5 ms** |

**2.7× of IRBEM's cost is the missing optimization flag alone.** Benchmarks therefore quote the
`-O2` build: a speedup measured against the as-shipped binary would be inflated 2.7× before we
wrote a single line of physics, and that is not a number worth reporting.

For scale: 15.5 ms per point is ~26 minutes for a 100 000-point ephemeris, single-threaded, and
IRBEM has no threading. That is the opportunity this module exists to take.

## 6. What stays strict regardless

- Never `-ffast-math`. `-ffp-contract=off` globally, so the reference lane reproduces across
  optimization levels.
- **Reductions are ordered** — Φ and I accumulate in a fixed sequence in fp64, never a tree
  reduction whose order varies with workgroup size.
- The CPU lanes carry bit-exact self-goldens, so a refactor that perturbs the last bit is caught
  even though agreement with IRBEM is a tolerance.
- The GPU lane cannot carry committed bit-goldens (FMA contraction is at the driver's discretion),
  so it is held to this budget against the CPU reference plus run-to-run self-consistency.

---

*Regenerate §2 and §5 with:* `tools/oracle/run.sh`
