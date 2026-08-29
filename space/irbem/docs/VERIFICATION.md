# What is implemented, what is verified, and how the domain was sampled

@brief The honest status of every `space.irbem` model — implemented vs verified vs neither, the
exact domain sampled against the Fortran reference, and what each number does and does not prove.

This file exists because "implemented" and "verified" are different claims, and a reimplementation
of a physics library is worth nothing if a reader cannot tell which one applies to the routine they
are about to trust. It is the first thing to read before extending this module, and it is kept
current with the code rather than written once.

Three things are recorded for every model: whether the code exists, what was measured against the
vendored Fortran oracle and **over which part of the domain**, and whether an independent
adversarial pass reproduced those numbers. A model can be green on all its own tests and still be
wrong; only the third column speaks to that.

## The verification ladder

Not all agreement is equal, so the words below are used precisely and always mean the same thing.

| term | what it means |
|---|---|
| **oracle parity** | Matches IRBEM to the internal-field budget (1e-6 relative on B or tighter). The model's published form *and* its fitted coefficients both exist, so any disagreement is our bug. |
| **published-form, measured gap** | The paper's form is implemented faithfully, but IRBEM evaluates something else — an unpublished revision, or coefficients that were never printed. The disagreement is **measured, quoted, and attributed**; it is never averaged away or hidden behind a loose tolerance. |
| **div B = 0** | A second-order stencil whose residual falls as h². This needs no oracle at all and is the strongest correctness check here: it cannot be satisfied by agreeing with a wrong reference. |
| **adversarially verified** | A second agent, told to *refute*, reproduced every number, perturbed the code three ways (a coefficient, a sign, a frame), and confirmed a test failed each time. Without this, a green suite may only prove the tests agree with the implementation. |

## Status, model by model

| kext | model | implemented | agreement with IRBEM | div B | adversarially verified |
|---|---|---|---|---|---|
| — | IGRF-14 (internal) | **yes** | **parity**, 8.8e-07 device-vs-host | n/a (potential field by construction) | yes |
| 1 | Mead & Fairfield 1975 | **yes** | **parity**, 2.1e-9 RMS rel | exact (degree-2 polynomial; residual h-independent) | **no — owed** |
| 4 | Tsyganenko 1989 | **yes** | **published-form, measured gap** | h² confirmed | yes |
| 5 | Olson–Pfitzer quiet 1977 | **yes** | see the header's own measurement | h² confirmed | **no — owed** |
| 6 | Olson–Pfitzer dynamic 1988 | **yes** | **published-form, measured gap: ~50 % RMS rel** | h² confirmed | **no — owed** |
| 8 | Ostapenko & Maltsev 1997 | **yes** | see the header's own measurement | h² confirmed | **no — owed** |
| 7 | T96 | **parked on `wave4-wip`** | not established | not established | no |
| 9 | T01 | **parked on `wave4-wip`** | not established | not established | no |
| 10 | T01-storm | **parked on `wave4-wip`** | not established | not established | no |
| 2/3 | T87 short/long | **not started** | — | — | — |
| 11–14 | T04/TS05, TS07D | **not started** | — | — | — |

Everything on `wave4-wip` compiles but is **not in the build and not on `main`**: T96's differential
caps were never measured, and T01/T01-storm have no test suite at all. That branch's commit message
lists exactly what each file is missing, including three deliberate perturbations that were never
reverted. Treat none of it as working code.

### The two gaps, stated rather than softened

**T89 (`kext=4`) is not the published T89.** IRBEM evaluates the 1992 *T89c* revision — two extra
tilt-modulation tail terms and an ISEE-1/2 refit — which was distributed only as source and whose
equations were never published. A free refit floating all 28 parameters of the published form
cannot close it (0.44 nT floor), which is evidence of a structural difference rather than a fitting
failure. Consequence on L\*: |dL\*| 0.028–0.088 across Kp bins 1–6 and 0.123 at bin 7, against
0.0066 for the internal field alone.

**OP-dynamic (`kext=6`) was never published at all.** Its 1988 citation is a conference abstract
with no equations; the OP77 parent is a McDonnell Douglas/AFOSR report. Neither the functional form
nor the coefficients are in the literature. What ships is the documented *structure* —
`B = s³·B_q(s·r) + ΔDst*·R(r)` with `s = (P/P₀)^(1/6)` — with every constant taken from a citable
source (CODATA; Chapman–Ferraro/Mead 1964; O'Brien & McPherron 2000; Dessler–Parker–Sckopke). It
sits **~50 % RMS-relative** from IRBEM's `kext=6` in the belts, with a measured structure floor of
67.8 %. The oracle's *tail* does follow the published pressure law (fitted `s` = 1.110/1.260/1.420/
1.475 against a predicted 1.122/1.260/1.414/1.468), but inside 6.5 Rₑ it behaves as an amplitude
scaling ≈ P^0.65 that no published relation predicts. **This model is not IRBEM-compatible and must
never be described as such.** Closing it requires reading LGPL source, which this clean room does
not do.

## How the domain was sampled

A tolerance is meaningless without the domain it was measured over: agreement at L = 4 on a quiet
day proves almost nothing about a storm-time tail. Every differential in this module is therefore
run over the corpus in [`tests/irbem_domain_corpus.hpp`](../../../tests/irbem_domain_corpus.hpp),
which samples two ways at once, because they fail differently.

**Four synthetic regimes**, chosen to sit in each regime's *interior* so a failure is not a boundary
artifact (boundaries are exercised separately by the validity tests):

| regime | Kp | Dst | what it is for |
|---|---|---|---|
| Quiet | < 2 | \|Dst\| < 20 nT | the easy case, and the one that proves least |
| Moderate | 3–4 | −30…−50 nT | typical disturbed conditions |
| Storm | 5–7 | −100…−200 nT | main phase, Bz strongly south |
| Extreme | 8–9 | −250…−400 nT | at or past most models' fitted envelopes |

**Four real storm events**, because synthetic sweeps find corners the solar wind rarely visits while
real events cluster where it actually goes: Halloween 2003, March 1989 (Hydro-Québec, the largest of
the space age, to Dst −589), St Patrick's Day 2015, and February 2022 (the storm that de-orbited 38
Starlink satellites). Their drivers are representative published *peak* values, not reconstructed
time series — enough to place a model in the right regime, not enough to reproduce a specific
minute. That limitation is stated here rather than implied by silence.

**Southward Bz is sampled densely on purpose.** It is what loads the tail and drives the storm-time
reconfiguration these models were fitted to reproduce; a corpus that samples Bz symmetrically spends
half its points where the physics is least interesting.

Sampled along with the drivers, in every differential:

- **Position** — all eight local times (noon and midnight matter most: they are the least
  dipole-like), on and off the equator, and radially out to each model's published validity limit.
- **Epoch** — multiple IGRF epochs including interpolated years, so the secular variation and the
  epoch-dependent rotations are exercised rather than assumed.
- **Dipole tilt** — at minimum three values spanning both signs, since every external model is
  tilt-modulated and a sign error there is invisible at zero tilt.
- **Validity boundaries from both sides** — just inside returns `Ok` with a value; just outside
  returns `OutOfValidityRange` **with the value still returned**, never a silent clamp.
- **Driver continuity** — for the continuously-driven models the drivers are swept continuously and
  *smoothness is asserted*: nearby drivers must give nearby fields. For the Kp-binned models the
  opposite is asserted — bin identity — because for them a smoothness test would be a bug.

Concretely, Mead's parity number (2.1e-9 RMS relative) is 870 grid points × 4 Kp bins × 3 tilts
in-suite, plus 1800 scattered points per bin at 1.2–16.7 Rₑ over six epochs in
[`tools/oracle/mead_diff.cpp`](../../../tools/oracle/mead_diff.cpp); worst single point 1.2e-7.

## What this does *not* establish

- **No adversarial pass exists for Mead, OP-quiet, OP-dynamic or Ostapenko.** Their verifiers were
  killed by a session limit. Their own suites are green and their numbers are self-reported. This is
  the single largest verification debt in the module.
- **`@test` tags are not yet enforced.** [`scripts/check_doc_tags.sh`](../../../scripts/check_doc_tags.sh)
  is written but wired into nothing, and running it reports 102 public functions under
  `space/irbem/` carrying no `@test` tag and 23 tags naming tests that do not exist. Until it is a
  gate stage, a `@test` tag is a claim, not a guarantee.
- **Fewer than 1 % of public functions carry a compiling `@code{.purr}` example**, so the
  purr-facing surface is far less exercised than the C++ one.
- **The oracle is a reference, not truth.** IRBEM's own discretization error at its recommended
  settings is 1.2e-3–4.0e-3 relative (0.010–0.017 absolute in L\* at L ≈ 6) — larger than the 0.01
  L\* target this module aims at — so a differential is only meaningful at matched `options`, and
  those are recorded per run. See [ERROR_BUDGET.md](ERROR_BUDGET.md).
- Nothing here is a statement about **physical** accuracy. These are empirical models fitted to
  spacecraft data and carry ~10 % errors of their own; matching IRBEM means reproducing the model,
  not the magnetosphere.
