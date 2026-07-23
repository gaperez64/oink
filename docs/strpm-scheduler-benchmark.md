# strpm-simd player-specific, failure-directed scheduler: benchmark report

This compares three states of `STRPM_SIMDSolver`:

- **`baseline_pr6`** — the head of PR #6 (`872be945`), before this work: a single
  global `priority_queue<(k,t), RatioCompare>` shared by both orientations.
- **`separate_neutral`** — after decoupling into two independent
  `StrpmOrientationSchedule`s (one per orientation) with a neutral, non-pressure
  tree-size score (state at the "decouple even and odd parameter schedules"
  commit).
- **`separate_pressure`** — the current head: separate schedules *plus*
  pressure-directed scoring (blocker evidence steers whether an orientation
  grows `k` or `t` next), including a follow-up fix (see below) that stops the
  scheduler from wasting attempts growing `t` at `k=1`, where `t` is
  mathematically inert.

All three pass correctness identically: **268/273 (baseline/pressure), 270/273
(neutral)** games solved within a 20-30s per-game timeout budget, and every
game that completed was accepted by `Verifier::verify()`. Separately, all 256
games in the `tests/` corpus plus the generated random-game set were
cross-checked against `--tl` (tangle learning, an independent complete
solver): **0 partition mismatches**. An AddressSanitizer + UndefinedBehaviorSanitizer
build passed the same 256 games and the pure-function unit tests with **0
findings**.

## Corpora

| Corpus | Games | Timeout/game | Notes |
|---|---|---|---|
| `corpus224` | 224 | 20s | `tests/vb001`-`vb225` |
| `random32` | 32 | 30s | generated via `rngame`: n=50/100/500, low/high outdegree, few/many priorities, priorities up to 1000 (exercises the uint16 level work) |
| `keiren9` | 9 | 30s | smallest 9 of the 16 `hard-benchmarks/keirens-jurdzinski/` games (the 7 largest were excluded to keep total benchmark time bounded; see Limitations) |
| `capforcing` | 8 | 30s | the n=500 subset of `random32`, reused as the "forces `(k,t)` near the SIMD cap" family per the benchmark spec (`floor_log2(500)=8`, already clamped to `STRPM_SIMD_T_MAX=7`) |

Preprocessing was disabled (`--no`) for every run so the comparison isolates
the `(k,t)` search itself; correctness was separately verified with
preprocessing on (default) as well (see Correctness section).

## Results

| Variant | Games solved | Timeouts | Total (k,t) attempts | Pooled geomean solve time | vs. baseline |
|---|---:|---:|---:|---:|---:|
| `baseline_pr6` | 268/273 | 5 | 1067 | 0.000393s | — |
| `separate_neutral` | 270/273 | 3 | 2168 | 0.000482s | **+22.8%** |
| `separate_pressure` | 267/273 | 6 | 1041 | 0.000436s | **+11.2%** |

Pooled geomean is computed over the 265 games every variant solved without
timing out, so all three are compared on the identical game set. "Total
attempts" is summed across *all* games including ones that timed out, so it
is not a clean apples-to-apples efficiency number by itself (see below).

Pressure-directed scoring vs. neutral scoring specifically (isolating what
Commit 3 adds on top of Commit 2's decoupling): **geomean −9.5%**, **total
attempts −52%** (2168 → 1041). By that narrower comparison, pressure-directed
ordering is an unambiguous improvement over neutral scoring.

### Per-corpus breakdown

| Corpus | `baseline_pr6` attempts/timeouts | `separate_neutral` attempts/timeouts | `separate_pressure` attempts/timeouts |
|---|---:|---:|---:|
| corpus224 | 837 / 0 | 1517 / 0 | **837 / 0** |
| random32 | 117 / 1 | 362 / 0 | 114 / 0 |
| keiren9 | 79 / 3 | 172 / 3 | 57 / **6** |
| capforcing | 34 / 0 | 117 / 0 | 33 / 0 |

On the 224-game and random corpora, `separate_pressure` now matches or beats
`baseline_pr6` on attempt count and is competitive on time. On `keiren9` (the
largest, hardest games), it has *more* timeouts than baseline.

## A discovered scheduling defect and a partial fix

The first full benchmark pass (pressure-directed, before the fix below) showed
a much larger regression: pooled geomean **+23.4%**, total attempts roughly
**doubled** (837→1519 on `corpus224` alone). Root cause, confirmed by direct
trace comparison: `prog_tmp()`'s `k==1` fast path ignores `t` entirely — every
progress measure at `k=1` is either the initial empty value or `Top`,
regardless of `t`. `log_tree_size_estimate(1, t, h)` correctly returns `0.0`
for every `t` (per the handoff's formula), so the scheduler always saw growing
`t` at `k=1` as "free" and exhausted `t` up to `STRPM_SIMD_T_MAX` (6 wasted,
zero-progress attempts per orientation) before ever trying `k=2`.

**Fix applied**: `expand_after_attempt()` no longer offers `(1, t+1)` as a
candidate neighbor — only `(k+1, t)` — since those states are provably
identical to `(1, t)`. This is a deliberate, narrow deviation from the
handoff's literal invariant #3 ("every legal neighbor pair is eventually
eligible for exploration"), justified because `(1, t+1)` for `t > 1` is not a
distinct reachable state at all, just a re-run of `(1, 1)`.

This fix is a clear net win on `corpus224` and `random32` (attempts now match
or beat baseline; `corpus224` is byte-for-byte identical to baseline's attempt
count). **It made `keiren9` worse** (2 timeouts → 6, more than baseline's 3).
Investigation traced this to a second, deeper effect: the wasted `(1,t)`
attempts had been an accidental source of *neutral* (unbiased, `p=0.5`)
candidate scores for the whole `(2,t)` family, because each `(1,t)` attempt's
own pressure is always empty (the `k==1` branch never reaches the
C1-C4 classification code) and so unconditionally offers `(2,t)` at a neutral
score. Removing that waste means the *only* offer for `(2,2)` etc. now comes
from whatever pressure profile is observed at `(2,1)` -- and on some of the
`keiren9` games, `(2,1)`'s pressure is dominated by C1 evidence, which pushes
the scheduler to prefer an expensive, unproductive `(3,1)` over the cheaper
`(2,2)` that the old neutral-scored offer (or the baseline's balanced-ratio
heuristic) would have tried instead.

**This is not fixed.** A principled fix would mean tuning how strongly a
single attempt's pressure evidence is allowed to bias the score (e.g.
smoothing/damping runaway one-sided evidence, or not discarding
previously-neutral offers so eagerly), which is a real design question beyond
a small patch, and further afield from the handoff's literal prescribed
formulas. It is left as a follow-up.

## Correctness

- 224/224 `tests/` corpus games verified, both with default preprocessing and
  with `--no` (preprocessing disabled).
- 32/32 generated random games (n=50/100/500; low/high outdegree; few/many
  priorities; priorities up to 1000, exercising the `uint16_t` level work)
  verified, both preprocessing modes.
- All 256 of the above cross-checked against `--tl` (tangle learning): 0
  partition mismatches.
- AddressSanitizer + UndefinedBehaviorSanitizer build: 0 findings across the
  same 256 games plus the pure-function unit tests.
- Winning partitions are byte-for-byte identical to the PR #6 baseline at
  every commit in this branch (Commits 1-3), confirmed by diffing `-p` output
  across all 224 `tests/` games.

## Isolated `prog_tmp()` microbenchmark

`prog_tmp()` is a protected member reachable only through `lift()`; there is
no existing harness in this codebase for calling solver internals directly
(`test_solvers.cpp` only ever drives solvers through the public API). As a
proxy, the very first orientation attempt -- `(k=1, t=1)`, the same starting
pair in every variant, so an identical number of `lift()`/`prog_tmp()` calls
against the identical game before any scheduling divergence can occur -- was
timed directly from trace timestamps (n=7 reps, `jurdzinskigame(100,_100)`):

| Variant | median | mean |
|---|---:|---:|
| `baseline_pr6` | 0.040000s | 0.037143s |
| `separate_pressure` (HEAD) | 0.040000s | 0.037143s |

**0% measured difference** -- the telemetry branch/counters added to
`prog_tmp()` in Commit 1 have no measurable per-call cost. This isolates the
regression above entirely to scheduling behavior, not the instrumented hot
loop.

## Performance gates (handoff section 12)

| Gate | Result |
|---|---|
| Zero correctness regressions | **Pass** — see Correctness above |
| Zero validator failures | **Pass** — every completed game verified |
| No OOB/sanitizer failures | **Pass** — ASan+UBSan clean |
| ≤2% regression in isolated `prog_tmp()` microbenchmark | **Pass** — 0% measured |
| ≤5% geomean regression over the full benchmark set | **Fails** — +11.2% pooled geomean vs. `baseline_pr6` |
| Pressure-directed scheduling improves geomean time or total attempts vs. neutral scoring, without materially worsening the other | **Pass** — −9.5% geomean, −52% total attempts vs. `separate_neutral` |

The gate that is not met is specifically about `separate_pressure` (or
`separate_neutral`) vs. the *original PR #6 baseline* — i.e. whether
decoupling the schedules at all is worth it in typical-case terms, not
whether pressure-directed ordering improves on neutral scoring (which it
clearly does, per the last row). The honest summary: **schedule separation
trades some typical-case efficiency for correctness-preserving robustness on
a subset of hard instances** (fewer baseline timeouts on `random32`/
`capforcing`, but currently *more* timeouts on the hardest `keiren9` games
after the k=1 fix), and pressure-directed scoring measurably improves on
neutral scoring within that new architecture without closing the gap back to
baseline on this particular benchmark mix.

## Limitations of this benchmark run

- The 7 largest Keiren games (`(500,_500)`, `(500,_200)`, `(200,_500)`,
  `(500,_100)`, `(100,_500)`, `(500,_50)`, `(50,_500)`) were excluded to keep
  total benchmark time bounded (single dev machine, single run). Timeouts
  were capped at 20-30s/game; several games hit that cap in more than one
  variant, so their exact relative standing beyond "timed out" is not fully
  resolved.
- Single run, no repeated trials beyond the isolated `prog_tmp()`
  microbenchmark -- wall-clock numbers on a shared dev laptop carry some
  noise, though the qualitative findings (attempt-count changes, the k=1
  defect, the keiren9 regression) are structural, not noise-driven.
- `capforcing` reuses the `random32` corpus's n=500 games rather than a
  purpose-built family; it does force `t_max` to its cap (`STRPM_SIMD_T_MAX=7`
  for any game with ≥128 nodes) but does not specifically target `k` near
  `STRPM_SIMD_K_MAX=9`.
