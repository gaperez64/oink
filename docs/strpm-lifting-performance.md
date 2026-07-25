# Speeding up a single `strpm-simd` lifting attempt

Baseline: `872be94` (`strpm: small-buffer-optimize measure storage`), the head of the
open PR. All work here is orthogonal to the `(k,t)` scheduler; it targets the cost of
one lifting attempt rather than the choice of attempts.

Machine: Apple Silicon (NEON, `-march=armv8.7-a`), GCC 16 (`/opt/homebrew/bin/g++-16`),
`CMAKE_BUILD_TYPE=Release`.

Corpora:

- **VB** — all 224 `tests/vb*`
- **JURD** — all 16 `hard-benchmarks/keirens-jurdzinski/*.gm` (7.6k … 751k vertices)
- **RND** — 19 generated games (`rngame`), 1k/10k/100k vertices, priorities 100–500, degree 1–4

Timing is `solving took` summed per corpus, best of N repetitions, `--no` (preprocessing
off) so the solver is measured in isolation.

## Result

| corpus | baseline | after | delta |
|---|---|---|---|
| JURD | 64.36 s | 55.82 s | **−13.3%** |
| VB | 0.0064 s | 0.0061 s | **−4.7%** |
| RND | 0.1252 s | 0.1213 s | **−3.1%** |

Correctness: `oink -v` (`Verifier::verify`) passes on all 224 VB games with and without
preprocessing (448/448), and on all 16 JURD plus all 19 RND games.

## Where the time actually goes

Establishing this first changed the plan twice, so it is worth recording.

**`lift` and `compare` are everything; `prog_tmp` is nothing.** `sample` on
`jurdzinskigame(500,500)`, top-of-stack: `lift` 4413, `compare` 1977, inner `run` 244,
**`prog_tmp` 35**. The most intricate and most parameterised routine in the file is
~0.5% of runtime.

**The solver lives at `k <= 3`.** Final `(k,t)` across all 224 VB games: `k=2,t=1` in
175, `k=1,t=1` in 22, `k=2,t=2` in 15, `k=3,t=2` in 12 — nothing above `k=3`. All 19 RND
games finish at `k<=2`; JURD finishes at `k=2,t=3` / `k=3,t=3`. Since `nlanes == k-1`,
the 8-lane registers carry **1–2 useful lanes in 100% of observed workloads**.
`compare()` already had a hand-written *scalar* `k==2` path, which is the same
observation made in code.

**One attempt dominates.** On `jurdzinskigame(100,100)` the `k=2,t=1` attempt is
1,628,061 of 1,863,787 lifts (87%) and ~91% of wall time. Every attempt made partial
progress, so scheduling wastes almost nothing here — independent confirmation that the
scheduler is not the lever.

**The worklist dominates, and its calls mostly scan.** Per-attempt counters added for
this work (`-t 1`), on the dominant `k=2,t=1` attempt of `jurdzinskigame(100,100)`:

```
lift-attempts=8978894  top=517265  fast=1761673  scan=6718320
succ-scanned=19574539  sweep=54701  worklist=8924193  strategy=18364
```

99.4% of lift attempts come from the worklist phase, not the initial sweep. Of the calls
that do work, 79% take the scan path (min/max over the whole successor list, mean 2.91
successors) rather than the single-target fast path. VB and RND agree qualitatively
(scan is 61% and 45% of calls respectively).

## What worked: measures as a value type (−13.3% JURD)

`tmp_*`/`best_*` were solver members, so every `to_tmp()`/`compare()` pair round-tripped
through memory: the store side wrote through `this`, and out-of-line `compare()` reloaded
it. They are now a local `Measure` value type passed by reference.

**On its own this made things 9.4% *worse*.** The win only appears once `compare` is
actually inlined — as a plain out-of-line call taking two references, the compiler must
spill both measures to the stack around every call, which costs more than the member
version. With `__attribute__((always_inline))` on `compare`, the successor loop stays in
registers and the same change becomes −13.3%. Worth remembering: this refactor is not
independently beneficial, it is beneficial *conditional on inlining*.

## What did not work: compact `k`-specialised storage (reverted)

`NodePM` is 56 bytes per vertex regardless of `k`, while at `k<=2` a measure has at most
one string. A 6-byte (then 8-byte) `NodePMSmall`, selected per attempt via
`template<bool SMALL>`, cut the working set on the 751k-vertex game from 42 MB to 4.5 MB
— a 9–14x reduction — with all measure logic unchanged (records are widened to the full
8-lane `Measure` for `prog`/`compare`).

Gate was ">=10% on JURD, no loss on VB". Measured:

| corpus | delta vs. value-type build |
|---|---|
| JURD | **+1.0%** |
| VB | +1.7% |
| RND | −6.7% |

Profiling showed exactly where it went: the init/marking work dropped from 979 to 462
samples (the memory saving is real), and `lift` rose from 4045 to 4528 samples (the
widening cost). **Writing a single lane of a `fixed_size_simd` round-trips through
memory**, and three such writes per load cost about as much as a 14x memory reduction
saves. Two attempted fixes — broadcast-and-mask instead of per-lane inserts, and padding
the record to 8-byte alignment — each left JURD still ~6.5% slower than without the
change.

The conclusion is that **Jurdzinski lifting is compute-bound, not memory-bound**. RND
improves because random games have poor locality; JURD does not because its access
pattern is already cache-resident. Reverted per the gate.

## What the evidence says to do next

Batching *vertices* (one SIMD lane per in-flight vertex, fed from the worklist) is the
remaining idea, and it is the one axis nobody has tried. Two prior findings frame it:

- `gaperez64/claude/parallelize-simd-successors-gofLh` (516 lines, never merged) batched
  *successors* with interleaved AVX2 packing. Its history is a losing fight ending in a
  `batch_can_skip` scalar pre-filter added to claw back time. The degree distribution
  explains why: `jurdzinskigame(100,100)` has mean in- and out-degree 2.66, so an 8-wide
  batch over either edge list runs at **29.3% lane utilisation** and forfeits the scalar
  loop's early exits. Batching the worklist instead gives full lanes by construction.
- Phase 2 warned that per-lane SIMD element access is expensive, which is a direct threat
  to a design built on gather/scatter/blend.

A microbenchmark settles it (`scratchpad/batchtest.cpp`): 8 sequential branchy scalar
`k=2` comparisons versus one 8-wide branch-free comparison *including* the cost of
gathering 8 scattered vertices into transposed lanes.

| working set | batched cost vs. scalar |
|---|---|
| 0.4 MB | 0.09x |
| 3.5 MB | 0.10x |
| 28 MB | 0.24x |
| 224 MB | 0.19x |

Batching is **4–11x cheaper on the kernel at every working-set size**, including fully
cache-resident ones — so the win is branch elimination and SIMD width, not memory-level
parallelism, and it does not depend on the memory hypothesis that Phase 2 disproved. The
gather pays for itself comfortably.

Design consequences for the implementation:

- Lane = vertex, fed from the worklist with **retire-and-refill** (a lane whose successor
  list ends finalises and pulls the next vertex from `Q`), which keeps utilisation near
  100% despite the ragged degree distribution.
- `pindex`, `want_max` and `do_prog` are constant for the whole of a lane's lift, so the
  control parameters never diverge within a lane.
- `prog` must be **vectorised across vertices**, not run scalar per lane with
  extract/insert — that round-trip is precisely what sank Phase 2. At `k<=2` this is
  tractable because a measure has one string, so `prog`'s cross-lane operations (prefix
  sums, `find_last_set`) collapse to scalar ops on a single value, making it a
  per-vertex-independent branch-free function.
- That `k=2` `prog` is the historically buggy code — a `k=2,t=1` specialisation was
  reverted for a real correctness bug (see the note at `prog()`'s head). It needs a
  lockstep differential mode asserting measure equality against the general path after
  every lift, not just end-to-end verification.

Not pursued, and why:

- **MLIR.** It does not attack the measured hot spots: `prog_tmp` — the one component
  whose complexity would justify a compiler — is 0.5% of runtime, and MLIR does not
  change memory layout. Its real value would be generating the kernel family across
  `(k,t)` and vector widths, which `template<int K>` covers on the `(k,t)` axis for free.
  Revisit only if cross-ISA breadth becomes a requirement, or for GPU offload of the
  million-vertex instances.
- **Cutting wasted lifts.** 9.6M lift attempts yield 1.86M successful lifts (5.2:1) on
  `jurdzinskigame(100,100)`. For a min-player vertex whose measure is the min over
  successors, if the cached argmin is not the successor that just changed, the min cannot
  have changed and the whole scan is skippable — measures only increase. Caching the
  argmin (`str` is already computed) could remove a large share of the 6.7M scan calls
  and 19.6M successor loads. Orthogonal to batching and likely cheaper; worth its own
  investigation.

## Reproducing

```bash
CXX=/opt/homebrew/bin/g++-16 CC=/opt/homebrew/bin/gcc-16 \
  cmake -S . -B build-release-gcc -DCMAKE_BUILD_TYPE=Release
cmake --build build-release-gcc -j8 --target solve

# per-attempt counters and timings
./build-release-gcc/oink --strpm-simd <game> --no -t 1

# correctness
./build-release-gcc/oink --strpm-simd <game> --no -v
```

`boost::filesystem::exists`/`directory_iterator` are broken on this machine (unrelated,
pre-existing), which silently no-ops `test_solvers`' directory scan, so correctness was
driven through the `oink` CLI in shell loops instead. `test/test_solvers.cpp` also still
needs an unrelated `boost::process` v1/v2 fix before it will compile here.
