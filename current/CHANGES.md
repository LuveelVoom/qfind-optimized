# qfind 2.4b — optimized version: what changed and why

This document explains every change made to `qfind.c` and `common.h` relative to the original qfind 2.4b (Matthias Merzenich). It covers why each change was made, how it works, what it measured, and how correctness was checked. The full line-by-line changes are in `qfind.c.diff` and `common.h.diff`.

**How to read the certainty labels**
- **[measured]**: timed or counted on the test machine.
- **[verified]**: checked by a test, such as identical output or a sanitizer run.
- **[inferred]**: reasoning only, not tested.

**About the timings**
- The test machine was a shared 2-core Linux VM, single thread (`-t 1`), noisy by about ±5–8%.
- Unless a run says otherwise, times are the minimum of several interleaved runs.
- Your machine will give different absolute numbers. On your Mac the speed gain was about 10%, versus about 25–30% here.

---

## 1. Summary

| Area | Result | Search results |
|---|---|---|
| Speed (speeds up to c/5) | about 25–30% faster, e.g. c/5 w8 odd: 5.5 s → 4.0 s [measured] | identical output [verified] |
| Speed (speeds above c/5) | a further 0–36% from turning the lookahead cache off by default [measured] | identical output [verified] |
| Table memory | about 40% less by default, about 57% less with `-DPACKLISTS` [measured] | identical table contents [verified] |
| Table memory growth (width ≥ 10) | **factored tables**: growth per extra width drops from about 7× to about 2.8×. Width 10: 1.48 GB → 6 MB; widths 11–12 now run in 16–26 MB where the old layout ran out of memory (section 4.4) [measured] | exhaustive; identical search when lists are sorted (`-DFACTORED_SORTED`) [verified] |
| Lookahead cost growth (width ≥ 12) | **split lookahead**: each lookahead test is solved on half rows joined on their 2 overlap bits, so its cost no longer grows with width. Equal-work time vs the previous lookahead: 2.5× slower at widths 7–8, about even at 12, 2.4× faster at 13 (first deepening). Used automatically from width 12 (section 3.7) [measured] | same answer on every test, so the same search [verified] |
| Split lookahead cache (width ≥ 12) | **shared by all threads**, 2^22 entries per half by default (1152 MB total), with prefetching: about 1.35× faster than the shared 2^20 cache and about 1.6× faster than the original per-thread cache at width 12 (sections 3.10, 3.12) [measured, M4 Pro] | same answers [verified] |
| Lookahead membership tests (widths 4–9, per-pair tables) | the innermost row33 test reads the small factored half-masks instead of a per-pair table opened only for that test: c/5 w9 odd 1.13× faster (1 and 4 threads), c/3 w9 1.04–1.06×, widths 7–8 1.02–1.04× (section 3.13) [measured, M4 Pro] | identical output [verified] |
| Bugs | 22 fixed; 19 also exist in the original (section 6) | no change on valid input [verified] |
| Readability | cryptic names renamed, an explanatory overview and comments added, dead code removed (section 7) | no change [verified] |
| New features | progress readouts with time and memory; `--lookahead-depth`; `--tables`; `--lookahead-method` / `--lookahead-budget`; `--help` performance notes | — |

**Exhaustiveness.** The program still finds a ship whenever one exists within the search limits.
- Every speed change re-implements the same test or the same data in a faster way.
- The default build was checked against the original in two ways:
  - It produces identical output on 11 test searches.
  - With `-DLA_VERIFY`, it runs both the fast and the original lookahead on every candidate and aborts on any disagreement. It never disagreed.

---

## 2. Building and options

**Building** works as before:

```
gcc -std=c11 -fopenmp -march=native -O3 -Wall -Wextra -o qfind qfind.c
```

| Compile flag | Meaning |
|---|---|
| `-DPACKLISTS` | Bit-packed lookup tables: about 24% less table memory than the default, about 25% slower. Use it only when tables would not otherwise fit in RAM. |
| `-DLA_VERIFY` | Check every lookahead against the original algorithm and abort on any difference. Slow; for auditing only. |
| `-DPERIOD=p -DOFFSET=o` | Fixed-speed build (unchanged). The lookahead cache is now compiled out unless `-DFORCECACHE` is also given. |
| `-DINNER_SPLIT` | Experimental: faster innermost lookahead test with factored tables and the cone lookahead (section 3.8). |
| `-DFACTORED_SORTED` | Testing only: sorts factored lists into the per-pair order, so `--tables factored` reproduces the original search exactly. |

| Runtime option | Change |
|---|---|
| `-c N` (lookahead cache) | **Default changed:** the cache is now off at all speeds (it used to be on above c/5). `-c N` still turns it on. |
| `--lookahead-depth 1\|2\|3` | **New.** 3 is the full test (the default). 1 or 2 prune less per candidate but cost less (section 5). |
| `--lookahead-method auto\|cone\|split\|hybrid` | **New.** How the full lookahead test is computed. auto (the default) uses split from width 12 and cone below (section 3.7). All methods give identical answers. |
| `--lookahead-budget N` | **New.** For `hybrid`: work units before a test switches from cone to split (default 256). |
| `--split-cache-bits N` | **New.** Split and hybrid: 2^N shared cache entries per half (default 22 = 1152 MB; 8–30). Larger is faster from width 11 up (section 3.10). |
| `--tables auto\|pair\|factored` | **New.** The lookup table representation. auto (the default) uses per-pair tables below width 10 and factored tables from width 10 (section 4.4). |
| `-w` | Now rejects widths above 14, which used to crash. |
| `-t` | Now rejects values ≤ 0. The default is at least 1 (the old default was 0 on a 1-CPU machine). |
| `-q` | Now requires at least 4 (smaller values crashed). |
| `--help` | Now ends with "Performance notes" describing this build. |

### 2.1 How program options affect speed and memory

These are qfind's own options (compile-time `-D` macros and runtime flags), including the ones that were already in the original. Compiler flags are not covered.

**Measurement setup.**
- **Search A:** `-r B3/S23 -v c/5 -w 8 -s odd -f 1`, stopping at the first ship (none exists at this width, so it runs to completion).
- **Search B:** `-r B3/S23 -v c/4 -w 8 -s even`, a complete search (2 ships).
- Runs were single-threaded (`-t 1`) and each figure is the best of 2 runs. Peak memory is from the new "Total time … peak memory" line.
- The ship counts were the same in every run, except where noted.
- These are two small searches on a noisy VM, so treat differences under about 10% as noise. Larger searches can behave differently.

**Compile-time options (`-D...`)**

| Option | In original? | Speed (A / B) | Memory | Notes |
|---|---|---|---|---|
| (default build) | — | 3.50 s / 2.52 s | 30 / 32 MB | |
| `-DPACKLISTS` | new | 3.76 s / 2.70 s (+7%) | 22 / 22 MB | Packs table entries into `width` bits. About 24–31% less table memory; the slowdown is about 7% here and up to about 25% on other searches [measured]. Use only when tables don't fit in RAM. |
| `-DLA_VERIFY` | new | 5.59 s / 3.56 s (+41–60%) | same | Runs both lookaheads and compares them. Auditing only. |
| `--lookahead-method split` (runtime option; auto from width 12) | new | — | 1152 MB shared by all threads (`--split-cache-bits`, 3.10) | See 3.7: 2.5× slower at widths 7–8, about even at 12, 2.4× faster at 13 at equal work. |
| `-DPERIOD=p -DOFFSET=o` (fixed speed) | yes | c/5: 3.25 s vs runtime-speed build 3.17 s; c/4: 2.27 s vs 2.39 s | same | No longer measurably faster than giving `-v` at run time. The original's comment says it may be "slightly faster"; here the difference was within noise. |
| `-DFACTORED_SORTED` | new | slower than unsorted factored (sorting each list) | same as factored | Testing only (exactness check). |
| `-DFORCECACHE` (with a fixed speed) | yes | c/4: 2.79 s vs 2.27 s without (+23%) | +1 MB | Compiles the lookahead cache in. Now never recommended. |
| `-DNOCACHE` | yes | — | — | Compiles the runtime cache code out. For fixed-speed builds this is now the default. |

**Runtime options**

| Option (default) | In original? | Speed A | Speed B | Memory | Notes |
|---|---|---|---|---|---|
| `-q N` queue bits (15) | yes | 12: 3.22 s, 18: 4.58 s, 20: 5.40 s | 12: 2.05 s, 18: 3.77 s, 20: 5.22 s | about 12 bytes × 2^N: N=15 0.4 MB, 20 12 MB, 24 200 MB, 31 about 25 GB [calculated from the allocations] | A smaller queue deepens more often and was faster here (−8%/−19% at 12). A bigger queue was slower. Don't go too small at wide widths: one node can add up to 2^width children, and if the queue overflows the search stops with "Exceeded … node limit". Minimum 4. |
| `-h N` hash bits (15) | yes | 0: 3.31 s, 12: 3.44 s, 20: 3.72 s, 24: 4.19 s | 0: **over 240 s (killed)**, 12: 2.55 s, 20: 2.76 s, 24: 3.62 s | 4 bytes × 2^N: 24 adds 64 MB (measured 30→94 MB) | `-h 0` turns off duplicate removal. That was harmless for A but made B more than 100× slower, so never use it. Values above the default were slower. |
| `-b N` base bits (4) | yes | 2: 3.42 s, 8: 3.19 s | 2: 2.56 s, 8: 2.41 s | 4 bytes × 2^(q−N), negligible | Within noise. |
| `-c N` lookahead cache MB per thread (now 0) | yes (default was 32 above c/5) | 32: 5.71 s (+63%) | 32: 3.48 s (+38%) | about N MB per thread (measured +32 MB) | Now off by default; see 3.6. |
| `-i N` minimum deepening (3) | yes | 1: 3.78 s, 8: 3.29 s | 1: 3.15 s, 8: 2.26 s | none | 1 was slower (+8%/+25%). 8 was slightly faster (−6%/−10%), but earlier tests at 3–12 were within noise. |
| `-n N` first deepening amount | yes | 10: 3.35 s | 10: 2.60 s | none | Within noise. |
| `--fixed-depth N` | yes | 5: 3.52 s | 5: 2.46 s | none | Within noise. |
| `-g N` minimum extension length to save (0 = save all) | yes | 100000 (never save): 5.06 s (+45%) | 3.05 s (+21%) | saved extensions use a little memory | Saving and reusing depth-first extensions is a big win; keep the default. |
| `--lookahead-depth N` (3) | new | 2: 3.02 s, 1: 3.02 s (default 3.22 s in that batch) | 2: 1.91 s, 1: 1.83 s (default 2.32 s) | none | Faster on these two, but slower at c/6 and c/7 (section 5). Time it per search. |
| `--tables pair\|factored` (auto) | new | — | — | factored: a few MB instead of 4^width-scale tables | See 4.4. Below width 10, factored was 0.86–1.9× the time of per-pair on complete searches. At width 10 it reached about the same depth in 90 s with the same search order, and on first-ship searches it can be much faster because no tables have to be built (c/3 w10: 9.2 s → 0.56 s). |
| `-t N` threads (CPUs − 1) | yes | 2: 1.83 s (1.9×) | 2: 1.85 s (1.4×) | the cache is per thread; other per-thread buffers are small | The test VM has only 2 cores. |
| `--disable-early-exit` | yes | 1.90 s vs 1.83 s (`-t 2`) | 1.77 s vs 1.85 s | none | Within noise. |
| `-m N` memory limit | yes | — | — | caps table memory | Only aborts when exceeded; no speed effect. |
| `-f`, `-e`, `-l`, `-j`, `-a`, `-d`, `-p`, output flags | yes | — | — | — | Control what is searched, saved or printed. No direct effect on speed per node. |

**Biggest memory users**
- **Lookup tables:** they grow roughly as 4^width. With every table built: w=9 208 MB, w=10 1391 MB (section 4.3). The pair index alone is 8 bytes × 4^width (w=12: 128 MB, w=14: 2 GB), but only the parts touched are committed to RAM.
- **Others:** after the tables come the queue (`-q`), the hash (`-h`) and the cache (`-c` × threads).

**What the measurements suggest for speed**
- Keep the defaults for `-h`, `-g` and `-c`.
- Try `-q 12`–`14` and `--lookahead-depth 1|2` on your own long searches, and compare deepening times using the progress readout.
- Use as many threads as you have performance cores.

---

## 3. Speed changes

### 3.1 Per-level lookahead (`lookAheadFast`, `LookaheadLevel`) — the main speed gain

**What the lookahead does.** Before a candidate row is accepted, `lookAhead()` checks that a small "cone" of future rows could exist. If no cone exists, the candidate cannot be extended, so dropping it never loses a ship.

The notation: `rowLG` is the row L periods further along the ship, G generations ahead, and `listLG` is the list of values it may take, read from the lookup tables. The cone exists if there are
- row11 in list11, row12 in list12 and row13 in list13,
- row22 in list22(row12,row11) and row23 in list23(row13,row12),
- such that the table for (row13,row23) can produce row22.

Profiling showed this function took about 70% of all instructions [measured].

**What was slow.** The depth-first search tries many sibling candidates at the same position, one after another. Only list11 depends on the candidate. list12, list13 and the tables behind list22/list23 depend only on rows already fixed above it. The original recomputed all of them for every candidate.

**What it does now.** Everything that doesn't depend on the candidate is kept per search-tree level in a `LookaheadLevel` record, built lazily in two stages:
- **Stage 0** (`laBuildStage0`): the rows the tables are looked up by, plus the table that gives list11. About 70% of candidates have no row11 at all and are rejected by a single bit test [measured].
- **Stage 1** (`laBuildStage1`): list12, list13, and slots for the per-row12 and per-row13 tables, which are fetched on first use.

A level's record is discarded whenever the search enters that level with a new parent (`laForgetLevel`), and for all levels at the start of each depth-first search (`laForgetAll`).

**Why it gives the same answers.** It asks exactly the same existence question; only the loop order differs (row12 outermost). The original is kept as `lookAheadReference()`. It is still used when the lookahead cache is on and when TRIPLEOFF == PERIOD, and `-DLA_VERIFY` compares the two on every call.

### 3.2 A bitmask in front of every table

Each table now starts with a bitmask of its non-empty buckets ("which result rows have at least one candidate"). The lookahead's existence tests become single bit tests instead of reading two offsets (`BITSET`). This is part of the compact layout described in 4.1.

### 3.3 Faster table construction (`buildTable`, was `makeRow`)

Lookup tables are built on demand, one per row pair (row1,row2). The building code now:
- evaluates the rule from precomputed "low half" and "high half" partial results (as before), but fused into a single pass;
- builds the table in a private per-thread buffer outside the lock, and only copies it into shared memory under the lock;
- hashes tables for deduplication (identical tables from different row pairs are stored once) with a 4-lane 64-bit hash (`hashTableWords`) instead of a byte-serial chain. Only the hash's distribution matters, not its values.

### 3.4 Start-up and dedupe structures

- **`pairTable` (was `gInd3`)**, one pointer per row pair: no longer filled with NULL in a loop. `calloc` already returns zeroed memory, and untouched pages are never committed. At width 10 that is 8 MB that no longer has to be written at start-up.
- **`tableDedupe` (was `rowHash`)**, the set of distinct tables: it starts at 1024 entries and doubles as needed. The original allocated it at full size (2·4^width entries) up front.

### 3.5 Huge pages for table memory (Linux only)

Table memory is taken from chunks of at least 8 MB, 2 MB-aligned via `mmap`, with transparent huge pages requested (`madvise(MADV_HUGEPAGE)`). Before this, page faults were about half of table-build time at wide widths [measured]. It only takes effect if `/sys/kernel/mm/transparent_hugepage/enabled` is `always` or `madvise`. macOS falls back to `calloc`. This is probably why the gain on your Mac was smaller [inferred].

### 3.6 Lookahead cache off by default

The original turned its lookahead cache on for speeds above c/5. The cache stores results keyed by the input lists, but it uses the original (slow) lookahead. Compared with the fast lookahead (3.1), turning the cache off was never slower in 19 test searches [measured]:

| search | cache on (old default) | cache off (new default) |
|---|---|---|
| c/4 w9 even, first ship | 14.0 s | 9.4 s |
| c/3 w10 even, first ship | 9.9 s | 6.8 s |
| 2c/5 w8 odd, first ship | 0.71 s | 0.50 s |
| 2c/7 w7 odd | 7.0 s | 6.1 s |
| c/3 w9 odd | 94.8 s | 94.1 s |

The output is identical either way, because both paths give exact answers. In the fixed-speed build at c/4: 3.7 s → 2.4 s.

### 3.7 Split lookahead (`--lookahead-method split`, automatic from width 12)

**The problem: the lookahead gets more expensive per test as width grows.** The lookahead (3.1) searches for a "cone" of six unknown rows: row11, row12, row13, row22, row23 and row33. The candidate lists it combines get longer with width, and they multiply. The cost of a test became increasingly lopsided [measured, c/5 odd searches]:

| width | tests costing more than 256 work units | their share of all lookahead work |
|---|---|---|
| 8 | about 0.15% of tests | about 32% |
| 10 | about 0.4% of tests | about 68% |
| 12 | about 0.6% of tests | about 92% (tests over 4096 units alone: about 63%; the largest over 250,000) |

So at width 12 the lookahead's total cost is dominated by a growing tail of very expensive tests.

**The idea.** The factored tables (4.4) rely on one fact: each rule evaluation splits into a low half and a high half that share only 2 "overlap" bits of each input row. The lookahead cone is just six such evaluations chained together. So:

> a full cone exists **if and only if** there is a choice of the 2 overlap bits for each unknown row such that a *low-half* cone and a *high-half* cone both exist with those overlap bits.

This holds because the low and high halves of every row share nothing except their overlap bits.

**How it is computed** (`lookAheadSplit`, `splitHalfSetT`, `splitMemo` in `qfind.c`).
- A half cone involves only half rows (`width/2 + 1` bits). Its candidate lists hold about 2 entries on average **at every width**.
- For each half, the program computes the set of overlap vectors for row12, row13, row22, row23 and row33 (10 bits, so a 1024-bit set) for which a half cone exists. These sets depend only on the half rows already fixed plus the half of row11, so they are kept in a cache shared by every level, node and thread (section 3.10; it was a per-thread cache of 2 × 65536 entries, about 18 MB per thread, before that change).
- A test then passes if, for some overlap value of row11, the union of the low-half sets and the union of the high-half sets have a vector in common.

**Why the cost stops growing.** The work per test fell as width grew, because more of the half-cone sets are already in the cache at larger widths [measured, c/5 even]:

| width | cache misses per lookup | work per miss | work per test |
|---|---|---|---|
| 8 | 34% | 98 | 38 |
| 10 | 17% | 153 | 28 |
| 12 | 6% | 300 | 17 |

The previous lookahead's work per test grows instead.

**Measured** (c/5 even, single thread). Both methods make identical decisions, so the k-th deepening is identical work:

| width | deepening | previous lookahead (cone) | split | split ÷ cone |
|---|---|---|---|---|
| 7 | 12 | 0.40 s | 1.00 s | 2.5 |
| 8 | 24 | 3.7 s | 9.2 s | 2.5 |
| 9 | 16 | 30.2 s | 55.5 s | 1.84 |
| 10 | 13 | 61.3 s | 109.2 s | 1.78 |
| 11 | 7 | 61.8 s | 83.2 s | 1.35 |
| 12 | 4 | 144.8 s | 137.3 s | 0.95 |
| 13 | 1 | 52.4 s | 22.1 s | 0.42 |
| 13 | 2 | over 270 s | 193.7 s | under 0.72 |

The ratio falls with every extra width, so this changes how cost grows rather than giving a fixed saving. At widths up to 11, split is slower: it must compute a full 1024-vector set, where the previous method stops at the first cone it finds. So `auto` uses split only from width 12 (`SPLIT_AUTO_WIDTH` in `common.h`).

**Hybrid (`--lookahead-method hybrid`).** This runs the previous lookahead with a work budget (`--lookahead-budget`, default 256 units). When a test exceeds the budget, it gives up and lets split answer that test. Measured at the same deepening:
- width 8: no gain;
- width 10: +1.6% (budget 1024) and +9% (256) vs the previous lookahead;
- width 11: +10–14%;
- width 12: −8%, about the same as split;
- width 13: −7% or better, but slower than pure split (251 s vs 194 s at deepening 2).

It is never clearly the best choice, so it is not used automatically. It is kept as an option for experiments.

**Notes.**
- `--lookahead-depth 1` and `2` always use the previous lookahead; split applies only to the full test.
- The lookahead cache (`-c`) path is unchanged.
- Correctness: `-DLA_VERIFY` compares every test against the original lookahead. It passed with split and with hybrid (budgets 0, 8, 16, 64, 256), and the regression searches give output identical to the original (section 8) [verified].
- The split cache (1152 MB by default, see 3.10) is not included in the "tables+index" readout; the start-up summary prints its size.

### 3.8 Inner split (`-DINNER_SPLIT`, experimental)

**What it does.** With factored tables (4.4) and the cone lookahead (3.1), the innermost lookahead check asks: is there a row23 in list23 and a row22 in list22 for which the table of (row13, row23) has a bucket for row22, i.e. some row33 exists? The old code joins both lists into full rows and tests every pair, costing |list23| × |list22|.

Both lists are joins of a low-half group and a high-half group that match on the 2 overlap bits. A row33 exists exactly when, for some overlap value u of row23 and v of row22, the low and high halves agree on row33's overlap. So for each (u, v) the new code ORs the 4-bit row33-overlap masks over all low-half pairs, does the same over all high-half pairs, and checks whether the two results share a bit. That costs |L23|·|L22| + |H23|·|H22| (half-lists) instead of the full product, and it still stops at the first hit. Row11s whose list22 is empty are filtered out once per row12, as before.

Code: `factoredInner()` and the `#ifdef INNER_SPLIT` block in `lookAheadFastT()` (`qfind.c`). It is used only when compiled with `-DINNER_SPLIT`, with `--tables factored` (the default from width 10) and the cone lookahead (the default below width 12; at 12 and above add `--lookahead-method cone`).

**Measured** (single thread, factored tables, cone lookahead; both versions make identical decisions, so more deepenings in the same time means faster) [measured]:

| search | without | with `-DINNER_SPLIT` |
|---|---|---|
| c/5 w8 even, complete | 2.32 s | 2.35 s (same) |
| 3c/7 w10 odd, 60 s | 28 deepenings | 41 deepenings |
| c/5 w12 even, 60 s | 3 deepenings | 4 deepenings (same as `--lookahead-method split`) |

**Checking (quick only).** The `-DLA_VERIFY` build, which compares every lookahead with the original, found 0 mismatches on 3 short searches (c/5 w8 odd, c/4 w7 odd, 3c/7 w7 odd) [verified]. The full regression was not rerun.

### 3.9 Parallel breadth-first lookaheads (`bfsPrecompute`, 4+ threads)

Only the deepening step used threads; the breadth-first step (`process()`) ran single-threaded, and almost all of its time is spent in `lookAhead()`. Now, before processing, the lookaheads of the next 256 queue nodes are computed on all threads. `process()` then uses the stored results and does the rest (duplicate check, enqueue, output) single-threaded in the original order, so the search is unchanged; the batch is dropped whenever the queue is compacted.

- Output identical to the previous build with `-t 2 --disable-early-exit` on 2 searches [verified]. With early exit on, multi-threaded runs differ between runs anyway.
- On this 2-core VM: 2–5% slower (c/4 w8 even, c/6 w7 even) [measured], so it switches on only with 4 or more threads. At c/5 w9 the breadth-first step was about 3% of wall time on 2 threads; on 16 cores it would be roughly a third without this change, so the expected gain there is about 1.3–1.5× in wall time [inferred, untested].

### 3.10 Shared split-lookahead cache (`--split-cache-bits`)

**What changed.** The split lookahead (3.7) caches half-cone sets keyed by (half-row context, row11 half). That cache was per thread, with 2^16 entries per half (about 18 MB per thread). It is now one cache shared by all threads, set with `--split-cache-bits N`. The default was first 2^20 entries per half (288 MB) and is now 2^22 (1152 MB in total, whatever the thread count; see the second table below).

**Why.** At widths 12–13 about 65–85% of the time was in the split lookahead, mostly in computing half-cone sets on cache misses [measured, `sample` profiles]. A larger cache misses less, and a shared one lets every thread reuse sets computed by the others. Making the private caches larger instead would multiply memory by the thread count (2^22 entries cost 11 GB at 10 threads).

**How it stays correct with threads.** Entries are written without locks, so each has a sequence counter (a seqlock; `splitMemoOr()` in `qfind.c`):
- 0 means never written; odd means a thread is writing it.
- A writer makes the counter odd with a compare-and-swap, writes the key and set, then makes it even again. A thread that loses that race just doesn't store its set.
- A reader copies the entry and checks that the counter was even and unchanged; otherwise it treats the entry as missing and computes the set itself. So a half-written entry is never used.
- The old code marked empty entries with the key ~0. An empty entry is now marked by counter 0, which also removes a theoretical case where a real key equal to ~0 would have read an unwritten set.

**Measured** (Apple M4 Pro, clang build, c/5 even; both versions make identical decisions, so reaching a deepening sooner means faster) [measured]:

| search | old (per-thread, 2^16) | shared, 2^20 (default) | shared, larger |
|---|---|---|---|
| w12, 1 thread, 3rd deepening | 59.9 s | 48.2 s (1.24×) | 40.8 s with 2^22 (1.47×) |
| w13, 1 thread, 1st deepening | 14.3 s | 13.0 s (1.10×) | 12.4 s with 2^22 (1.15×) |
| w13, 10 threads, 4th deepening | 43.6 s | 37.0 s (1.18×) | 36.4 s with 2^23 |
| w14, 10 threads, 2nd deepening | over 60 s | 59.7 s | 45.4 s with 2^23 |

**Default raised to 2^22** (single thread, c/5 w12 even, time to the 4th deepening, two rounds): 2^20 66.3 / 68.5 s, 2^22 57.9 / 58.6 s (about 1.15×), 2^23 57.1 / 60.4 s, 2^24 62.1 / 59.7 s. Beyond 2^22 there was no further gain, so the default is 2^22 (1152 MB) [measured]. Note: these and the 3.12 timings ran while another 12-thread search was using most of the machine, so they are single-thread comparisons only.

At widths 10–11 the cone lookahead is still as fast or faster than split with the larger cache (c/4 w11 odd: 16 vs 14 deepenings in 60 s; 3c/7 w10 odd: 30 vs 27), so `auto` still switches to split at width 12.

**Checking.**
- `-DLA_VERIFY` (every split answer compared with the original lookahead) with 8 threads and a deliberately tiny cache (`--split-cache-bits 10`, to force collisions and write races): 0 mismatches on c/5 w9 odd, c/4 w8 even, 3c/7 w8 odd, c/6 w7 even and 2c/5 w9 even [verified].
- The 12 regression searches with `--lookahead-method split` (clang with `--split-cache-bits 12`, and gcc with the default) and the standard `tools/regress.sh`: output identical to the original [verified].
- ThreadSanitizer (clang, 4 threads, c/4 w7 even, split, `--split-cache-bits 8`): no reports in the split cache. It did report races in the older per-pair table builder (`tableAlloc`/`buildTable`); those predate this change and were not investigated [verified].

### 3.11 Building on a Mac

gcc on macOS has no native thread-local storage: every `_Thread_local` access calls `__emutls_get_address`, which took 8–10% of the time in profiles [measured]. Apple clang has native thread-local storage, and Homebrew's `libomp` supplies OpenMP:

```
brew install libomp
clang -std=c11 -O3 -mcpu=native -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include \
      -L$(brew --prefix libomp)/lib -lomp -DINNER_SPLIT -o qfind qfind.c
```

Measured against Homebrew gcc-16 on an M4 Pro, single thread, same searches: c/5 w8 odd 5.6 → 5.1 s, c/4 w9 even 11.5 → 11.0 s, c/6 w7 even 12.0 → 10.9 s, c/5 w12 even 35.0 → 29.3 s to the 3rd deepening (4–16% faster) [measured]. The output is the same; the regression passes with both compilers [verified].

Two other build options were measured (single thread, 3 interleaved runs each) [measured]:
- **`-flto`:** no reliable gain. qfind is one translation unit (`qfind.c` includes `common.h`), so there is nothing extra to optimize across files. gcc: c/5 w8 odd 5.7 → 5.6 s, c/6 w7 even 10.9 → 12.0 s; clang: 5.1 → 4.8 s and 10.7 → 10.7 s.
- **Profile-guided optimization with clang** (`-fprofile-instr-generate`, run a few short searches, `xcrun llvm-profdata merge`, rebuild with `-fprofile-instr-use`): 0–9% faster (c/4 w9 even 1.09×, 2c/7 w7 odd 1.05×, c/6 w7 even 1.03×, c/5 w8 odd 1.00×). Optional; it adds a build step for a small gain.

### 3.12 Prefetching split-cache entries

With a 1 GB cache almost every lookup misses the CPU caches, and `lookAheadSplit()` did its lookups one after another, each waiting for memory. It now first computes the address of every entry the test may read (one per half of row11, for each allowed overlap value) and asks the CPU to start loading them (`splitPrefetch()`, two cache lines per 144-byte entry), then does the test as before. The memory accesses overlap.

Measured (single thread, c/5 w12 even, time to the 4th deepening, two rounds): 59.7 / 59.5 s without, 49.8 / 51.6 s with prefetching, about 1.15–1.2× faster [measured]. Prefetching only changes timing, so the answers are unchanged (same checks as 3.10, rerun with this change: 0 `-DLA_VERIFY` mismatches with 4 threads and `--split-cache-bits 10`; the regression searches identical to the original, standard and with `--lookahead-method split`, clang and gcc) [verified].

### 3.13 Membership tests from the half-masks (`useFactoredInner`, widths 4–9)

The idea came from someone else's edit of this program (shared by the user on 20 October 2026), which reported about 14% at c/5 w9.

**What it does.** In the full lookahead, the innermost test asks, for each row23, whether any row22 is a row reachable from the pair (row13, row23). With per-pair tables that meant opening table (row13, row23), a fresh pair almost every time, just to read its mask. The tables form a pool of hundreds of MB, so that read was usually a cache miss. The factored half-masks (section 4.4) hold exactly the same set (the low and high halves of a row are independent apart from their two overlap bits), take a few hundred kB and stay in cache. Per-pair searches now also build the half-tables and answer that test from them (`factoredBases()` + `factoredHasRowBase()`, the index bases computed once per pair). The candidate lists still come from the per-pair tables, so the search visits rows in the same order.

Applied in both places with that loop: `lookAheadFastT()` (the `half` argument; `lookAheadFastPairHalf()`) and `lookAheadReference()` (used for every c/3 lookahead and with the lookahead cache). Also from the same edit: the masks of the row12 tables are prefetched before the loop that reads them, the two lookahead stage builders are inlined with the table representation as a constant, and `lookaheadDepth` is read once per call.

**Where it does not help.** The other membership tests in the lookahead (stage 0's row11 test, the row23 test, skipping empty row22 lists) were tried the same way and were 0–3.5% slower: each of those tables is read right afterwards for a list, so the cache miss happens anyway and the extra test only adds work [measured].

Measured (single thread unless stated, two runs at a time on an M4 Pro, best of two):

| Search | before | after | |
|---|---|---|---|
| c/5 w9 odd | 36.1 s | 31.9 s | 1.13× |
| c/5 w9 odd, `-t 4` | 9.8 s | 8.7 s | 1.13× |
| c/3 w9 even / odd | 4.37 / 38.2 s | 4.12 / 36.6 s | 1.06× / 1.04× |
| c/4 w8 even, 2c/5 w8 odd, c/6 w7 even | | | 1.02–1.04× |
| c/4 w10 even (factored tables), c/5 w8 with `-c 20` | | | about 1.00× |

The gain grows with width (the per-pair pool grows) and with other programs competing for the CPU caches: with five searches running at once, c/5 w9 odd was 1.20× faster.

Verified: `tools/regress.sh` all OK against the original 2.4b (default, `-DPACKLISTS`, `-DLA_VERIFY`, `-DFACTORED_SORTED` builds, gcc); identical output to the previous version on the 8 searches above and on c/5 w9 odd with `-t 4 --disable-early-exit` [verified].

---

## 4. Memory changes

### 4.1 Compact table layout (default)

**Original layout:** one 16-bit offset for every possible result row (2^width offsets per table), followed by the lists. Only about 15% of buckets are non-empty [measured].

**New layout** (16-bit words, documented above `bucketRank()` in `common.h`):

```
mask[MASKWORDS]   bit r set iff bucket r is non-empty
pre[PREWORDS]     number of non-empty buckets before each 64-bit chunk of the mask
total             table size in words
dbase             where the list data starts
off[ne+1]         start of each NON-EMPTY bucket (plus the end)
data              the lists, in bucket order (same contents and order as before)
```

A bucket is located by counting the set mask bits before it (`bucketRank`, a popcount). The access functions are `bucketFind`, `bucketRead` and `bucketEntries`.

### 4.2 Optional bit-packed lists (`-DPACKLISTS`)

Each list entry is stored in `width` bits instead of 16, and offsets count entries rather than words. Lists are decoded on demand into buffers: one per search level (`levelDecodeBuf`) plus a few inside the lookahead.

The lookahead cache used to identify a list by its memory address. A decoded list has no stable address, so lists are now identified by a 64-bit value (table address × 2^16 + position).

The cost is about 51M small decodes in a 5-second run; the time goes into the chain of dependent lookups, not the bit arithmetic [measured]. Several faster decoders were tried and did not help. It is therefore off by default.

### 4.3 Results

Table memory with every table built (the most a search at that width can need), as counted by the program itself [measured]:

| width | original | compact (default) | packed (`-DPACKLISTS`) |
|---|---|---|---|
| 9 | 333 MB | 208 MB | 144 MB |
| 10 | 2434 MB | 1391 MB | 1052 MB |

### 4.4 Factored lookup tables (`--tables factored`, automatic from width 10)

**The problem.** Real searches build nearly every per-pair table. For c/5 even, peak memory was 206 MB at width 9 and 1.48 GB at width 10, growing about 7× per extra width and heading toward 8×. On this 8 GB machine, widths 11 and 12 were killed for running out of memory before the first deepening finished [measured].

**The idea.** The program already evaluated the rule in two halves that overlap by 2 bits (the low/high split inside `buildTable`):
- the low result bits of row4 depend only on the low `width/2+1` bits of row1, row2 and row3;
- the high result bits depend only on the upper bits of the three rows.

So the per-pair table is just a product of two small pieces that don't depend on the whole pair. The factored representation stores each piece once, globally:
- **low half-table:** for every (row1-low, row2-low), the row3-low values grouped by (row4-low, overlap value);
- **high half-table:** the same for the high halves.

**How it is used.**
- A successor list for (row1, row2, row4) is built on demand by joining the two groups that share the same overlap value.
- The "is this list non-empty" test (used constantly by the lookahead) is one AND of two 4-bit masks, as cheap as the per-pair bit test.

Memory is O(2^(1.5·width)) instead of O(4^width) tables, about 2.8× per extra width instead of about 7×. The core search (BFS + deepening + the same 3-generation lookahead) is unchanged; only where the successor rows come from changes.

**Zero cost when not used.** Each search picks one representation at startup. The hot lookahead is compiled once per representation (`lookAheadFastT` with a constant `factored` argument, through the `TABLE_REF` / `TABLE_HAS_ROW` / `TABLE_ROWS` macros), so per-pair mode runs the same code as before. Measured: 6.45 s vs 6.49 s and 2.24 s vs 2.23 s [measured]. A hybrid that mixed both representations per table pair cost 11–22% even when all tables fitted, so it was dropped.

**Search order.** Factored lists contain the same rows but in a different order. Like changing `-i` or `-q`, this changes timing and which variants of a ship are printed, never whether one is found. With `-DFACTORED_SORTED`, lists are sorted into the per-pair order. The factored search then reproduced the original's output exactly on 13 test searches, including the packed, cache and `LA_VERIFY` builds [verified].

**Measurements** (single thread):

| search | per-pair | factored |
|---|---|---|
| c/5 w8 even, complete | 2.4 s, 32 MB | 2.8 s (+17%), 3 MB |
| c/6 w7 even, complete | 6.9 s, 7 MB | 11.6 s (+68%), 3 MB |
| c/4 w8 even, complete | 2.7 s, 32 MB | 2.6 s, 3 MB |
| 2c/5 w8 even, complete | 0.50 s, 30 MB | 0.43 s, 3 MB |
| c/4 w9 even, first ship (same order, sorted) | 9.1 s, 203 MB | 10.6 s, 4 MB |
| c/3 w10 even, first ship (same order, sorted) | 9.2 s, 1059 MB | 0.56 s, 5 MB (no table building) |
| c/5 w10 even, 90 s, same order | 15 deepenings, 1.48 GB | 14 deepenings, 6 MB |
| c/5 w11 even, 90 s | killed: out of memory (8 GB machine) | 10 deepenings, 16 MB |
| c/5 w12 even, 90 s | killed: out of memory | 4 deepenings, 26 MB |

Half-table memory by width [calculated from the allocation sizes; matches the measured figures above]:
- width 10: 3 MB;
- width 12: about 23 MB;
- width 14: about 180 MB.

The per-pair index alone would be 128 MB at width 12 and 2 GB at width 14, before any tables. Widths up to 14 run, with no sanitizer errors at 13 or 14 [verified].

---

## 5. New features

**Progress readouts.**
- Every "Queue full"/"Depth" line now ends with `[elapsed s, mem X MB (peak Y), tables+index Z MB]`. Current memory use (RSS) is shown on Linux only; other systems show only the peak.
- A "Total time / peak memory" line prints at the end.
- Z is the program's own accounting (tables, the pair index and the hash), not RSS.
- `qfind-original-with-readouts/` holds the original program with the same readouts, for side-by-side comparison.

**`--lookahead-depth 1|2|3`.** This option controls how much of the lookahead test to apply:
- **3** (the default): the full test described in 3.1.
- **2**: require row11, row12, row22, row13 and row23 to exist, but skip the final row33 test.
- **1**: require only row11, row12 and row22 to exist.

(The first version of this option used a weaker depth 1, "some row11 exists", which did not match these measurements. It was much slower, e.g. 26.8 s instead of 7.3 s on c/6 w7 even [measured], and has been corrected.)

Every depth is a necessary condition, so the search stays exhaustive. A weaker depth prunes less but costs less per candidate. Because the queue then evolves differently, which variants of a ship get printed can change. The original behaves the same way when `-i` or `-q` change: c/4 w7 odd prints 12, 13 or 8 ships with the default, `-i 5` or `-q 18` [measured].

Measured times (single runs) [measured]:

| search | depth 3 | depth 2 | depth 1 |
|---|---|---|---|
| c/4 w8 even | 2.45 s | 2.24 s | 1.91 s |
| c/5 w8 odd, first ship | 3.75 s | 3.41 s | 3.38 s |
| c/4 w7 odd | 1.26 s | 1.46 s | 1.44 s |
| 2c/7 w7 odd | 6.09 s | 7.39 s | 7.86 s |
| c/6 w7 even | 7.27 s | 8.88 s | 9.11 s |
| c/7 w6 odd | 1.41 s | 2.02 s | 2.36 s |

Depth 3 is the right default. Depths 1–2 can be worth timing on long c/4 or c/5 searches.

---

## 6. Bug fixes

These were found by reading the code, by two independent review passes, and by running the program under AddressSanitizer and UndefinedBehaviorSanitizer. The sanitizer runs covered dozens of configurations plus dump, load, split and preview.

"Orig" says whether the original 2.4b has the bug. Of the 23, 20 are in the original, 1 partly (#20), and 2 came from the optimizations (#18, #19).

| # | Bug | Effect | Orig | Fix |
|---|---|---|---|---|
| 1 | Table pointers shared between threads with "relaxed" atomics | On weakly ordered CPUs (ARM, e.g. Apple Silicon), a thread could see a new table's pointer before its contents and read garbage. Not possible in practice on x86 [inferred] | yes | release store / acquire load |
| 2 | Default thread count `cpus−1` is 0 on a 1-CPU machine; `-t 0` accepted | zero-sized per-thread buffers, heap overflow | yes | default ≥ 1; `-t` ≤ 0 rejected |
| 3 | `-w` above 14 accepted | header fields overflow, crash [verified] | yes | rejected |
| 4 | `-q` below 4 accepted | negative shifts, crash in compaction [verified] | yes | `-q` ≥ 4 required |
| 5 | `-v "(0,-1)c/4"` accepted | negative offset, endless recursion in `gcd()` [verified] | yes | rejected |
| 6 | `enqueue()` when the queue is full | writes past the end of the queue arrays | yes | stops at the end; the search ends with the normal "Exceeded … node limit" message |
| 7 | `extensionIndex` allocated with `sizeof(pointer)` | accidentally twice the needed size, which hid a real need for 2·QSIZE+1 entries during compaction (the original is off by one at the extreme) | yes | explicit `EXTINDEX_SIZE` with a comment |
| 8 | `-e` initial-rows file read into `char[14]` with unbounded `%s` | stack overflow at width 14; short rows read uninitialised memory | yes | bounded read plus a length check |
| 9 | `parseDumpRoot()` copied fixed lengths from shorter strings and built a pointer before an array | out-of-bounds reads; with a dump root loaded from a file, overflow of a 60-byte buffer [verified by UBSan] | yes | bounded `snprintf` |
| 10 | Index read before its bounds check (`dumpState`, `doCompactPart2`) | read one past the end | yes | condition order swapped |
| 11 | Duplicate-state hash used signed overflow | undefined behaviour on every run [verified by UBSan] | yes | unsigned arithmetic giving identical hash values |
| 12 | Cache key `513*rowsKey` overflows `int` at width ≥ 11 | undefined behaviour [verified] | yes | 64-bit unsigned |
| 13 | `memlimit = -1 << 20` when no `-m` is given | undefined behaviour on every run [verified] | yes | computed only when ≥ 0 |
| 14 | Rules like `B3-/S23` or `B3~4-/S23` (minus sign with no letters) | silently removed B3 or forbade all births | yes | rejected as malformed |
| 15 | `loadState()` trusted sizes in the dump file | a corrupt dump writes out of bounds | yes | bounds checks; load fails cleanly |
| 16 | Saved extension that doesn't match its node | slot leaked (could eventually exhaust slots); unbounded search for the saved row | yes | checked first, slot freed, searches bounded |
| 17 | Lookahead cache keyed an empty list by an address shared with the next non-empty list | could in principle return another query's cached "no" and wrongly prune. **0 errors in 38M checked cache hits** on the original [measured], so not seen in practice | yes | empty list12/list13 answer "no" before the cache is used |
| 18 | `-DPACKLISTS` at width 1 | build buffer too small, heap overflow [verified] | new | larger buffer |
| 19 | possible `memset(NULL, 0, 0)` in the fast lookahead | undefined behaviour (never hit) | new | guarded |
| 20 | table chunk allocation not NULL-checked; memory counted before rounding | crash on out-of-memory; slightly low memory figure | partly | checked; counted after rounding |
| 21 | Minor issues (listed below) | — | yes | fixed or removed |
| 23 | "Longest partial" and "Maximum depth" only tracked during the breadth-first phase | when a deepening kills every node (`-> 0/0`), the search ends reporting a shallow depth and a short partial, though deepening went much further (e.g. c/5 w7 even `-i 40`: depth 9 and a 1-row partial reported, instead of 75 and 14 rows) [measured] | yes | `depthFirst()` also records partials deeper than the current longest (under the existing output lock) |
| 22 | Row-ordering statistics (`genStatCounts`) overflow a signed `int` at width 14 | undefined behaviour [verified by UBSan]; only affects the heuristic row order | yes | unsigned 32-bit counts, giving the same wrapped values the original produced in practice |

The minor issues in #21:
- the base-bits clamp used `>` instead of `>=` (and its message said `-q`);
- huge `-i`/`-n` values were silently cut to 16 bits;
- per-thread buffers were not NULL-checked;
- a latent `1<<31` overflow (unreachable at width ≤ 14);
- subperiod tables were built twice.

On valid input none of these fixes changes which ships are found. The one qualification is #17: it would change results only if the collision it prevents ever occurred, and none occurred in the 38M checked cache hits.

---

## 7. Readability ("unobfuscation")

**Renames.** These are whole-word renames applied to code only, never inside comments or strings, and verified by identical output. Old → new:

| Area | Renames |
|---|---|
| Lookup tables | `gInd3`→`pairTable`, `getoffset`→`tableForPair`, `getoffset2`→`tableFor`, `getoffsetcount`→`getSuccessors`, `makeRow`→`buildTable`, `bmalloc`/`bbuf`/`bbuf_left`→`tableAlloc`/`tableChunk`/`tableChunkLeft`, `rowHash`→`tableDedupe`, `hashRow`→`hashTableWords`, `gWorkConcat`→`buildScratch`, `gcount`→`rowOrderWeight`, `valorder`→`rowOrder`, `nttable`/`nttable2`→`ruleTable`/`ruleTable3x3`, `fasterTable`→`makeRuleTable3x3` |
| Table access | `tabRank`/`tabFind`/`tabGet`/`tabEntriesB`→`bucketRank`/`bucketFind`/`bucketRead`/`bucketEntries`, `LISTCAP`→`DECODE_CAP` |
| BFS queue | `rows`→`queueRows`, `base`→`queueBase`, `ROFFSET`→`PARENT_OFFSET`, `MAXOFFSET`→`MAX_PARENT_OFFSET`, `putnum`→`printCompactNumber` |
| Duplicate detection | `hash`→`seenHash`, `hashFunction`→`seenHashIndex`, `same`/`sameFlipped`→`sameState`/`sameStateFlipped`, `flip`→`flipRow` |
| Search state | `pRows`→`pathRows`, `pInd`/`pRemain`→`candEnd`/`candLeft` (end of a level's candidate list / candidates left), `pPhase`→`rowPhase`, `riStart`/`numRows` (in `process`)→`candidates`/`numCandidates` |
| Saved extensions | `deepRows`→`extensionRows`, `deepRowIndices`→`extensionIndex`, `deepQHead`/`deepQTail`/`oldDeepQHead`→`extQHead`/`extQTail`/`prevExtQHead`, `saveDepthFirst`/`reloadDepthFirst`→`saveExtension`/`reloadExtension` |
| Lookahead | `lookAhead`'s internals → `lookAheadReference` (original algorithm) and `lookAheadFast`; `riStartXY`/`numRowsXY`→`listXY`/`lenXY`; candidate index `a`→`candRow`; cache functions `getCacheKey`/`setCacheKey`→`cacheLookup`/`cacheStore` (fields `p1..p3`, `abn` → `id1..id3`, `rowsKey`) |
| Table building | locals `gWork`/`gWork2`/`gWork4`/`cnt`/`tb`/`good` → `row4List`/`row3List`/`partialRow4`/`bucketCount`/`tbl`/`numEntries` |
| Output | `srows`/`ssrows`→`rowsLow`/`rowsHigh`, `sxsAlloc*`/`oldsrows`/`oldnrows`→`outRows*`/`prevOutRows*`/`prevOutRowCount`, `swidth`/`nrows`→`patternWidth`/`numRows` |
| Rule and evolution | `parseRule`/`checkConditions` locals `bs`/`dig`/`neg`/`tempTab`→`ruleOffset`/`digit`/`negated`/`letterHits`; `evolveRow*` locals `row1_s`…/`s`/`t`→`row1Shifted`…/`symShift`/`boundaryShift` |
| Factored tables (new code) | half-table data `lowHalfRows`/`lowHalfOffsets`/`lowHalfMasks` (and `highHalf*`), sizes `lowHalfBits`/`lowResultBits`/`highHalfShift`/`highHalfBits`/`highResultBits`, functions `factoredTableRef`/`factoredHasRow`/`factoredRows`/`buildHalfTable`/`buildFactoredTables`, macros `TABLE_REF`/`TABLE_HAS_ROW`/`TABLE_ROWS` (with `PAIR_*` and `FACTORED_*` variants) |

**Comments added:**
- an overview of how the search works, at the top of `qfind.c`;
- what the phase offsets mean (`backOff`, `fwdOff`, `doubleOff`, `tripleOff`);
- a notation block for the lookahead;
- the table layout;
- the queue representation;
- the saved-extension format;
- what each major global is for.

**Dead code removed:**
- the unused globals `gRows`, `phase` and `hashPhase`;
- the unused functions `getcount` and `unbmalloc`;
- the unused macro `TABLEMASK`;
- the duplicate subperiod-table setup in `main()`.

---

## 8. How correctness was verified

1. **Identical output.** 11 searches covering speeds c/3–c/7, 2c/5, 2c/6, 2c/7; odd, even, asymmetric and gutter symmetry; a wave search; and rule B36/S23. Each was run at `-t 1` and compared with the original on progress lines, ships and final counts. Four builds were compared this way: default, `-DPACKLISTS`, `-DLA_VERIFY` and `-DLA_VERIFY -DPACKLISTS`. The fixed-speed build was compared too.
2. **Lookahead cross-check.** `-DLA_VERIFY` compares the fast and original lookahead on every call. There were no mismatches in any test.
3. **Table contents.** A layout-independent checksum over every logical (row pair, bucket, entry) of every table matches the original at widths 8 and 9.
4. **Sanitizers.** AddressSanitizer and UndefinedBehaviorSanitizer reported nothing on the modified code. The runs covered threads, small queues (to force compaction and saved extensions), the cache on, packed tables, dump, load, split and preview.
5. **Split and hybrid lookahead.** The 13 regression searches were run with `--lookahead-method split`, `hybrid --lookahead-budget 0` and `hybrid --lookahead-budget 64`. Each ran both with every lookahead test cross-checked (`-DLA_VERIFY`) and without it, and each gave output identical to the original. Sanitizer runs were clean at widths 4, 5, 13 and 14, with threads, packed tables and small queues.
6. **Factored tables.** With `-DFACTORED_SORTED`, `--tables factored` reproduced the original's output exactly on 13 searches, in the plain, `-DPACKLISTS`, `-DLA_VERIFY` and cache (`-c 8`) builds. With unsorted factored lists, `-DLA_VERIFY` found no mismatch on any of them. Sanitizer runs were clean at widths 4, 5, 13 and 14 and with threads, small queues and packed tables.
7. **Caveat.** With `-t 2` or more and unlimited ship counts, the order in which ships are printed varies from run to run, even in the original. Compare sets of ships, not order.

---

## 9. Ideas tried and rejected (so they don't need repeating)

| Idea | Result |
|---|---|
| 4-generation lookahead | Sound, and prunes harder (accepted nodes −45%, calls −35%), but 2.3× slower overall: about 2.5 µs per check against about 1.3 µs saved per pruned node [measured] |
| "No next row can pass stage 0" pruning rule | Sound; pruned about 30% of accepted candidates, but time changed only −12%/+5%/−2% (noise). It catches dead ends that are already cheap [measured] |
| Learned candidate ordering | Best case (trained on the same search) only 5–6% faster. Doing it during a search would need versioned orderings in saved extensions to stay exhaustive [measured] |
| Caching repeated states, per-level "Q" sets, decoded-table cache, flattened loops, move-to-front ordering, inlining stage 0, prefetching, branchless code, 64-byte alignment, PGO | no gain or slower [measured] |
| Inner split extended to row11 (`-DINNER_SPLIT2`, with `-DINNER_SPLIT`): ORs the masks over all low halves and all high halves of row11 per overlap value, so cost per (row12,row13) is a sum of half-lists instead of a product | Exact (0 `-DLA_VERIFY` mismatches on 2 searches), but slower: c/5 w8 even 2.23 s → 2.93 s; 3c/7 w10 odd 41 → 36 deepenings in 60 s; c/5 w12 even 4 → 4 [measured]. It loses the early exit at the first working row11. Code left in, off unless `-DINNER_SPLIT2` is given |
| Bitset AND for the pair-mode row22 test (union of row22 lists as a 2^width-bit set, ANDed with table33's mask) | Exact, but 10–15% slower (c/6 w7 even 9.5→11.3 s; c/4 w9 even 11.7→12.7 s): lists are short and the first entry usually hits [measured] |
| Half-relaxed deepening test: before a full depth-first search, check that the low half and the high half can each be extended on their own (overlap bits free), memoized by half-state across the whole search | Motivated by a measurement: in dead depth-first searches, full nodes ÷ (distinct low-half + distinct high-half row sequences) grows with width (c/5 even: 0.65 at w7, 1.1 at w9, 2.3 at w11, 3.5 at w12, 6.7 at w13), so the search repeats much of each half's work for every variant of the other half. The test is sound (it never failed on a live node), but with free overlap bits each half is too loose: it caught 31% of dead nodes at w8, 11% at w10 and 0% at w11–13 [measured]. The dead ends come from the halves interacting, so exploiting the redundancy needs the overlap sequences matched exactly, which amounts to a column-based (LLSSS-style) search. Not built |
| Output-sensitive candidate generation (produce only rows that pass stage 0, column by column) | Candidates per node grow about 1.45× per width (c/5 even: 5.4 at w7, 13 at w10, 50 at w13) while stage-0 survivors stay at 2.4–4.6, so the saving would grow with width; but at w12–13 reading candidates is only 2–4% of the time [measured]. Not built |
| Exact split search (low-half search recording overlap-bit sequences, then high half restricted to them), pruned only by half tables and half stage 0 | Joining is cheap (high half 12–4,700 nodes), but without the full lookahead the low half explodes (6–4,600× the nodes of the real search) and proves only 0–38% of dead nodes dead (0% on large dead searches) [measured, c/5, c/4, 3c/7 at w10–13]. Would need the lookahead inside each half, i.e. a column-decomposed rewrite of `depthFirst()` |
| Decision-diagram (BDD) queue, columns interleaved | On real queues (`-q 20`, about 131k states) the diagram is smaller than the list only for shallow queues (depth 4–5: 1.15, 1.81, 3.19 states per node at w9, w11, w13) and larger for deeper ones (w9 depth 10: 0.42; w11 depth 5: 0.85) [measured]. Not built |
| Merging horizontally shifted duplicate states in asymmetric searches | 0–12% of the queue at w8–9, 0–0.5% at w10 [measured]. Not built |
| Per-level "goal sets" in the per-pair cone lookahead (for each row12, the union of the table masks over all row13/row23, built once per level as a bitset; each candidate then needs only bit tests) | Exact (0 `-DLA_VERIFY` mismatches), but 12–17% slower; built only on a row12's second use, still 5–18% slower. The early exit usually finds a row33 within the first few tests [measured] |
| Prefetching per-pair table masks in the cone lookahead | No measurable change (c/5 w9 odd 118.5 vs 118.7 s, c/4 w9 even 11.2 vs 11.0 s) [measured]. Removed |
| Per-pair tables at width 10 (about 1.4 GB) instead of factored | Slower on the M4 Pro: c/5 w10 even 14 deepenings in 86.5 s vs 16 in 80.2 s factored; c/4 w10 odd 30 in 89.5 s vs 35 in 83.9 s (building the tables costs several seconds up front) [measured]. `auto` unchanged |
| Persistent failed-state memo across deepenings (state = last 2·period rows + phase, kept between passes) | Hit rate only 0.02–2.6% of lookups (c/5 w8 odd 997 / 4.0M; c/6 w7 even 13k / 16M; c/4 w8 even 43k / 1.6M) [measured]; not worth building |
| Hybrid tables (per-pair until a memory budget, then factored) | 11–22% slower even when all tables fit (a branch in the hot loop), and slower than pure factored at widths 10–12 [measured]. Replaced by choosing one representation per search |
| Split-lookahead speed-ups: per-half code specialization (25% fewer instructions, same time), a larger split cache (2^20–2^22 entries: about 12% faster at width 8, not adopted), a small per-level cache of half-cone sets | no gain or within noise [measured] |
| Larger `-q`/`-h`, different `-i` | larger `-q`/`-h` slower; `-i` changes were within noise [measured] |

**Where the remaining time goes [measured].** About 57–90% of the expensive lookahead work is spent in branches that turn out to be dead ends. 39–61% is in dead ends that die within one or two levels.

A perfect, free rule for spotting those would therefore save a lot. But any real check has to look one or two levels ahead, which is the same work the search already does. Only a test several times cheaper than exploring those levels would pay off, and none is known [inferred].

---

## 10. Files

- `qfind-mem2/qfind.c`, `common.h`: the optimized program.
- `qfind-mem2/qfind.c.diff`, `common.h.diff`: unified diffs against the original 2.4b.
- `qfind-mem2/CHANGES.md`: this document (it replaces the earlier `OPTIMIZATIONS.md`).
- `qfind-original-with-readouts/`: the original 2.4b with only the progress readouts added, for side-by-side timing.
