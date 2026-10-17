# qfind optimization project: handoff

Everything needed to continue this work in a new chat or account. Start a new chat, attach this zip, and paste the prompt in section 9.

---

## 1. What this project is

- **Program:** qfind 2.4b, a spaceship search program for Conway's Game of Life and other cellular automata, by Matthias Merzenich (forum name **Sokwe**). It's C with OpenMP: `qfind.c` plus `common.h`.
- **User:** posts on the ConwayLife.com forum as **LuveelVoom**, sharing results in the "qfind - a spaceship search program" thread. Sokwe reads the posts and replies; the thread is in `context/conwaylife_forum_thread.txt`.
- **Goal:** make qfind faster and use less memory **while keeping the search exhaustive**: it must still find every ship within the search limits.
- **Platform:**
  - The user runs a Mac.
  - All measurements here were made on a shared, noisy 2-core Linux VM, single thread unless stated. Expect ±5–8% noise.
  - The user's Mac showed about a 10% gain where the VM showed 25–30%, partly because macOS has no huge-page support.

## 2. Files in this zip

| Path | What |
|---|---|
| `current/qfind.c`, `current/common.h` | **Latest version.** Includes everything in section 4. |
| `current/CHANGES.md` | Full user-facing documentation of every change: options, measurements, bugs fixed, renames, rejected ideas. Kept up to date with each change. |
| `current/*.diff` | Unified diffs against the original 2.4b. |
| `original-2.4b/` | Untouched original qfind 2.4b. |
| `original-with-readouts/` | Original 2.4b plus progress lines showing time and memory; the reference for regression checks and timing comparisons. |
| `tools/regress.sh` | Regression test: builds four variants and compares their output with the original on 12 searches. |
| `experiments/*.c` | Full program copies with measurement hooks for the rejected experiments (section 6). Not part of the program. |
| `context/conwaylife_forum_thread.txt` | Text of the forum thread, including Sokwe's feedback. |

**Building:**
```
gcc -std=c11 -fopenmp -march=native -O3 -DINNER_SPLIT -o qfind current/qfind.c
```
- **Mac:** use Homebrew `gcc-14`; on Apple Silicon use `-mcpu=native` if `-march=native` fails. Apple clang has no OpenMP, so drop `-fopenmp` there; the program then runs single-threaded.
- **The last shipped build used `-DINNER_SPLIT`.**

**Testing:**
```
gcc -std=c11 -fopenmp -O3 -o origr original-with-readouts/qfind.c
tools/regress.sh current ./origr
```
Every line should say `OK`. Use `-t 1` for exact comparisons: with more threads the deepening's early exit makes the output differ from run to run, even in the original.

## 3. How qfind works (short)

- **Rows:** a partial ship is a sequence of rows. Row index a in "phase" ph must satisfy evolve(row a−2P, row a−P, row a) = row a−P+backOff[ph], where P = period.
- **Lookup tables:** these give, for rows (row1, row2) and a required result row, the list of every row3 that works.
- **Breadth-first queue:** `process()` takes a node and adds each candidate row that passes `lookAhead()` and isn't a duplicate (`seenHash`).
- **Deepening:** when the queue fills up, `deepen()` runs `depthFirst()` on every queue node in parallel. Nodes that can't be extended by the deepening amount are removed, and the queue is compacted.
- **Lookahead:** a 3-generation "cone" test that rejects rows which can't be extended. It is the hot spot, around 60–70% of the time.

## 4. What has been done (all in `current/`; details in CHANGES.md)

**Speed:**
1. **Per-level lookahead** (`lookAheadFast`, `LookaheadLevel`): work that doesn't depend on the candidate is cached per search-tree level. Main speed gain: about 25–30% on the VM, about 10% on the user's Mac.
2. **Lookahead cache off by default:** the original turned it on above c/5. Turning it off was 0–36% faster at high speeds.
3. **Split lookahead** (`--lookahead-method split`, automatic from width 12): the cone test is done on half-rows joined on 2 overlap bits, so its cost per test stops growing with width. Slower below width 11, about even at 12, 2.4× faster at 13. A hybrid mode exists too.
4. **Inner split** (`-DINNER_SPLIT`, factored tables with the cone lookahead): the innermost lookahead test ORs bit masks over the half-tables instead of trying every pair of full rows. Measured: 3c/7 w10 odd went from 28 to 41 deepenings in 60 s; c/5 w12 even from 3 to 4; no change at width 8.
5. **Parallel breadth-first lookaheads** (`bfsPrecompute`, only on at 4+ threads): lookaheads for the next 256 queue nodes run in parallel. 2–5% slower on 2 threads, untested on 4+.
6. **`--lookahead-depth 1|2|3`:** a weaker but cheaper lookahead. Still exhaustive, but which ship variants get printed changes. Ranged from 25% faster to 60% slower depending on the search.

**Memory:**
- **Compact table layout:** about 40% less table memory.
- **`-DPACKLISTS`:** bit-packed lists, about 24% less table memory again but about 25% slower.
- **Factored tables** (`--tables factored`, automatic from width 10):
  - Memory growth per width drops from about 7× to about 2.8×. At width 10, 1.48 GB becomes 6 MB; widths 11–14 become feasible.
  - **But** Sokwe measured them at about half the speed of per-pair tables at 3c/7 w10, so on machines with enough RAM, `--tables pair` is faster at widths 10–11.
- **Huge pages** on Linux for the tables.

**Other:**
- **Progress readouts:** each progress line now ends with elapsed time and memory.
- **`--help`:** notes on all of the above.
- **Bugs:** 23 fixed, 20 of them also in the original. The most notable:
  - ARM memory ordering for the shared table pointers;
  - the crash on a 1-CPU machine (default thread count was 0);
  - widths above 14 and `-q` below 4, both of which crashed;
  - #23: "longest partial" and "Maximum depth" were wrong when a search ends with `-> 0/0`. Now fixed: `depthFirst` records partials too.
- **Renaming:** cryptic identifiers were renamed throughout. **Sokwe complained this also renamed his own original names; he only wanted new variables renamed.** Not yet reverted; see section 7.

## 5. Verification approach used

- **Same output:** compare with the original-with-readouts build using `tools/regress.sh` at `-t 1`.
- **`-DLA_VERIFY`:** runs the fast and the reference lookahead on every call and aborts on any disagreement.
- **`-DFACTORED_SORTED` with `--tables factored`:** must reproduce the original's output exactly.
- **Sanitizers:** AddressSanitizer and UndefinedBehaviorSanitizer were used for bug sweeps.
- **Not re-verified:** the most recent rounds (inner split, longest-partial fix, `bfsPrecompute`) had only quick checks, at the user's request. The full regression wasn't rerun after them. **Run `tools/regress.sh` first in the new chat.**

## 6. Tried and rejected (don't repeat)

All measured; the numbers are in CHANGES.md section 9.

**Lookahead and pruning:**
- 4-generation lookahead: 2.3× slower.
- "No child can pass stage 0" pruning rule: no gain.
- 2-level dead-end analysis: only a perfect, free check would pay.
- Pair pruning (`experiments/pair_prune.c`): 0 dead row pairs at widths 6–8.

**Caching and memoization:**
- Failed-state memo within one deepening: 16 hits in 4M lookups.
- Failed-state memo across deepenings (`experiments/memo_across_deepenings.c`): 0.02–2.6% hit rate, and falling with width (c/4: 2.6% at w8, 1.3% at w9).
- Decoded-table cache and other caches: no gain.

**Ordering:**
- Learned candidate ordering: 5–6% best case.

**Table layouts:**
- Hybrid per-pair/factored tables: 11–22% slower.
- Low half kept as full rows (#5 in the chat): memory near 5.7^w; not built.

**Low-level:**
- Bitset AND for the row22 test in pair mode: 10–15% slower.
- Inner split extended to row11 (`-DINNER_SPLIT2`, code left in, off by default): slower.
- SIMD / forced vectorization: no gain.
- Bit-parallel masks, prefetching, branchless code, 64-byte alignment, PGO: no gain.

**Other:**
- Larger `-q`/`-h`: slower.
- **Meet in the middle** (searching backward from the tail and joining), Phase 0 feasibility test in `experiments/phase0_meet_in_middle.c` (env `PHASE0=<depth>`), on c/4 w6 even with no lookahead:
  - Forward frontier stays around 10k–60k states; backward explodes: 1.8M at depth 4, 17M by depth 7 (cap hit).
  - Verdict: **no-go**. The rule barely constrains rows from the tail end.
- **Where the remaining time goes:** 57–90% of the expensive lookahead work is in branches that die, 39–61% in branches that die within 1–2 levels. Catching those would need a check several times cheaper than exploring those levels; none is known.

## 7. Open items and ideas not yet done

**Sokwe's stated priorities (from the forum):**
1. Fix the gcd(period, offset) > 1 bug in duplicate detection without simply requiring matching phase (that requirement makes gcd > 1 searches slow).
2. Keep duplicate nodes linked to the existing copy instead of deleting them, so tagalongs and pushalongs aren't lost. He expects this to need a tree and queue rewrite.
3. Limit renaming to new variables only. Possible approach: rename back the identifiers that existed in 2.4b, using CHANGES.md section 7's table as the map.
4. Optionally, switch rows to 32 bits (`uint32_t`, via a flag like `-DWIDE`) for widths above 14. He says this isn't urgent.
5. He wants changes kept limited so he can follow the differences.

**Ideas suggested but not done:**
- **`--tables auto` choosing by available RAM:** use per-pair tables at widths 10–11 when they fit in memory. Easy, and recovers about 2× at width 10.
- **Running several settings at once and keeping the first ship** (different `-i`, `-q`, `--lookahead-depth`, `--tables` on separate cores or machines): a script only; probably 1.5–3× faster time to the first ship.
- **Distributing `-j` split pieces across machines:** a script only.
- **Packing hot tables together in memory:** maybe 1.1–1.3×.
- **Tail-completion check:** report ships a few rows earlier; no pruning.
- **Testing `bfsPrecompute`** on a machine with 4+ cores.

## 8. User preferences (follow these)

- **Usage:** the user often asks to **use as few tokens / as little usage as possible**. Keep testing minimal unless asked, try at most one to three ideas, and stop and report if something gets hard.
- **Code changes:** make minimal changes, and don't start fixing unrelated bugs unless asked.
- **CHANGES.md:** always document changes there.
- **Confidence:** label how certain claims are (measured, verified or inferred).
- **Style:** be concise; don't agree just to agree; no unfounded "deep" statements. The user's general preference is to start each chat fresh.
- **Exhaustiveness:** must be preserved; say clearly if anything risks it.
- **Memory limit for new speed ideas:** growth of roughly 5.5^w to 6^w per width is acceptable.

## 9. Prompt to paste into the new chat

> I'm continuing a project optimizing qfind 2.4b (Game of Life spaceship search, C/OpenMP). The attached zip has everything. Read HANDOFF.md first, then current/CHANGES.md. The latest source is in current/ (build with `-DINNER_SPLIT`). Before changing anything, run tools/regress.sh against the original-with-readouts build to confirm the current state. Follow the preferences in HANDOFF.md section 8: minimal usage, minimal code changes, document changes in CHANGES.md, label certainty, and preserve exhaustiveness.
