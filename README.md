# qfind optimized

A faster, leaner version of **qfind 2.4b**, the spaceship search program for Conway's Game of Life and related
cellular automata by Matthias Merzenich (Sokwe). It takes the same options as qfind and stays exhaustive: it
still finds every ship within the search limits. Discussed in the ConwayLife.com thread "qfind - a spaceship
search program".

## What's different

- **Faster:** about 25–30% at speeds up to c/5 on Linux (about 10% on Apple silicon), and up to
  36% more above c/5.
  From width 12 a new "split" lookahead takes over whose cost doesn't grow with width (2.4× faster at width 13).
- **Much less memory at large widths:** from width 10, lookup tables are stored factored. Width 10 drops from
  1.48 GB to 6 MB, and widths 11–12 now run in tens of MB where the original ran out of memory.
- **Progress readouts:** every depth line shows elapsed time and memory; a total is printed at the end.
- **22 bugs fixed**, 19 of which are also in the original (for example, `-w` above 14 or `-q`
  below 4 used to crash).
- **More readable source:** renamed variables, an overview and comments.

The defaults pick the fastest methods automatically. New options, for tuning:

| Option | What it does |
|---|---|
| `--lookahead-depth 1\|2\|3` | How much of the lookahead test to apply (default 3, the full test). 1–2 can be faster on long c/4–c/5 searches; time it. |
| `--lookahead-method auto\|cone\|split\|hybrid` | How the lookahead is computed. auto uses split from width 12. All give identical answers. |
| `--split-cache-bits N` | Split lookahead cache size, 2^N entries per half (default 22 = 1152 MB, shared by all threads) |
| `--tables auto\|pair\|factored` | Lookup table layout. auto uses factored tables from width 10. |
| `-c N` | The lookahead cache is now **off by default**; `-c N` turns it back on |

`./qfind --help` ends with performance notes for this build. Everything changed, with measurements:
[current/CHANGES.md](current/CHANGES.md).

## Build

```
gcc -std=c11 -fopenmp -march=native -O3 -DINNER_SPLIT -o qfind current/qfind.c
```

On a Mac, Apple clang with Homebrew's `libomp` is the fastest build (4–16% faster than gcc, which emulates
thread-local storage on macOS):

```
brew install libomp
clang -std=c11 -O3 -mcpu=native -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include \
      -L$(brew --prefix libomp)/lib -lomp -DINNER_SPLIT -o qfind current/qfind.c
```

Then run it exactly like qfind, e.g. `./qfind -r B3/S23 -v c/4 -w 8 -s even`.

## Files

| Path | What it is |
|---|---|
| `current/` | The optimized `qfind.c` and `common.h`, plus diffs against 2.4b |
| `original-2.4b/` | Untouched original qfind 2.4b |
| `original-with-readouts/` | The original with the same time/memory readouts, for side-by-side comparison |
| `tools/regress.sh` | Checks the optimized build against the original on 12 searches |

To run the regression check:

```
gcc -std=c11 -fopenmp -O3 -o origr original-with-readouts/qfind.c
tools/regress.sh current ./origr      # every line should say OK
```

Developer notes (how qfind works internally, open items, rejected experiments in `experiments/`) are in
`HANDOFF.md`.
