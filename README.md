# qfind optimized

Speed and memory optimizations of **qfind 2.4b**, the spaceship search program for Conway's Game of Life
and related cellular automata by Matthias Merzenich (Sokwe). The search stays exhaustive: it still finds
every ship within the search limits. Discussed in the ConwayLife.com thread "qfind - a spaceship search program".

| Path | What it is |
|---|---|
| `current/` | Latest version (`qfind.c`, `common.h`) and diffs against 2.4b |
| `current/CHANGES.md` | Every change: options, measurements, bugs fixed, rejected ideas |
| `original-2.4b/` | Untouched original qfind 2.4b |
| `original-with-readouts/` | Original plus time and memory progress lines; the regression reference |
| `tools/regress.sh` | Regression test against the original on 12 searches |
| `experiments/` | Program copies with measurement hooks for rejected experiments |
| `HANDOFF.md` | Project state, how qfind works, open items |

## Build

```
gcc -std=c11 -fopenmp -march=native -O3 -DINNER_SPLIT -o qfind current/qfind.c
```

On a Mac use Homebrew `gcc-14` (Apple clang has no OpenMP; without `-fopenmp` it runs single-threaded).
On Apple Silicon use `-mcpu=native` if `-march=native` fails.

## Test

```
gcc -std=c11 -fopenmp -O3 -o origr original-with-readouts/qfind.c
tools/regress.sh current ./origr      # every line should say OK
```
