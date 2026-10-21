# Verification

Everything claimed in the README traces to an executed check. Three layers:

1. **C++ unit gate** — `tests/test_quadrature.cpp`, built and run by CTest
   (`cmake --build build && ctest --test-dir build`).
2. **Bounded-output regression** — `tests/test_bench_regression.py`
   (pytest, needs `g++` and Python ≥ 3.9).
3. **Reproducible results** — `examples/run_verification.py` regenerates
   every figure, the GIF and `run_summary.json` from source.

## Layer 1 — unit gate (groups 1–7)

| # | Group | Pinned property |
| --- | --- | --- |
| 1 | Cubic exactness | `∫₀²(x³-2x+1)` = 2 to < 1e-14; evaluations = n+1 |
| 2 | Known integrals | sin on [0,π] → 2, exp on [0,1] → e−1, both < 1e-9 at tol 1e-10, status OK |
| 3 | Efficiency (honest) | narrow Gaussian σ²=1e-6: adaptive 81 evals beats fixed (needs n=2048); smooth Runge case reported without dominance claim |
| 4 | Tolerance hierarchy | error decreases monotonically across the ladder; tightest < 1e-9, loosest stays > 1e-4 |
| 5 | Exhaustion honesty | BUDGET_EXHAUSTED / MAX_DEPTH return partial-but-finite values, never OK |
| 6 | Input validation | null f, reversed bounds, odd n, zero/negative tolerances, insufficient budget, NaN bound → BAD_INPUT |
| 7 | Ladder regression | tolerance ladder finite, strictly decreasing, within [1e-11, 1e-3], exactly 5 rungs |

Group 7 pins the semantics whose violation once produced a runaway
benchmark (below).

## Layer 2 — bounded-output regression (why it exists)

On 2026-10-18 the benchmark driver's tolerance loop multiplied its
tolerance **up** (`1e-3, 1e-1, 1e+1, …`) while its continuation condition
required a lower bound (`tol ≥ 1e-11`). The loop could never exit; the
CSV grew to **25,525,479,948 bytes** (~25.5 GB) and filled the disk. The
ladder was repaired (`tol *= 1e-2`, shrinking 1e-3 → 1e-11), and the
failure mode was turned into a permanent regression:

`test_bench_regression.py` compiles the real driver and asserts, under a
**hard 120 s subprocess timeout**:

- exit 0 and stdout reports `rows=<n>`;
- exactly **51 data rows** (finite expected row count);
- every tolerance parsed from the CSV is finite and ∈ {1e-3 … 1e-11};
  every fixed parameter ∈ {2 … 4096} powers of two (finite valid values);
- every error/evaluation field finite, errors ≥ 0, evaluations > 0, no
  duplicate grid points, full lattice covered;
- file size < **64 KiB** (real output ≈ 2.4 KB);
- a nonexistent output directory → **exit 1** with `cannot open …
  efficiency.csv` on stderr (the driver checks `fopen`/`fprintf`/
  `fflush`/`fclose` and never writes through a null/erroring `FILE*`).

The driver additionally enforces a hard row cap (512) at runtime, so even
a future ladder edit that breaks finiteness aborts within one row instead
of filling a disk.

## Layer 3 — reproducible results

```bash
python examples/run_verification.py
```

compiles the test and benchmark with `g++ -std=c++17 -O2 -Wall -Wextra
-Wpedantic` (300 s timeout per subprocess, 64 KiB CSV cap before
plotting), runs both, then renders `efficiency.png/.svg`,
`refinement_progression.gif` (38 frames from the recorded partition
geometry) and `run_summary.json`, and copies the raw grid to
`docs/results/efficiency.csv`. Measured full runtime ≈ 24 s (Python
3.13.5, GCC 16.2.0, Windows 11; CI runs Ubuntu 24.04 / GCC 13).

## Measured reference numbers

From `docs/results/efficiency.csv` (case `runge_1_25x2`, exact
`atan(5)/5`):

| Setting | Error | Evaluations |
| --- | --- | --- |
| adaptive tol=1e-9 | 2.55e-09 | 129 |
| adaptive tol=1e-11 | 8.93e-11 | 285 |
| fixed n=32 | 4.55e-09 | 33 |
| fixed n=64 | 2.61e-10 | 65 |

Read honestly: on this *smooth* integrand fixed Simpson is competitive
(33 evals vs 129 for comparable error) — the adaptive advantage is
demonstrated on the *localized* narrow-Gaussian case (81 vs 2049
evaluations, [model.md §4](model.md#4-when-adaptive-wins-honestly)).
