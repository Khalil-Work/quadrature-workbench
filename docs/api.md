# API reference

Header: [`src/quadrature.hpp`](../src/quadrature.hpp) — include and link
`qw` (CMake target) or compile `quadrature.cpp` alongside your file.
C++17, no external dependencies.

```cpp
#include "quadrature.hpp"

using qw::composite_simpson, qw::adaptive_simpson, qw::Status;
```

## Types

```cpp
using Fcn = double (*)(double);   // plain function pointer, no capture

enum class Status {
    OK,               // converged within tolerance and budget
    BUDGET_EXHAUSTED, // hit max evaluations; value is partial
    MAX_DEPTH,        // hit recursion depth; value is partial
    BAD_INPUT         // invalid arguments; value/evaluations are zeroed
};

struct Result {
    double value;           // computed integral
    double error_estimate;  // sum of accepted local |S2-S1|/15 (heuristic)
    long   evaluations;     // function evaluations actually used
    int    subintervals;    // accepted subintervals at exit
    int    max_depth_used;  // deepest recursion level reached
    Status status;
};

const char* status_message(Status s);  // human-readable status line
```

## `composite_simpson`

```cpp
Result composite_simpson(Fcn f, double a, double b, int n);
```

Uniform-grid composite Simpson, `n` **even and ≥ 2**, `a < b`, `f` finite
at all `n+1` nodes. Returns `evaluations = n+1`, `error_estimate = 0.0`
(the adaptive estimator is the primary one). `BAD_INPUT` for null `f`,
non-finite/reversed bounds, odd or too-small `n`, or any non-finite sample
(the partially accumulated value is discarded, never returned as NaN).

## `adaptive_simpson`

```cpp
Result adaptive_simpson(Fcn f, double a, double b,
                        double atol, double rtol,
                        long max_evaluations, int max_depth);
```

Recursive bisection with the `|S2-S1|/15` local estimator and the
per-interval criterion `|S2-S1|/15 ≤ atol + rtol·|S2|` (see
[model.md](model.md#2-adaptive-simpson-with-the-s2---s1--15-estimator)).

Parameter contract (violations return `BAD_INPUT` with zeroed fields):

| Parameter           | Requirement                     | Notes |
| --- | --- | --- |
| `f`                 | non-null, finite on sampled points | non-finite samples end refinement honestly |
| `a`, `b`            | finite, `a < b`                 | reversed/degenerate bounds rejected |
| `atol`, `rtol`      | finite, `> 0`                   | both criteria combined additively |
| `max_evaluations`   | `≥ 7`                           | 3 initial + ≥2 per bisection step |
| `max_depth`         | `≥ 1`                           | bounds recursion; 50–60 is ample for double precision |

`status` is `OK` only when every interval was accepted inside both limits;
`BUDGET_EXHAUSTED` / `MAX_DEPTH` return the partial value with the
evaluation/interval counters reached — check `status` before trusting
`value` at tight tolerances.

## Benchmark driver `qw_bench`

`src/bench.cpp` builds `qw_bench` (CMake target `qw_bench`, or compile
`bench.cpp + quadrature.cpp`). Usage: `qw_bench [output_dir]` (default
`build`). Writes `output_dir/efficiency.csv` with the fixed 51-row lattice
(documented in [model.md §3](model.md#3-benchmark-grid-and-bounded-output))
and prints `rows=<n>` on stdout. Failure contract:

| Failure | Behaviour |
| --- | --- |
| output file cannot be opened | message with `errno`/`strerror` on stderr, exit **1** |
| non-finite computed value | message naming the case/parameter, exit **2** |
| row cap (512) exceeded | message naming the offending grid point, exit **3** |
| `fprintf` write failure | message, exit **4** (5/6 for `fflush`/`fclose`) |

## Worked example

```cpp
#include <cmath>
#include <cstdio>
#include "quadrature.hpp"

static double runge(double x) { return 1.0 / (1.0 + 25.0 * x * x); }

int main() {
    const auto r = qw::adaptive_simpson(runge, 0.0, 1.0, 1e-10, 1e-10,
                                        100000, 50);
    std::printf("value=%.12f (exact %.12f)\nstatus=%s\nevals=%ld\n",
                r.value, std::atan(5.0) / 5.0,
                qw::status_message(r.status), r.evaluations);
}
```

Output (measured, GCC 16.2.0 `-O2`, 2026-10-19):

```
value=0.274680153689 (exact 0.274680153389)
status=OK (converged within tolerance and budget)
evals=201
```
