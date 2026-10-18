# Model and method

This document states exactly what Quadrature Workbench computes, the
derivations its error estimator relies on, and the validity limits of every
claim. Symbols: `f` the integrand, `[a, b]` the integration interval with
`a < b`, `n` the (even) number of uniform subintervals, `h = (b-a)/n`.

## 1. Composite Simpson rule

The elementary Simpson rule on one panel `[x_{i}, x_{i+2}]` with midpoint
`x_{i+1}` is the exact (degree-≤3) Newton–Cotes interpolant integral:

```
S_i = (h/3) * ( f(x_i) + 4 f(x_{i+1}) + f(x_{i+2}) ),   h = (x_{i+2}-x_i)/2
```

Summing over `n/2` disjoint panels gives the composite rule implemented in
`qw::composite_simpson`:

```
S(n) = (h/3) [ f(a) + 4 Σ_odd f + 2 Σ_even f + f(b) ],  evaluations = n+1
```

For `f ∈ C⁴[a,b]` the leading error term is

```
E(n) = -(b-a) h⁴ f⁗(ξ) / 180,   ξ ∈ (a,b)
```

so the observed convergence order is 4 while `h` is large enough that
roundoff is negligible; below `h⁴·|f⁗| ~ u·|S|` (`u ≈ 1.1e-16`, double
precision) the error floor is dominated by summation roundoff. The test
suite pins cubic exactness (error < 1e-14 on `∫₀² (x³-2x+1)dx = 2`) —
Simpson is exact for cubics because `f⁗ ≡ 0` there — and the measured
order-4 regime in `efficiency.csv`.

## 2. Adaptive Simpson with the |S2 - S1| / 15 estimator

On one interval `[a, b]` let `S1` be the single Simpson estimate and `S2`
the sum of the two half-interval estimates. A classical Richardson
argument for a fourth-order rule gives `E(S2) ≈ (S2 - S1)/15`: if
`I - S1 = C h⁴` then `I - S2 ≈ C (h/2)⁴ = (I - S1)/16`, hence
`S2 - S1 = (I - S1) - (I - S2) ≈ 15 (I - S2)` and

```
|I - S2| ≈ |S2 - S1| / 15
```

`qw::adaptive_simpson` recurses: on each interval it accepts `S2` when

```
|S2 - S1| / 15  ≤  atol + rtol·|S2|
```

(`max(atol, rtol·|S|)`-style local criterion) and otherwise bisects. The
reported `error_estimate` is the **sum of accepted local estimates** — a
heuristic aggregate, **not a guaranteed bound**; the tests measure the
actual attained error `|value - exact|` against closed-form references
instead of trusting the estimator.

**Termination contract.** Three conditions end recursion, each visible in
the returned status:
- all intervals accepted → `OK`;
- evaluation budget `max_evaluations` reached → `BUDGET_EXHAUSTED`
  (partial result, honest status, never silently OK);
- recursion depth `max_depth` reached → `MAX_DEPTH` (partial result).

The per-interval tolerance means a *tighter requested tolerance forces
refinement*, but roundoff sets a floor: on smooth integrands the attainable
error saturates near `1e-13…1e-16`, so tolerances below that are met only
in the estimator's sense. The recorded tolerance ladder stops at `1e-11`
for exactly this reason.

## 3. Benchmark grid and bounded output

The driver sweeps the documented lattice:

- cases: `sin` on `[0,π]` (exact 2), Runge `1/(1+25x²)` on `[0,1]`
  (exact `atan(5)/5`), narrow Gaussian `exp(-x²/0.02)` on `[0,1]`
  (exact `√(0.02π)/2 · erf(1/√0.02)`);
- adaptive: tolerances `1e-3, 1e-5, 1e-7, 1e-9, 1e-11` (shrinking by
  a factor `1e-2`, five rungs);
- fixed: `n = 2 … 4096` (twelve powers of two).

That is `3 × (5 + 12) = 51` data rows — finite by construction. The
driver enforces this with a hard row cap (512) and checked file I/O; the
regression suite pins the exact count, the finiteness/validity of every
tolerance and `n`, a 64 KiB size cap and a 120 s subprocess timeout. See
[verification.md](verification.md) for why this contract exists.

## 4. When adaptive wins, honestly

Adaptive refinement concentrates evaluations where `|S2 - S1|` is large —
i.e. where `f⁗` is large. Two honest regimes are demonstrated:

- **Localized feature (dominance):** a narrow Gaussian
  (`exp(-x²/(2·1e-6))`, width ≪ interval) needs refinement only near the
  peak. Measured: adaptive reaches ~`4e-11` error with **81 evaluations**
  while fixed Simpson needs `n = 2048` (**2049 evaluations**) — ~25× fewer.
- **Smooth Runge-type case (no dominance claimed):** fixed Simpson at
  `n = 64` is already competitive (2.6e-10 vs 3.0e-10 at 201 adaptive
  evaluations). The README reports this as measured, without a general
  superiority claim.

## 5. Assumptions and validity limits

- `f` must be finite at every sampled point; a non-finite sample aborts
  refinement of that interval with a partial status (never a silent NaN).
- Singular integrands (integrable or not) are outside the guaranteed
  scope: the estimator only sees finite samples and cannot detect an
  upstream singularity between them.
- The `|S2-S1|/15` estimate assumes `f ∈ C⁴` locally and can
  under-estimate on functions with cancellation or high-frequency
  oscillation relative to `h`.
- Double precision throughout; no extended precision, no Gauss–Kronrod
  pairing, no oscillatory transforms — those are deliberate scope cuts.

## References

- J. Stoer, R. Bulirsch, *Introduction to Numerical Analysis*, Ch. 3
  (interpolatory quadrature, Richardson extrapolation).
- P. Davis, P. Rabinowitz, *Methods of Numerical Integration*, 2nd ed.
  (adaptive strategies, error estimation practice).
- W. Gander, W. Gautschi, "Adaptive Quadrature — Revisited",
  *BIT Numerical Mathematics* 40 (2000) 84–101 (why naive adaptive
  criteria fail; termination design).
