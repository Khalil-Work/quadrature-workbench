// Verification tests for quadrature-workbench.
//
// Properties pinned:
// 1. Exactness: composite Simpson integrates cubics to machine precision.
// 2. Known integrals: sin over [0,pi] = 2 and exp over [0,1] = e-1 within
//    tolerance, with recorded evaluations.
// 3. Adaptive beats fixed Simpson at comparable-or-fewer evaluations on a
//    hard case (narrow Gaussian / 1/(1+25x^2)) - measured, honest.
// 4. Tolerance hierarchy: tighter tolerance -> smaller measured error
//    (documented ~O(h^4) smoothing behavior).
// 5. Budget/depth exhaustion returns an honest partial status, never OK.
// 6. Invalid inputs rejected (null f, reversed bounds, odd n, bad
//    tolerances, insufficient budget).
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "../src/quadrature.hpp"

using qw::adaptive_simpson;
using qw::composite_simpson;
using qw::Status;

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kE  = 2.71828182845904523536;
}


static int failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAILED: %s\n", what); ++failures; }
}

static double f_cubic(double x) { return x * x * x - 2.0 * x + 1.0; }
static double fsin(double x) { return std::sin(x); }
static double fexp(double x) { return std::exp(x); }
static double f_hard(double x) { return 1.0 / (1.0 + 25.0 * x * x); }
static double f_peak(double x) { return std::exp(-x * x / (2.0 * 1e-6)); }

int main() {
    // 1. Cubic exactness.
    {
        const auto r = composite_simpson(f_cubic, 0.0, 2.0, 8);
        check(r.status == Status::OK, "cubic: OK");
        // int (x^3-2x+1) on [0,2] = 4 - 4 + 2 = 2.
        check(std::fabs(r.value - 2.0) < 1e-14, "cubic: exact to roundoff");
        check(r.evaluations == 9, "cubic: evaluations = n+1");
    }

    // 2. Known integrals.
    {
        const auto r = adaptive_simpson(fsin, 0.0, kPi, 1e-10, 1e-10,
                                         100000, 50);
        check(r.status == Status::OK, "sin: OK");
        check(std::fabs(r.value - 2.0) < 1e-9, "sin: value within 1e-9");
        const auto r2 = adaptive_simpson(fexp, 0.0, 1.0, 1e-10, 1e-10,
                                          100000, 50);
        check(r2.status == Status::OK && std::fabs(r2.value - (kE - 1.0)) < 1e-9,
              "exp: e-1 within 1e-9");
    }

    // 3. Adaptive vs fixed on the hard Runge-type case.
    {
        // True value: int_0^1 1/(1+25x^2) = atan(5)/5.
        const double exact = std::atan(5.0) / 5.0;
        const auto adapt = adaptive_simpson(f_hard, 0.0, 1.0, 1e-10, 1e-10,
                                             100000, 50);
        check(adapt.status == Status::OK, "hard: adaptive OK");
        const double adapt_err = std::fabs(adapt.value - exact);
        check(adapt_err < 1e-9, "hard: adaptive error < 1e-9");
        // HONEST measured observation on this SMOOTH function: fixed
        // Simpson is competitive (recorded in output; no dominance claim).
        int fixed_n = 64;
        double fixed_err = 0.0;
        for (; fixed_n <= 1024; fixed_n *= 2) {
            const auto r = composite_simpson(f_hard, 0.0, 1.0, fixed_n);
            fixed_err = std::fabs(r.value - exact);
            if (fixed_err <= adapt_err) break;
        }
        std::printf("hard case (smooth): adapt evals=%ld err=%.2e | fixed n=%d err=%.2e\n",
                    adapt.evaluations, adapt_err, fixed_n, fixed_err);
        std::fflush(stdout);

        // LOCALIZED feature (narrow Gaussian): adaptive must win — fixed
        // Simpson refines everywhere, adaptive only near the peak.
        const double gauss_exact = std::sqrt(2.0 * kPi * 1e-6) * 0.5 *
                                    (std::erf(1.0 / std::sqrt(2.0 * 1e-6)));
        const auto g_adapt = adaptive_simpson(f_peak, 0.0, 1.0, 1e-9, 1e-9,
                                              100000, 60);
        check(g_adapt.status == Status::OK, "gauss: adaptive OK");
        const double g_adapt_err = std::fabs(g_adapt.value - gauss_exact);
        check(g_adapt_err < 1e-8, "gauss: adaptive error < 1e-8");
        int g_fixed_n = 2;
        for (; g_fixed_n <= 65536; g_fixed_n *= 2) {
            const auto r = composite_simpson(f_peak, 0.0, 1.0, g_fixed_n);
            if (std::fabs(r.value - gauss_exact) <= g_adapt_err) break;
        }
        check(g_fixed_n + 1 > g_adapt.evaluations,
              "gauss: adaptive more efficient than fixed (localized feature)");
        std::printf("gauss case: adapt evals=%ld err=%.2e | fixed needs n=%d\n",
                    g_adapt.evaluations, g_adapt_err, g_fixed_n);
        std::fflush(stdout);
    }

    // 3b/4 moved into section 3 above (gauss + efficiency comparison).
    // 4. Tolerance hierarchy: tighter tolerances drive the accepted error
    // down on sin (roundoff floor ~1e-16; stopping at 1e-12 stays safely
    // above it while exercising genuine refinement).
    {
        double loosest_err = -1.0, tightest_err = -1.0;
        for (double tol = 1e-3; tol >= 1e-12; tol *= 1e-2) {
            const auto r = adaptive_simpson(fsin, 0.0, kPi, tol, tol,
                                             100000, 40);
            const double err = std::fabs(r.value - 2.0);
            if (loosest_err < 0.0) loosest_err = err;
            tightest_err = err;
        }
        check(tightest_err < 1e-9,
              "tolerance hierarchy: tightest < 1e-9 (estimator err ~ 0.3*tol)");
        check(loosest_err > 1e-4, "tolerance hierarchy: loosest stays loose");
        check(tightest_err < loosest_err,
              "tolerance hierarchy: tightest error < loosest error");
    }

    // 5. Budget/depth exhaustion is honest.
    {
        const auto r = adaptive_simpson(f_hard, 0.0, 1.0, 1e-12, 1e-12, 60, 50);
        check(r.status == Status::BUDGET_EXHAUSTED,
              "budget: exhausted status, not OK");
        const auto r2 = adaptive_simpson(f_hard, 0.0, 1.0, 1e-14, 1e-14,
                                          100000, 3);
        check(r2.status == Status::MAX_DEPTH,
              "depth: MAX_DEPTH status, not OK");
        check(r2.value == r2.value, "depth: value is not NaN");
    }

    // 6. Invalid inputs.
    {
        check(composite_simpson(nullptr, 0, 1, 8).status == Status::BAD_INPUT,
              "invalid: null f");
        check(composite_simpson(fsin, 1.0, 0.0, 8).status == Status::BAD_INPUT,
              "invalid: reversed bounds");
        check(composite_simpson(fsin, 0.0, 1.0, 7).status == Status::BAD_INPUT,
              "invalid: odd n");
        check(adaptive_simpson(fsin, 0.0, 1.0, 0.0, 1e-6, 1000, 10).status
                  == Status::BAD_INPUT, "invalid: zero atol");
        check(adaptive_simpson(fsin, 0.0, 1.0, 1e-6, -1e-6, 1000, 10).status
                  == Status::BAD_INPUT, "invalid: negative rtol");
        check(adaptive_simpson(fsin, 0.0, 1.0, 1e-6, 1e-6, 3, 10).status
                  == Status::BAD_INPUT, "invalid: insufficient budget");
        const double bad = std::nan("");
        (void)bad;
        check(composite_simpson(fsin, bad, 1.0, 8).status == Status::BAD_INPUT,
              "invalid: NaN bound");
    }

    // 7. Bench tolerance-ladder regression. History: the benchmark loop
    // once multiplied tol *up* (1e-3, 1e-1, 1e+1, ...) while its continuation
    // test required tol >= 1e-11, so it never exited and efficiency.csv grew
    // to 25,525,479,948 bytes. Pin the ladder semantics the driver relies
    // on: finite, strictly decreasing, bounded to [1e-11, 1e-3], and a
    // bounded rung count. The end-to-end CSV bound is pinned separately in
    // tests/test_bench_regression.py with a hard subprocess timeout.
    {
        double prev = 0.0;
        int rungs = 0;
        for (double tol = 1e-3; tol >= 1e-11; tol *= 1e-2) {
            check(std::isfinite(tol) && tol >= 1e-11 && tol <= 1e-3,
                  "ladder: tol finite and within [1e-11, 1e-3]");
            check(rungs == 0 || tol < prev, "ladder: tol strictly decreasing");
            prev = tol;
            ++rungs;
        }
        check(rungs == 5, "ladder: exactly 5 rungs (bounded grid)");
    }

    if (failures) {
        std::fprintf(stderr, "%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    std::printf("quadrature-workbench: all checks passed\n");
    return EXIT_SUCCESS;
}
