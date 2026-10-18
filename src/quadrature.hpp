// quadrature-workbench: adaptive and composite quadrature with explicit
// error control and evaluation budgets.
//
// Composite Simpson on a uniform grid (n even): O(h^4) for smooth
// integrands. Adaptive Simpson: recursive interval bisection using the
// classic local error estimate |S2 - S1| / 15 with the standard split
// into two half-interval Simpson rules. The tolerance contract is
// explicit: the reported error_estimate is the sum of the accepted local
// |S2-S1|/15 estimates (a heuristic estimate, NOT a bound); termination
// is by absolute OR relative criterion max(atol, rtol*|S|) on each
// subinterval, subject to a recursion budget and a maximum depth.
//
// All termination paths carry a reason; a budget exhaustion returns
// CONVERGED=false with partial diagnostics rather than silently wrong
// answers. Every entry point validates its inputs.
#ifndef QUADRATURE_WORKBENCH_HPP
#define QUADRATURE_WORKBENCH_HPP

namespace qw {

using Fcn = double (*)(double);

enum class Status {
    OK = 0,             // converged within tolerance and budget
    BUDGET_EXHAUSTED,   // hit max evaluations; result is partial
    MAX_DEPTH,          // hit recursion depth; result is partial
    BAD_INPUT           // invalid arguments (null f, bad bounds/tolerances)
};

struct Result {
    double value;            // computed integral
    double error_estimate;   // sum of accepted local estimates (heuristic)
    long evaluations;        // function evaluations actually used
    int subintervals;        // accepted subintervals at exit
    int max_depth_used;
    Status status;
};

// Composite Simpson on a uniform grid of n intervals (n must be even).
// Returns evaluations = n+1. BAD_INPUT on invalid arguments.
Result composite_simpson(Fcn f, double a, double b, int n);

// Adaptive Simpson with the |S2-S1|/15 local estimator.
// atol, rtol > 0; max_evaluations > 0; max_depth > 0.
Result adaptive_simpson(Fcn f, double a, double b, double atol,
                        double rtol, long max_evaluations, int max_depth);

const char* status_message(Status s);

}  // namespace qw

#endif  // QUADRATURE_WORKBENCH_HPP
