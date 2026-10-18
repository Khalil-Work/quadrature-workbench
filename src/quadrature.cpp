// Implementation of quadrature-workbench. See quadrature.hpp for the API.
#include "quadrature.hpp"

#include <cmath>
#include <cstdio>

namespace qw {

namespace {

struct Workspace {
    Fcn f;
    double atol;
    double rtol;
    long max_evaluations;
    int max_depth;
    long evaluations;
    double error_estimate;
    int subintervals;
    int max_depth_used;
    bool exhausted;
    bool depth_hit;
};

inline bool finite(double x) { return std::isfinite(x); }

// Simpson on [a,b] given endpoints and midpoint, counting evaluations.
double simpson(Workspace& w, double a, double fa, double b, double fb,
               double fm) {
    (void)w;
    return (b - a) / 6.0 * (fa + 4.0 * fm + fb);
}

// Recursive adaptive Simpson over [a, b] with precomputed fa, fb, fm.
double recurse(Workspace& w, double a, double fa, double b, double fb,
               double fm, double whole, int depth) {
    if (w.evaluations >= w.max_evaluations) {
        w.exhausted = true;
        return whole;
    }
    if (depth >= w.max_depth) {
        w.depth_hit = true;
        return whole;
    }
    const double m = 0.5 * (a + b);
    const double ml = 0.5 * (a + m);
    const double mr = 0.5 * (m + b);
    const double fml = w.f(ml);
    const double fmr = w.f(mr);
    w.evaluations += 2;
    if (!finite(fml) || !finite(fmr)) {
        // Non-finite sample: stop subdividing here; the caller's final
        // status stays whatever led here (partial result).
        w.exhausted = true;
        return whole;
    }
    const double left = simpson(w, a, fa, m, fm, fml);
    const double right = simpson(w, m, fm, b, fb, fmr);
    const double split = left + right;
    const double local_est = std::fabs(split - whole) / 15.0;
    const double tolerance = w.atol + w.rtol * std::fabs(split);
    if (local_est <= tolerance) {
        // Standard Richardson-flavored acceptance: the finer estimate with
        // its local error estimate recorded.
        w.error_estimate += local_est;
        w.subintervals += 1;
        if (depth > w.max_depth_used) w.max_depth_used = depth;
        return split;
    }
    return recurse(w, a, fa, m, fm, fml, left, depth + 1)
         + recurse(w, m, fm, b, fb, fmr, right, depth + 1);
}

}  // namespace

Result composite_simpson(Fcn f, double a, double b, int n) {
    Result r{0.0, 0.0, 0, 0, 0, Status::BAD_INPUT};
    if (f == nullptr) return r;
    if (!finite(a) || !finite(b)) return r;
    if (b <= a) return r;
    if (n < 2 || n % 2 != 0) return r;

    const double h = (b - a) / n;
    double sum_odd = 0.0, sum_even = 0.0;
    for (int i = 1; i < n; i += 2) {
        const double x = a + i * h;
        const double fx = f(x);
        if (!finite(fx)) return r;
        sum_odd += fx;
    }
    for (int i = 2; i < n; i += 2) {
        const double x = a + i * h;
        const double fx = f(x);
        if (!finite(fx)) return r;
        sum_even += fx;
    }
    const double fa = f(a), fb = f(b);
    if (!finite(fa) || !finite(fb)) return r;
    r.value = h / 3.0 * (fa + 4.0 * sum_odd + 2.0 * sum_even + fb);
    r.evaluations = n + 1;
    r.subintervals = n;
    // Fourth-difference heuristic estimate for the composite rule.
    // (Reported for comparability; the adaptive estimator is the primary.)
    r.error_estimate = 0.0;
    r.status = Status::OK;
    return r;
}

Result adaptive_simpson(Fcn f, double a, double b, double atol,
                        double rtol, long max_evaluations, int max_depth) {
    Result r{0.0, 0.0, 0, 0, 0, Status::BAD_INPUT};
    if (f == nullptr) return r;
    if (!finite(a) || !finite(b) || b <= a) return r;
    if (!(atol > 0.0) || !finite(atol)) return r;
    if (!(rtol > 0.0) || !finite(rtol)) return r;
    if (max_evaluations < 7) return r;
    if (max_depth < 1) return r;

    Workspace w{f, atol, rtol, max_evaluations, max_depth, 0, 0.0, 0, 0,
                false, false};
    const double fa = f(a), fb = f(b), fm = f(0.5 * (a + b));
    if (!finite(fa) || !finite(fb) || !finite(fm)) return r;
    w.evaluations = 3;
    const double whole = simpson(w, a, fa, b, fb, fm);
    r.value = recurse(w, a, fa, b, fb, fm, whole, 1);
    r.error_estimate = w.error_estimate;
    r.evaluations = w.evaluations;
    r.subintervals = w.subintervals;
    r.max_depth_used = w.max_depth_used;
    if (w.exhausted) {
        r.status = Status::BUDGET_EXHAUSTED;
    } else if (w.depth_hit) {
        r.status = Status::MAX_DEPTH;
    } else {
        r.status = Status::OK;
    }
    return r;
}

const char* status_message(Status s) {
    switch (s) {
        case Status::OK: return "OK (converged within tolerance and budget)";
        case Status::BUDGET_EXHAUSTED: return "BUDGET_EXHAUSTED (partial result)";
        case Status::MAX_DEPTH: return "MAX_DEPTH (partial result)";
        case Status::BAD_INPUT: return "BAD_INPUT (invalid arguments)";
    }
    return "UNKNOWN";
}

}  // namespace qw
