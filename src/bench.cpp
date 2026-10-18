// Benchmark driver: runs both methods on the documented case set across a
// tolerance/n grid and writes efficiency.csv for plotting.
//
// Output is bounded by construction: a fixed case list, a finite tolerance
// ladder (1e-3 shrinking by a factor 1e-2 down to 1e-11) and a finite n
// ladder (2..4096). Two guards make that invariant explicit:
//   * a hard row cap aborts the run if a future edit ever breaks the
//     ladders (this is the regression against the 2026-10 runaway that
//     produced a 25.5 GB CSV when tolerance *grew* while the loop required
//     a lower bound);
//   * every file operation is checked; a failed open/write closes cleanly
//     and exits nonzero instead of writing through a null/erroring stream.
// Any non-finite computed value aborts with the offending case name rather
// than polluting the CSV with nan/inf fields.
#define _CRT_SECURE_NO_WARNINGS 1
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "quadrature.hpp"

namespace {
constexpr double kPi = 3.14159265358979323846;
// Hard cap on data rows: the documented grid is 3 cases x (5 tolerances +
// 12 n values) = 51 rows; the cap leaves headroom for grid edits but stops
// any unbounded loop within one row.
constexpr long kMaxRows = 512;
}  // namespace

static double fsin(double x) { return std::sin(x); }
static double f_hard(double x) { return 1.0 / (1.0 + 25.0 * x * x); }
static double f_gauss(double x) {
    return std::exp(-x * x / (2.0 * 0.01));
}

// Returns 0 on success; nonzero after printing a diagnostic to stderr.
static int write_row(std::FILE* eff, long& rows, const char* case_name,
                     const char* method, const char* parameter, double error,
                     long evaluations, int status) {
    if (!std::isfinite(error) || evaluations < 0) {
        std::fprintf(stderr,
                     "qw_bench: non-finite result for %s/%s/%s (error=%g, "
                     "evaluations=%ld); aborting to keep output valid\n",
                     case_name, method, parameter, error, evaluations);
        return 2;
    }
    if (++rows > kMaxRows) {
        std::fprintf(stderr,
                     "qw_bench: row cap %ld exceeded at %s/%s/%s; the case/"
                     "tolerance/n ladders are no longer bounded — aborting\n",
                     kMaxRows, case_name, method, parameter);
        return 3;
    }
    if (std::fprintf(eff, "%s,%s,%s,%.9e,%ld,%d\n", case_name, method,
                     parameter, error, evaluations, status) < 0) {
        std::fprintf(stderr, "qw_bench: fprintf failed at row %ld\n", rows);
        return 4;
    }
    return 0;
}

int main(int argc, char** argv) {
    const char* out = argc > 1 ? argv[1] : "build";

    struct Case {
        const char* name;
        double (*f)(double);
        double a, b, exact;
    } cases[] = {
        {"sin_0_pi", fsin, 0.0, kPi, 2.0},
        {"runge_1_25x2", f_hard, 0.0, 1.0, std::atan(5.0) / 5.0},
        {"gauss_s0.1", f_gauss, 0.0, 1.0,
         std::sqrt(2.0 * kPi * 0.01) * 0.5 *
             (std::erf(1.0 / std::sqrt(2.0 * 0.01)))},
    };

    char path[512];
    std::snprintf(path, sizeof(path), "%s/efficiency.csv", out);
    std::FILE* eff = std::fopen(path, "w");
    if (eff == nullptr) {
        std::fprintf(stderr, "qw_bench: cannot open '%s' for writing "
                             "(errno=%d, %s)\n",
                     path, errno, std::strerror(errno));
        return 1;
    }
    long rows = 0;
    if (std::fprintf(eff, "case,method,parameter,error,evaluations,status\n") < 0) {
        std::fprintf(stderr, "qw_bench: failed to write CSV header to '%s'\n",
                     path);
        std::fclose(eff);
        return 4;
    }

    for (const auto& c : cases) {
        // Adaptive across tolerances: 1e-3 *shrinking* by 1e-2 to 1e-11.
        // (History: this ladder once multiplied tolerance *up* while the
        // continuation test required tol >= 1e-11, so it never exited and
        // wrote a 25.5 GB CSV. It is pinned by test_bench_regression.py.)
        for (double tol = 1e-3; tol >= 1e-11; tol *= 1e-2) {
            const auto r = qw::adaptive_simpson(c.f, c.a, c.b, tol, tol,
                                                5000000, 60);
            const double err = std::fabs(r.value - c.exact);
            char param[32];
            std::snprintf(param, sizeof(param), "tol=%.0e", tol);
            const int rc = write_row(eff, rows, c.name, "adaptive", param,
                                     err, r.evaluations, (int)r.status);
            if (rc != 0) { std::fclose(eff); return rc; }
        }
        // Fixed Simpson across n (2..4096, even).
        for (int n = 2; n <= 4096; n *= 2) {
            const auto r = qw::composite_simpson(c.f, c.a, c.b, n);
            const double err = std::fabs(r.value - c.exact);
            char param[32];
            std::snprintf(param, sizeof(param), "n=%d", n);
            const int rc = write_row(eff, rows, c.name, "fixed", param,
                                     err, r.evaluations, (int)r.status);
            if (rc != 0) { std::fclose(eff); return rc; }
        }
    }

    if (std::fflush(eff) != 0) {
        std::fprintf(stderr, "qw_bench: fflush failed on '%s' (errno=%d)\n",
                     path, errno);
        std::fclose(eff);
        return 5;
    }
    if (std::fclose(eff) != 0) {
        std::fprintf(stderr, "qw_bench: fclose failed on '%s' (errno=%d)\n",
                     path, errno);
        return 6;
    }
    std::printf("wrote %s/efficiency.csv rows=%ld\n", out, rows);
    return 0;
}
