"""End-to-end bounded-output regression for the Quadrature Workbench driver.

History (2026-10): the benchmark's tolerance ladder once multiplied tol
*up* while its continuation test required tol >= 1e-11, so the loop never
exited and efficiency.csv grew to 25,525,479,948 bytes (~25.5 GB), filling
the disk. These checks pin the exact bounded contract of ``qw_bench``:

* it terminates well inside a hard subprocess timeout;
* the CSV has the exact expected row count (finite by construction);
* every tolerance/parameter value is finite and drawn from the documented
  ladders (tol in {1e-3 .. 1e-11, factor 1e-2}; n in {2 .. 4096, powers
  of two);
* every numeric field is finite and the file stays far below a hard size
  cap;
* a failed output-file open exits nonzero with a useful message instead of
  dereferencing a null FILE*.

Run:  python -m pytest tests/test_bench_regression.py -q
(requires ``g++`` on PATH; C++17).
"""

from __future__ import annotations

import pathlib
import re
import shutil
import subprocess
import tempfile

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]

EXPECTED_CASES = {"sin_0_pi", "runge_1_25x2", "gauss_s0.1"}
EXPECTED_TOLS = {10.0**-e for e in range(3, 12, 2)}      # 1e-3 .. 1e-11
EXPECTED_NS = {2**k for k in range(1, 13)}               # 2 .. 4096
EXPECTED_ROWS = len(EXPECTED_CASES) * (len(EXPECTED_TOLS) + len(EXPECTED_NS))
SIZE_CAP_BYTES = 64 * 1024   # real output is ~4 KB; cap leaves huge margin
HARD_TIMEOUT_S = 120         # real runtime is well under 1 s


def _gxx() -> str:
    gxx = shutil.which("g++") or shutil.which("g++-13") or shutil.which("g++-12")
    if gxx is None:
        pytest.skip("g++ not available on PATH")
    return gxx


@pytest.fixture(scope="module")
def bench_exe(tmp_path_factory: pathlib.Path) -> pathlib.Path:
    out = tmp_path_factory.mktemp("qw-build") / "qw_bench"
    suffix = ".exe" if subprocess.run(["python", "-c",
                                       "import sys;print(sys.platform=='win32')"],
                                      capture_output=True, text=True,
                                      check=True).stdout.strip() == "True" else ""
    out = out.with_name(out.name + suffix)
    subprocess.run([_gxx(), "-std=c++17", "-O2", "-Wall", "-Wextra",
                    "-Wpedantic", str(ROOT / "src" / "bench.cpp"),
                    str(ROOT / "src" / "quadrature.cpp"), "-o", str(out)],
                   check=True)
    return out


def test_bounded_csv_under_hard_timeout(bench_exe: pathlib.Path,
                                         tmp_path: pathlib.Path) -> None:
    outdir = tmp_path / "out"
    outdir.mkdir()
    proc = subprocess.run([str(bench_exe), str(outdir)], capture_output=True,
                          text=True, timeout=HARD_TIMEOUT_S)
    assert proc.returncode == 0, proc.stderr

    csv = outdir / "efficiency.csv"
    assert csv.exists(), "driver reported success but wrote no CSV"
    size = csv.stat().st_size
    assert size < SIZE_CAP_BYTES, f"CSV size {size} exceeds cap {SIZE_CAP_BYTES}"

    lines = csv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "case,method,parameter,error,evaluations,status"
    rows = lines[1:]
    assert len(rows) == EXPECTED_ROWS, (
        f"expected exactly {EXPECTED_ROWS} data rows, got {len(rows)}")

    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        case, method, parameter, error, evaluations, status = row.split(",")
        assert case in EXPECTED_CASES
        err = float(error)
        evals = int(evaluations)
        assert err == err and err not in (float("inf"),) and err >= 0.0, (
            f"non-finite or negative error in row: {row}")
        assert evals > 0, f"non-positive evaluation count in row: {row}"
        int(status)  # parseable integer status code
        if method == "adaptive":
            match = re.fullmatch(r"tol=([0-9.e+-]+)", parameter)
            assert match is not None, f"malformed adaptive parameter: {parameter}"
            tol = float(match.group(1))
            assert tol == tol and abs(tol) != float("inf")
            assert tol in EXPECTED_TOLS, f"tolerance outside ladder: {tol}"
        elif method == "fixed":
            match = re.fullmatch(r"n=(\d+)", parameter)
            assert match is not None, f"malformed fixed parameter: {parameter}"
            assert int(match.group(1)) in EXPECTED_NS
        else:
            pytest.fail(f"unknown method: {method}")
        key = (case, method, parameter)
        assert key not in seen, f"duplicate grid point: {key}"
        seen.add(key)

    assert len(seen) == EXPECTED_ROWS, "grid is not the full documented lattice"
    # Driver stdout reports the bounded row count for the audit trail.
    assert re.search(r"\brows=\d+\b", proc.stdout)


def test_failed_output_open_exits_nonzero_with_message(
        bench_exe: pathlib.Path, tmp_path: pathlib.Path) -> None:
    # A nonexistent output directory makes fopen fail; the driver must
    # report it and exit nonzero, never write through a null FILE*.
    missing = tmp_path / "no" / "such" / "dir"
    proc = subprocess.run([str(bench_exe), str(missing)], capture_output=True,
                          text=True, timeout=HARD_TIMEOUT_S)
    assert proc.returncode == 1
    assert "cannot open" in proc.stderr.lower()
    assert "efficiency.csv" in proc.stderr
    assert not list(tmp_path.rglob("efficiency.csv"))
