#!/usr/bin/python3

import os
import sys
import argparse
import subprocess
import shutil
import tempfile
import glob

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MICROMECHANICS_DIR = os.path.dirname(SCRIPT_DIR)
STANDARD_INPUTS_DIR = os.path.join(SCRIPT_DIR, "standard_inputs")
STANDARD_RESULTS_DIR = os.path.join(SCRIPT_DIR, "standard_results")

EXPECTED_OUTPUT_FILES = ["str_str.out", "vm.out", "err.out", "conv.out"]

PATH_PLACEHOLDER = "/absolute/path/to/"

DEFAULT_TOLERANCE = 1e-6

# ---------------------------------------------------------------------------
# Test definitions
# ---------------------------------------------------------------------------
# Each entry maps a test name to its configuration.  The "input_dir" tells
# the runner where the supporting parameter / microstructure files live so
# the placeholder paths inside the main input file can be resolved.
# ---------------------------------------------------------------------------
TESTS = [
    {
        "name": "EVPFFT_standalone",
        "solver": "evpfft",
        "input_file": "example_evpfft_standalone_inputfile.txt",
        "input_dir": os.path.join(
            STANDARD_INPUTS_DIR, "EVPFFT", "example_input_files"
        ),
        "nprocs": 1,
        "description": "EVPFFT standalone Cu 8x8x8 RVE (30 steps)",
        "comparison_file": "str_str.out",
    },
    {
        "name": "LS-EVPFFT_standalone",
        "solver": "ls-evpfft",
        "input_file": "example_evpfft_standalone_inputfile.txt",
        "input_dir": os.path.join(
            STANDARD_INPUTS_DIR, "EVPFFT", "example_input_files"
        ),
        "nprocs": 1,
        "description": "LS-EVPFFT standalone Cu 8x8x8 RVE (30 steps)",
        "comparison_file": "str_str.out",
    },
]

# ---------------------------------------------------------------------------
# Binary discovery
# ---------------------------------------------------------------------------

def find_binary(solver, user_binary=None):
    """Locate the solver executable, either from an explicit path or by
    scanning known build-directory patterns under the micromechanics tree."""

    if user_binary:
        binary = os.path.abspath(user_binary)
        if os.path.isfile(binary) and os.access(binary, os.X_OK):
            return binary
        raise FileNotFoundError(
            f"Binary not found or not executable: {binary}"
        )

    if solver == "evpfft":
        search_dir = os.path.join(MICROMECHANICS_DIR, "EVPFFT")
        patterns = [
            os.path.join(search_dir, "evpfft_*", "evpfft"),
            os.path.join(search_dir, "build*", "evpfft"),
        ]
    elif solver == "ls-evpfft":
        search_dir = os.path.join(MICROMECHANICS_DIR, "LS-EVPFFT")
        patterns = [
            os.path.join(search_dir, "ls-evpfft_*", "ls-evpfft"),
            os.path.join(search_dir, "build*", "ls-evpfft"),
        ]
    else:
        raise ValueError(f"Unknown solver: {solver}")

    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        for match in matches:
            if os.access(match, os.X_OK):
                return os.path.abspath(match)

    raise FileNotFoundError(
        f"Could not auto-detect {solver} binary under {search_dir}.\n"
        f"  Use --evpfft-binary or --ls-evpfft-binary to specify the path."
    )

# ---------------------------------------------------------------------------
# Input-file preparation
# ---------------------------------------------------------------------------

def prepare_input_file(input_file_path, input_dir, work_dir):
    """Copy the master input file into *work_dir*, replacing every
    occurrence of the placeholder prefix with the real absolute path to
    *input_dir* so that EVPFFT can locate its supporting files."""

    abs_input_dir = os.path.abspath(input_dir)

    with open(input_file_path, "r") as fh:
        content = fh.read()

    content = content.replace(PATH_PLACEHOLDER, abs_input_dir + "/")

    dest = os.path.join(work_dir, os.path.basename(input_file_path))
    with open(dest, "w") as fh:
        fh.write(content)

    return dest

# ---------------------------------------------------------------------------
# Output parsing & comparison
# ---------------------------------------------------------------------------

def parse_csv_output(filepath):
    """Return (headers, data_rows) from a comma-separated output file
    such as str_str.out or vm.out produced by EVPFFT."""

    with open(filepath, "r") as fh:
        header = fh.readline().strip().split(",")
        data = []
        for line in fh:
            line = line.strip()
            if line:
                data.append([float(v) for v in line.split(",")])
    return header, data


def compare_results(result_file, standard_file, tolerance=DEFAULT_TOLERANCE):
    """Compare a result CSV against its standard reference.
    Returns (passed: bool, errors: list[str])."""

    errors = []
    res_hdr, res_data = parse_csv_output(result_file)
    std_hdr, std_data = parse_csv_output(standard_file)

    if res_hdr != std_hdr:
        errors.append(
            f"Header mismatch:\n"
            f"  Result:   {res_hdr}\n"
            f"  Standard: {std_hdr}"
        )
        return False, errors

    if len(res_data) != len(std_data):
        errors.append(
            f"Row count mismatch: result has {len(res_data)} rows, "
            f"standard has {len(std_data)} rows"
        )
        return False, errors

    for row_idx, (r_row, s_row) in enumerate(zip(res_data, std_data)):
        if len(r_row) != len(s_row):
            errors.append(
                f"Step {row_idx + 1}: column count mismatch "
                f"({len(r_row)} vs {len(s_row)})"
            )
            continue
        for col_idx, (r_val, s_val) in enumerate(zip(r_row, s_row)):
            diff = abs(r_val - s_val)
            if diff > tolerance:
                col_name = (
                    res_hdr[col_idx]
                    if col_idx < len(res_hdr)
                    else f"col{col_idx}"
                )
                errors.append(
                    f"Step {row_idx + 1}, {col_name}: "
                    f"result={r_val:.10e}, standard={s_val:.10e}, "
                    f"diff={diff:.10e}"
                )

    return len(errors) == 0, errors

# ---------------------------------------------------------------------------
# MPI helpers
# ---------------------------------------------------------------------------

def mpi_cmd_prefix(nprocs):
    """Build the mpirun prefix, adding --allow-run-as-root when
    the process is running as UID 0 (common in CI containers)."""

    cmd = ["mpirun"]
    if os.getuid() == 0:
        cmd.append("--allow-run-as-root")
    cmd.extend(["-np", str(nprocs)])
    return cmd

# ---------------------------------------------------------------------------
# Test runner
# ---------------------------------------------------------------------------

def save_standard_results(test_name, work_dir):
    """Copy solver output files from *work_dir* into the standard-results
    tree so they can be used as the baseline for future regression runs."""

    dest_dir = os.path.join(STANDARD_RESULTS_DIR, test_name)
    os.makedirs(dest_dir, exist_ok=True)

    saved = []
    for fname in EXPECTED_OUTPUT_FILES:
        src = os.path.join(work_dir, fname)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(dest_dir, fname))
            saved.append(fname)

    return dest_dir, saved


def run_test(test, binary, nprocs_override=None, timeout=600,
             save_results=False):
    """Execute a single regression test.
    Returns (passed: bool, message: str)."""

    test_name = test["name"]
    input_dir = test["input_dir"]
    input_file = test["input_file"]
    nprocs = nprocs_override or test["nprocs"]

    input_file_path = os.path.join(input_dir, input_file)
    if not os.path.isfile(input_file_path):
        return False, f"Input file not found: {input_file_path}"

    work_dir = tempfile.mkdtemp(prefix=f"evpfft_reg_{test_name}_")

    try:
        prepared_input = prepare_input_file(
            input_file_path, input_dir, work_dir
        )

        cmd = mpi_cmd_prefix(nprocs) + [
            binary, f"--infile={prepared_input}"
        ]

        print(f"  Command: {' '.join(cmd)}")
        print(f"  Working dir: {work_dir}")

        result = subprocess.run(
            cmd,
            cwd=work_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )

        if result.returncode != 0:
            stdout_tail = result.stdout.decode("utf-8", errors="replace")[-500:]
            stderr_tail = result.stderr.decode("utf-8", errors="replace")[-500:]
            return False, (
                f"Solver exited with code {result.returncode}\n"
                f"  stdout (last 500 chars): ...{stdout_tail}\n"
                f"  stderr (last 500 chars): ...{stderr_tail}"
            )

        missing = [
            f for f in EXPECTED_OUTPUT_FILES
            if not os.path.isfile(os.path.join(work_dir, f))
        ]
        if missing:
            return False, f"Missing expected output files: {', '.join(missing)}"

        if save_results:
            dest_dir, saved = save_standard_results(test_name, work_dir)
            return True, (
                f"Standard results saved to {dest_dir}\n"
                f"    Files: {', '.join(saved)}"
            )

        comparison_file = test.get("comparison_file")
        if comparison_file:
            std_path = os.path.join(
                STANDARD_RESULTS_DIR, test_name, comparison_file
            )
            if os.path.isfile(std_path):
                res_path = os.path.join(work_dir, comparison_file)
                passed, errs = compare_results(res_path, std_path)
                if not passed:
                    summary = "\n    ".join(errs[:10])
                    extra = (
                        f"\n    ... and {len(errs) - 10} more errors"
                        if len(errs) > 10 else ""
                    )
                    return False, (
                        f"Result comparison failed:\n    {summary}{extra}"
                    )
            else:
                print(
                    f"  Note: No standard results at {std_path}\n"
                    f"        Running as smoke test only."
                )

        return True, "Solver completed successfully, output files verified"

    except subprocess.TimeoutExpired:
        return False, f"Solver timed out after {timeout} seconds"
    except Exception as exc:
        return False, f"Unexpected error: {exc}"
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Regression test suite for EVPFFT and LS-EVPFFT solvers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  python3 regression_tests.py\n"
            "  python3 regression_tests.py EVPFFT_standalone\n"
            "  python3 regression_tests.py --solver evpfft\n"
            "  python3 regression_tests.py --evpfft-binary /path/to/evpfft\n"
            "  python3 regression_tests.py --save-standard-results\n"
            "  python3 regression_tests.py --save-standard-results EVPFFT_standalone\n"
            "  python3 regression_tests.py --list\n"
        ),
    )
    parser.add_argument(
        "test_name", nargs="?",
        help="Run only this test (default: run all)",
    )
    parser.add_argument(
        "--evpfft-binary",
        help="Path to the EVPFFT executable (auto-detected if omitted)",
    )
    parser.add_argument(
        "--ls-evpfft-binary",
        help="Path to the LS-EVPFFT executable (auto-detected if omitted)",
    )
    parser.add_argument(
        "--solver", choices=["evpfft", "ls-evpfft"],
        help="Run only tests for a specific solver",
    )
    parser.add_argument(
        "--nprocs", type=int, default=None,
        help="Override MPI process count for every test",
    )
    parser.add_argument(
        "--timeout", type=int, default=600,
        help="Per-test timeout in seconds (default: 600)",
    )
    parser.add_argument(
        "--tolerance", type=float, default=DEFAULT_TOLERANCE,
        help=f"Numerical comparison tolerance (default: {DEFAULT_TOLERANCE})",
    )
    parser.add_argument(
        "--save-standard-results", action="store_true",
        help="Run solver(s) and save output as the new standard results "
             "baseline (use after valid code changes)",
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List available tests and exit",
    )
    args = parser.parse_args()

    if args.list:
        print("Available tests:")
        for t in TESTS:
            print(f"  {t['name']:30s} [{t['solver']}] {t['description']}")
        sys.exit(0)

    tests_to_run = list(TESTS)

    if args.test_name:
        tests_to_run = [t for t in tests_to_run if t["name"] == args.test_name]
        if not tests_to_run:
            available = ", ".join(t["name"] for t in TESTS)
            print(f"Error: Test '{args.test_name}' not found. "
                  f"Available: {available}")
            sys.exit(1)

    if args.solver:
        tests_to_run = [t for t in tests_to_run if t["solver"] == args.solver]
        if not tests_to_run:
            print(f"Error: No tests found for solver '{args.solver}'")
            sys.exit(1)

    user_binaries = {
        "evpfft": args.evpfft_binary,
        "ls-evpfft": args.ls_evpfft_binary,
    }

    solvers_needed = {t["solver"] for t in tests_to_run}
    binaries = {}
    for solver in solvers_needed:
        try:
            binaries[solver] = find_binary(solver, user_binaries.get(solver))
            print(f"Found {solver} binary: {binaries[solver]}")
        except FileNotFoundError as exc:
            print(f"Warning: {exc}")
            print(f"  Tests requiring '{solver}' will be skipped.\n")

    passed_tests = []
    failed_tests = []
    skipped_tests = []

    save_mode = args.save_standard_results

    print(f"\n{'=' * 60}")
    print("  EVPFFT / LS-EVPFFT Regression Test Suite")
    print(f"  Running {len(tests_to_run)} test(s)")
    if save_mode:
        print("  Mode: SAVE STANDARD RESULTS")
    print(f"{'=' * 60}\n")

    for test in tests_to_run:
        name = test["name"]
        solver = test["solver"]

        print(f"--- {name} ---")
        print(f"  {test['description']}")

        if solver not in binaries:
            print(f"  SKIPPED: {solver} binary not available\n")
            skipped_tests.append(name)
            continue

        passed, msg = run_test(
            test,
            binaries[solver],
            nprocs_override=args.nprocs,
            timeout=args.timeout,
            save_results=save_mode,
        )

        if passed:
            print(f"\n  *** PASSED: {msg} ***\n")
            passed_tests.append(name)
        else:
            print(f"\n  *** FAILED: {msg} ***\n")
            failed_tests.append(name)

    print(f"\n{'=' * 60}")
    print("  TEST SUMMARY")
    print(f"{'=' * 60}")
    total = len(passed_tests) + len(failed_tests) + len(skipped_tests)
    print(f"  Total:   {total}")
    print(f"  Passed:  {len(passed_tests)}")
    print(f"  Failed:  {len(failed_tests)}")
    print(f"  Skipped: {len(skipped_tests)}")

    if passed_tests:
        print("\n  Passed:")
        for t in passed_tests:
            print(f"    + {t}")

    if failed_tests:
        print("\n  Failed:")
        for t in failed_tests:
            print(f"    - {t}")

    if skipped_tests:
        print("\n  Skipped:")
        for t in skipped_tests:
            print(f"    ~ {t}")

    print()

    sys.exit(1 if failed_tests else 0)


if __name__ == "__main__":
    main()
