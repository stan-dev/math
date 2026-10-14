#!/usr/bin/env python3
"""Build and run Stan Math tests with CMake."""

import os
from pathlib import Path
import re
import subprocess
import sys
import time
from argparse import ArgumentParser, RawTextHelpFormatter

maxChangedTests = 20
testsfx = "_test.cpp"

def processCLIArgs():
    """
    Define and process the command line interface to the runTests.py script.
    """
    cli_description = "Generate and run stan math library tests."
    cli_epilog = "See more information at: https://github.com/stan-dev/math"

    parser = ArgumentParser(
        description=cli_description,
        epilog=cli_epilog,
        formatter_class=RawTextHelpFormatter,
    )

    # Now define all the rules of the command line args and opts
    parser.add_argument(
        "-j", metavar="N", type=int, default=1, help="parallel build jobs (default: 1)"
    )

    parser.add_argument(
        "-e",
        metavar="M",
        type=int,
        default=-1,
        help="number of files to split expressions tests in",
    )
    f_help_msg = "Only tests with file names matching these will be executed.\n"
    f_help_msg += "Example: '-f chol', '-f opencl', '-f prim'"
    parser.add_argument("-f", type=str, default=[], action="append", help=f_help_msg)
    changed_help = (
        "Use git to determine which tests may have changed, and run only those"
    )
    parser.add_argument("--changed", action="store_true", help=changed_help)
    parser.add_argument(
        "-d",
        "--debug",
        dest="debug",
        action="store_true",
        help="request additional script debugging output.",
    )
    parser.add_argument(
        "-m",
        "--make-only",
        dest="make_only",
        action="store_true",
        help="Don't run tests, just try to make them.",
    )
    parser.add_argument(
        "-t",
        "--test-only",
        dest="test_only",
        action="store_true",
        help="Don't make tests, just try to run them (must already be made).",
    )
    parser.add_argument(
        "--run-all",
        dest="run_all",
        action="store_true",
        help="Don't stop at the first test failure, run all of them.",
    )
    parser.add_argument(
        "--only-functions",
        nargs="+",
        type=str,
        default=[],
        help="Function names to run expression tests for. Default: all functions",
    )
    parser.add_argument(
        "--jumbo-test",
        dest="do_jumbo",
        action="store_true",
        help="accepted for compatibility, directory tests already use one executable",
    )

    tests_help_msg = "The path(s) to the test case(s) to run.\n"
    tests_help_msg += "Example: 'test/unit', 'test/prob', and/or\n"
    tests_help_msg += "         'test/unit/math/prim/fun/abs_test.cpp'"
    parser.add_argument("tests", nargs="*", type=str, help=tests_help_msg)
    parser.add_argument("--build-dir", default="build", help="CMake build directory")
    parser.add_argument("--fresh", action="store_true",
                        help="reset the CMake cache before configuring (requires CMake 3.24 or newer)")
    parser.add_argument("--cmake-arg", action="append", default=[],
                        help="extra configure argument, e.g. --cmake-arg=-DSTAN_THREADS=ON")
    parser.add_argument("--config", default="Release", help="CMake build configuration")
    # And parse the command line against those rules
    inputs = parser.parse_args()
    if inputs.j < 1:
        parser.error("-j must be positive")
    if inputs.make_only and inputs.test_only:
        parser.error("--make-only and --test-only cannot be combined")
    if inputs.fresh and inputs.test_only:
        parser.error("--fresh cannot be combined with --test-only")
    return inputs


def stopErr(msg, returncode):
    """Report an error message to stderr and exit with a given code."""
    sys.stderr.write("%s\n" % msg)
    sys.stderr.write("exit now (%s)\n" % time.strftime("%x %X %Z"))
    sys.exit(returncode)


def findChangedTests(debug):
    changed_files = subprocess.run(
        ["git", "diff", "--name-only", "--diff-filter=d", "origin/develop...HEAD"],
        text=True, capture_output=True, check=True, cwd=Path(__file__).resolve().parent
    ).stdout.splitlines()
    if debug:
        print("Changed files:", changed_files)

    # changes in prim should also test other non-prim signatures
    test_deps = {"prim/": ["prim/", "rev/", "mix/", "fwd/"], "rev/": ["rev/", "mix/"], "fwd/": ["fwd/", "mix/"]}

    changed_tests = set()
    for f in changed_files:
        if 'stan/' in f and f.endswith('.hpp') and 'opencl' not in f: # changed fn
            maybe_test = f.replace('stan/', 'test/unit/').replace('.hpp', testsfx)
            for (path, replacements) in test_deps.items():
                if path in maybe_test:
                    for rep in replacements:
                        maybe_test = maybe_test.replace(path, rep)
                        if os.path.exists(maybe_test):
                            changed_tests.add(maybe_test)

        if f.endswith(testsfx) and 'opencl' not in f:
            changed_tests.add(f)

    if len(changed_tests) > maxChangedTests:
        stopErr("Number of changed tests excluded maximum, not running priority tests", 0)

    return sorted(changed_tests)


def main():
    inputs = processCLIArgs()
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "test"))
    from cmake import run
    try:
        sys.exit(run(inputs))
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        stopErr(str(error), 1)



if __name__ == "__main__":
    main()
