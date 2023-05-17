#!/usr/bin/python

"""
Replacement for runtest target in Makefile.

Call script with '-h' as an option to see a helpful message.
"""

from __future__ import print_function

import glob
import os
import os.path
import platform
import re
import subprocess
import sys
import time
from argparse import ArgumentParser, RawTextHelpFormatter

winsfx = ".exe"
testsfx = "_test.cpp"

allowed_paths_with_jumbo = [
    "test/unit/math/prim/",
    "test/unit/math/",
    "test/unit/math/rev/",
    "test/unit/math/fwd/",
    "test/unit/math/mix/",
    "test/unit/math/mix/fun/",
    "test/unit/math/opencl/",
    "test/unit/",
]

jumbo_folders = [
    "test/unit/math/prim/core",
    "test/unit/math/prim/err",
    "test/unit/math/prim/fun",
    "test/unit/math/prim/functor",
    "test/unit/math/prim/meta",
    "test/unit/math/prim/prob",
    "test/unit/math/rev/core",
    "test/unit/math/rev/err",
    "test/unit/math/rev/fun",
    "test/unit/math/rev/functor",
    "test/unit/math/rev/meta",
    "test/unit/math/rev/prob",
    "test/unit/math/fwd/core",
    "test/unit/math/fwd/fun",
    "test/unit/math/fwd/functor",
    "test/unit/math/fwd/meta",
    "test/unit/math/fwd/prob",
    "test/unit/math/mix/core",
    "test/unit/math/mix/fun",
    "test/unit/math/mix/functor",
    "test/unit/math/mix/meta",
    "test/unit/math/mix/prob",
    "test/unit/math/opencl/device_functions",
    "test/unit/math/opencl/kernel_generator",
    "test/unit/math/opencl/prim",
    "test/unit/math/opencl/rev",
]


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
        "-j", metavar="N", type=int, default=1, help="number of cores for make to use"
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
        help="Build/run jumbo tests.",
    )

    tests_help_msg = "The path(s) to the test case(s) to run.\n"
    tests_help_msg += "Example: 'test/unit', 'test/prob', and/or\n"
    tests_help_msg += "         'test/unit/math/prim/fun/abs_test.cpp'"
    parser.add_argument("tests", nargs="*", type=str, help=tests_help_msg)
    # And parse the command line against those rules
    return parser.parse_args()


def stopErr(msg, returncode):
    """Report an error message to stderr and exit with a given code."""
    sys.stderr.write("%s\n" % msg)
    sys.stderr.write("exit now (%s)\n" % time.strftime("%x %X %Z"))
    sys.exit(returncode)


def isWin():
    return platform.system().lower().startswith(
        "windows"
    ) or os.name.lower().startswith("windows")


batchSize = 20 if isWin() else 200
jumboSize = 5 if isWin() else 12
maxChangedTests = 20


def mungeName(name):
    """Set up the makefile target name"""
    if name.endswith(testsfx):
        name = name.replace(testsfx, "_test")

        if isWin():
            name += winsfx
            name = name.replace("\\", "/")
    return name


def doCommand(command, exit_on_failure=True):
    """Run command as a shell command and report/exit on errors."""
    print("------------------------------------------------------------")
    print("%s" % command)
    p1 = subprocess.Popen(command, shell=True)
    p1.wait()
    if exit_on_failure and (not (p1.returncode is None) and not (p1.returncode == 0)):
        stopErr("%s failed" % command, p1.returncode)


def generateTests(j):
    """Generate all tests and pass along the j parameter to make."""
    if isWin():
        doCommand("mingw32-make -j%d generate-tests -s" % (j or 1))
    else:
        doCommand("make -j%d generate-tests -s" % (j or 1))


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i : i + n]

include_pattern = re.compile(r'^\s*#include\s+<.*>$')

def generateJumboTests(paths, debug=False):
    jumbo_files_to_create = []
    jumbo_files = []
    for p in paths:
        if not p.endswith(testsfx) and not p.endswith("/"):
            p = p + "/"
        if p in allowed_paths_with_jumbo:
            jumbo_files_to_create.extend([x for x in jumbo_folders if x.startswith(p)])
        else:
            stopErr("The --jumbo flag is only allowed with top level folders.", 10)
    for jf in jumbo_files_to_create:
        tests_in_subfolder = sorted(os.path.join(jf, x) for x in os.listdir(jf) if x.endswith(testsfx))
        for i, tests in enumerate(divide_chunks(tests_in_subfolder, jumboSize)):

            includes = []
            jumbo_contents = []
            for t in tests:
                with open(t, "r") as test_file:
                    contents_raw = test_file.read()

                test_contents = []
                for l in contents_raw.splitlines():
                    if include_pattern.fullmatch(l.strip()):
                        includes.append(l.strip())
                        continue
                    test_contents.append(l)

                jumbo_contents.append('\n'.join(test_contents))

            jumbo_file_path = jf + "_" + str(i) + testsfx
            jumbo_files.append(jumbo_file_path)
            if debug:
                print("Generating jumbo test file '{}' with tests: {}".format( jumbo_file_path, ', '.join(tests)))
            with open(jumbo_file_path, "w") as jumbo_file:
                jumbo_file.write('\n'.join(list(dict.fromkeys(includes))))
                jumbo_file.write('\n')
                jumbo_file.write('\n'.join(jumbo_contents))
    return jumbo_files


def cleanupJumboTests(paths):
    for f in paths:
        if os.path.exists(f):
            os.remove(f)


def makeTest(name, j):
    """Run the make command for a given single test."""
    if isWin():
        doCommand("mingw32-make -j%d %s" % (j or 1, name))
    else:
        doCommand("make -j%d %s" % (j or 1, name))


def commandExists(command):
    p = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    p.wait()
    return p.returncode != 127


def runTest(name, run_all=False, mpi=False, j=1):
    executable = mungeName(name).replace("/", os.sep)
    xml = mungeName(name).replace(winsfx, "")
    command = '%s --gtest_output="xml:%s.xml"' % (executable, xml)
    if mpi:
        if not commandExists("mpirun"):
            stopErr(
                "Error: need to have mpi (and mpirun) installed to run mpi tests",
                -1,
            )
        if "mpi_" in name:
            j = j > 2 and j or 2
        else:
            j = 1
        command = "mpirun -np {} {}".format(j, command)
    doCommand(command, not run_all)


def test_files_in_folder(folder):
    """Returns a list of test files (*_test.cpp) in the folder and all
    its subfolders recursively. The folder can be written with
    wildcards as with the Unix find command.
    """
    files = []
    for f in glob.glob(folder):
        if os.path.isdir(f):
            files.extend(test_files_in_folder(f + os.sep + "**"))
        else:
            if f.endswith(testsfx):
                files.append(f)
    return files


def findTests(base_path, filter_names, do_jumbo=False):
    tests = []
    for path in base_path:
        if (not os.path.isdir(path)) and path.endswith("_test"):
            tests.append(path + ".cpp")
        else:
            tests.extend(test_files_in_folder(path))
    tests = map(mungeName, tests)
    tests = [
        test
        for test in tests
        if all(filter_name in test for filter_name in filter_names)
    ]
    if do_jumbo:
        filtered_jumbo_tests = []
        for t in tests:
            add = True
            for k in jumbo_folders:
                k = k + "/"
                if t.startswith(k):
                    add = False
                    break
            if add:
                filtered_jumbo_tests.append(t)
        tests = filtered_jumbo_tests
    return tests


def findChangedTests(debug):
    import subprocess

    changed_files = subprocess.run(
        ["git", "diff", "--name-only", "origin/develop...HEAD"], text=True, capture_output=True
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

    return list(map(mungeName, changed_tests))


def batched(tests):
    return [tests[i : i + batchSize] for i in range(0, len(tests), batchSize)]


def handleExpressionTests(tests, only_functions, n_test_files):
    expression_tests = False
    for n, i in list(enumerate(tests))[::-1]:
        if "test/expressions" in i or "test\\expressions" in i:
            del tests[n]
            expression_tests = True
    if expression_tests:
        HERE = os.path.dirname(os.path.realpath(__file__))
        sys.path.append(os.path.join(HERE, "test"))
        sys.path.append(os.path.join(HERE, "test/expressions"))
        import generate_expression_tests

        generate_expression_tests.main(only_functions, n_test_files)
        for i in range(n_test_files):
            tests.append("test/expressions/tests%d_test.cpp" % i)
    elif only_functions:
        stopErr(
            "--only-functions can only be specified if running expression tests (test/expressions)",
            -1,
        )


def checkToolchainPathWindows():
    if isWin():
        p1 = subprocess.Popen(
            "where.exe mingw32-make",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        out, err = p1.communicate()
        if re.search(" |\(|\)", out):
            stopErr(
                "The RTools toolchain is installed in a path with spaces or bracket. Please reinstall to a valid path.",
                -1,
            )


def main():
    inputs = processCLIArgs()
    checkToolchainPathWindows()

    try:
        with open("make/local") as f:
            stan_mpi = "STAN_MPI" in f.read()
    except IOError:
        stan_mpi = False

    jumboFiles = []

    if inputs.changed:
        tests = findChangedTests(inputs.debug)
    else:
        # pass 0: generate all auto-generated tests
        if any("test/prob" in arg for arg in inputs.tests):
            generateTests(inputs.j)

        if inputs.do_jumbo:
            jumboFiles = generateJumboTests(inputs.tests)
        if inputs.e == -1:
            if inputs.j == 1:
                num_expr_test_files = 1
            else:
                num_expr_test_files = inputs.j * 4
        else:
            num_expr_test_files = inputs.e
        handleExpressionTests(inputs.tests, inputs.only_functions, num_expr_test_files)

        tests = findTests(inputs.tests, inputs.f, inputs.do_jumbo)

    if not tests:
        if inputs.changed:
            stopErr("No changed tests were found!", 0)
        stopErr(
            "No matching tests found. Check the path passed on the command line", -1
        )

    try:
        # pass 1: make test executables
        for batch in batched(tests):
            if inputs.debug:
                print("Test batch: ", batch)
            makeTest(" ".join(batch), inputs.j)
        if not inputs.make_only:
            # pass 2: run test targets
            for t in tests:
                if inputs.debug:
                    print("run single test: %s" % t)
                runTest(t, inputs.run_all, mpi=stan_mpi, j=inputs.j)
    except BaseException as e:
       print(e, file=sys.stderr)
       sys.exit(1)
    finally:
        cleanupJumboTests(jumboFiles)
        pass


if __name__ == "__main__":
    main()
