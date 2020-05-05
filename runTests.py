#!/usr/bin/python

"""
Replacement for runtest target in Makefile.

Call script with '-h' as an option to see a helpful message.
"""

from __future__ import print_function
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import os.path
import platform
import subprocess
import sys
import time

winsfx = ".exe"
testsfx = "_test.cpp"

allowed_paths_with_jumbo = [
    "test/unit/math/prim/",
    "test/unit/math/rev/",
    "test/unit/math/fwd/",
    "test/unit/math/mix/",
    "test/unit/math/rev/",
    "test/unit/math/rev/",
    "test/unit/math/",
    "test/unit/",
]

jumbo_folders = [
    #"test/unit/math/prim/core",
    #"test/unit/math/prim/err",
    "test/unit/math/prim/fun",
    #"test/unit/math/prim/functor",
    "test/unit/math/prim/meta",
    "test/unit/math/prim/prob",
    #"test/unit/math/rev/core",
    #"test/unit/math/rev/err",
    "test/unit/math/rev/fun",
    #"test/unit/math/rev/functor",
    "test/unit/math/rev/meta",
    #"test/unit/math/rev/prob",
    "test/unit/math/fwd/core",
    "test/unit/math/fwd/fun",
    "test/unit/math/fwd/functor",
    "test/unit/math/fwd/meta",
    "test/unit/math/fwd/prob",
    "test/unit/math/mix/core",
    #"test/unit/math/mix/fun",
    #"test/unit/math/mix/functor",
    #"test/unit/math/mix/meta",
    "test/unit/math/mix/prob"
]

def processCLIArgs():
    """
    Define and process the command line interface to the runTests.py script.
    """
    cli_description = "Generate and run stan math library tests."
    cli_epilog = "See more information at: https://github.com/stan-dev/math"

    parser = ArgumentParser(description=cli_description,
                            epilog=cli_epilog,
                            formatter_class=RawTextHelpFormatter)

    # Now define all the rules of the command line args and opts
    parser.add_argument("-j", metavar="N", type=int, default=1,
                        help="number of cores for make to use")

    tests_help_msg = "The path(s) to the test case(s) to run.\n"
    tests_help_msg += "Example: 'test/unit', 'test/prob', and/or\n"
    tests_help_msg += "         'test/unit/math/prim/fun/abs_test.cpp'"
    parser.add_argument("tests", nargs="+", type=str,
                        help=tests_help_msg)
    f_help_msg = "Only tests with file names matching these will be executed.\n"
    f_help_msg += "Example: '-f chol', '-f opencl', '-f prim'"
    parser.add_argument("-f", type=str, default = [], action="append",
                        help=f_help_msg)
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
                        help="request additional script debugging output.")
    parser.add_argument("-m", "--make-only", dest="make_only",
                        action="store_true", help="Don't run tests, just try to make them.")
    parser.add_argument("--run-all", dest="run_all", action="store_true",
                        help="Don't stop at the first test failure, run all of them.")
    parser.add_argument("--jumbo-test", dest="do_jumbo", action="store_true",
                        help="Build/run jumbo tests.")

    # And parse the command line against those rules
    return parser.parse_args()


def stopErr(msg, returncode):
    """Report an error message to stderr and exit with a given code."""
    sys.stderr.write('%s\n' % msg)
    sys.stderr.write('exit now (%s)\n' % time.strftime('%x %X %Z'))
    sys.exit(returncode)


def isWin():
    return (platform.system().lower().startswith("windows")
            or os.name.lower().startswith("windows"))


batchSize = 20 if isWin() else 200


def mungeName(name):
    """Set up the makefile target name"""
    if (name.endswith(testsfx)):
        name = name.replace(testsfx, "_test")

        if (isWin()):
            name += winsfx
            name = name.replace("\\", "/")

    return name

def doCommand(command, exit_on_failure=True):
    """Run command as a shell command and report/exit on errors."""
    print("------------------------------------------------------------")
    print("%s" % command)
    p1 = subprocess.Popen(command, shell=True)
    p1.wait()
    if exit_on_failure and (not(p1.returncode is None) and not(p1.returncode == 0)):
        stopErr('%s failed' % command, p1.returncode)


def generateTests(j):
    """Generate all tests and pass along the j parameter to make."""
    if isWin():
        doCommand('mingw32-make -j%d generate-tests -s' % (j or 1))
    else:
        doCommand('make -j%d generate-tests -s' % (j or 1))

def divide_chunks(l, n): 
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n] 

def generateJumboTests(paths):
    jumbo_files_to_create = []
    jumbo_files = []
    for p in paths:
        if not p.endswith(testsfx) and not p.endswith("/"):
            p = p + "/"
        if p in allowed_paths_with_jumbo:
            jumbo_files_to_create.extend(
                [x for x in jumbo_folders if x.startswith(p)]
            )
        else:
            stopErr('The --jumbo flag is only allowed with top level folders.', 10)
    for jf in jumbo_files_to_create:
        tests_in_subfolder = sorted([x for x in os.listdir(jf) if x.endswith(testsfx)])
        chunked_tests = divide_chunks(tests_in_subfolder, 30)
        i = 0
        for tests in chunked_tests:
            i = i + 1
            jumbo_file_path = jf + "_" + str(i) + testsfx
            jumbo_files.append(jumbo_file_path)
            f = open(jumbo_file_path, "w")
            for t in tests:
                f.write("#include <"+jf+"/"+t+">\n")
            f.close()
    return jumbo_files

def cleanupJumboTests(paths):
    for f in paths:
        if os.path.exists(f):
            os.remove(f)

def makeTest(name, j):
    """Run the make command for a given single test."""
    if isWin():
        doCommand('mingw32-make -j%d %s' % (j or 1, name))
    else:
        doCommand('make -j%d %s' % (j or 1, name))

def commandExists(command):
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p.wait()
    return p.returncode != 127

def runTest(name, run_all=False, mpi=False, j=1):
    executable = mungeName(name).replace("/", os.sep)
    xml = mungeName(name).replace(winsfx, "")
    command = '%s --gtest_output="xml:%s.xml"' % (executable, xml)
    if mpi:
        if not commandExists("mpirun"):
            stopErr("Error: need to have mpi (and mpirun) installed to run mpi tests"
                    + "\nCheck https://github.com/stan-dev/stan/wiki/Parallelism-using-MPI-in-Stan for more details."
                    , -1)
        if "mpi_" in name:
            j = j > 2 and j or 2
        else:
            j = 1
        command = "mpirun -np {} {}".format(j, command)
    doCommand(command, not run_all)

def findTests(base_path, filter_names, do_jumbo = False):
    folders = filter(os.path.isdir, base_path)
    nonfolders = list(set(base_path) - set(folders))
    tests = nonfolders + [os.path.join(root, n)
            for f in folders
            for root, _, names in os.walk(f)
            for n in names
            if n.endswith(testsfx)]
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
    tests = map(mungeName, tests)
    tests = [test
            for test in tests
            if all(filter_name in test
                   for filter_name in filter_names)]
    return tests

def batched(tests):
    return [tests[i:i + batchSize] for i in range(0, len(tests), batchSize)]


def main():
    inputs = processCLIArgs()

    try:
        with open("make/local") as f:
            stan_mpi =  "STAN_MPI" in f.read()
    except IOError:
        stan_mpi = False

    # pass 0: generate all auto-generated tests
    if any(['test/prob' in arg for arg in inputs.tests]):
        generateTests(inputs.j)
    jumboFiles = []
    if inputs.do_jumbo:
        jumboFiles = generateJumboTests(inputs.tests)
    
    tests = findTests(inputs.tests, inputs.f, inputs.do_jumbo)
    if not tests:
        stopErr("No matching tests found.", -1)
    if inputs.debug:
        print("Collected the following tests:\n", tests)
    #pass 1: make test executables
    for batch in batched(tests):
        if inputs.debug:
            print("Test batch: ", batch)
        makeTest(" ".join(batch), inputs.j)

    if not inputs.make_only:
        # pass 2: run test targets
        for t in tests:
            if inputs.debug:
                print("run single test: %s" % testname)
            runTest(t, inputs.run_all, mpi = stan_mpi, j = inputs.j)

    cleanupJumboTests(jumboFiles)

if __name__ == "__main__":
    main()
