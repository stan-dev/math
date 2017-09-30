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
batchSize = 25


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
    tests_help_msg += "         'test/unit/math/prim/scal/fun/abs_test.cpp'"
    parser.add_argument("tests", nargs="+", type=str,
                        help=tests_help_msg)

    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
                        help="request additional script debugging output.")

    # And parse the command line against those rules
    return parser.parse_args()


def stopErr(msg, returncode):
    """Report an error message to stderr and exit with a given code."""
    sys.stderr.write('%s\n' % msg)
    sys.stderr.write('exit now (%s)\n' % time.strftime('%x %X %Z'))
    sys.exit(returncode)


def isWin():
    if (platform.system().lower().startswith("windows")
            or os.name.lower().startswith("windows")):
        return True
    return False


def mungeName(name):
    """Set up the makefile target name"""
    if (name.endswith(testsfx)):
        name = name.replace(testsfx, "_test")

        if (isWin()):
            name += winsfx
            name = name.replace("\\", "/")

    return name


def doCommand(command):
    """Run command as a shell command and report/exit on errors."""
    print("------------------------------------------------------------")
    print("%s" % command)
    p1 = subprocess.Popen(command, shell=True)
    p1.wait()
    if (not(p1.returncode is None) and not(p1.returncode == 0)):
        stopErr('%s failed' % command, p1.returncode)


def generateTests(j):
    """Generate all tests and pass along the j parameter to make."""
    doCommand('make -j%d generate-tests -s' % (j or 1))

def makeTest(name, j):
    """Run the make command for a given single test."""
    doCommand('make -j%d %s' % (j or 1, mungeName(name)))

def runTest(name):
    executable = mungeName(name).replace("/", os.sep)
    xml = mungeName(name).replace(winsfx, "")
    command = '%s --gtest_output="xml:%s.xml"' % (executable, xml)
    doCommand(command)

def findTests(args):
    folders = filter(os.path.isdir, args)
    nonfolders = list(set(args) - set(folders))
    tests = nonfolders + [os.path.join(root, n)
            for f in folders
            for root, _, names in os.walk(f)
            for n in names]
    return filter(lambda n: n.endswith(testsfx), tests)

def batched(tests):
    return [tests[i:i + batchSize] for i in range(0, len(tests), batchSize)]


def main():
    inputs = processCLIArgs()

    # pass 0: generate all auto-generated tests
    if any(['test/prob' in arg for arg in inputs.tests]):
        generateTests(inputs.j)

    tests = findTests(inputs.tests)
    if not tests:
        stopErr("No matching tests found.", -1)

    # pass 1: make test executables
    for batch in batched(tests):
        if inputs.debug:
            print("Test batch: ", batch)
        makeTest(" ".join(batch), inputs.j)

    # pass 2: run test targets
    for t in tests:
        if inputs.debug:
            print("run single test: %s" % testname)
        runTest(t)


if __name__ == "__main__":
    main()
