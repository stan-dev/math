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
    if j is None:
        command = 'make generate-tests -s'
    else:
        command = 'make -j%d generate-tests -s' % j
    doCommand(command)


def makeTest(name, j):
    """Run the make command for a given single test."""
    target = mungeName(name)
    if j is None:
        command = 'make %s' % target
    else:
        command = 'make -j%d %s' % (j, target)
    doCommand(command)


def makeTests(dirname, filenames, j, debug):
    targets = list()
    for name in filenames:
        if (not name.endswith(testsfx)):
            continue

        target = "/".join([dirname, name])
        target = mungeName(target)
        targets.append(target)

    if (len(targets) > 0):
        if (debug):
            print('# targets: %d' % len(targets))
        startIdx = 0
        endIdx = batchSize

        while (startIdx < len(targets)):
            if j is None:
                command = 'make %s' % ' '.join(targets[startIdx:endIdx])
            else:
                command = 'make -j%d %s' % (j, ' '.join(targets[startIdx:endIdx]))

            if (debug):
                print('start %d, end %d' % (startIdx, endIdx))
                print(command)

            doCommand(command)
            startIdx = endIdx
            endIdx = startIdx + batchSize

            if (endIdx > len(targets)):
                endIdx = len(targets)


def runTest(name):
    executable = mungeName(name).replace("/", os.sep)
    xml = mungeName(name).replace(winsfx, "")
    command = '%s --gtest_output="xml:%s.xml"' % (executable, xml)
    doCommand(command)


def main():
    inputs = processCLIArgs()

    # pass 0: generate all auto-generated tests
    if any(['test/prob' in arg for arg in inputs.tests]):
        generateTests(inputs.j)

    # pass 1: call make to compile test targets
    for testname in inputs.tests:
        # Ensure that the test actually exists at all
        if (not(os.path.exists(testname))):
            stopErr('%s: no such file or directory' % testname, -1)

        # If it exists and is a not a directory, run that one test
        if (not(os.path.isdir(testname))):
            if (not(testname.endswith(testsfx))):
                stopErr('%s: not a testfile' % testname, -1)

            if (inputs.debug):
                print("make single test: %s" % testname)

            makeTest(testname, inputs.j)

        # If it exists and is a directory, run all tests inside
        else:
            for root, dirs, files in os.walk(testname):
                if (inputs.debug):
                    print("make root: %s" % root)

                makeTests(root, files, inputs.j, inputs.debug)

    # pass 2: run test targets
    for testname in inputs.tests:
        if (not(os.path.isdir(testname))):
            if (inputs.debug):
                print("run single test: %s" % testname)
            runTest(testname)
        else:
            for root, dirs, files in os.walk(testname):
                for name in files:
                    if (name.endswith(testsfx)):
                        if (inputs.debug):
                            print("run dir,test: %s,%s" % (root, name))
                        runTest(os.sep.join([root, name]))


if __name__ == "__main__":
    main()
