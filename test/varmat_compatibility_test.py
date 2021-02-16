#!/usr/bin/python
from __future__ import print_function

import contextlib
import concurrent.futures
import itertools
import json
import numbers
import os
import subprocess
import sys
import tempfile
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

HERE = os.path.dirname(os.path.realpath(__file__))
TEST_FOLDER = os.path.abspath(os.path.join(HERE, "..", "test"))
sys.path.append(TEST_FOLDER)
WORKING_FOLDER = "test/varmat-compatibility"

from sig_utils import *

BENCHMARK_TEMPLATE = """
static void {benchmark_name}() {{
{setup}
{code}
}}
"""

overload_scalar = {
    "Prim": "double",
    "Rev": "stan::math::var",
    "Fwd": "stan::math::fvar<double>",
    "Mix": "stan::math::fvar<stan::math::var>",
}


def run_command(command):
    """
    Runs given command and waits until it finishes executing.
    :param command: command to execute
    """
    proc = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = proc.communicate()
    
    if proc.poll() == 0:
        return (True, stdout, stderr)
    else:
        return (False, stdout, stderr)

def build(exe_filepath):
    """
    Builds a file using make.
    :param exe_filepath: File to build
    """
    return run_command([make, exe_filepath])

def main(functions_or_sigs, results_file, cores):
    """
    yada yada
    :param functions_or_sigs: List of function names and/or signatures to benchmark
    :param results_file: File to use as a results cache
    :param cores: Number of cores to use for compiling
    """
    all_signatures = get_signatures()
    functions, signatures = handle_function_list(functions_or_sigs)

    requested_functions = set(functions)

    compatible_signatures = set()
    incompatible_signatures = set()
    impossible_signatures = set()
    if results_file and os.path.exists(results_file):
        with open(results_file, "r") as f:
            results_file_contents = json.load(f)

            if "compatible_signatures" in results_file_contents:
                compatible_signatures = set(results_file_contents["compatible_signatures"])
            
            if "incompatible_signatures" in results_file_contents:
                incompatible_signatures = set(results_file_contents["incompatible_signatures"])

            if "impossible_signatures" in results_file_contents:
                impossible_signatures = set(results_file_contents["impossible_signatures"])

    skip_signatures = set(compatible_signatures)

    signatures_to_check = set()
    for signature in all_signatures:
        return_type, function_name, stan_args = parse_signature(signature)

        if len(requested_functions) > 0:
            if function_name not in requested_functions:
                continue
        elif signature in skip_signatures:
            continue
        
        signatures_to_check.add(signature)

    max_size = 2

    def chunk(myiterable, chunk_size):
        mylist = list(myiterable)
        for start in range(0, len(mylist), chunk_size):
            yield (start, mylist[start:min(start + chunk_size, len(mylist))])

    for start, signatures_to_check_chunk in chunk(signatures_to_check, cores):
        test_files_to_compile = {}

        for signature in signatures_to_check_chunk:
            fg = FunctionGenerator(signature)

            if fg.is_high_order():
                test_files_to_compile[signature] = None
                continue

            result = ""
            any_overload_uses_varmat = False
            for n, arg_overloads in enumerate(fg.overloads()):
                # generate one benchmark
                setup, code, uses_varmat = fg.cpp(arg_overloads, max_size)

                result += BENCHMARK_TEMPLATE.format(
                    benchmark_name=f"{fg.function_name}_{n}",
                    setup=setup,
                    code=code,
                )
                any_overload_uses_varmat |= uses_varmat

            if any_overload_uses_varmat:
                f = tempfile.NamedTemporaryFile("w", dir = WORKING_FOLDER, prefix = f"{fg.function_name}_", suffix = "_test.cpp", delete = False)

                f.write("#include <test/expressions/expression_test_helpers.hpp>\n\n")
                f.write(result)

                f.close()

                test_files_to_compile[signature] = os.path.join(WORKING_FOLDER, os.path.basename(f.name))
            else:
                test_files_to_compile[signature] = None

        def do_build(signature, cpp_path):
            attempted = False
            successful = False

            if cpp_path:
                attempted = True

                object_path = cpp_path.replace(".cpp", ".o")
                dependency_path = cpp_path.replace(".cpp", ".d")
                stdout_path = cpp_path.replace(".cpp", ".stdout")
                stderr_path = cpp_path.replace(".cpp", ".stderr")
                    
                successful, stdout, stderr = build(object_path)

                # Only clean up files if check succesful
                if successful:
                    with contextlib.suppress(FileNotFoundError):
                        os.remove(cpp_path)
                        os.remove(dependency_path)
                        os.remove(object_path)
                else:
                    with open(stdout_path, "w") as stdout_f:
                        stdout_f.write(stdout.decode("utf-8"))

                    with open(stderr_path, "w") as stderr_f:
                        stderr_f.write(stderr.decode("utf-8"))
                
            return signature, successful, attempted

        with concurrent.futures.ThreadPoolExecutor(cores) as e:
            futures = []
            for signature, cpp_path in test_files_to_compile.items():
                futures.append(e.submit(do_build, signature, cpp_path))

            for n, finished in enumerate(concurrent.futures.as_completed(futures)):
                signature, successful, attempted = finished.result()
        
                print("Check results of test {0} / {1}, {2} ... ".format(n + start, len(signatures_to_check), signature.strip()), end = '')

                if signature in compatible_signatures:
                    compatible_signatures.remove(signature)

                if signature in incompatible_signatures:
                    incompatible_signatures.remove(signature)

                if signature in impossible_signatures:
                    impossible_signatures.remove(signature)

                if attempted:
                    if successful:
                        print("Success!")
                        compatible_signatures.add(signature)
                    else:
                        print("Fail!")
                        incompatible_signatures.add(signature)
                else:
                    print("Impossible!")
                    impossible_signatures.add(signature)

                if results_file:
                    with open(results_file, "w") as f:
                        json.dump({
                            "compatible_signatures" : list(compatible_signatures),
                            "incompatible_signatures" : list(incompatible_signatures),
                            "impossible_signatures" : list(impossible_signatures)
                        }, f, indent = 4, sort_keys = True)


class FullErrorMsgParser(ArgumentParser):
    """
    Modified ArgumentParser that prints full error message on any error.
    """

    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


def processCLIArgs():
    """
    Define and process the command line interface to the benchmark.py script.
    """
    parser = FullErrorMsgParser(
        description="Generate and run_command benchmarks.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--functions",
        nargs="+",
        type=str,
        default=[],
        help="Signatures and/or function names to benchmark. Ignores any finished checks in results file (if given).",
    )
    parser.add_argument(
        "-j",
        type=int,
        default=1,
        help="Number of parallel cores to use.",
    )
    parser.add_argument(
        "--cpp",
        metavar="filename",
        type=str,
        default="benchmark.cpp",
        help="Filename of the cpp file to generate.",
    )
    parser.add_argument(
        "--results_file",
        type=str,
        default=None,
        help="File to save results in. If it already exists, it is used to skip already finished checks.",
    )
    args = parser.parse_args()

    main(functions_or_sigs=args.functions, results_file = args.results_file, cores = args.j)


if __name__ == "__main__":
    processCLIArgs()
