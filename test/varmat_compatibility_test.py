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

def build_signature(signature, cpp_path):
    attempted = False
    successful = False

    if cpp_path:
        attempted = True

        object_path = cpp_path.replace(".cpp", ".o")
        dependency_path = cpp_path.replace(".cpp", ".d")
        stdout_path = cpp_path.replace(".cpp", ".stdout")
        stderr_path = cpp_path.replace(".cpp", ".stderr")
            
        successful, stdout, stderr = run_command([make, object_path])

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

def main(functions_or_sigs, results_file, cores, chunk_size):
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

    signatures_to_check = set()
    for signature in all_signatures:
        return_type, function_name, stan_args = parse_signature(signature)

        if len(requested_functions) > 0 and function_name not in requested_functions:
            continue
        
        signatures_to_check.add(signature)

    max_size = 2

    def chunk(myiterable, chunk_size):
        mylist = list(myiterable)
        for start in range(0, len(mylist), chunk_size):
            yield (start, mylist[start:min(start + chunk_size, len(mylist))])

    for start, signatures_to_check_chunk in chunk(signatures_to_check, chunk_size * cores):
        test_files_to_compile = {}

        for signature in signatures_to_check_chunk:
            fg = FunctionGenerator(signature)

            if fg.is_high_order():
                test_files_to_compile[signature] = None
                continue

            result = ""
            any_overload_uses_varmat = False
            
            for n, overloads in enumerate(itertools.product(("Prim", "Rev", "RevSOA"), repeat = fg.number_arguments())):
                fg.reset()

                arg_list_base = fg.build_arguments(overloads, max_size)

                arg_list = []
                for overload, arg in zip(overloads, arg_list_base):
                    if arg.is_reverse_mode() and arg.is_matrix_like() and overload.endswith("SOA"):
                        any_overload_uses_varmat = True
                        arg = fg.add(FunctionCallAssign("stan::math::to_var_value", arg.name + "_varmat", arg))
                    
                    arg_list.append(arg)

                fg.add(FunctionCallAssign(f"stan::math::{fg.function_name}", "result", *arg_list))

                result += BENCHMARK_TEMPLATE.format(
                    benchmark_name=f"{fg.function_name}_{n}",
                    code=fg.cpp(),
                )

            if any_overload_uses_varmat:
                f = tempfile.NamedTemporaryFile("w", dir = WORKING_FOLDER, prefix = f"{fg.function_name}_", suffix = "_test.cpp", delete = False)

                f.write("#include <test/expressions/expression_test_helpers.hpp>\n\n")
                f.write(result)

                f.close()

                test_files_to_compile[signature] = os.path.join(WORKING_FOLDER, os.path.basename(f.name))
            else:
                test_files_to_compile[signature] = None

        with concurrent.futures.ThreadPoolExecutor(cores) as e:
            futures = []
            for signature, cpp_path in test_files_to_compile.items():
                futures.append(e.submit(build_signature, signature, cpp_path))

            for n, finished in enumerate(concurrent.futures.as_completed(futures)):
                signature, successful, attempted = finished.result()

                print("Check results of test {0} / {1}, {2} ... ".format(n + start, len(signatures_to_check), signature.strip()), end = '')

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
        "-c",
        type=int,
        default=10,
        help="Work in tests in chunks of size c * j.",
    )
    parser.add_argument(
        "results_file",
        type=str,
        default=None,
        help="File to save results in.",
    )
    args = parser.parse_args()

    main(functions_or_sigs=args.functions, results_file = args.results_file, cores = args.j, chunk_size = args.c)


if __name__ == "__main__":
    processCLIArgs()
