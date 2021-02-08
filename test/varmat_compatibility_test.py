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

def benchmark(
        functions_or_sigs,
        results_file="results.json",
        overloads=("Prim", "Rev"),
        max_size_param=None,
        cores=1,
):
    """
    Generates benchmark code, compiles it and runs the benchmark.
    :param functions_or_sigs: List of function names and/or signatures to benchmark
    :param cpp_filename: filename of cpp file to use
    :param overloads: Which overloads to benchmark
    :param multiplier_param: Multiplyer, by which to increase argument size.
    :param max_size_param: Maximum argument size.
    :param max_dim: Maximum number of argument dimensions to benchmark. Signatures with any argument with
             larger number of dimensions are skipped."
    :param n_repeats: Number of times to repeat each benchmark.
    :param skip_similar_signatures: Whether to skip similar signatures. Two signatures are similar if they
             difffer only in similar vector types, which are vector, row_vector and real[].
    :param csv_out_file: Filename of the csv file to store benchmark results in.
    """
    all_signatures = get_signatures()
    functions, signatures = handle_function_list(functions_or_sigs)

    requested_functions = set(functions)

    parsed_signatures = {}

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

    ref_signatures = set()
    for signature in all_signatures:
        return_type, function_name, stan_args = parse_signature(signature)

        if len(requested_functions) > 0:
            if function_name not in requested_functions:
                continue
        elif signature in skip_signatures:
            continue

        reference_args = tuple(reference_vector_argument(i) for i in stan_args)
        
        parsed_signatures[signature] = [return_type, function_name, stan_args]
        ref_signatures.add((function_name, reference_args))

    max_args_with_max_dimm = 0
    default_max_size = 2

    parsed_signatures_list = list(parsed_signatures.items())

    def chunk(mylist, chunk_size):
        for start in range(0, len(mylist), chunk_size):
            yield (start, mylist[start:min(start + chunk_size, len(mylist))])

    for start, parsed_signatures_chunk in chunk(parsed_signatures_list, cores):
        test_files_to_compile = {}

        for signature, (return_type, function_name, stan_args) in parsed_signatures_chunk:
            if any(arg.find("ode") >= 0 for arg in stan_args):
                test_files_to_compile[signature] = None
                continue

            result = ""
            dimm = 0
            args_with_max_dimm = 0
            for arg in stan_args:
                arg_dimm = 0
                if "vector" in arg:
                    arg_dimm = 1
                if "matrix" in arg:
                    arg_dimm = 2
                if "[" in arg:
                    arg_dimm += len(arg.split("[")[1])
                if arg_dimm == dimm:
                    args_with_max_dimm += 1
                elif arg_dimm > dimm:
                    dimm = arg_dimm
                    args_with_max_dimm = 1
            max_args_with_max_dimm = max(max_args_with_max_dimm, args_with_max_dimm)
            if max_size_param is None:
                if dimm == 0:  # signature with only scalar arguments
                    max_size = 1
                else:
                    max_size = default_max_size
            else:
                max_size = max_size_param

            cpp_arg_templates = []
            overload_opts = []

            for n, stan_arg in enumerate(stan_args):
                cpp_arg_template = get_cpp_type(stan_arg)
                arg_overload_opts = ["Prim"]
                if "SCALAR" in cpp_arg_template and not (
                        function_name in non_differentiable_args
                        and n in non_differentiable_args[function_name]
                ):
                    arg_overload_opts = ["Prim", "Rev", "Rev_SOA"]
                cpp_arg_templates.append(cpp_arg_template)
                overload_opts.append(arg_overload_opts)
                
            any_overload_uses_varmat = False
            for n, arg_overloads in enumerate(itertools.product(*overload_opts)):
                # generate one benchmark
                fg = FunctionGenerator(function_name, return_type, stan_args, arg_overloads, max_size)
                setup, code = fg.cpp()

                result += BENCHMARK_TEMPLATE.format(
                    benchmark_name=f"{function_name}_{n}",
                    setup=setup,
                    code=code,
                    max_size=max_size,
                )

                any_overload_uses_varmat |= fg.uses_varmat()

            if any_overload_uses_varmat:
                f = tempfile.NamedTemporaryFile("w", dir = WORKING_FOLDER, prefix = f"{function_name}_", suffix = "_test.cpp", delete = False)

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
        
                print("Check results of test {0} / {1}, {2} ... ".format(n + start, len(parsed_signatures), signature.strip()), end = '')

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


def main(
        functions_or_sigs, results_file, cores
):
    """
    Generates benchmark code, compiles it and runs the benchmark. Optionally plots the results.
    :param functions_or_sigs: List of function names and/or signatures to benchmark
    """
    benchmark(functions_or_sigs, results_file = results_file, cores = cores)


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
