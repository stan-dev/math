#!/usr/bin/python

import itertools
import json
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import os
import queue as Queue
import subprocess
import re
import sys
import tempfile
import threading

from sig_utils import make, handle_function_list, get_signatures
from signature_parser import SignatureParser
from code_generator import CodeGenerator

HERE = os.path.dirname(os.path.realpath(__file__))
TEST_FOLDER = os.path.abspath(os.path.join(HERE, "..", "test"))
sys.path.append(TEST_FOLDER)
WORKING_FOLDER = "test/varmat-compatibility"

TEST_TEMPLATE = """
static void {test_name}() {{
{code}
}}
"""

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

def build_signature(prefix, cpp_code, debug):
    """
    Try to build the given cpp code

    Return true if the code was successfully built

    :param prefix: Prefix to give file names so easier to debug
    :param cpp_code: Code to build
    :param debug: If true, don't delete temporary files
    """
    f = tempfile.NamedTemporaryFile("w", dir = WORKING_FOLDER, prefix = prefix + "_", suffix = "_test.cpp", delete = False)
    f.write("#include <test/expressions/expression_test_helpers.hpp>\n\n")
    f.write(cpp_code)
    f.close()

    cpp_path = os.path.join(WORKING_FOLDER, os.path.basename(f.name))

    object_path = cpp_path.replace(".cpp", ".o")
    dependency_path = cpp_path.replace(".cpp", ".d")
    stdout_path = cpp_path.replace(".cpp", ".stdout")
    stderr_path = cpp_path.replace(".cpp", ".stderr")

    successful, stdout, stderr = run_command([make, object_path])

    if successful or not debug:
        try:
            os.remove(cpp_path)
        except OSError:
            pass

        try:
            os.remove(dependency_path)
        except OSError:
            pass

        try:
            os.remove(object_path)
        except OSError:
            pass
    else:
        if debug:
            with open(stdout_path, "w") as stdout_f:
                stdout_f.write(stdout.decode("utf-8"))

            with open(stderr_path, "w") as stderr_f:
                stderr_f.write(stderr.decode("utf-8"))
        
    return successful

def main(functions_or_sigs, results_file, cores, debug):
    """
    Attempt to build all the signatures in functions_or_sigs, or all the signatures
    associated with all the functions in functions_or_sigs, or if functions_or_sigs
    is empty every signature the stanc3 compiler exposes.

    Results are written to a results json file. Individual signatures are classified
    as either compatible, incompatible, or irrelevant.

    Compatible signatures can be compiled with varmat types in every argument that
    could possibly be a varmat (the matrix-like ones).

    Incompatible signatures cannot all be built, and for irrelevant signatures it does
    not make sense to try to build them (there are no matrix arguments, or the function
    does not support reverse mode autodiff, etc).

    Compilation is done in parallel using the number of specified cores.

    :param functions_or_sigs: List of function names and/or signatures to benchmark
    :param results_file: File to use as a results cache
    :param cores: Number of cores to use for compiling
    :param debug: If true, don't delete temporary files
    """
    all_signatures = get_signatures()
    functions, signatures = handle_function_list(functions_or_sigs)

    requested_functions = set(functions)

    compatible_signatures = set()
    incompatible_signatures = set()
    irrelevant_signatures = set()

    # Read the arguments and figure out the exact list of signatures to test
    signatures_to_check = set()
    for signature in all_signatures:
        sp = SignatureParser(signature)

        if len(requested_functions) > 0 and sp.function_name not in requested_functions:
            continue
        
        signatures_to_check.add(signature)

    work_queue = Queue.Queue()

    # For each signature, generate cpp code to test
    for signature in signatures_to_check:
        sp = SignatureParser(signature)

        if sp.is_high_order():
            work_queue.put((n, signature, None))
            continue

        cpp_code = ""
        any_overload_uses_varmat = False
        
        for m, overloads in enumerate(itertools.product(("Prim", "Rev", "RevVarmat"), repeat = sp.number_arguments())):
            cg = CodeGenerator()

            arg_list_base = cg.build_arguments(sp, overloads, size = 1)

            arg_list = []
            for overload, arg in zip(overloads, arg_list_base):
                if arg.is_reverse_mode() and arg.is_varmat_compatible() and overload.endswith("Varmat"):
                    any_overload_uses_varmat = True
                    arg = cg.to_var_value(arg)
                
                arg_list.append(arg)

            cg.function_call_assign("stan::math::" + sp.function_name, *arg_list)

            cpp_code += TEST_TEMPLATE.format(
                test_name = sp.function_name + repr(m),
                code=cg.cpp(),
            )

        if any_overload_uses_varmat:
            work_queue.put((work_queue.qsize(), signature, cpp_code))
        else:
            print("{0} ... Irrelevant".format(signature.strip()))
            irrelevant_signatures.add(signature)

    output_lock = threading.Lock()

    if not os.path.exists(WORKING_FOLDER):
        os.mkdir(WORKING_FOLDER)

    work_queue_original_length = work_queue.qsize()

    # Test if each cpp file builds and update the output file
    # This part is done in parallel
    def worker():
        while True:
            try:
                n, signature, cpp_code = work_queue.get(False)
            except Queue.Empty:
                return # If queue is empty, worker quits

            # Use signature as filename prefix to make it easier to find
            prefix = re.sub('[^0-9a-zA-Z]+', '_', signature.strip())

            # Test the signature
            successful = build_signature(prefix, cpp_code, debug)

            # Acquire a lock to do I/O
            with output_lock:
                if successful:
                    result_string = "Success"
                    compatible_signatures.add(signature)
                else:
                    result_string = "Fail"
                    incompatible_signatures.add(signature)

                print("Results of test {0} / {1}, {2} ... ".format(n, work_queue_original_length, signature.strip()) + result_string)

            work_queue.task_done()

    for i in range(cores):
        threading.Thread(target = worker).start()

    work_queue.join()

    with open(results_file, "w") as f:
        json.dump({ "compatible_signatures" : list(compatible_signatures),
                    "incompatible_signatures" : list(incompatible_signatures),
                    "irrelevant_signatures" : list(irrelevant_signatures)
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
        help="Signatures and/or function names to benchmark.",
    )
    parser.add_argument(
        "-j",
        type=int,
        default=1,
        help="Number of parallel cores to use.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Keep cpp, stdout, and stderr for incompatible functions.",
    )
    parser.add_argument(
        "results_file",
        type=str,
        default=None,
        help="File to save results in.",
    )
    args = parser.parse_args()

    main(functions_or_sigs=args.functions, results_file = args.results_file, cores = args.j, debug = args.debug)

if __name__ == "__main__":
    processCLIArgs()
