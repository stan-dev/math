import argparse
import json
import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import sig_utils

def main(results_file, functions, print_which, print_names):
    """
    Generates benchmark code, compiles it and runs the benchmark. Optionally plots the results.
    :param functions_or_sigs: List of function names and/or signatures to benchmark
    """
    with open(results_file, "r") as f:
        results = json.load(f)

    if print_which == "compatible":
        signatures = results["compatible_signatures"]
    elif print_which == "incompatible":
        signatures = results["incompatible_signatures"]
    else:
        signatures = results["impossible_signatures"]

    requested_functions = set(functions)

    names_to_print = set()

    for signature in signatures:
        return_type, function_name, stan_args = sig_utils.parse_signature(signature)

        if len(requested_functions) > 0 and function_name not in requested_functions:
            continue

        if print_names:
            names_to_print.add(function_name)
        else:
            names_to_print.add(signature)
    
    for name in sorted(names_to_print):
        print(name.strip())

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
        "results_file",
        type=str,
        default=[],
        help="File with results of varmat compatibility test.",
    )
    parser.add_argument(
        "--functions",
        nargs="+",
        type=str,
        default=[],
        help="Function names to summarize. By default summarize everything.",
    )
    parser.add_argument(
        "--which",
        default = "compatible",
        choices = ["compatible", "incompatible", "impossible"],
        help = "Print which signatures."
    )
    parser.add_argument(
        "--print_names",
        default=False,
        help="Print function names, not signatures.",
        action="store_true",
    )
    args = parser.parse_args()

    main(args.results_file, args.functions, args.which, args.print_names)


if __name__ == "__main__":
    processCLIArgs()
