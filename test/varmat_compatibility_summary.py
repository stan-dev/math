import argparse
import json
import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import sig_utils

def main(results_file, functions, print_which, print_fully, print_names):
    """
    Generates benchmark code, compiles it and runs the benchmark. Optionally plots the results.
    :param results_file: File containing results of varmat compatibility calculation
    :param functions: List of function names to check compatibility for
    :param print_which: For which type of compatibility should functions be printed
    :param print_fully: Only print fully compatible/incompatible signatures (no partial varmat support)
    :param print_names: Print function names, not signatures
    """
    with open(results_file, "r") as f:
        results = json.load(f)

    compatible_signatures = set()
    incompatible_signatures = set()
    impossible_signatures = set()
    if "compatible_signatures" in results:
        compatible_signatures = set(results["compatible_signatures"])
    if "incompatible_signatures" in results:
        incompatible_signatures = set(results["incompatible_signatures"])
    if "impossible_signatures" in results:
        impossible_signatures = set(results["impossible_signatures"])

    requested_functions = set(functions)

    names_to_print = set()

    def convert_signatures_list_to_functions(signatures):
        functions = set()
        for signature in signatures:
            return_type, function_name, stan_args = sig_utils.parse_signature(signature)

            if len(requested_functions) > 0 and function_name not in requested_functions:
                continue

            functions.add(function_name)
        return functions

    compatible_functions = convert_signatures_list_to_functions(compatible_signatures)
    incompatible_functions = convert_signatures_list_to_functions(incompatible_signatures)
    impossible_functions = convert_signatures_list_to_functions(impossible_signatures)

    if print_fully:
        fully_compatible_functions = compatible_functions - incompatible_functions
        fully_incompatible_functions = incompatible_functions - compatible_functions
        fully_impossible_functions = impossible_functions - compatible_functions - incompatible_functions
        compatible_functions = fully_compatible_functions
        incompatible_functions = fully_incompatible_functions
        impossible_functions = fully_impossible_functions
    
    if print_names:
        if print_which == "compatible":
            names_to_print = compatible_functions
        elif print_which == "incompatible":
            names_to_print = incompatible_functions
        else:
            names_to_print = impossible_functions
    else:
        if print_which == "compatible":
            names_to_print = compatible_signatures
        elif print_which == "incompatible":
            names_to_print = incompatible_signatures
        else:
            names_to_print = impossible_signatures
    
    for name in sorted(names_to_print):
        if print_names:
            print(name.strip())
        else:
            print(name.replace("scalar_return_t", "real").strip())

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
    Define and process the command line interface to the varmat compatibility summarize script.
    """
    parser = FullErrorMsgParser(
        description="Process results of varmat compatibility test.",
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
        help = "Print by compatibility."
    )
    parser.add_argument(
        "--fully",
        default = False,
        help = "When printing compatible or incompatible function names, print those that are fully compatible or incompatible, ignoring impossible functions. When printing impossible functions, print only those with no compatible or incompatible versions.",
        action = "store_true"
    )
    parser.add_argument(
        "--names",
        default=False,
        help="Print function names, not signatures.",
        action="store_true",
    )
    args = parser.parse_args()

    main(
        args.results_file, args.functions,
        print_which = args.which,
        print_fully = args.fully,
        print_names = args.names
    )


if __name__ == "__main__":
    processCLIArgs()
