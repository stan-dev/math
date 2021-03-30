#!/usr/bin/python

import argparse
import json
import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from signature_parser import SignatureParser

def convert_signatures_list_to_functions(signatures):
    """
    Take a list of signatures and extract the unique set of functions they represent

    :param signatures: List of signatures
    """
    functions = set()
    for signature in signatures:
        sp = SignatureParser(signature)
        functions.add(sp.function_name)
    return functions

def select_signatures_matching_functions(signatures, functions):
    """
    Select from a collection of signatures those which have function names in the given function list

    :param signatures: List of signatures
    :param functions: Function names to match against
    """
    signatures_out = set()
    for signature in signatures:
        sp = SignatureParser(signature)

        if sp.function_name not in functions:
            continue

        signatures_out.add(signature)
    return signatures_out

def remove_signatures_matching_functions(signatures, functions):
    """
    Filter from a collection of signatures those which have function names in the given function list

    :param signatures: List of signatures
    :param functions: Function names to filter against
    """
    signatures_out = set()
    for signature in signatures:
        sp = SignatureParser(signature)

        if sp.function_name in functions:
            continue

        signatures_out.add(signature)
    return signatures_out

def process_results(results, functions, which, fully, names):
    """
    Processes results of varmat compatiblity results file and returns list of functions/signatures that
    match the given criteria

    :param results: Results of varmat compatibility calculation
    :param functions: List of function names to check compatibility for
    :param which: For which type of compatibility should functions be printed (possible values: compatible, incompatible, irrelevant)
    :param fully: Only print functions/signatures for which all signatures with the same function name are fully compatible/incompatible signatures (no partial varmat support)
    :param names: Print function names, not signatures
    """
    compatible_signatures = set()
    incompatible_signatures = set()
    irrelevant_signatures = set()
    if "compatible_signatures" in results:
        compatible_signatures = set(results["compatible_signatures"])
    if "incompatible_signatures" in results:
        incompatible_signatures = set(results["incompatible_signatures"])
    if "irrelevant_signatures" in results:
        irrelevant_signatures = set(results["irrelevant_signatures"])

    requested_functions = set(functions)

    if len(requested_functions) > 0:
        compatible_signatures = select_signatures_matching_functions(compatible_signatures, requested_functions)
        incompatible_signatures = select_signatures_matching_functions(incompatible_signatures, requested_functions)
        irrelevant_signatures = select_signatures_matching_functions(irrelevant_signatures, requested_functions)

    compatible_functions = convert_signatures_list_to_functions(compatible_signatures)
    incompatible_functions = convert_signatures_list_to_functions(incompatible_signatures)
    irrelevant_functions = convert_signatures_list_to_functions(irrelevant_signatures)

    if not names:
        if which == "compatible":
            if fully:
                return remove_signatures_matching_functions(compatible_signatures, incompatible_functions)
            else:
                return compatible_signatures
        elif which == "incompatible":
            if fully:
                return remove_signatures_matching_functions(incompatible_signatures, compatible_functions)
            else:
                return incompatible_signatures
        else:
            if fully:
                return remove_signatures_matching_functions(irrelevant_signatures, set(compatible_functions) | set(incompatible_functions))
            else:
                return irrelevant_signatures
    else:
        if which == "compatible":
            if fully:
                return compatible_functions - incompatible_functions
            else:
                return compatible_functions
        elif which == "incompatible":
            if fully:
                return incompatible_functions - compatible_functions
            else:
                return incompatible_functions
        else:
            if fully:
                return irrelevant_functions - compatible_functions - incompatible_functions
            else:
                return irrelevant_functions

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
        choices = ["compatible", "incompatible", "irrelevant"],
        help = "Print by compatibility."
    )
    parser.add_argument(
        "--fully",
        default = False,
        help = "When printing compatible or incompatible function names, print those that are fully compatible or incompatible, ignoring irrelevant functions. When printing irrelevant functions, print only those with no compatible or incompatible versions.",
        action = "store_true"
    )
    parser.add_argument(
        "--names",
        default=False,
        help="Print function names, not signatures.",
        action="store_true",
    )
    args = parser.parse_args()

    with open(args.results_file, "r") as f:
        results = json.load(f)

    names_to_print = process_results(
        results, args.functions,
        print_which = args.which,
        print_fully = args.fully,
        print_names = args.names
    )

    for name in sorted(names_to_print):
        print(name.strip())

if __name__ == "__main__":
    processCLIArgs()