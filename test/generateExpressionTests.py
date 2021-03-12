import numbers
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from sig_utils import *
from function_generator import *

src_folder = "./test/expressions/"
build_folder = "./test/expressions/"


test_code_template = """
TEST(ExpressionTest{overload}, {test_name}) {{
{code}
}}
"""

class FullErrorMsgParser(ArgumentParser):
    """
    Modified ArgumentParser that prints full error message on any error.
    """

    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)

def save_tests_in_files(N_files, tests):
    """
    Saves tests in files
    :param N_files: number of files to distribute tests into
    :param tests: list of test sources
    """
    for i in range(N_files):
        start = i * len(tests) // N_files
        end = (i + 1) * len(tests) // N_files
        with open(src_folder + "tests%d_test.cpp" % i, "w") as out:
            out.write("#include <test/expressions/expression_test_helpers.hpp>\n\n")
            for test in tests[start:end]:
                out.write(test)


def main(functions=(), j=1):
    """
    Generates expression tests.

    For every signature prim, rev and fwd instantiations are tested (with all scalars
    of type double/var/fvar<double>). Tests check the following:
     - signatures can be compiled with expressions arguments
     - results when given expressions are same as when given plain matrices
       (including derivatives)
     - functions evaluate expressions at most once

    :param functions: functions to generate tests for. This can contain names of functions
    already supported by stanc3, full function signatures or file names of files containing
    any of the previous two. Default: all signatures supported by stanc3
    :param j: number of files to split tests in
    """

    tests = []
    functions, signatures = handle_function_list(functions)

    signatures = set(signatures)

    if not functions and not signatures:
        default_checks = True
        signatures |= get_signatures()
    else:
        for signature in get_signatures():
            return_type, function_name, function_args = parse_signature(signature)
            if function_name in functions:
                signatures.add(signature)
        default_checks = False

    for n, signature in enumerate(signatures):
        for overload in ("Prim", "Rev", "Fwd"):
            fg = FunctionGenerator(signature)

            # skip ignored signatures if running the default checks
            if default_checks and fg.function_name in ignored:
                continue

            # skip signatures without inputs that can be eigen types
            if not fg.has_vector_arg():
                continue

            if fg.is_rng() and not overload == "Prim":
                function_args.append("rng")
            if fg.no_fwd_overload() and overload == "Fwd":
                continue
            if fg.no_rev_overload() and overload == "Rev":
                continue

            arg_list_no_expression = fg.build_arguments(len(fg.stan_args) * [overload], 2)
            arg_list_expression_base = fg.build_arguments(len(fg.stan_args) * [overload], 2)

            is_reverse_mode = not fg.returns_int() and overload == "Rev"

            arg_list_expression = []
            for arg in arg_list_expression_base:
                if arg.is_matrix_like():
                    arg = fg.add(ExpressionArgument(overload, arg.name + "_expr", arg))
                
                arg_list_expression.append(arg)

            # Check results with expressions and without are the same
            result = fg.add(FunctionCallAssign(f"stan::math::{fg.function_name}", "result", *arg_list_expression))
            result_no_expression = fg.add(FunctionCallAssign(f"stan::math::{fg.function_name}", "result_no_expression", *arg_list_no_expression))
            fg.add(FunctionCall("EXPECT_STAN_EQ", result, result_no_expression))

            # Check that expressions evaluated only once
            for arg in arg_list_expression:
                if arg.is_expression():
                    fg.add(FunctionCall("EXPECT_LEQ_ONE", arg.counter))

            # If reverse mode, check the adjoints            
            if is_reverse_mode:
                summed_result = fg.add(FunctionCallAssign("stan::test::recursive_sum", "summed_result", result))
                summed_result_no_expression = fg.add(FunctionCallAssign("stan::test::recursive_sum", "summed_result_no_expression", result_no_expression))
                sum_of_sums = fg.add(FunctionCallAssign("stan::math::add", "sum_of_sums", summed_result, summed_result_no_expression))
                fg.add(FunctionCall("stan::test::grad", sum_of_sums))

                for arg, arg_no_expression in zip(arg_list_no_expression, arg_list_expression):
                    fg.add(FunctionCall("stan::test::expect_adj_eq", arg, arg_no_expression))

                fg.add(FunctionCall("stan::math::recover_memory"))

            tests.append(
                test_code_template.format(
                    overload=overload,
                    test_name=f"{fg.function_name}_{n}",
                    code = fg.cpp(),
                )
            )
    
    save_tests_in_files(j, tests)

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
    args = parser.parse_args()

    main(functions=args.functions, j = args.j)

if __name__ == "__main__":
    processCLIArgs()

