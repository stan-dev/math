#!/usr/bin/python

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import numbers
import sys

from sig_utils import handle_function_list, get_signatures
from signature_parser import SignatureParser
from code_generator import CodeGenerator

src_folder = "./test/expressions/"
build_folder = "./test/expressions/"

test_code_template = """
TEST(ExpressionTest{overload}, {test_name}) {{
/*
{comment}
 */
{code}
}}
"""

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
        signatures |= set(get_signatures())
    else:
        for signature in get_signatures():
            sp = SignatureParser(signature)
            if sp.function_name in functions:
                signatures.add(signature)
        default_checks = False

    for n, signature in enumerate(signatures):
        for overload in ("Prim", "Rev", "Fwd"):
            sp = SignatureParser(signature)

            # skip ignored signatures if running the default checks
            if default_checks and sp.is_ignored():
                continue

            # skip signatures without inputs that can be eigen types
            if not sp.has_eigen_compatible_arg():
                continue

            # skip functions in forward mode with no forward mode overload, same for reverse mode
            if overload == "Fwd" and not sp.is_fwd_compatible():
                continue
            if overload == "Rev" and not sp.is_rev_compatible():
                continue

            is_reverse_mode = overload == "Rev" and not sp.returns_int()

            cg = CodeGenerator()

            # Generate two sets of arguments, one will be expressions one will not
            arg_list_expression_base = cg.build_arguments(sp, sp.number_arguments() * [overload], size = 1)
            arg_list_expression = [cg.convert_to_expression(arg, size = 1) if arg.is_eigen_compatible() else arg for arg in arg_list_expression_base]
            arg_list_no_expression = cg.build_arguments(sp, sp.number_arguments() * [overload], size = 1)

            # Check results with expressions and without are the same
            cpp_function_name = "stan::math::" + sp.function_name
            result = cg.function_call_assign(cpp_function_name, *arg_list_expression)
            result_no_expression = cg.function_call_assign(cpp_function_name, *arg_list_no_expression)
            cg.expect_eq(result, result_no_expression)

            # Check that expressions evaluated only once
            for arg in arg_list_expression:
                if arg.is_expression():
                    cg.expect_leq_one(arg.counter)

            # If reverse mode, check the adjoints
            if is_reverse_mode:
                summed_result = cg.recursive_sum(result)
                summed_result_no_expression = cg.recursive_sum(result_no_expression)
                sum_of_sums = cg.add(summed_result, summed_result_no_expression)
                cg.grad(sum_of_sums)

                for arg, arg_no_expression in zip(arg_list_no_expression, arg_list_expression):
                    cg.expect_adj_eq(arg, arg_no_expression)

                cg.recover_memory()

            tests.append(
                test_code_template.format(
                    overload=overload,
                    comment=signature.strip(),
                    test_name=sp.function_name + repr(n),
                    code = cg.cpp(),
                )
            )
    
    save_tests_in_files(j, tests)
