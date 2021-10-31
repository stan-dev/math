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
TEST(PrimSigsTests, {test_name}) {{
/*
{comment}
 */    
{code}
{check}
}}
"""


def main():
    """
    Generates tests to check primitive scalar functions.

    The tests check whether the primitive scalar functions are defined in the
    stan::math namespace and that no stan::math::var objects are created when
    they are used. If a var object is created, some of the input arguments
    were cast into a var and the reverse mode function was called.

    """

    tests = []

    default_checks = True
    signatures = set(get_signatures())

    for n, signature in enumerate(signatures):
        sp = SignatureParser(signature)

        # skip ignored signatures if running the default checks
        if default_checks and sp.is_ignored():
            continue

        # skip signatures without inputs that can be eigen types
        if (
            not sp.has_scalar_only_args()
            or sp.returns_complex()
            or sp.is_nullary()
            or sp.is_rng()
        ):
            continue

        cg = CodeGenerator()

        arg_list_no_expression = cg.build_arguments(
            sp, sp.number_arguments() * ["Prim"], size=1
        )

        # Check results with expressions and without are the same
        cpp_function_name = "stan::math::" + sp.function_name
        result_1 = cg.function_call_assign(cpp_function_name, *arg_list_no_expression)
        result_2 = cg.function_call_assign(cpp_function_name, *arg_list_no_expression)
        cg.expect_eq(result_1, result_2)

        var_check = (
            "EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);\n"
            + "EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(), 0);\n"
            + "EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);\n"
        )

        tests.append(
            test_code_template.format(
                overload="Prim",
                comment=signature.strip(),
                test_name=sp.function_name + repr(n),
                code=cg.cpp(),
                check=var_check,
            )
        )

    with open(src_folder + "prim_sig_test.cpp", "w") as out:
        out.write(
            "#include <stan/math.hpp>\n\n#include <test/expressions/expression_test_helpers.hpp>\n\n"
        )
        for test in tests:
            out.write(test)
