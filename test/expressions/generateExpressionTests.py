import numbers

from sig_utils import *

src_folder = "./test/expressions/"
build_folder = "./test/expressions/"

eigen_types = ["matrix", "vector", "row_vector"]

test_code_template = """
TEST(ExpressionTest{overload}, {function_name}{signature_number}) {{
{matrix_argument_declarations}
  auto res_mat = stan::math::{function_name}({matrix_argument_list});

{expression_argument_declarations}
  auto res_expr = stan::math::{function_name}({expression_argument_list});

  EXPECT_STAN_EQ(res_expr, res_mat);

{checks}
}}
"""


def make_arg_code(arg, scalar, var_name, var_number, function_name):
    """
    Makes code for declaration and initialization of an argument to function.

    Default argument range (between 0 and 1) works for most function, but some need
    values outside this range - these require special handling. Specific lambda is
    also hardcoded - if we use more functor arguments in future this may need to be
    extended.
    :param arg: stan lang type of the argument
    :param scalar: scalar type used in argument
    :param var_name: name of the variable to create
    :param var_number: number of variable in the function call
    :param function_name: name of the function that will be tested using this argument
    :return: code for declaration and initialization of an argument
    """
    if (
        function_name in non_differentiable_args
        and var_number in non_differentiable_args[function_name]
    ):
        scalar = "double"
    if arg == "(vector, vector, data real[], data int[]) => vector":
        return "  stan::test::simple_eq_functor " + var_name
    elif "=>" in arg:  # functors
        return_type = arg_types[arg.split(" => ")[1]].replace("SCALAR", scalar)
        return "  stan::test::test_functor " + var_name
    arg_type = arg_types[arg].replace("SCALAR", scalar)
    if (function_name in special_arg_values) and (arg != "rng"):
        if isinstance(special_arg_values[function_name][var_number], numbers.Number):
            return "  {} {} = stan::test::make_arg<{}>({})".format(
                arg_type,
                var_name,
                arg_type,
                special_arg_values[function_name][var_number],
            )
        elif isinstance(special_arg_values[function_name][var_number], str):
            return "  {} {} = stan::test::{}<{}>()".format(
                arg_type,
                var_name,
                special_arg_values[function_name][var_number],
                arg_type,
            )
    if (
        function_name in ("csr_to_dense_matrix", "csr_matrix_times_vector")
        and var_number == 4
    ):
        return "  {} {}{{1, 2}}".format(
            arg_type,
            var_name,
        )
    else:
        return "  %s %s = stan::test::make_arg<%s>()" % (
            arg_type,
            var_name,
            arg_type,
        )


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

    test_n = {}
    tests = []
    signatures = get_signatures()
    functions, extra_signatures = handle_function_list(functions)
    signatures += extra_signatures
    remaining_functions = set(functions)
    for signature in signatures:
        return_type, function_name, function_args = parse_signature(signature)
        # skip ignored signatures
        if (
            function_name in ignored
            and not functions
            and signature not in extra_signatures
        ):
            continue
        # skip default if we have list of function names/signatures to test
        if (
            (functions or extra_signatures)
            and function_name not in functions
            and signature not in extra_signatures
        ):
            continue
        # skip signatures without eigen inputs
        for arg2test in eigen_types:
            if arg2test in function_args:
                break
        else:
            continue
        if function_name in remaining_functions:
            remaining_functions.remove(function_name)
        func_test_n = test_n.get(function_name, 0)
        test_n[function_name] = func_test_n + 1

        if function_name.endswith("_rng"):
            function_args.append("rng")
        for overload, scalar in (
            ("Prim", "double"),
            ("Rev", "stan::math::var"),
            ("Fwd", "stan::math::fvar<double>"),
        ):
            if function_name.endswith("_rng") and overload != "Prim":
                continue
            if function_name in no_fwd_overload and overload == "Fwd":
                continue
            if function_name in no_rev_overload and overload == "Rev":
                continue
            mat_declarations = ""
            for n, arg in enumerate(function_args):
                mat_declarations += (
                    make_arg_code(arg, scalar, "arg_mat%d" % n, n, function_name)
                    + ";\n"
                )
            mat_arg_list = ", ".join("arg_mat%d" % n for n in range(len(function_args)))

            expression_declarations = ""
            for n, arg in enumerate(function_args):
                expression_declarations += (
                    make_arg_code(arg, scalar, "arg_expr%d" % n, n, function_name)
                    + ";\n"
                )
                if arg in eigen_types:
                    expression_declarations += "  int counter%d = 0;\n" % n
                    expression_declarations += (
                        "  stan::test::counterOp<%s> counter_op%d(&counter%d);\n"
                        % (scalar, n, n)
                    )
            expression_arg_list = ""
            for n, arg in enumerate(function_args[:-1]):
                if arg in eigen_types:
                    if arg == "matrix":
                        expression_arg_list += (
                            "arg_expr%d.block(0,0,1,1).unaryExpr(counter_op%d), "
                            % (
                                n,
                                n,
                            )
                        )
                    else:
                        expression_arg_list += (
                            "arg_expr%d.segment(0,1).unaryExpr(counter_op%d), "
                            % (
                                n,
                                n,
                            )
                        )
                else:
                    expression_arg_list += "arg_expr%d, " % n
            if function_args[-1] in eigen_types:
                if function_args[-1] == "matrix":
                    expression_arg_list += (
                        "arg_expr%d.block(0,0,1,1).unaryExpr(counter_op%d)"
                        % (
                            len(function_args) - 1,
                            len(function_args) - 1,
                        )
                    )
                else:
                    expression_arg_list += (
                        "arg_expr%d.segment(0,1).unaryExpr(counter_op%d)"
                        % (
                            len(function_args) - 1,
                            len(function_args) - 1,
                        )
                    )
            else:
                expression_arg_list += "arg_expr%d" % (len(function_args) - 1)
            checks = ""
            for n, arg in enumerate(function_args):
                if arg in eigen_types:
                    # besides evaluating its input rank also accesses one of the elements,
                    # resulting in counter being incremented twice.
                    if function_name == "rank":
                        checks += "  EXPECT_LE(counter%d, 2);\n" % n
                    else:
                        checks += "  EXPECT_LE(counter%d, 1);\n" % n
            if overload == "Rev" and (
                return_type.startswith("real")
                or return_type.startswith("vector")
                or return_type.startswith("row_vector")
                or return_type.startswith("matrix")
            ):
                checks += "  (stan::test::recursive_sum(res_mat) + stan::test::recursive_sum(res_expr)).grad();\n"
                for n, arg in enumerate(function_args):
                    # functors don't have adjoints to check
                    if "=>" in arg:
                        continue
                    checks += "  EXPECT_STAN_ADJ_EQ(arg_expr%d,arg_mat%d);\n" % (
                        n,
                        n,
                    )
            tests.append(
                test_code_template.format(
                    overload=overload,
                    function_name=function_name,
                    signature_number=func_test_n,
                    matrix_argument_declarations=mat_declarations,
                    matrix_argument_list=mat_arg_list,
                    expression_argument_declarations=expression_declarations,
                    expression_argument_list=expression_arg_list,
                    checks=checks,
                )
            )
    if remaining_functions:
        raise NameError("Functions not found: " + ", ".join(remaining_functions))
    save_tests_in_files(j, tests)
