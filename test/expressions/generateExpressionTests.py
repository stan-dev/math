import re
import os
import sys
import subprocess
from argparse import ArgumentParser, RawTextHelpFormatter

if os.name == "nt":  # Windows
    make = "mingw32-make"
else:
    make = "make"

src_folder = "./test/expressions/"
build_folder = "./test/expressions/"
exceptions_list_location = "./test/expressions/stan_math_sigs_exceptions.expected"

eigen_types = ["matrix", "vector", "row_vector"]
arg_types = {
    "int": "int",
    "int[]": "std::vector<int>",
    "int[,]": "std::vector<std::vector<int>>",
    "real": "SCALAR",
    "real[]": "std::vector<SCALAR>",
    "real[,]": "std::vector<std::vector<SCALAR>>",
    "vector": "Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>",
    "vector[]": "std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>>",
    "row_vector": "Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>",
    "row_vector[]": "std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>",
    "matrix": "Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>",
    "(vector, vector, data real[], data int[]) => vector": "auto",
    "rng": "std::minstd_rand",
}

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


def get_ignored_signatures():
    """
    Loads list of ignored signatures from the file listing the exceptions.
    :return: set of ignored signatures
    """
    part_sig = ""
    ignored = set()
    for signature in open(exceptions_list_location):
        signature = part_sig + signature
        part_sig = ""
        if not signature.endswith(")\n"):
            part_sig = signature
            continue
        ignored.add(signature)
    return ignored


def parse_signature_file(sig_file):
    """
    Parses signatures from a file of signatures
    :param sig_file: file-like object to pares
    :return: list of signatures
    """
    res = []
    part_sig = ""
    for signature in sig_file:
        signature = part_sig + signature
        part_sig = ""
        if not signature.endswith(")\n"):
            part_sig = signature
            continue
        res.append(signature)
    return res

def add_extra_signatures(res):
    """
    Adds signatures not defined in the stan language that can still accept
     expression types.
    """
    res.extend(["vector unit_vector_constrain(vector)",
        "vector unit_vector_constrain(row_vector)",
        "vector unit_vector_constrain(vector, real)",
        "vector unit_vector_constrain(row_vector, real)",
        "vector unit_vector_free(vector)",
        "vector unit_vector_free(row_vector)",
        "int is_cholesky_factor(matrix)",
        "int is_cholesky_factor_corr(matrix)",
        "int is_column_index(matrix, int)",
        "int is_column_index(vector, int)",
        "int is_corr_matrix(matrix)",
        "int is_cholesky_factor(matrix)",
        "int is_lower_triangular(matrix)",
        "int is_mat_finite(matrix)",
        "int is_mat_finite(vector)",
        "int is_matching_dims(matrix, matrix)",
        "int is_matching_dims(vector, matrix)",
        "int is_matching_dims(matrix, vector)",
        "int is_matching_dims(row_vector, matrix)",
        "int is_matching_dims(matrix, row_vector)",
        "int is_matching_dims(matrix, matrix)",
        "int is_matching_dims(row_vector, row_vector)",
        "int is_matching_dims(vector, row_vector)",
        "int is_matching_dims(row_vector, vector)",
        "int is_matching_dims(vector, vector)",
        "int is_pos_definite(matrix)",
        "int is_square(matrix)",
        "int is_square(vector)",
        "int is_square(row_vector)",
        "int is_symmetric(matrix)",
        "int is_unit_vector(vector)",
        "int is_unit_vector(row_matrix)"])
    return res

def get_signatures():
    """
    Retrieves function signatures from stanc3
    :return: list of signatures
    """
    if os.name == "nt":
        stanc3 = ".\\test\\expressions\\stanc.exe"
    else:
        stanc3 = "./test/expressions/stanc"
    p = subprocess.Popen((make, stanc3))
    if p.wait() != 0:
        sys.stderr.write("Error in making stanc3!")
        sys.exit(-1)

    p = subprocess.Popen(
        (stanc3 + " --dump-stan-math-signatures"),
        stdout=subprocess.PIPE,
        universal_newlines=True,
        shell=True,
    )
    res = parse_signature_file(p.stdout)
    if p.wait() != 0:
        sys.stderr.write("Error in getting signatures from stanc3!\n")
        sys.exit(-1)

    res_full = add_extra_signatures(res)
    return res_full


def parse_signature(signature):
    """
    Parses one signature
    :param signature: stanc3 function signature
    :return: return type, fucntion name and list of function argument types
    """
    return_type, rest = signature.split(" ", 1)
    function_name, rest = rest.split("(", 1)
    args = re.findall(r"(?:[(][^()]+[)][^,()]+)|(?:[^,()]+(?:,*[]])?)", rest)
    args = [i.strip() for i in args if i.strip()]
    return return_type, function_name, args


# list of function arguments that need special scalar values.
# None means to use the default argument value.
special_arg_values = {
    "acosh": [1.4],
    "log1m_exp": [-0.6],
    "categorical_log": [None, 1],
    "categorical_rng": [1, None],
    "categorical_lpmf": [None, 1],
    "dirichlet_log" : [1, None],
    "dirichlet_lpdf" : [1, None],
    "hmm_hidden_state_prob": [None, 1, 1],
    "hmm_latent_rng": [None, 1, 1, None],
    "hmm_marginal": [None, 1, 1],
    "lkj_corr_lpdf": [1, None],
    "lkj_corr_log": [1, None],
    "log_diff_exp": [3, None],
    "log_inv_logit_diff": [1.2, 0.4],
    "multinomial_log" : [None, 1],
    "multinomial_lpmf" : [None, 1],
    "multinomial_rng" : [1, None, None],
    "unit_vector_free" : [1.0],

}


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
    arg_type = arg_types[arg].replace("SCALAR", scalar)
    if arg == "(vector, vector, data real[], data int[]) => vector":
        return (
            "  %s %s = [](const auto& a, const auto&, const auto&, const auto&){return a;}"
            % (arg_type, var_name)
        )
    elif (
        function_name in special_arg_values
        and special_arg_values[function_name][var_number] is not None
    ):
        return "  %s %s = stan::test::make_arg<%s>(%f)" % (
            arg_type,
            var_name,
            arg_type,
            special_arg_values[function_name][var_number],
        )
    else:
        return "  %s %s = stan::test::make_arg<%s>()" % (arg_type, var_name, arg_type,)


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


def handle_function_list(functions_input, signatures):
    """
    Handles list of functions, splitting elements between functions and signatures.
    :param functions_input: This can contain names of functions
    already supported by stanc3, full function signatures or file names of files containing
    any of the previous two.
    :param signatures:
    :return:
    """
    function_names = []
    function_signatures = []
    for f in functions_input:
        if "." in f or "/" in f or "\\" in f:
            functions_input.extend(parse_signature_file(open(f)))
        elif " " in f:
            function_signatures.append(f)
            signatures.append(f)
        else:
            function_names.append(f)
    return function_names, function_signatures


# lists of functions that do not support fwd or rev autodiff
no_rev_overload = ["hmm_hidden_state_prob"]
no_fwd_overload = ["hmm_hidden_state_prob"]


def main(functions=(), j=1):
    """
    Generates expression tests. Functions that do not support expressions yet are listed
    in stan_math/tests/expressions/stan_math_sigs_exceptions.expected

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
    ignored = get_ignored_signatures()

    test_n = {}
    tests = []
    signatures = get_signatures()
    functions, extra_signatures = handle_function_list(functions, signatures)
    remaining_functions = set(functions)
    for signature in signatures:
        return_type, function_name, function_args = parse_signature(signature)
        # skip ignored signatures
        if signature in ignored and not functions and signature not in extra_signatures:
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
                    expression_arg_list += "arg_expr%d.unaryExpr(counter_op%d), " % (
                        n,
                        n,
                    )
                else:
                    expression_arg_list += "arg_expr%d, " % n
            if function_args[-1] in eigen_types:
                expression_arg_list += "arg_expr%d.unaryExpr(counter_op%d)" % (
                    len(function_args) - 1,
                    len(function_args) - 1,
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
                    if arg == "(vector, vector, data real[], data int[]) => vector":
                        continue
                    checks += "  EXPECT_STAN_ADJ_EQ(arg_expr%d,arg_mat%d);\n" % (n, n,)
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


def processCLIArgs():
    """
    Define and process the command line interface to the generateExpressionTests
    .py script.
    """
    parser = ArgumentParser(
        description="Generate and run stan math expression tests.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-j",
        metavar="N",
        type=int,
        default=1,
        help="Number of cores for make to use. Also number of files tests are split in.",
    )
    parser.add_argument(
        "functions",
        nargs="+",
        type=str,
        default=[],
        help="Names of the functions to test. By default tests all functions.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    cli_args = processCLIArgs()
    main(cli_args.functions, cli_args.j)
