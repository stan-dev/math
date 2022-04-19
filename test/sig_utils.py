import itertools
import numbers
import os
import re
import subprocess
import sys

if os.name == "nt":  # Windows
    make = "mingw32-make"
    exe_extension = ".exe"
else:
    make = "make"
    exe_extension = ""

arg_types = {
    "int": "int",
    "array[] int": "std::vector<int>",
    "array[,] int": "std::vector<std::vector<int>>",
    "real": "SCALAR",
    "complex": "std::complex<SCALAR>",
    "array[] real": "std::vector<SCALAR>",
    "array[,] real": "std::vector<std::vector<SCALAR>>",
    "vector": "Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>",
    "array[] vector": "std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>>",
    "row_vector": "Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>",
    "array[] row_vector": "std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>",
    "matrix": "Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>",
    "complex_vector": "Eigen::Matrix<std::complex<SCALAR>, Eigen::Dynamic, 1>",
    "array[] complex_vector": "std::vector<Eigen::Matrix<std::complex<SCALAR>, Eigen::Dynamic, 1>>",
    "complex_row_vector": "Eigen::Matrix<std::complex<SCALAR>, 1, Eigen::Dynamic>",
    "array[] complex_row_vector": "std::vector<Eigen::Matrix<std::complex<SCALAR>, 1, Eigen::Dynamic>>",
    "complex_matrix": "Eigen::Matrix<std::complex<SCALAR>, Eigen::Dynamic, Eigen::Dynamic>",
}

scalar_stan_types = ("int", "real", "rng", "ostream_ptr")


def parse_array(stan_arg):
    """
    parses stan array type
    :param stan_arg: stan type, possibly an array
    :return: number of nested arrays, inner type
    """
    if stan_arg.startswith("array["):
        commas, inner_type = stan_arg.lstrip("array[").split("]")
        return len(commas) + 1, inner_type.strip()
    return 0, stan_arg.strip()


def get_cpp_type(stan_type):
    """
    Determines cpp type that implements given stan type.
    :param stan_type: stan type
    :return: cpp type
    """
    n_vec, inner_type = parse_array(stan_type)
    res = arg_types[inner_type]
    for i in range(n_vec):
        res = "std::vector<{}>".format(res)
    return res


simplex = "simplex"
pos_definite = "positive_definite_matrix"
scalar_return_type = "scalar_return_type"

make_special_arg_values = {
    simplex: "make_simplex",
    pos_definite: "make_pos_definite_matrix",
}

# list of function arguments that need special scalar values.
# None means to use the default argument value.
special_arg_values = {
    "acosh": [1.4],
    "algebra_solver": [None, None, None, None, None, None, None, 10],
    "algebra_solver_newton": [None, None, None, None, None, None, None, 10],
    "log1m_exp": [-0.6],
    "categorical_log": [None, simplex],
    "categorical_rng": [simplex, None],
    "categorical_lpmf": [None, simplex],
    "cholesky_decompose": [pos_definite, None],
    "csr_to_dense_matrix": [1, 1, None, [1], [1, 2]],
    "csr_matrix_times_vector": [1, 1, None, [1], [1, 2], None],
    "distance": [0.6, 0.4],
    "dirichlet_log": [simplex, None],
    "dirichlet_lpdf": [simplex, None],
    "hmm_hidden_state_prob": [None, simplex, simplex],
    "hmm_latent_rng": [None, simplex, simplex, None],
    "hmm_marginal": [None, simplex, simplex],
    "lambert_wm1": [-0.3],
    "lkj_corr_lpdf": [1, None],
    "lkj_corr_log": [1, None],
    "log_diff_exp": [3, None],
    "log_inv_logit_diff": [1.2, 0.4],
    "mdivide_left_spd": [pos_definite, None],
    "multinomial_log": [None, simplex],
    "multinomial_lpmf": [None, simplex],
    "multinomial_rng": [simplex, None, None],
    "multi_normal_lpdf": [None, None, pos_definite],
    "multi_normal_log": [None, None, pos_definite],
    "multi_normal_rng": [None, pos_definite],
    "multi_normal_prec_lpdf": [None, None, pos_definite],
    "multi_normal_prec_log": [None, None, pos_definite],
    "multi_normal_prec_rng": [None, pos_definite],
    "multi_student_t_lpdf": [None, None, None, pos_definite],
    "multi_student_t_log": [None, None, None, pos_definite],
    "multi_student_t_rng": [None, None, pos_definite],
    "ode_adams_tol": [None, None, 0.2, 0.4, None, None, 10, None, None, None],
    "ode_adams": [None, None, 0.2, 0.4, None, None, None],
    "ode_bdf_tol": [None, None, 0.2, 0.4, None, None, 10, None, None, None],
    "ode_bdf": [None, None, 0.2, 0.4, None, None, None],
    "ode_rk45_tol": [None, None, 0.2, 0.4, None, None, 10, None, None, None],
    "ode_rk45": [None, None, 0.2, 0.4, None, None, None],
    "pareto_cdf": [1.5, 0.7, None],
    "pareto_cdf_log": [1.5, 0.7, None],
    "pareto_lcdf": [1.5, 0.7, None],
    "pareto_type_2_cdf": [1.5, 0.7, None, None],
    "pareto_type_2_cdf_log": [1.5, 0.7, None, None],
    "pareto_type_2_lcdf": [1.5, 0.7, None, None],
    "positive_ordered_constrain": [None, scalar_return_type],
    "positive_ordered_free": [1.0],
    "ordered_constrain": [None, scalar_return_type],
    "ordered_free": [1.0],
    "simplex_constrain": [None, scalar_return_type],
    "simplex_free": [simplex],
    "student_t_cdf": [0.8, None, 0.4, None],
    "student_t_cdf_log": [0.8, None, 0.4, None],
    "student_t_ccdf_log": [0.8, None, 0.4, None],
    "student_t_lccdf": [0.8, None, 0.4, None],
    "student_t_lcdf": [0.8, None, 0.4, None],
    "unit_vector_constrain": [None, scalar_return_type],
    "unit_vector_free": [simplex],
    "uniform_cdf": [None, 0.2, 0.9],
    "uniform_ccdf_log": [None, 0.2, 0.9],
    "uniform_cdf_log": [None, 0.2, 0.9],
    "uniform_lccdf": [None, 0.2, 0.9],
    "uniform_lcdf": [None, 0.2, 0.9],
    "uniform_log": [None, 0.2, 0.9],
    "uniform_lpdf": [None, 0.2, 0.9],
    "uniform_rng": [0.2, 1.9, None],
    "wiener_log": [0.8, None, 0.4, None, None],
    "wiener_lpdf": [0.8, None, 0.4, None, None],
}

# list of functions we do not test. These are mainly functions implemented in compiler
# (not in Stan Math).
ignored = [
    "lchoose",
    "lmultiply",
    "assign_add",
    "assign_divide",
    "assign_elt_divide",
    "assign_elt_times",
    "assign_multiply",
    "assign_subtract",
    "if_else",
]

# list of function argument indices, for which real valued arguments are not differentiable
# - they need to be double even in autodiff overloads
non_differentiable_args = {
    "algebra_solver": [3],
    "algebra_solver_newton": [3],
    "ode_adams_tol": [4, 5, 6],
    "ode_bdf_tol": [4, 5, 6],
    "ode_rk45_tol": [4, 5, 6],
}

# lists of functions that do not support fwd or rev autodiff
no_rev_overload = ["hmm_hidden_state_prob", "quantile"]
no_fwd_overload = [
    "algebra_solver",
    "algebra_solver_newton",
    "hmm_hidden_state_prob",
    "map_rect",
    "ode_adams",
    "ode_adams_tol",
    "ode_bdf",
    "ode_bdf_tol",
    "ode_rk45",
    "ode_rk45_tol",
    "quantile",
]

internal_signatures = [
    "unit_vector_constrain(vector) => vector",
    "unit_vector_constrain(vector, real) => vector",
    "unit_vector_free(vector) => vector",
    "positive_ordered_constrain(vector) => vector",
    "positive_ordered_constrain(vector, real) => vector",
    "positive_ordered_free(vector) => vector",
    "ordered_constrain(vector) => vector",
    "ordered_constrain(vector, real) => vector",
    "ordered_free(vector) => vector",
    "simplex_constrain(vector) => vector",
    "simplex_constrain(vector, real) => vector",
    "simplex_free(vector) => vector",
    "is_cholesky_factor(matrix) => int",
    "is_cholesky_factor_corr(matrix) => int",
    "is_column_index(matrix, int) => int",
    "is_column_index(vector, int) => int",
    "is_corr_matrix(matrix) => int",
    "is_cholesky_factor(matrix) => int",
    "is_lower_triangular(matrix) => int",
    "is_mat_finite(matrix) => int",
    "is_mat_finite(vector) => int",
    "is_matching_dims(matrix, matrix) => int",
    "is_matching_dims(vector, matrix) => int",
    "is_matching_dims(matrix, vector) => int",
    "is_matching_dims(row_vector, matrix) => int",
    "is_matching_dims(matrix, row_vector) => int",
    "is_matching_dims(matrix, matrix) => int",
    "is_matching_dims(row_vector, row_vector) => int",
    "is_matching_dims(vector, row_vector) => int",
    "is_matching_dims(row_vector, vector) => int",
    "is_matching_dims(vector, vector) => int",
    "is_pos_definite(matrix) => int",
    "is_square(matrix) => int",
    "is_square(vector) => int",
    "is_square(row_vector) => int",
    "is_symmetric(matrix) => int",
    "is_unit_vector(vector) => int",
    # variadic functions: these are tested with one vector for variadic args
    "ode_adams((real, vector, ostream_ptr, vector) => vector, vector, real, array[] real, ostream_ptr, vector) => array[] vector",
    "ode_adams_tol((real, vector, ostream_ptr, vector) => vector, vector, real, array[] real, real, real, real, ostream_ptr, vector) => array[] vector",
    "ode_bdf((real, vector, ostream_ptr, vector) => vector, vector, real, array[] real, ostream_ptr, vector) => array[] vector",
    "ode_bdf_tol((real, vector, ostream_ptr, vector) => vector, vector, real, array[] real, real, real, real, ostream_ptr, vector) => array[] vector",
    "ode_rk45((real, vector, ostream_ptr, vector) => vector, vector, real, array[] real, ostream_ptr, vector) => array[] vector",
    "ode_rk45_tol((real, vector, ostream_ptr, vector) => vector, vector, real, array[] real, real, real, real, ostream_ptr, vector) => array[] vector",
    "reduce_sum(array[] real, int, vector) => real",
]


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
        if signature.endswith(",\n"):
            part_sig = signature
            continue
        res.append(signature)
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
    p = subprocess.Popen((make, stanc3.replace("\\", "/")))
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

    return res + internal_signatures


def parse_signature(signature):
    """
    Parses one signature
    :param signature: stanc3 function signature
    :return: return type, fucntion name and list of function argument types
    """
    rest, return_type = signature.rsplit("=>", 1)
    function_name, rest = rest.split("(", 1)
    args = re.findall(r"(?:[(][^()]+[)][^,()]+)|(?:[^,()]+(?:,*[]][^,()]+)?)", rest)
    #  regex parts:        ^^^^^^functor^^^^^^     ^^^^any other arg^^^^^^^
    args = [
        i.lstrip("data").strip() if "data" in i else i.strip()
        for i in args
        if i.strip()
    ]
    return return_type.strip(), function_name.strip(), args


def handle_function_list(functions_input):
    """
    Handles list of functions, splitting elements between functions and signatures.
    :param functions_input: This can contain names of functions
    already supported by stanc3, full function signatures or file names of files containing
    any of the previous two.
    :return:
    """
    function_names = []
    function_signatures = []
    for f in functions_input:
        if ("." in f) or ("/" in f) or ("\\" in f):
            with open(f) as fh:
                functions_input.extend(parse_signature_file(fh))
        elif "(" in f:
            function_signatures.append(f)
        else:
            function_names.append(f)
    return function_names, function_signatures


def reference_vector_argument(arg):
    """
    Determines a reference argument, so as not to duplicate arrays of reals, vectors and row vectors,
    which usually have the same implementation.
    :param arg: argument
    :return: reference argument
    """
    if arg in ("array[] real", "row_vector"):
        return "vector"
    return arg


overload_scalar = {
    "Prim": "double",
    "Rev": "stan::math::var",
    "RevVarmat": "stan::math::var",
    "Fwd": "stan::math::fvar<double>",
    "Mix": "stan::math::fvar<stan::math::var>",
}
