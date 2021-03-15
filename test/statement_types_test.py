from statement_types import *
import pytest

@pytest.fixture
def int_var():
    return IntVariable("myint", 5)

@pytest.fixture
def real_var():
    return RealVariable("Rev", "myreal", 0.5)

@pytest.fixture
def matrix_var():
    return MatrixVariable("Rev", "mymatrix", "matrix", 0.5, 2)

def test_int(int_var):
    assert int_var.cpp() == "int myint = 5;"
    assert int_var.is_expression() == False
    assert int_var.is_matrix_like() == False
    assert int_var.is_reverse_mode() == False
    assert int_var.is_varmat_compatible() == False

def test_real(real_var):
    assert real_var.cpp() == "stan::math::var myreal = 0.5;"
    assert real_var.is_expression() == False
    assert real_var.is_matrix_like() == False
    assert real_var.is_reverse_mode() == True
    assert real_var.is_varmat_compatible() == False

def test_matrix(matrix_var):
    assert matrix_var.cpp() == "auto mymatrix = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>>(0.5, 2);"
    assert matrix_var.is_expression() == False
    assert matrix_var.is_matrix_like() == True
    assert matrix_var.is_reverse_mode() == True
    assert matrix_var.is_varmat_compatible() == True

def test_simplex():
    simplex_var = SimplexVariable("Prim", "mysimplex", 0.1, 3)
    assert simplex_var.cpp() == "auto mysimplex = stan::test::make_simplex<Eigen::Matrix<double, Eigen::Dynamic, 1>>(0.1, 3);"
    assert simplex_var.is_expression() == False
    assert simplex_var.is_matrix_like() == True
    assert simplex_var.is_reverse_mode() == False
    assert simplex_var.is_varmat_compatible() == True

def test_positive_definite_matrix():
    positive_definite_matrix_var = PositiveDefiniteMatrixVariable("Fwd", "mymatrix", 0.7, 2)
    assert positive_definite_matrix_var.cpp() == "auto mymatrix = stan::test::make_pos_definite_matrix<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic>>(0.7, 2);"
    assert positive_definite_matrix_var.is_expression() == False
    assert positive_definite_matrix_var.is_matrix_like() == True
    assert positive_definite_matrix_var.is_reverse_mode() == False
    assert positive_definite_matrix_var.is_varmat_compatible() == True

def test_algebra_solver_functor_variable():
    algebra_solver_functor_variable = AlgebraSolverFunctorVariable("myfunc")
    assert algebra_solver_functor_variable.cpp() == "stan::test::simple_eq_functor myfunc;"

def test_ode_functor_variable():
    ode_functor_variable = OdeFunctorVariable("myfunc")
    assert ode_functor_variable.cpp() == "stan::test::test_functor myfunc;"

def test_return_type_t_variable(int_var, real_var):
    return_type_t_variable = ReturnTypeTVariable("ret_t", int_var, real_var)
    assert return_type_t_variable.cpp() == "stan::return_type_t<decltype(myint),decltype(myreal)> ret_t = 0;"
    assert return_type_t_variable.is_reverse_mode() == True
    return_type_t_prim_variable = ReturnTypeTVariable("ret_t", int_var)
    assert return_type_t_prim_variable.is_reverse_mode() == False

def test_rng_variable():
    rng_variable = RngVariable("myrng")
    assert rng_variable.cpp() == "std::minstd_rand myrng;"

def test_ostream_variable():
    ostream_variable = OStreamVariable("myostream")
    assert ostream_variable.cpp() == "std::ostream* myostream = &std::cout;"

def test_array_variable(matrix_var):
    array_variable = ArrayVariable("Rev", "myarray", 3, "real", matrix_var, 2)
    assert array_variable.cpp() == "std::vector<std::vector<std::vector<decltype(mymatrix)>>> myarray = {{{mymatrix,mymatrix},{mymatrix,mymatrix}},{{mymatrix,mymatrix},{mymatrix,mymatrix}}};"
    assert array_variable.is_expression() == False
    assert array_variable.is_matrix_like() == False
    assert array_variable.is_reverse_mode() == True
    assert array_variable.is_varmat_compatible() == True

def test_array_variable_list_initializer():
    array_variable = ArrayVariable("Rev", "myarray", 2, "real", [(1, 2), [2, 2]], 2)
    assert array_variable.cpp() == "std::vector<std::vector<stan::math::var>> myarray = {{1, 2}, {2, 2}};"

def test_function_call_assign(int_var, real_var):
    function_call_assign = FunctionCallAssign("stan::math::add", "result", int_var, real_var)
    assert function_call_assign.cpp() == "auto result = stan::math::eval(stan::math::add(myint,myreal));"

def test_function_call(real_var):
    function_call = FunctionCall("stan::math::subtract", real_var, real_var)
    assert function_call.cpp() == "stan::math::subtract(myreal,myreal);"

def test_expression_variable(matrix_var):
    expression_variable = ExpressionVariable("myexpr", matrix_var, 2)
    assert expression_variable.cpp() == """int myexpr_counter = 0;
stan::test::counterOp<stan::math::var> myexpr_counter_op(&myexpr_counter);
auto myexpr = mymatrix.block(0,0,2,2).unaryExpr(myexpr_counter_op);"""
    assert expression_variable.is_expression() == True
    assert expression_variable.is_matrix_like() == True
    assert expression_variable.is_reverse_mode() == True
    assert expression_variable.is_varmat_compatible() == False