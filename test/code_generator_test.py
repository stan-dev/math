from signature_parser import SignatureParser
from code_generator import CodeGenerator
from statement_types import IntVariable, RealVariable, MatrixVariable
import pytest

@pytest.fixture
def add():
    return SignatureParser("add(real, vector) => vector")

@pytest.fixture
def int_var():
    return IntVariable("myint", 5)

@pytest.fixture
def real_var1():
    return RealVariable("Rev", "myreal1", 0.5)

@pytest.fixture
def real_var2():
    return RealVariable("Rev", "myreal2", 0.5)

@pytest.fixture
def matrix_var():
    return MatrixVariable("Rev", "mymatrix", "matrix", 0.5, 2)

@pytest.fixture
def cg():
    return CodeGenerator()

def test_prim_prim(add, cg):
    cg.build_arguments(add, ["Prim", "Prim"])
    assert cg.cpp() == """double real0 = 0.4;
auto matrix1 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, 1>>(0.4, 1);"""

def test_prim_rev(add, cg):
    cg.build_arguments(add, ["Prim", "Rev"])
    assert cg.cpp() == """double real0 = 0.4;
auto matrix1 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 1);"""

def test_rev_rev(add, cg):
    cg.build_arguments(add, ["Rev", "Rev"])
    assert cg.cpp() == """stan::math::var real0 = 0.4;
auto matrix1 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 1);"""

def test_size(add, cg):
    cg.build_arguments(add, ["Rev", "Rev"], 2)
    assert cg.cpp() == """stan::math::var real0 = 0.4;
auto matrix1 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 2);"""

def test_add(cg, real_var1, real_var2):
    cg.add(real_var1, real_var2)
    assert cg.cpp() == "auto sum_of_sums0 = stan::math::eval(stan::math::add(myreal1,myreal2));"

def test_convert_to_expression(cg, matrix_var):
    cg.convert_to_expression(matrix_var)
    assert cg.cpp() == """int mymatrix_expr0_counter = 0;
stan::test::counterOp<stan::math::var> mymatrix_expr0_counter_op(&mymatrix_expr0_counter);
auto mymatrix_expr0 = mymatrix.block(0,0,2,2).unaryExpr(mymatrix_expr0_counter_op);"""

def test_expect_adj_eq(cg, real_var1, real_var2):
    cg.expect_adj_eq(real_var1, real_var2)
    assert cg.cpp() == "stan::test::expect_adj_eq(myreal1,myreal2);"

def test_expect_eq(cg, real_var1, real_var2):
    cg.expect_eq(real_var1, real_var2)
    assert cg.cpp() == "EXPECT_STAN_EQ(myreal1,myreal2);"

def test_expect_leq_one(cg, int_var):
    cg.expect_leq_one(int_var)
    assert cg.cpp() == "EXPECT_LEQ_ONE(myint);"

def test_function_call_assign(cg, real_var1, real_var2):
    cg.function_call_assign("stan::math::add", real_var1, real_var2)
    assert cg.cpp() == "auto result0 = stan::math::eval(stan::math::add(myreal1,myreal2));"

def test_grad(cg, real_var1):
    cg.grad(real_var1)
    assert cg.cpp() == "stan::test::grad(myreal1);"

def test_recover_memory(cg):
    cg.recover_memory()
    assert cg.cpp() == "stan::math::recover_memory();"

def test_recursive_sum(cg, real_var1):
    cg.recursive_sum(real_var1)
    assert cg.cpp() == "auto summed_result0 = stan::math::eval(stan::test::recursive_sum(myreal1));"

def test_to_var_value(cg, matrix_var):
    cg.to_var_value(matrix_var)
    assert cg.cpp() == "auto mymatrix_varmat0 = stan::math::eval(stan::math::to_var_value(mymatrix));"
