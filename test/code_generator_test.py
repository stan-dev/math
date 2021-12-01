from signature_parser import SignatureParser
from code_generator import CodeGenerator
from statement_types import IntVariable, RealVariable, MatrixVariable
import unittest

class CodeGeneratorTest(unittest.TestCase):
    def setUp(self):
        self.add = SignatureParser("add(real, vector) => vector")
        self.int_var = IntVariable("myint", 5)
        self.real_var1 = RealVariable("Rev", "myreal1", 0.5)
        self.real_var2 = RealVariable("Rev", "myreal2", 0.5)
        self.matrix_var = MatrixVariable("Rev", "mymatrix", "matrix", 2, 0.5)
        self.cg = CodeGenerator()

    def test_prim_prim(self):
        self.cg.build_arguments(self.add, ["Prim", "Prim"], 1)
        self.assertEqual(self.cg.cpp(), """double real0 = 0.4;
auto matrix1 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, 1>>(0.4, 1);""")

    def test_prim_rev(self):
        self.cg.build_arguments(self.add, ["Prim", "Rev"], 1)
        self.assertEqual(self.cg.cpp(), """double real0 = 0.4;
auto matrix1 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 1);""")

    def test_rev_rev(self):
        self.cg.build_arguments(self.add, ["Rev", "Rev"], 1)
        self.assertEqual(self.cg.cpp(), """stan::math::var real0 = 0.4;
auto matrix1 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 1);""")

    def test_size(self):
        self.cg.build_arguments(self.add, ["Rev", "Rev"], 2)
        self.assertEqual(self.cg.cpp(), """stan::math::var real0 = 0.4;
auto matrix1 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 2);""")

    def test_add(self):
        self.cg.add(self.real_var1, self.real_var2)
        self.assertEqual(self.cg.cpp(), "auto sum_of_sums0 = stan::math::eval(stan::math::add(myreal1,myreal2));")

    def test_convert_to_expression(self):
        self.cg.convert_to_expression(self.matrix_var)
        self.assertEqual(self.cg.cpp(), """int mymatrix_expr0_counter = 0;
stan::test::counterOp<stan::math::var> mymatrix_expr0_counter_op(&mymatrix_expr0_counter);
auto mymatrix_expr0 = mymatrix.block(0,0,2,2).unaryExpr(mymatrix_expr0_counter_op);""")

    def test_expect_adj_eq(self):
        self.cg.expect_adj_eq(self.real_var1, self.real_var2)
        self.assertEqual(self.cg.cpp(), "stan::test::expect_adj_eq(myreal1,myreal2);")

    def test_expect_eq(self):
        self.cg.expect_eq(self.real_var1, self.real_var2)
        self.assertEqual(self.cg.cpp(), "EXPECT_STAN_EQ(myreal1,myreal2);")

    def test_expect_leq_one(self):
        self.cg.expect_leq_one(self.int_var)
        self.assertEqual(self.cg.cpp(), """int int0 = 1;
EXPECT_LE(myint,int0);""")

    def test_function_call_assign(self):
        self.cg.function_call_assign("stan::math::add", self.real_var1, self.real_var2)
        self.assertEqual(self.cg.cpp(), "auto result0 = stan::math::eval(stan::math::add(myreal1,myreal2));")

    def test_grad(self):
        self.cg.grad(self.real_var1)
        self.assertEqual(self.cg.cpp(), "stan::test::grad(myreal1);")

    def test_recover_memory(self):
        self.cg.recover_memory()
        self.assertEqual(self.cg.cpp(), "stan::math::recover_memory();")

    def test_recursive_sum(self):
        self.cg.recursive_sum(self.real_var1)
        self.assertEqual(self.cg.cpp(), "auto summed_result0 = stan::math::eval(stan::test::recursive_sum(myreal1));")

    def test_to_var_value(self):
        self.cg.to_var_value(self.matrix_var)
        self.assertEqual(self.cg.cpp(), "auto mymatrix_varmat0 = stan::math::eval(stan::math::to_var_value(mymatrix));")

if __name__ == '__main__':
    unittest.main()