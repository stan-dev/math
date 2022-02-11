from statement_types import *
import unittest

class StatementTypesTest(unittest.TestCase):
    def setUp(self):
        self.int_var = IntVariable("myint", 5)
        self.real_var = RealVariable("Rev", "myreal", 0.5)
        self.matrix_var = MatrixVariable("Rev", "mymatrix", "matrix", 2, 0.5)

    def test_int(self):
        self.assertEquals(self.int_var.cpp(), "int myint = 5;")
        self.assertEquals(self.int_var.is_expression(), False)
        self.assertEquals(self.int_var.is_eigen_compatible(), False)
        self.assertEquals(self.int_var.is_reverse_mode(), False)
        self.assertEquals(self.int_var.is_varmat_compatible(), False)

    def test_real(self):
        self.assertEquals(self.real_var.cpp(), "stan::math::var myreal = 0.5;")
        self.assertEquals(self.real_var.is_expression(), False)
        self.assertEquals(self.real_var.is_eigen_compatible(), False)
        self.assertEquals(self.real_var.is_reverse_mode(), True)
        self.assertEquals(self.real_var.is_varmat_compatible(), False)

    def test_matrix(self):
        self.assertEquals(self.matrix_var.cpp(), "auto mymatrix = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>>(0.5, 2);")
        self.assertEquals(self.matrix_var.is_expression(), False)
        self.assertEquals(self.matrix_var.is_eigen_compatible(), True)
        self.assertEquals(self.matrix_var.is_reverse_mode(), True)
        self.assertEquals(self.matrix_var.is_varmat_compatible(), True)

    def test_simplex(self):
        simplex_var = SimplexVariable("Prim", "mysimplex", 3, 0.1)
        self.assertEquals(simplex_var.cpp(), "auto mysimplex = stan::test::make_simplex<Eigen::Matrix<double, Eigen::Dynamic, 1>>(0.1, 3);")
        self.assertEquals(simplex_var.is_expression(), False)
        self.assertEquals(simplex_var.is_eigen_compatible(), True)
        self.assertEquals(simplex_var.is_reverse_mode(), False)
        self.assertEquals(simplex_var.is_varmat_compatible(), True)

    def test_positive_definite_matrix(self):
        positive_definite_matrix_var = PositiveDefiniteMatrixVariable("Fwd", "mymatrix", 2, 0.7)
        self.assertEquals(positive_definite_matrix_var.cpp(), "auto mymatrix = stan::test::make_pos_definite_matrix<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic>>(0.7, 2);")
        self.assertEquals(positive_definite_matrix_var.is_expression(), False)
        self.assertEquals(positive_definite_matrix_var.is_eigen_compatible(), True)
        self.assertEquals(positive_definite_matrix_var.is_reverse_mode(), False)
        self.assertEquals(positive_definite_matrix_var.is_varmat_compatible(), True)

    def test_algebra_solver_functor_variable(self):
        algebra_solver_functor_variable = AlgebraSolverFunctorVariable("myfunc")
        self.assertEquals(algebra_solver_functor_variable.cpp(), "stan::test::simple_eq_functor myfunc;")

    def test_ode_functor_variable(self):
        ode_functor_variable = OdeFunctorVariable("myfunc")
        self.assertEquals(ode_functor_variable.cpp(), "stan::test::test_functor myfunc;")

    def test_return_type_t_variable(self):
        return_type_t_variable = ReturnTypeTVariable("ret_t", self.int_var, self.real_var)
        self.assertEquals(return_type_t_variable.cpp(), "stan::return_type_t<decltype(myint),decltype(myreal)> ret_t = 0;")
        self.assertEquals(return_type_t_variable.is_reverse_mode(), True)
        return_type_t_prim_variable = ReturnTypeTVariable("ret_t", self.int_var)
        self.assertEquals(return_type_t_prim_variable.is_reverse_mode(), False)

    def test_rng_variable(self):
        rng_variable = RngVariable("myrng")
        self.assertEquals(rng_variable.cpp(), "std::minstd_rand myrng;")

    def test_ostream_variable(self):
        ostream_variable = OStreamVariable("myostream")
        self.assertEquals(ostream_variable.cpp(), "std::ostream* myostream = &std::cout;")

    def test_array_variable(self):
        array_variable = ArrayVariable("Rev", "myarray", 3, "real", 2, self.matrix_var)
        self.assertEquals(array_variable.cpp(), "std::vector<std::vector<std::vector<decltype(mymatrix)>>> myarray = {{{mymatrix,mymatrix},{mymatrix,mymatrix}},{{mymatrix,mymatrix},{mymatrix,mymatrix}}};")
        self.assertEquals(array_variable.is_expression(), False)
        self.assertEquals(array_variable.is_eigen_compatible(), False)
        self.assertEquals(array_variable.is_reverse_mode(), True)
        self.assertEquals(array_variable.is_varmat_compatible(), True)

    def test_array_variable_list_initializer(self):
        array_variable = ArrayVariable("Rev", "myarray", 1, "real", 1, (1, 2))
        self.assertEquals(array_variable.cpp(), "std::vector<stan::math::var> myarray = {1,2};")

    def test_function_call_assign(self):
        function_call_assign = FunctionCall("stan::math::add", "result", self.int_var, self.real_var)
        self.assertEquals(function_call_assign.cpp(), "auto result = stan::math::eval(stan::math::add(myint,myreal));")

    def test_function_call(self):
        function_call = FunctionCall("stan::math::subtract", None, self.real_var, self.real_var)
        self.assertEquals(function_call.cpp(), "stan::math::subtract(myreal,myreal);")

    def test_expression_variable(self):
        expression_variable = ExpressionVariable("myexpr", self.matrix_var, 2)
        self.assertEquals(expression_variable.cpp(), """int myexpr_counter = 0;
stan::test::counterOp<stan::math::var> myexpr_counter_op(&myexpr_counter);
auto myexpr = mymatrix.block(0,0,2,2).unaryExpr(myexpr_counter_op);""")
        self.assertEquals(expression_variable.is_expression(), True)
        self.assertEquals(expression_variable.is_eigen_compatible(), True)
        self.assertEquals(expression_variable.is_reverse_mode(), True)
        self.assertEquals(expression_variable.is_varmat_compatible(), False)

if __name__ == '__main__':
    unittest.main()