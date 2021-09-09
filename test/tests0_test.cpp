#include <test/expressions/expression_test_helpers.hpp>

TEST(ExpressionTestFwd, generalized_inverse23215) {
/*
generalized_inverse(matrix) => matrix
 */
auto matrix0 = stan::test::make_arg<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
int matrix0_expr1_counter = 0;
stan::test::counterOp<stan::math::fvar<double>> matrix0_expr1_counter_op(&matrix0_expr1_counter);
auto matrix0_expr1 = matrix0.block(0,0,1,1).unaryExpr(matrix0_expr1_counter_op);
auto matrix2 = stan::test::make_arg<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
auto result3 = stan::math::eval(stan::math::generalized_inverse(matrix0_expr1));
auto result4 = stan::math::eval(stan::math::generalized_inverse(matrix2));
EXPECT_STAN_EQ(result3,result4);
int int5 = 1;
EXPECT_LE(matrix0_expr1_counter,int5);
}

TEST(ExpressionTestPrim, generalized_inverse23215) {
/*
generalized_inverse(matrix) => matrix
 */
auto matrix0 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
int matrix0_expr1_counter = 0;
stan::test::counterOp<double> matrix0_expr1_counter_op(&matrix0_expr1_counter);
auto matrix0_expr1 = matrix0.block(0,0,1,1).unaryExpr(matrix0_expr1_counter_op);
auto matrix2 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
auto result3 = stan::math::eval(stan::math::generalized_inverse(matrix0_expr1));
auto result4 = stan::math::eval(stan::math::generalized_inverse(matrix2));
EXPECT_STAN_EQ(result3,result4);
int int5 = 1;
EXPECT_LE(matrix0_expr1_counter,int5);
}