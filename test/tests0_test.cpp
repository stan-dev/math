#include <test/expressions/expression_test_helpers.hpp>


TEST(ExpressionTestPrim, multi_normal_log1093) {
/*
multi_normal_log(vector, array[] vector, matrix) => real
 */
auto matrix0 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, 1>>(0.4, 1);
auto matrix1 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, 1>>(0.4, 1);
std::vector<decltype(matrix1)> array1 = {matrix1};
auto positive_definite_matrix2 = stan::test::make_pos_definite_matrix<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
int matrix0_expr3_counter = 0;
stan::test::counterOp<double> matrix0_expr3_counter_op(&matrix0_expr3_counter);
auto matrix0_expr3 = matrix0.segment(0,1).unaryExpr(matrix0_expr3_counter_op);
int positive_definite_matrix2_expr4_counter = 0;
stan::test::counterOp<double> positive_definite_matrix2_expr4_counter_op(&positive_definite_matrix2_expr4_counter);
auto positive_definite_matrix2_expr4 = positive_definite_matrix2.block(0,0,1,1).unaryExpr(positive_definite_matrix2_expr4_counter_op);
auto matrix5 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, 1>>(0.4, 1);
auto matrix6 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, 1>>(0.4, 1);
std::vector<decltype(matrix6)> array6 = {matrix6};
auto positive_definite_matrix7 = stan::test::make_pos_definite_matrix<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
auto result8 = stan::math::eval(stan::math::multi_normal_log(matrix0_expr3,array1,positive_definite_matrix2_expr4));
auto result9 = stan::math::eval(stan::math::multi_normal_log(matrix5,array6,positive_definite_matrix7));
EXPECT_STAN_EQ(result8,result9);
int int10 = 1;
EXPECT_LE(matrix0_expr3_counter,int10);
int int11 = 1;
EXPECT_LE(positive_definite_matrix2_expr4_counter,int11);
}

TEST(ExpressionTestRev, multi_normal_log1093) {
/*
multi_normal_log(vector, array[] vector, matrix) => real
 */
auto matrix0 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 1);
auto matrix1 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 1);
std::vector<decltype(matrix1)> array1 = {matrix1};
auto positive_definite_matrix2 = stan::test::make_pos_definite_matrix<Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
int matrix0_expr3_counter = 0;
stan::test::counterOp<stan::math::var> matrix0_expr3_counter_op(&matrix0_expr3_counter);
auto matrix0_expr3 = matrix0.segment(0,1).unaryExpr(matrix0_expr3_counter_op);
int positive_definite_matrix2_expr4_counter = 0;
stan::test::counterOp<stan::math::var> positive_definite_matrix2_expr4_counter_op(&positive_definite_matrix2_expr4_counter);
auto positive_definite_matrix2_expr4 = positive_definite_matrix2.block(0,0,1,1).unaryExpr(positive_definite_matrix2_expr4_counter_op);
auto matrix5 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 1);
auto matrix6 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(0.4, 1);
std::vector<decltype(matrix6)> array6 = {matrix6};
auto positive_definite_matrix7 = stan::test::make_pos_definite_matrix<Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
auto result8 = stan::math::eval(stan::math::multi_normal_log(matrix0_expr3,array1,positive_definite_matrix2_expr4));
auto result9 = stan::math::eval(stan::math::multi_normal_log(matrix5,array6,positive_definite_matrix7));
EXPECT_STAN_EQ(result8,result9);
int int10 = 1;
EXPECT_LE(matrix0_expr3_counter,int10);
int int11 = 1;
EXPECT_LE(positive_definite_matrix2_expr4_counter,int11);
auto summed_result12 = stan::math::eval(stan::test::recursive_sum(result8));
auto summed_result13 = stan::math::eval(stan::test::recursive_sum(result9));
auto sum_of_sums14 = stan::math::eval(stan::math::add(summed_result12,summed_result13));
stan::test::grad(sum_of_sums14);
stan::test::expect_adj_eq(matrix5,matrix0_expr3);
stan::test::expect_adj_eq(array6,array1);
stan::test::expect_adj_eq(positive_definite_matrix7,positive_definite_matrix2_expr4);
stan::math::recover_memory();
}

TEST(ExpressionTestFwd, multi_normal_log1093) {
/*
multi_normal_log(vector, array[] vector, matrix) => real
 */
auto matrix0 = stan::test::make_arg<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, 1>>(0.4, 1);
auto matrix1 = stan::test::make_arg<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, 1>>(0.4, 1);
std::vector<decltype(matrix1)> array1 = {matrix1};
auto positive_definite_matrix2 = stan::test::make_pos_definite_matrix<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
int matrix0_expr3_counter = 0;
stan::test::counterOp<stan::math::fvar<double>> matrix0_expr3_counter_op(&matrix0_expr3_counter);
auto matrix0_expr3 = matrix0.segment(0,1).unaryExpr(matrix0_expr3_counter_op);
int positive_definite_matrix2_expr4_counter = 0;
stan::test::counterOp<stan::math::fvar<double>> positive_definite_matrix2_expr4_counter_op(&positive_definite_matrix2_expr4_counter);
auto positive_definite_matrix2_expr4 = positive_definite_matrix2.block(0,0,1,1).unaryExpr(positive_definite_matrix2_expr4_counter_op);
auto matrix5 = stan::test::make_arg<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, 1>>(0.4, 1);
auto matrix6 = stan::test::make_arg<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, 1>>(0.4, 1);
std::vector<decltype(matrix6)> array6 = {matrix6};
auto positive_definite_matrix7 = stan::test::make_pos_definite_matrix<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic>>(0.4, 1);
auto result8 = stan::math::eval(stan::math::multi_normal_log(matrix0_expr3,array1,positive_definite_matrix2_expr4));
auto result9 = stan::math::eval(stan::math::multi_normal_log(matrix5,array6,positive_definite_matrix7));
EXPECT_STAN_EQ(result8,result9);
int int10 = 1;
EXPECT_LE(matrix0_expr3_counter,int10);
int int11 = 1;
EXPECT_LE(positive_definite_matrix2_expr4_counter,int11);
}
