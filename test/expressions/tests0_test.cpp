#include <test/expressions/expression_test_helpers.hpp>


TEST(ExpressionTestPrim, bad_wrong_value0) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> arg_mat0 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>();

  auto res_mat = stan::math::bad_wrong_value(arg_mat0);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> arg_expr0 = stan::test::make_arg<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>();
  int counter0 = 0;
  stan::test::counterOp<double> counter_op0(&counter0);

  auto res_expr = stan::math::bad_wrong_value(arg_expr0.unaryExpr(counter_op0));

  EXPECT_STAN_EQ(res_expr, res_mat);

  EXPECT_LE(counter0, 1);

}

TEST(ExpressionTestRev, bad_wrong_value0) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> arg_mat0 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>>();

  auto res_mat = stan::math::bad_wrong_value(arg_mat0);

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> arg_expr0 = stan::test::make_arg<Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>>();
  int counter0 = 0;
  stan::test::counterOp<stan::math::var> counter_op0(&counter0);

  auto res_expr = stan::math::bad_wrong_value(arg_expr0.unaryExpr(counter_op0));

  EXPECT_STAN_EQ(res_expr, res_mat);

  EXPECT_LE(counter0, 1);
  (stan::test::recursive_sum(res_mat) + stan::test::recursive_sum(res_expr)).grad();
  EXPECT_STAN_ADJ_EQ(arg_expr0,arg_mat0);

}

TEST(ExpressionTestFwd, bad_wrong_value0) {
  Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic> arg_mat0 = stan::test::make_arg<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic>>();

  auto res_mat = stan::math::bad_wrong_value(arg_mat0);

  Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic> arg_expr0 = stan::test::make_arg<Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic>>();
  int counter0 = 0;
  stan::test::counterOp<stan::math::fvar<double>> counter_op0(&counter0);

  auto res_expr = stan::math::bad_wrong_value(arg_expr0.unaryExpr(counter_op0));

  EXPECT_STAN_EQ(res_expr, res_mat);

  EXPECT_LE(counter0, 1);

}
