#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform, row_stochastic_rt0) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> x(4, 4);
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    x(i) = static_cast<double>(i);
  }
  double lp = 0;
  Matrix<double, Dynamic, Dynamic> x_test
      = stan::math::stochastic_row_constrain<false>(x, lp);
  EXPECT_EQ(lp, 0.0);
  Matrix<double, Dynamic, Dynamic> x_res(4, 5);
  double lp_orig = 0.0;
  for (Eigen::Index i = 0; i < x.cols(); ++i) {
    x_res.row(i) = stan::math::simplex_constrain<false>(x.row(i), lp_orig);
  }
  EXPECT_EQ(lp_orig, 0.0);
  EXPECT_MATRIX_EQ(x_test, x_res);
  Matrix<double, Dynamic, Dynamic> x_lp_test
      = stan::math::stochastic_row_constrain<true>(x, lp);
  for (Eigen::Index i = 0; i < x.cols(); ++i) {
    x_res.row(i) = stan::math::simplex_constrain<true>(x.row(i), lp_orig);
  }
  EXPECT_MATRIX_EQ(x_lp_test, x_res);
}

TEST(prob_transform, row_stochastic_constrain_free) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> x(4, 4);
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    x(i) = static_cast<double>(i);
  }
  double lp = 0;
  Matrix<double, Dynamic, Dynamic> x_test = stan::math::stochastic_row_free(
      stan::math::stochastic_row_constrain<false>(x, lp));
  EXPECT_MATRIX_NEAR(x, x_test, 1e-9);

  Matrix<double, Dynamic, Dynamic> x_lp_test = stan::math::stochastic_row_free(
      stan::math::stochastic_row_constrain<true>(x, lp));
  EXPECT_MATRIX_NEAR(x, x_lp_test, 1e-9);
}
