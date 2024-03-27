#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform, stochastic_column_rt0) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> x(4, 4);
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    x(i) = static_cast<double>(i);
  }
  double lp = 0;
  Matrix<double, Dynamic, Dynamic> x_test
      = stan::math::stochastic_column_constrain<false>(x, lp);
  EXPECT_EQ(lp, 0.0);
  Matrix<double, Dynamic, Dynamic> x_res(5, 4);
  double lp_orig = 0.0;
  for (Eigen::Index i = 0; i < x.cols(); ++i) {
    x_res.col(i) = stan::math::simplex_constrain<false>(x.col(i), lp_orig);
  }
  EXPECT_EQ(lp_orig, 0.0);
  EXPECT_MATRIX_EQ(x_test, x_res);
  Matrix<double, Dynamic, Dynamic> x_lp_test
      = stan::math::stochastic_column_constrain<true>(x, lp);
  for (Eigen::Index i = 0; i < x.cols(); ++i) {
    x_res.col(i) = stan::math::simplex_constrain<true>(x.col(i), lp_orig);
  }
  EXPECT_MATRIX_EQ(x_lp_test, x_res);
}

TEST(prob_transform, stochastic_column_constrain_and_free) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> x(4, 4);
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    x(i) = static_cast<double>(i);
  }
  double lp = 0;
  Matrix<double, Dynamic, Dynamic> x_test = stan::math::stochastic_column_free(
      stan::math::stochastic_column_constrain<false>(x, lp));
  EXPECT_MATRIX_NEAR(x, x_test, 1e-9);

  Matrix<double, Dynamic, Dynamic> x_lp_test
      = stan::math::stochastic_column_free(
          stan::math::stochastic_column_constrain<true>(x, lp));
  EXPECT_MATRIX_NEAR(x, x_lp_test, 1e-9);
}
