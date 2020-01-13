#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(prob_transform, lkj_cov_matrix_rt) {
  unsigned int K = 4;
  unsigned int K_choose_2 = 6;
  Matrix<double, Dynamic, 1> x(K_choose_2 + K);
  x << -1.0, 2.0, 0.0, 1.0, 3.0, -1.5, 1.0, 2.0, -1.5, 2.5;
  Matrix<double, Dynamic, Dynamic> y
      = stan::math::cov_matrix_constrain_lkj(x, K);
  Matrix<double, Dynamic, 1> xrt = stan::math::cov_matrix_free_lkj(y);
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(x[i], xrt[i]);
  }
}

TEST(prob_transform, lkj_cov_matrix_free_exception) {
  Matrix<double, Dynamic, Dynamic> y(0, 0);

  EXPECT_THROW(stan::math::cov_matrix_free_lkj(y), std::invalid_argument);
  y.resize(0, 10);
  EXPECT_THROW(stan::math::cov_matrix_free_lkj(y), std::invalid_argument);
  y.resize(10, 0);
  EXPECT_THROW(stan::math::cov_matrix_free_lkj(y), std::invalid_argument);
  y.resize(1, 2);
  EXPECT_THROW(stan::math::cov_matrix_free_lkj(y), std::invalid_argument);

  y.resize(2, 2);
  y << 0, 0, 0, 0;
  EXPECT_THROW(stan::math::cov_matrix_free_lkj(y), std::domain_error);
}
