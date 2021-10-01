#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, log_sum_exp_signed) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::log_sum_exp_signed;
  using stan::math::log_sum_exp;
  using stan::math::dot_product;
  using stan::math::log;
  using stan::math::exp;

  Eigen::VectorXd m1(4);
  m1 << 1, 2, 3, 4;

  Eigen::VectorXi signs(4);
  signs << 1, 1, 1, 1;
  
  EXPECT_FLOAT_EQ(log_sum_exp(m1), log_sum_exp_signed(m1, signs));
}
