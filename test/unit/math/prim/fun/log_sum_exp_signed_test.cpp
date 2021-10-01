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

  signs << 1, -1, -1, 1;

  EXPECT_FLOAT_EQ(log(dot_product(exp(m1), signs)),
                  log_sum_exp_signed(m1, signs));
}

TEST(MathFunctions, log_sum_exp_signed_matrix) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::log_sum_exp_signed;
  using stan::math::log_sum_exp;
  using stan::math::dot_product;
  using stan::math::log;
  using stan::math::exp;
  using stan::math::sum;
  using stan::math::elt_multiply;

  Eigen::MatrixXd m1(2, 2);
  m1 << 4, 6, 10, 5;

  Eigen::MatrixXi signs(2, 2);
  signs << 1, 1, 1, 1;

  EXPECT_FLOAT_EQ(log_sum_exp(m1), log_sum_exp_signed(m1, signs));

  signs << -1, -1, 1, 1;

  EXPECT_FLOAT_EQ(log(sum(elt_multiply(exp(m1), signs))),
                  log_sum_exp_signed(m1, signs));
}
